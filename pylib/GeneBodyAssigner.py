"""Gene-body read counting: credit every read to a gene, intronic reads included.

Single-nucleus libraries are mostly unspliced pre-mRNA. Isoform quantification
assigns a read only when it is structurally compatible with an annotated isoform
(config["fraction_read_align_overlap"] of the read on the isoform's exons), so a read
lying in an intron is assigned to nothing and drops out of the gene counts too. This
module gives such a read to the single most appropriate gene, and labels every read
spliced / unspliced / ambiguous (S / U / A) against the gene it is credited to.

It generalizes the `--oversimplify` best-overlap assigner (LRAA
_run_oversimplify_best_overlap), which scores exonic overlap only, per contig, and
replaces EM. This module works per gene (exon union + gene span), scores intronic
overlap too, and is additive: reads that isoform quantification assigns keep that
assignment and its EM split, and only the reads it leaves unassigned are placed here.

Pure geometry, no pysam or splice graph, so it is unit-testable on its own. The
caller (StreamingQuant.stream_assign via GeneBodyRecorder) supplies read segments.

Coordinates are 1-based inclusive throughout, as in Transcript.get_exon_segments().

Decision order for a read with no isoform assignment (all candidates are on the
read's transcribed strand; this runs inside one contig-strand pass, so an antisense
gene is never a candidate):

  0. candidate: the gene span (first exon start .. last exon end over all isoforms)
     holds >= gene_body_min_read_frac_in_span of the read's aligned bases.
  1. most read junctions that match an annotated intron of the gene;
  2. most exonic bases (overlap with the gene's exon union); exonic overlap below
     gene_body_min_feature_bp counts as none, so a read end grazing an exon does not
     outrank a gene whose intron holds the read;
  3. the read lies wholly inside the gene span;
  4. the read's 3' end (strand-aware) lies inside the gene span;
  5. split among the tied genes in proportion to their isoform-assigned reads in this
     contig-strand unit (full library; config gene_body_tie_weighting
     "isoform_reads", or "none" to skip). A gene with no isoform-assigned read gets
     none. This is what keeps an annotated readthrough gene (MRPS28-TPD52,
     ENSG00000276418) from taking the intronic reads of the gene it overlaps: in the
     ALS test window it has 0 isoform reads against TPD52's 1,536;
  6. no tied gene has isoform reads: the smaller gene span (innermost gene);
  7. still tied: ambiguous. gene_body_ambiguous_policy "drop" records the read
     against every tied gene with frac 0 (counted nowhere); "split" gives 1/n each.

Tiers 5-7 need the unit's isoform counts, which exist only at the end of the
stream, so a tied read is spilled to a temp file and resolved by
GeneBodyRecorder.finish(); the sidecar is not in read order.

Splice class, against the gene the read is credited to:
  U  the read covers >= gene_body_min_feature_bp of that gene's intronic sequence
     (gene span minus exon union), spliced or not -- pre-mRNA or partially spliced;
  S  otherwise, and the read carries at least one junction;
  A  otherwise: unspliced and exonic, compatible with mature and pre-mRNA alike.
"""

import bisect
import tempfile
from collections import defaultdict

import LRAA_Globals

SPLICED = "S"
UNSPLICED = "U"
AMBIGUOUS = "A"

POLICY_DROP = "drop"
POLICY_SPLIT = "split"


def _merge_intervals(intervals):
    merged = []
    for lend, rend in sorted(intervals):
        if merged and lend <= merged[-1][1] + 1:
            if rend > merged[-1][1]:
                merged[-1][1] = rend
        else:
            merged.append([lend, rend])
    return [tuple(x) for x in merged]


def _overlap_bp(segments, intervals):
    """Bases of `segments` inside `intervals`; both sorted, non-overlapping, inclusive."""
    ov = 0
    i = j = 0
    while i < len(segments) and j < len(intervals):
        a1, a2 = segments[i]
        b1, b2 = intervals[j]
        if a2 < b1:
            i += 1
        elif b2 < a1:
            j += 1
        else:
            ov += min(a2, b2) - max(a1, b1) + 1
            if a2 < b2:
                i += 1
            else:
                j += 1
    return ov


def segments_from_cigar(reference_start, cigartuples):
    """Aligned blocks and introns of one alignment, 1-based inclusive.

    Blocks break only at 'N' (introns); deletions and insertions stay inside a block,
    unlike pysam's get_blocks(), which splits at both.
    """
    segments = []
    introns = []
    pos = reference_start  # 0-based
    block_start = None
    for op, length in cigartuples:
        if op in (0, 7, 8, 2):  # M = X D consume reference, inside a block
            if block_start is None:
                block_start = pos
            pos += length
        elif op == 3:  # N
            if block_start is not None:
                segments.append((block_start + 1, pos))
                block_start = None
            introns.append((pos + 1, pos + length))
            pos += length
        # I S H P consume no reference
    if block_start is not None:
        segments.append((block_start + 1, pos))
    # A leading or trailing N is not an intron between two aligned blocks.
    if segments:
        introns = [x for x in introns if x[0] > segments[0][0] and x[1] < segments[-1][1]]
    else:
        introns = []
    return segments, introns


class GeneModel:
    __slots__ = ("gene_id", "lend", "rend", "exons", "introns", "intronic")

    def __init__(self, gene_id, exon_segment_lists):
        self.gene_id = gene_id
        all_exons = [tuple(e) for segs in exon_segment_lists for e in segs]
        self.exons = _merge_intervals(all_exons)
        self.lend = self.exons[0][0]
        self.rend = self.exons[-1][1]
        introns = set()
        for segs in exon_segment_lists:
            segs = sorted(tuple(e) for e in segs)
            for k in range(1, len(segs)):
                introns.add((segs[k - 1][1] + 1, segs[k][0] - 1))
        self.introns = introns
        # gene span minus exon union
        self.intronic = [
            (self.exons[k - 1][1] + 1, self.exons[k][0] - 1)
            for k in range(1, len(self.exons))
        ]

    def span_len(self):
        return self.rend - self.lend + 1


class ReadGeneEvidence:
    __slots__ = (
        "gene",
        "span_bp",
        "exonic_bp",
        "intronic_bp",
        "shared_junctions",
        "contained",
        "three_prime_in_span",
    )


class GeneBodyCall:
    """Genes left after tiers 0-4: [(gene_id, splice_class, evidence)], sorted by id,
    and a parallel list flagging the innermost (smallest-span) ones."""

    __slots__ = ("genes", "innermost")

    def __init__(self, genes, innermost):
        self.genes = genes
        self.innermost = innermost


class GeneBodyIndex:
    """Gene models of one contig-strand, built from the annotation's transcripts."""

    def __init__(self, transcripts, contig_strand):
        by_gene = defaultdict(list)
        for t in transcripts or []:
            if contig_strand and t.get_strand() != contig_strand:
                continue
            by_gene[t.get_output_gene_id()].append(t.get_exon_segments())
        self.contig_strand = contig_strand
        self.genes = {gid: GeneModel(gid, segs) for gid, segs in by_gene.items()}
        self._sorted = sorted(self.genes.values(), key=lambda g: (g.lend, g.rend, g.gene_id))
        self._lends = [g.lend for g in self._sorted]
        # running max of rend lets an overlap query stop scanning left early
        self._max_rend_prefix = []
        m = 0
        for g in self._sorted:
            m = max(m, g.rend)
            self._max_rend_prefix.append(m)

    def __len__(self):
        return len(self.genes)

    def overlapping_genes(self, lend, rend):
        hits = []
        k = bisect.bisect_right(self._lends, rend) - 1
        while k >= 0:
            if self._max_rend_prefix[k] < lend:
                break
            g = self._sorted[k]
            if g.rend >= lend:
                hits.append(g)
            k -= 1
        return hits

    # ------------------------------------------------------------------ settings
    @staticmethod
    def _settings():
        cfg = LRAA_Globals.config
        return (
            float(cfg.get("gene_body_min_read_frac_in_span", 0.5)),
            int(cfg.get("gene_body_min_feature_bp", 20)),
            str(cfg.get("gene_body_ambiguous_policy", POLICY_DROP)),
        )

    # ------------------------------------------------------------------ evidence
    def evidence(self, gene, segments, read_introns, strand):
        ev = ReadGeneEvidence()
        ev.gene = gene
        ev.span_bp = _overlap_bp(segments, [(gene.lend, gene.rend)])
        ev.exonic_bp = _overlap_bp(segments, gene.exons)
        ev.intronic_bp = _overlap_bp(segments, gene.intronic)
        ev.shared_junctions = sum(1 for x in read_introns if x in gene.introns)
        r_lend, r_rend = segments[0][0], segments[-1][1]
        ev.contained = gene.lend <= r_lend and r_rend <= gene.rend
        three_prime = r_rend if strand != "-" else r_lend
        ev.three_prime_in_span = gene.lend <= three_prime <= gene.rend
        return ev

    @staticmethod
    def splice_class(ev, read_introns, min_feature_bp):
        if ev is not None and ev.intronic_bp >= min_feature_bp:
            return UNSPLICED
        if read_introns:
            return SPLICED
        return AMBIGUOUS

    def classify_for_gene(self, gene_id, segments, read_introns, strand=None):
        """Splice class of a read already credited to gene_id (isoform-assigned reads)."""
        _, min_bp, _ = self._settings()
        gene = self.genes.get(gene_id)
        ev = self.evidence(gene, segments, read_introns, strand or self.contig_strand) if gene else None
        return self.splice_class(ev, read_introns, min_bp), ev

    # ------------------------------------------------------------------ choice
    def assign(self, segments, read_introns, strand=None):
        """Candidate gene(s) for a read isoform quantification did not assign, or None.

        Applies tiers 0-4. More than one gene left is a tie for resolve() to settle;
        `innermost` marks the tied genes with the smallest span (tier 6).
        """
        if not segments:
            return None
        strand = strand or self.contig_strand
        min_frac, min_bp, _ = self._settings()
        aligned_bp = sum(r - l + 1 for l, r in segments)
        if aligned_bp <= 0:
            return None
        cands = []
        for gene in self.overlapping_genes(segments[0][0], segments[-1][1]):
            ev = self.evidence(gene, segments, read_introns, strand)
            if ev.span_bp / aligned_bp >= min_frac:
                cands.append(ev)
        if not cands:
            return None

        def key(ev):
            return (
                ev.shared_junctions,
                ev.exonic_bp if ev.exonic_bp >= min_bp else 0,
                1 if ev.contained else 0,
                1 if ev.three_prime_in_span else 0,
            )

        best = max(key(ev) for ev in cands)
        tied = sorted((ev for ev in cands if key(ev) == best), key=lambda e: e.gene.gene_id)
        smallest = min(ev.gene.span_len() for ev in tied)
        return GeneBodyCall(
            [
                (ev.gene.gene_id, self.splice_class(ev, read_introns, min_bp), ev)
                for ev in tied
            ],
            [ev.gene.span_len() == smallest for ev in tied],
        )

    @staticmethod
    def resolve(call, gene_weights=None):
        """Fractions for a call's genes: [(gene_id, frac, splice_class, evidence)].

        One gene: 1. Several (tier 5): split in proportion to gene_weights (the
        genes' isoform-assigned reads), genes with no weight getting none; with no
        weight on any (or no weights given), the innermost gene(s) (tier 6); still
        several: gene_body_ambiguous_policy, "drop" (0 each) or "split" (1/n).
        """
        genes = call.genes
        if len(genes) == 1:
            g, cls, ev = genes[0]
            return [(g, 1.0, cls, ev)], False
        if gene_weights is not None:
            w = [float(gene_weights.get(g, 0.0)) for g, _, _ in genes]
            tot = sum(w)
            if tot > 0:
                return [
                    (g, wi / tot, cls, ev)
                    for (g, cls, ev), wi in zip(genes, w)
                    if wi > 0
                ], False
        inner = [x for x, keep in zip(genes, call.innermost) if keep]
        if len(inner) == 1:
            g, cls, ev = inner[0]
            return [(g, 1.0, cls, ev)], False
        policy = str(LRAA_Globals.config.get("gene_body_ambiguous_policy", POLICY_DROP))
        frac = 1.0 / len(inner) if policy == POLICY_SPLIT else 0.0
        return [(g, frac, cls, ev) for g, cls, ev in inner], True


GENE_BODY_TRACKING_HEADER = [
    "gene_id",
    "read_name",
    "frac_assigned",
    "splice_class",
    "source",
    "n_genes",
    "aligned_bp",
    "exonic_bp",
    "intronic_bp",
]


def gene_body_tracking_suffix():
    return ".quant.gene_body.tracking.gz"


class _Ev:
    """Evidence fields a spilled read needs back: what the sidecar row prints."""

    __slots__ = ("exonic_bp", "intronic_bp")

    def __init__(self, exonic_bp, intronic_bp):
        self.exonic_bp = exonic_bp
        self.intronic_bp = intronic_bp


class GeneBodyRecorder:
    """Writes one sidecar row per (read, gene) for every retained read of a contig-strand.

    source = "isoform": the read's isoform assignment (EM split summed per gene).
    source = "gene_body:<why>": placed by GeneBodyIndex; <why> is the reason isoform
             assignment failed (no_path, spacer_path, unassignable).
    Reads with no candidate gene are counted, not written.

    A read tied between genes after tiers 0-4 is spilled to a temp file and resolved by
    finish(), once the whole stream has been seen: tier 5 weights the tied genes by
    their isoform-assigned reads over the FULL library of this unit, which only exists
    at the end. Memory stays bounded by the gene count, whatever the tie count.
    """

    def __init__(self, index, fh, spill_dir=None):
        self.index = index
        self.fh = fh
        self.counts = defaultdict(float)
        self.gene_isoform_reads = defaultdict(float)
        self._spill = None
        self._spill_dir = spill_dir

    def _write(self, gene_id, read_name, frac, cls, source, n, aligned_bp, ev):
        self.fh.write(
            "{}\t{}\t{:.6f}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                gene_id,
                read_name,
                frac,
                cls,
                source,
                n,
                aligned_bp,
                ev.exonic_bp if ev is not None else "",
                ev.intronic_bp if ev is not None else "",
            )
        )

    def record_isoform(self, read_name, segments, read_introns, rows):
        """rows: (gene_id, transcript_id, ..., frac, ...) as written to quant.tracking."""
        by_gene = defaultdict(float)
        for row in rows:
            by_gene[row[0]] += float(row[5])
        aligned_bp = sum(r - l + 1 for l, r in segments)
        n = len(by_gene)
        for gene_id in sorted(by_gene):
            cls, ev = self.index.classify_for_gene(gene_id, segments, read_introns)
            self._write(gene_id, read_name, by_gene[gene_id], cls, "isoform", n, aligned_bp, ev)
            self.counts["isoform_" + cls] += by_gene[gene_id]
            self.gene_isoform_reads[gene_id] += by_gene[gene_id]
        self.counts["reads_isoform"] += 1

    def _emit(self, read_name, source, aligned_bp, resolved, ambiguous):
        n = len(resolved)
        for gene_id, frac, cls, ev in resolved:
            self._write(gene_id, read_name, frac, cls, source, n, aligned_bp, ev)
            self.counts["gene_body_" + cls] += frac
        self.counts["reads_gene_body_ambiguous" if ambiguous else "reads_gene_body"] += 1

    def record_unassigned(self, read_name, segments, read_introns, reason):
        call = self.index.assign(segments, read_introns)
        if call is None:
            self.counts["reads_no_gene"] += 1
            return
        aligned_bp = sum(r - l + 1 for l, r in segments)
        source = "gene_body:" + reason
        if len(call.genes) == 1:
            resolved, ambiguous = GeneBodyIndex.resolve(call)
            self._emit(read_name, source, aligned_bp, resolved, ambiguous)
            return
        self.counts["reads_gene_tied"] += 1
        if self._spill is None:
            self._spill = tempfile.TemporaryFile(mode="w+t", dir=self._spill_dir)
        fields = [read_name, source, str(aligned_bp)]
        for (g, cls, ev), inner in zip(call.genes, call.innermost):
            fields += [g, cls, str(ev.exonic_bp), str(ev.intronic_bp), "1" if inner else "0"]
        self._spill.write("\t".join(fields) + "\n")

    def finish(self):
        """Resolve the spilled ties with this unit's isoform-assigned gene counts."""
        if self._spill is None:
            return
        use_weights = (
            str(LRAA_Globals.config.get("gene_body_tie_weighting", "isoform_reads"))
            == "isoform_reads"
        )
        weights = self.gene_isoform_reads if use_weights else None
        self._spill.seek(0)
        for line in self._spill:
            f = line.rstrip("\n").split("\t")
            read_name, source, aligned_bp = f[0], f[1], f[2]
            genes, innermost = [], []
            for k in range(3, len(f), 5):
                genes.append((f[k], f[k + 1], _Ev(f[k + 2], f[k + 3])))
                innermost.append(f[k + 4] == "1")
            resolved, ambiguous = GeneBodyIndex.resolve(GeneBodyCall(genes, innermost), weights)
            self._emit(read_name, source, aligned_bp, resolved, ambiguous)
        self._spill.close()
        self._spill = None
