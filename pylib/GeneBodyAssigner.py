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
callers (StreamingQuant.stream_assign and the --oversimplify path, via
GeneBodyRecorder) supply read segments.

Coordinates are 1-based inclusive throughout, as in Transcript.get_exon_segments().

Strand: every candidate is on the read's transcribed strand (each contig-strand pass
sees only its own reads and genes), so an antisense gene is never a candidate. An
unspliced read has no `ts` evidence and takes its aligned strand: inputs are expected
to be strand-specific (oriented 10x / Kinnex reads); an unstranded library would put
half its unspliced reads on the wrong strand. No special handling, but GeneBodyRecorder
logs a warning when the spliced reads show it (STRAND_CHECK_*).

Decision order for a read with no isoform assignment:

  0. candidate: the gene span (first exon start .. last exon end over all isoforms)
     holds >= gene_body_min_read_frac_in_span of the read's aligned bases. With
     gene_body_exclude_readthrough_genes (default on), a gene whose transcripts all
     carry the GTF tag "readthrough_transcript" (or whose gene line carries
     "readthrough_gene") is never a candidate; its isoform-assigned reads still count;
  1. most read junctions that match an annotated intron of the gene;
  2. most exonic bases (overlap with the gene's exon union); exonic overlap below
     gene_body_min_feature_bp counts as none, so a read end grazing an exon does not
     outrank a gene whose intron holds the read;
  3. the read lies wholly inside the gene span;
  4. weight by isoform-assigned reads in this contig-strand unit (full library;
     gene_body_tie_weighting "isoform_reads", or "none" to skip): a gene holding
     >= gene_body_dominance_frac (0.9) of the tied genes' weight takes the whole read;
     otherwise the read is split in proportion, genes with no weight getting none;
  5. no tied gene has isoform reads: the read's 3' end (strand-aware) inside the span.
     For internally primed pre-mRNA the 3' end is a priming site, so it comes after
     the expression evidence;
  6. the smaller gene span (innermost gene);
  7. still tied: ambiguous. gene_body_ambiguous_policy "drop" records the read
     against every tied gene with frac 0 (counted nowhere); "split" gives 1/n each.

Rules 4-7 need the unit's isoform counts, which exist only at the end of the
stream, so a tied read is spilled to a temp file and resolved by
GeneBodyRecorder.finish(); the sidecar is not in read order.

Splice class, against the gene the read is credited to:
  U  the read covers >= gene_body_min_feature_bp of that gene's intronic sequence
     (gene span minus exon union), spliced or not -- pre-mRNA or partially spliced;
  S  otherwise, and the read carries at least one junction;
  A  otherwise: unspliced and exonic, compatible with mature and pre-mRNA alike.
gene_full = S + U + A.
"""

import bisect
import logging
import tempfile
from collections import defaultdict

import LRAA_Globals

logger = logging.getLogger(__name__)

SPLICED = "S"
UNSPLICED = "U"
AMBIGUOUS = "A"

POLICY_DROP = "drop"
POLICY_SPLIT = "split"

READTHROUGH_TRANSCRIPT_TAG = "readthrough_transcript"
READTHROUGH_GENE_TAG = "readthrough_gene"

# Strandedness check. A spliced read's transcribed strand is known from its junctions,
# so in an oriented library its aligned strand agrees with the unit's strand nearly
# always. Unstranded cDNA gives about half antisense, and then every unspliced read,
# which can only take its aligned strand, lands on the wrong strand half the time.
# Each unit logs a warning and also writes its counts to the gene-body sidecar as one
# comment line (STRAND_CHECK_PREFIX, then unit, spliced reads, antisense reads); the
# merges keep these lines, and summarize_strand_check() totals them for the run.
STRAND_CHECK_MIN_SPLICED_READS = 200
STRAND_CHECK_MAX_ANTISENSE_FRAC = 0.2
STRAND_CHECK_PREFIX = "#strand_check"
STRAND_CHECK_TSV_COLUMNS = [
    "units",
    "spliced_reads_checked",
    "spliced_reads_antisense",
    "antisense_frac",
    "max_antisense_frac",
    "min_spliced_reads",
    "verdict",
]


def is_strand_check_line(line):
    return line.startswith(STRAND_CHECK_PREFIX + "\t")


def summarize_strand_check(lines):
    """Run-level strand check from the sidecar's STRAND_CHECK_PREFIX lines.

    Returns a dict keyed by STRAND_CHECK_TSV_COLUMNS. verdict: "stranded",
    "unstranded_or_reverse" (antisense fraction above STRAND_CHECK_MAX_ANTISENSE_FRAC),
    or "too_few_spliced_reads" (fewer than STRAND_CHECK_MIN_SPLICED_READS in all).
    """
    units = 0
    n = 0
    anti = 0
    for line in lines:
        if not is_strand_check_line(line):
            continue
        f = line.rstrip("\n").split("\t")
        units += 1
        n += int(float(f[2]))
        anti += int(float(f[3]))
    frac = anti / n if n else float("nan")
    if n < STRAND_CHECK_MIN_SPLICED_READS:
        verdict = "too_few_spliced_reads"
    elif frac > STRAND_CHECK_MAX_ANTISENSE_FRAC:
        verdict = "unstranded_or_reverse"
    else:
        verdict = "stranded"
    return {
        "units": units,
        "spliced_reads_checked": n,
        "spliced_reads_antisense": anti,
        "antisense_frac": frac,
        "max_antisense_frac": STRAND_CHECK_MAX_ANTISENSE_FRAC,
        "min_spliced_reads": STRAND_CHECK_MIN_SPLICED_READS,
        "verdict": verdict,
    }


def write_strand_check_tsv(path, summary):
    with open(path, "wt") as ofh:
        ofh.write("\t".join(STRAND_CHECK_TSV_COLUMNS) + "\n")
        vals = []
        for col in STRAND_CHECK_TSV_COLUMNS:
            v = summary[col]
            vals.append("{:.6f}".format(v) if isinstance(v, float) else str(v))
        ofh.write("\t".join(vals) + "\n")

# how a read tied after rules 0-3 was settled (GeneBodyRecorder.counts keys)
RULE_DOMINANCE = "dominance"
RULE_WEIGHT_SPLIT = "weight_split"
RULE_THREE_PRIME = "three_prime"
RULE_INNERMOST = "innermost"
RULE_AMBIGUOUS = "ambiguous"


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


def introns_from_segments(segments):
    """Gaps between consecutive aligned blocks (Pretty_alignment segments)."""
    segments = sorted(segments)
    return [
        (segments[k - 1][1] + 1, segments[k][0] - 1)
        for k in range(1, len(segments))
        if segments[k][0] - segments[k - 1][1] > 1
    ]


class GeneModel:
    __slots__ = ("gene_id", "lend", "rend", "exons", "introns", "intronic", "readthrough")

    def __init__(self, gene_id, exon_segment_lists, readthrough=False):
        self.gene_id = gene_id
        self.readthrough = readthrough
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


def _transcript_is_readthrough(t):
    tags = t.get_gtf_tags() if hasattr(t, "get_gtf_tags") else ()
    return READTHROUGH_TRANSCRIPT_TAG in tags or READTHROUGH_GENE_TAG in tags


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


class TiedGene:
    """What resolve() needs of one gene left after rules 0-3 (also the spill record)."""

    __slots__ = ("gene_id", "splice_class", "exonic_bp", "intronic_bp", "three_prime", "span_len")

    def __init__(self, gene_id, splice_class, exonic_bp, intronic_bp, three_prime, span_len):
        self.gene_id = gene_id
        self.splice_class = splice_class
        self.exonic_bp = exonic_bp
        self.intronic_bp = intronic_bp
        self.three_prime = three_prime
        self.span_len = span_len

    def to_fields(self):
        return [
            self.gene_id,
            self.splice_class,
            str(self.exonic_bp),
            str(self.intronic_bp),
            "1" if self.three_prime else "0",
            str(self.span_len),
        ]

    @classmethod
    def from_fields(cls, f):
        return cls(f[0], f[1], f[2], f[3], f[4] == "1", int(f[5]))


TIED_GENE_NFIELDS = 6


class GeneBodyIndex:
    """Gene models of one contig-strand, built from the annotation's transcripts."""

    def __init__(self, transcripts, contig_strand):
        by_gene = defaultdict(list)
        all_readthrough = {}
        for t in transcripts or []:
            if contig_strand and t.get_strand() != contig_strand:
                continue
            gid = t.get_output_gene_id()
            by_gene[gid].append(t.get_exon_segments())
            all_readthrough[gid] = all_readthrough.get(gid, True) and _transcript_is_readthrough(t)
        self.contig_strand = contig_strand
        self.genes = {
            gid: GeneModel(gid, segs, readthrough=all_readthrough[gid])
            for gid, segs in by_gene.items()
        }
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
            bool(cfg.get("gene_body_exclude_readthrough_genes", True)),
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
        """Genes left after rules 0-3 for a read isoform quant did not assign, or None.

        Returns a list of TiedGene sorted by gene id. More than one is a tie for
        resolve() to settle.
        """
        if not segments:
            return None
        strand = strand or self.contig_strand
        min_frac, min_bp, exclude_readthrough = self._settings()
        aligned_bp = sum(r - l + 1 for l, r in segments)
        if aligned_bp <= 0:
            return None
        cands = []
        for gene in self.overlapping_genes(segments[0][0], segments[-1][1]):
            if exclude_readthrough and gene.readthrough:
                continue
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
            )

        best = max(key(ev) for ev in cands)
        tied = sorted((ev for ev in cands if key(ev) == best), key=lambda e: e.gene.gene_id)
        return [
            TiedGene(
                ev.gene.gene_id,
                self.splice_class(ev, read_introns, min_bp),
                ev.exonic_bp,
                ev.intronic_bp,
                ev.three_prime_in_span,
                ev.gene.span_len(),
            )
            for ev in tied
        ]

    @staticmethod
    def resolve(tied, gene_weights=None):
        """Settle rules 4-7: ([(TiedGene, frac)], rule) for a list of TiedGene.

        rule is None when only one gene was left after rules 0-3.
        """
        if len(tied) == 1:
            return [(tied[0], 1.0)], None
        cfg = LRAA_Globals.config
        if gene_weights is not None:
            w = [float(gene_weights.get(t.gene_id, 0.0)) for t in tied]
            tot = sum(w)
            if tot > 0:
                dominance = float(cfg.get("gene_body_dominance_frac", 0.9))
                top = max(range(len(tied)), key=lambda k: w[k])
                if w[top] / tot >= dominance:
                    return [(tied[top], 1.0)], RULE_DOMINANCE
                return [(t, wi / tot) for t, wi in zip(tied, w) if wi > 0], RULE_WEIGHT_SPLIT
        three = [t for t in tied if t.three_prime]
        if len(three) == 1:
            return [(three[0], 1.0)], RULE_THREE_PRIME
        if three:
            tied = three
        smallest = min(t.span_len for t in tied)
        inner = [t for t in tied if t.span_len == smallest]
        if len(inner) == 1:
            return [(inner[0], 1.0)], RULE_INNERMOST
        policy = str(cfg.get("gene_body_ambiguous_policy", POLICY_DROP))
        frac = 1.0 / len(inner) if policy == POLICY_SPLIT else 0.0
        return [(t, frac) for t in inner], RULE_AMBIGUOUS


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


class GeneBodyRecorder:
    """Writes one sidecar row per (read, gene) for every retained read of a contig-strand.

    source = "isoform": the read's isoform assignment (EM split summed per gene).
    source = "gene_body:<why>": placed by GeneBodyIndex; <why> is the reason isoform
             assignment failed (no_path, spacer_path, unassignable, no_quant_unit,
             oversimplify_no_exon_overlap).
    Reads with no candidate gene are counted, not written.

    A read tied between genes after rules 0-3 is spilled to a temp file and resolved by
    finish(), once the whole stream has been seen: rule 4 weights the tied genes by
    their isoform-assigned reads over the FULL library of this unit, which only exists
    at the end. Memory stays bounded by the gene count, whatever the tie count.
    """

    def __init__(self, index, fh, spill_dir=None, label=None):
        self.index = index
        self.fh = fh
        self.counts = defaultdict(float)
        self.gene_isoform_reads = defaultdict(float)
        self._spill = None
        self._spill_dir = spill_dir
        self.label = label if label is not None else str(index.contig_strand)

    def note_aligned_strand(self, aligned_strand, spliced):
        """Count a spliced read whose aligned strand disagrees with this unit's strand."""
        if not spliced:
            return
        self.counts["strand_check_spliced"] += 1
        if aligned_strand != self.index.contig_strand:
            self.counts["strand_check_spliced_antisense"] += 1

    def _warn_if_unstranded(self):
        n = self.counts.get("strand_check_spliced", 0)
        if n < STRAND_CHECK_MIN_SPLICED_READS:
            return
        frac = self.counts.get("strand_check_spliced_antisense", 0) / n
        if frac > STRAND_CHECK_MAX_ANTISENSE_FRAC:
            logger.warning(
                "[%s] gene-body counts: %.1f%% of %d spliced reads align to the strand "
                "opposite the one their junctions are transcribed from. The library looks "
                "unstranded (or reverse-stranded); unspliced reads take their aligned "
                "strand, so their gene-body counts will be wrong.",
                self.label, 100.0 * frac, n,
            )

    def _write(self, gene_id, read_name, frac, cls, source, n, aligned_bp, exonic_bp, intronic_bp):
        self.fh.write(
            "{}\t{}\t{:.6f}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                gene_id, read_name, frac, cls, source, n, aligned_bp, exonic_bp, intronic_bp
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
            self._write(
                gene_id, read_name, by_gene[gene_id], cls, "isoform", n, aligned_bp,
                ev.exonic_bp if ev is not None else "",
                ev.intronic_bp if ev is not None else "",
            )
            self.counts["isoform_" + cls] += by_gene[gene_id]
            self.gene_isoform_reads[gene_id] += by_gene[gene_id]
        self.counts["reads_isoform"] += 1

    def _emit(self, read_name, source, aligned_bp, resolved, rule):
        n = len(resolved)
        for t, frac in resolved:
            self._write(
                t.gene_id, read_name, frac, t.splice_class, source, n, aligned_bp,
                t.exonic_bp, t.intronic_bp,
            )
            self.counts["gene_body_" + t.splice_class] += frac
        self.counts["reads_gene_body"] += 1
        if rule is not None:
            self.counts["tie_" + rule] += 1

    def record_unassigned(self, read_name, segments, read_introns, reason):
        tied = self.index.assign(segments, read_introns)
        if tied is None:
            self.counts["reads_no_gene"] += 1
            return
        aligned_bp = sum(r - l + 1 for l, r in segments)
        source = "gene_body:" + reason
        if len(tied) == 1:
            self._emit(read_name, source, aligned_bp, [(tied[0], 1.0)], None)
            return
        self.counts["reads_gene_tied"] += 1
        if self._spill is None:
            self._spill = tempfile.TemporaryFile(mode="w+t", dir=self._spill_dir)
        fields = [read_name, source, str(aligned_bp)]
        for t in tied:
            fields += t.to_fields()
        self._spill.write("\t".join(fields) + "\n")

    def _write_strand_check(self):
        n = int(self.counts.get("strand_check_spliced", 0))
        if n == 0:
            return
        anti = int(self.counts.get("strand_check_spliced_antisense", 0))
        self.fh.write("{}\t{}\t{}\t{}\n".format(STRAND_CHECK_PREFIX, self.label, n, anti))

    def finish(self):
        """Resolve the spilled ties with this unit's isoform-assigned gene counts."""
        self._warn_if_unstranded()
        self._write_strand_check()
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
            tied = [
                TiedGene.from_fields(f[k : k + TIED_GENE_NFIELDS])
                for k in range(3, len(f), TIED_GENE_NFIELDS)
            ]
            resolved, rule = GeneBodyIndex.resolve(tied, weights)
            self._emit(read_name, source, aligned_bp, resolved, rule)
        self._spill.close()
        self._spill = None
