#!/usr/bin/env python3

"""Extract a mini contig, BAM and GTF for one CHUNK of a contig-strand.

Chunks are strictly disjoint coordinate intervals that tile the contig-strand.
There is no halo and no widening: a chunk emits exactly the records that lie
wholly inside it. Two different rules get it there, and the difference between
them is the whole design.

Annotated loci: a boundary may not cut one, ever
------------------------------------------------
A record spanning ``[s, e]`` BLOCKS the cut positions ``[s - margin, e + margin - 1]``.
With ``margin == 0`` that is exactly the set of cuts the record spans; a positive
margin additionally demands clearance on both sides. Every same-strand GENE locus
blocks, and blocking is absolute: ``extract_partition`` refuses the region.

It has to be absolute because a locus cannot be dropped and recovered. LRAA's GTF
reader filters features by containment, so a locus straddling a boundary is
contained by NEITHER neighbour and both omit it -- its isoforms would vanish from
the run with nothing anywhere to attribute the loss to. Blocking is on the GENE
span rather than each isoform because ``genes_contained`` emits a gene whole or
not at all; that makes it the strictly stronger form of "no cut inside an
annotated isoform" and the form that matches emission.

Merging the blocked intervals yields ISLANDS: maximal runs of coordinate that
cannot be divided. The gaps between them are the admissible CUT ZONES, and this
is the hard constraint the cut selector in ``select_contig_cut_points.py``
respects when it searches a window for a position.

Alignments: an overhanging read is DROPPED, counted and named
------------------------------------------------------------
A retained primary alignment crossing a boundary does NOT forbid it. The
alignment is dropped: not truncated, which would invent an alignment the aligner
never made, and not widened around, which would break disjointness. Dropping is
what buys equal-span chunks -- demanding that nothing cross a boundary makes
granularity a function of where the alignments happen to leave gaps, and on real
data the largest indivisible island then sets a floor no chunk count can beat.

The cost is paid in reads and must therefore be visible. Every drop is counted in
the manifest (``counts.alignments_dropped_overhang``) and every dropped name is
written to ``{output_prefix}.dropped_reads.txt``, because the comparison
methodology builds a pruned baseline BAM from exactly that set and a count alone
would not let anyone reproduce it. The selector minimises the count in advance by
choosing, within each wiggle window, the position the fewest alignments span.

What disjointness proves, and what it does not
---------------------------------------------
Quantification splits a read's mass across the transcripts it may be assigned to,
so a chunk that sees only one of two competing transcripts awards it the whole
read; two such chunks then merge to 2.0 where an unpartitioned run gives 1.0.
No duplicate transcript id and no repeated read name is involved, only the missing
competitor, which is why the argument has to be made rather than assumed.

PROVEN, by construction: no duplicate processing. Each emitted alignment is
physically present in exactly one chunk's BAM, so no assignment path can award
its mass twice. Each retained primary alignment is emitted exactly once or
dropped, never both, and the dropped set is enumerated. That excludes duplication
only -- it does NOT bound merged mass against the unpartitioned answer in either
direction, because a chunk's candidate set differs from the whole contig's: a
read left unassigned in an unpartitioned run can become assigned in a chunk, and
vice versa. Equivalence, and the sign of any delta, are empirical questions. This
does not depend on what the assignment path does, which matters, because it is
broader than coordinate overlap: with the shipped default
``rescue_unassigned_reads_via_transcriptome_alignment`` an unassigned read is
re-aligned against transcript SEQUENCES and reprojected onto the hit's genomic
path (IsoformReadRescue.py) with no test that the read's original blocks overlap
that transcript. So a read can be assigned to a transcript it does not overlap at
all, and no overlap-based argument would be safe.

NOT PROVEN, and empirical only: equivalence to an unpartitioned run. That same
rescue path is why. A read's baseline target need not overlap it, so under hard
cuts that target may sit in another chunk where the read is absent, and the read
is then assigned elsewhere or not at all -- merged mass short of baseline, never
over. Measured in ``__LRAA_local_runs/SegmentedQuant.md``: exact on all 82,884
read-tracking keys at 100, 50 and 25 annotated genes per segment, and one read of
82,884 short at 10 genes per segment. Dropped overhanging reads are a second,
separate deficit on top of that, which is why the baseline they are compared
against has to be pruned by the same set.

``overlap_offenders`` verifies the coordinate half of this empirically over a
whole partitioning: no EMITTED alignment overlaps a transcript outside its own
chunk. Overlap, not splice compatibility, because overlap is the weaker and
therefore safer notion.

What the margin buys
--------------------
Because a locus blocks every cut within ``margin`` of itself, no annotated model
sits within ``margin`` of either edge: every emitted GTF feature has that much
clean flanking sequence inside the mini contig. At the default 200 bp that covers
every sequence-context window LRAA evaluates near a model's ends -- the 40 bp
polyA-signal window (``polyA_signal_window`` [-40, -10]), the 20 bp
internal-priming window (``Util_funcs.check_internal_priming``, check_length=20)
and the 50 bp read-enrichment window (``TSS_window_read_enrich_len``) -- as well
as the 50 bp TSS/polyA site-snapping distances.

Stated rather than implied: the margin does NOT extend to alignments. A read may
sit flush against a boundary with no flank at all. Extending the margin to reads
would mean dropping reads that are wholly contained, which would make the dropped
count larger than the spanning count the selector minimised and computed in
advance -- and an accounting identity that holds is worth more here than a
flanking guarantee for reads, since a read's ends are evidence rather than a model
whose sequence context is being judged.

Cut alignment to the depth-window grid
--------------------------------------
Boundaries are multiples of ``depth_window`` (100 bp) so that no depth window
spans a cut. A window straddling a boundary would draw aligned bases from two
chunks, giving a different median, a different acceptance probability, and
different reads kept near the threshold -- per-chunk normalization would not
reproduce the whole-contig result. A cut at 1-based ``b`` leaves ``b`` as the last
base of a window and ``b + 1`` as the first base of the next, so with the grid
origin pinned at absolute 0 the condition is ``b % depth_window == 0``.

Each chunk's manifest carries ``window_origin``, the value to pass
``normalize_bam_by_strand.py --window_origin``: the rebase offset in the GENOME
frame, unreduced, since absolute 0-based == chunk-local 0-based + offset. Because
every boundary is grid-aligned, so is every offset. The whole-contig control run
must be given ``--window_origin 0`` explicitly; the unset default anchors on the
first aligned base per contig and no chunk grid can match it.

Normalization scope otherwise: a mini contig is a different input, so whatever a
downstream run derives from its input as a whole is derived over the chunk
instead -- notably the mapped-primary library count. A chunked run should pass
the global count with ``--num_total_reads``, or each chunk scales
``RPM_total_reads`` to its own slice. TPM is unaffected: it is normalised over
the run's own assigned reads, not the library. Secondary alignments need no such
care: LRAA
processes primary non-supplementary alignments only, so there is no rescue
grouping left to span a boundary.

Strandless chunks, and the one ordering that must not be reversed
-----------------------------------------------------------------
A region with no strand suffix (``chr1:1-10000000``) is a STRANDLESS chunk: one
mini contig, one GTF carrying both strands' features, one BAM carrying both
orientations, and the strand split deferred to a per-chunk step downstream. That
is the whole point -- the split is the largest serial phase in the pipeline, and
cutting first moves it into the parallel phase. It also halves the mini-FASTA:
two strand chunks over the same interval each wrote the same sequence, and one
strandless chunk writes it once.

Nothing in the emission path needed widening for it. ``_strand_matches`` returns
True for both orientations on a falsy strand, ``find_islands`` unions both
strands' loci, and ``admissibility_offenders`` therefore refuses a boundary that
cuts a locus on EITHER strand. What is new is the accounting -- per-orientation
emitted counts, ``"strand": null`` and ``"strand_split_required"`` -- and a
refusal.

The refusal is the ordering constraint. ``separate_bam_by_strand.py`` REWRITES
``is_reverse`` when a read's inferred transcribed strand disagrees with the
aligner, and the strand filter here reads the RAW flag. So a chunk must be
extracted strandlessly and split AFTERWARDS. Extracting ``chr1+:...`` from a raw
bam would silently assign every flipped read to the wrong strand and still
produce output that looks fine, which is why the CLI refuses a strand-suffixed
region whose source bam still holds the other orientation rather than trusting
the caller to remember the order.
"""

import argparse
import bisect
import collections
import contextlib
import json
import logging
import os
import re
import subprocess
import sys

import pysam

sys.path.insert(
    0, os.path.sep.join([os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"])
)

import LRAA_Globals
import Util_funcs

logger = logging.getLogger(__name__)

# Retention filter shared by cut selection and emission: a spanning count is only
# worth anything if it is computed over exactly the alignments that will be
# written out.
NONPRIMARY_FLAGS = 0x100 | 0x800

# CIGAR ops, in pysam's integer encoding.
_CIGAR_REF_CONSUMING = frozenset((0, 2, 7, 8))  # M D = X
_CIGAR_INTRON = 3  # N

# Bases of clearance a cut must leave on both sides of an annotated locus. The
# rationale for the value lives with it in LRAA_Globals: 4x the largest
# boundary-snapping distance there, so no snapping can reach across a cut.
DEFAULT_MARGIN = LRAA_Globals.config["chunk_margin"]


class ExtractionError(RuntimeError):
    """Malformed input, an inadmissible cut, or a record that would not fit."""


Region = collections.namedtuple("Region", "chrom strand lend rend")


class TranscriptModel:

    __slots__ = ("transcript_id", "gene_id", "strand", "lend", "rend", "exons", "lines")

    def __init__(self, transcript_id, gene_id, strand):
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.strand = strand
        self.lend = None
        self.rend = None
        self.exons = []
        self.lines = []

    def add_line(self, line, feature, lend, rend):
        self.lines.append(line)
        if feature == "exon":
            self.exons.append((lend, rend))
        self.lend = lend if self.lend is None else min(self.lend, lend)
        self.rend = rend if self.rend is None else max(self.rend, rend)

    def finalize(self):
        if self.exons:
            self.exons.sort()
        else:
            # a transcript line with no exon lines is a single-block model
            self.exons = [(self.lend, self.rend)]


class GeneModel:

    __slots__ = ("gene_id", "strand", "lend", "rend", "lines", "transcript_ids")

    def __init__(self, gene_id, strand):
        self.gene_id = gene_id
        self.strand = strand
        self.lend = None
        self.rend = None
        self.lines = []  # gene-level lines, i.e. lines carrying no transcript_id
        self.transcript_ids = []

    def note_span(self, lend, rend):
        self.lend = lend if self.lend is None else min(self.lend, lend)
        self.rend = rend if self.rend is None else max(self.rend, rend)

    @property
    def label(self):
        return "gene locus {} ({}{}:{}-{})".format(
            self.gene_id, "", self.strand or ".", self.lend, self.rend
        )


class _IntervalIndex:
    """Static interval set supporting overlap queries, sorted by start."""

    def __init__(self, items):
        # items: iterable of (lend, rend, payload)
        self._items = sorted(items, key=lambda it: (it[0], it[1]))
        self._lends = [it[0] for it in self._items]
        self._max_rend = []
        running = 0
        for item in self._items:
            running = max(running, item[1])
            self._max_rend.append(running)

    def overlapping(self, lend, rend):
        idx = bisect.bisect_right(self._lends, rend)
        for i in range(idx - 1, -1, -1):
            if self._max_rend[i] < lend:
                break  # nothing at or left of i can reach lend
            item = self._items[i]
            if item[1] >= lend:
                yield item[2]


class Annotation:
    """Gene and transcript loci of one contig, optionally one strand."""

    def __init__(self, transcripts=None, genes=None):
        self.transcripts = transcripts if transcripts is not None else {}
        self.genes = genes if genes is not None else {}
        self._gene_index = _IntervalIndex(
            (g.lend, g.rend, g) for g in self.genes.values()
        )
        self._transcript_index = _IntervalIndex(
            (t.lend, t.rend, t) for t in self.transcripts.values()
        )

    def genes_overlapping(self, lend, rend):
        return self._gene_index.overlapping(lend, rend)

    def transcripts_overlapping(self, lend, rend):
        return self._transcript_index.overlapping(lend, rend)

    def genes_contained(self, lend, rend):
        return [g for g in self.genes_overlapping(lend, rend) if g.lend >= lend and g.rend <= rend]


def parse_region(region_str):
    """``chr20-:100-200`` -> Region(chrom='chr20', strand='-', lend=100, rend=200).

    The chrom group is non-greedy rather than a negated class. It used to be
    ``[^:+-]+``, which excluded the hyphen and so refused every contig with one in
    its name: ``LRAA --chunk`` died at stage 3 with "Cannot parse region string" on
    GRCh38 analysis-set HLA contigs (``HLA-A*``), on 1 of the 50 contigs in
    ``testing/sep_contigs`` (``TRIM39-RPP21^ENSG...``), and on all 6 of
    ``testing/ont_sep_contigs`` (``minig-ENSG...``), which could therefore not be
    chunked at all.

    Non-greedy plus the anchored digits tail resolves the shapes by backtracking:
    ``HLA-A:1-100`` gives chrom ``HLA-A``; ``chr20-:100-200`` still gives chrom
    ``chr20`` strand ``-``; and a name containing colons, ``HLA-A:01:01:1-100``,
    gives chrom ``HLA-A:01:01``, because only the final digits-dash-digits can
    satisfy the tail.

    One ambiguity is irreducible in this encoding and is unchanged: a contig whose
    name genuinely ends in ``+`` or ``-``, addressed without a strand suffix, cannot
    be told from the same name carrying one. The suffix wins.
    """

    region_str = region_str.replace(",", "")
    m = re.match(r"^(.+?)([+-]?):(\d+)-(\d+)$", region_str)
    if m is None:
        raise ExtractionError(
            "Cannot parse region string: {} (expected chrom[+-]:lend-rend)".format(
                region_str
            )
        )
    lend = int(m.group(3))
    rend = int(m.group(4))
    if lend < 1 or rend < lend:
        raise ExtractionError("Nonsensical region bounds: {}".format(region_str))
    return Region(m.group(1), m.group(2), lend, rend)


def _attribute(attributes, key):
    m = re.search(r'{}[\s=]+"?([^";]+)"?'.format(re.escape(key)), attributes)
    return m.group(1).strip() if m else None


class _GtfIngest:
    """Accumulate GTF lines into gene and transcript models for one contig-strand.

    One parser serves both the full scan and the indexed region fetch. Two copies
    would be free to disagree about which lines belong to which gene, and the
    symptom of a disagreement is a silently dropped locus -- the exact failure
    ``admissibility_offenders`` exists to refuse -- rather than an error.
    """

    __slots__ = ("transcripts", "genes", "chrom", "strand")

    def __init__(self, chrom, strand=""):
        self.transcripts = {}
        self.genes = {}
        self.chrom = chrom
        self.strand = strand

    def ingest(self, line, where):
        if line.startswith("#") or not line.strip():
            return
        vals = line.rstrip("\n").split("\t")
        if len(vals) < 9:
            raise ExtractionError(
                "{}: GTF line has {} fields, need 9".format(where, len(vals))
            )
        if vals[0] != self.chrom:
            return
        feature_strand = vals[6]
        if self.strand and feature_strand != self.strand:
            return
        lend, rend = int(vals[3]), int(vals[4])
        if rend < lend:
            raise ExtractionError(
                "{}: GTF feature end {} precedes start {}".format(where, rend, lend)
            )
        attributes = vals[8]
        transcript_id = _attribute(attributes, "transcript_id")
        gene_id = _attribute(attributes, "gene_id")
        if transcript_id is None and gene_id is None:
            raise ExtractionError(
                "{}: GTF line carries neither gene_id nor transcript_id; "
                "it cannot be assigned to a partition".format(where)
            )
        if gene_id is None:
            gene_id = transcript_id

        gene = self.genes.get(gene_id)
        if gene is None:
            gene = self.genes[gene_id] = GeneModel(gene_id, feature_strand)
        elif gene.strand != feature_strand:
            raise ExtractionError(
                "{}: gene {} appears on both strands ({} and {})".format(
                    where, gene_id, gene.strand, feature_strand
                )
            )
        gene.note_span(lend, rend)

        if transcript_id is None:
            gene.lines.append(line)
            return

        transcript = self.transcripts.get(transcript_id)
        if transcript is None:
            transcript = self.transcripts[transcript_id] = TranscriptModel(
                transcript_id, gene_id, feature_strand
            )
            gene.transcript_ids.append(transcript_id)
        elif transcript.gene_id != gene_id:
            raise ExtractionError(
                "{}: transcript {} claims gene {} and gene {}".format(
                    where, transcript_id, transcript.gene_id, gene_id
                )
            )
        transcript.add_line(line, vals[2], lend, rend)

    def finish(self):
        for transcript in self.transcripts.values():
            transcript.finalize()
        return Annotation(self.transcripts, self.genes)


def load_gtf(gtf_filename, chrom, strand=""):
    """Index a GTF for one contig-strand, grouping every line under its gene.

    Reads the file end to end. ``load_gtf_for_region`` is the indexed equivalent
    and is what the per-chunk path uses; this stays as the fallback for a GTF
    that cannot be indexed, and as the reference the indexed path is tested
    against.
    """

    ingest = _GtfIngest(chrom, strand)
    with open(gtf_filename, "rt") as fh:
        for lineno, line in enumerate(fh, start=1):
            ingest.ingest(line, "{}:{}".format(gtf_filename, lineno))
    return ingest.finish()


# A sorted, bgzipped, tabix-indexed copy of a GTF, beside the original. The
# original is never modified: reference GTFs are shared, and often read-only.
GTF_INDEX_SUFFIX = ".lraa_tabix.gtf.gz"
GTF_INDEX_STAMP_SUFFIX = ".lraa_tabix.json"

# A tabix-indexed GTF is just its path; `ensure_gtf_index` returns one or None.


def _gtf_index_key(gtf_filename):
    """What makes a cached index stale.

    Size and mtime rather than a content hash: digesting a 1.5 GB GTF on every
    chunk would cost more than the scan the index exists to replace.
    """

    st = os.stat(gtf_filename)
    return {
        "path": os.path.realpath(gtf_filename),
        "size": st.st_size,
        "mtime_ns": st.st_mtime_ns,
    }


def _gtf_index_basename(gtf_filename):
    base = os.path.basename(gtf_filename)
    for suffix in (".gtf.gz", ".gff3.gz", ".gff.gz", ".gtf", ".gff3", ".gff"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def _gtf_index_homes(gtf_filename, cache_dir):
    """Where an index may live, best first.

    Beside the GTF is preferred, so that repeat runs against the same reference
    reuse one index instead of rebuilding per run. Reference directories are
    often read-only, hence the fallback.
    """

    homes = [os.path.dirname(os.path.realpath(gtf_filename)) or "."]
    if cache_dir:
        homes.append(cache_dir)
    return homes


def _read_gtf_stamp(stamp_path):
    try:
        with open(stamp_path, "rt") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _build_gtf_index(gtf_filename, gz_path):
    """Sort, compress and index a copy of ``gtf_filename``.

    Sorting is delegated to ``sort``, which spills to disk: a whole-genome GTF
    does not fit in the memory chunking exists to bound.

    Written to temporaries and moved into place, because chunk workers race to
    build the same index and a half-written one must never be visible.
    """

    tmp_sorted = "{}.{}.sorting".format(gz_path, os.getpid())
    tmp_gz = "{}.{}.tmp".format(gz_path, os.getpid())

    try:
        env = dict(os.environ, LC_ALL="C")
        with open(tmp_sorted, "wt") as ofh:
            sorter = subprocess.Popen(
                ["sort", "-k1,1", "-k4,4n"],
                stdin=subprocess.PIPE,
                stdout=ofh,
                env=env,
                text=True,
            )
            with open(gtf_filename, "rt") as fh:
                for line in fh:
                    # Comments and blanks are dropped rather than sorted: the
                    # readers ignore them, and a blank line is not a record
                    # tabix can index.
                    if line.startswith("#") or not line.strip():
                        continue
                    sorter.stdin.write(line)
            sorter.stdin.close()
            if sorter.wait() != 0:
                raise ExtractionError(
                    "sort failed with exit {} while indexing {}".format(
                        sorter.returncode, gtf_filename
                    )
                )



        pysam.tabix_compress(tmp_sorted, tmp_gz, force=True)
        pysam.tabix_index(tmp_gz, preset="gff", force=True)
        os.replace(tmp_gz + ".tbi", gz_path + ".tbi")
        os.replace(tmp_gz, gz_path)
    finally:
        for path in (tmp_sorted, tmp_gz, tmp_gz + ".tbi"):
            try:
                os.unlink(path)
            except OSError:
                pass

    return None


def ensure_gtf_index(gtf_filename, cache_dir=None):
    """Path to a tabix-indexed copy of ``gtf_filename``, building one if needed.

    Returns None when none could be built or reused, which tells the caller to
    fall back to ``load_gtf``. A missing index is a performance problem, not a
    correctness one, so refusing the run over it would be the worse failure.
    """

    key = _gtf_index_key(gtf_filename)
    base = _gtf_index_basename(gtf_filename)
    homes = _gtf_index_homes(gtf_filename, cache_dir)

    for home in homes:
        gz = os.path.join(home, base + GTF_INDEX_SUFFIX)
        cached = _read_gtf_stamp(os.path.join(home, base + GTF_INDEX_STAMP_SUFFIX))
        if (
            cached is not None
            and cached.get("key") == key
            and os.path.exists(gz)
            and os.path.exists(gz + ".tbi")
        ):
            return gz

    for home in homes:
        gz = os.path.join(home, base + GTF_INDEX_SUFFIX)
        stamp = os.path.join(home, base + GTF_INDEX_STAMP_SUFFIX)
        try:
            os.makedirs(home, exist_ok=True)
            _build_gtf_index(gtf_filename, gz)
            tmp_stamp = "{}.{}.tmp".format(stamp, os.getpid())
            with open(tmp_stamp, "wt") as ofh:
                json.dump({"key": key}, ofh)
            os.replace(tmp_stamp, stamp)
        except (OSError, ValueError, ExtractionError) as err:
            logger.warning(
                "could not build a tabix index for %s under %s: %s",
                gtf_filename,
                home,
                err,
            )
            continue
        logger.info("indexed %s at %s", gtf_filename, gz)
        return gz

    logger.warning(
        "no writable location for a tabix index of %s; falling back to a full "
        "scan per region, which is what makes extraction cost scale with chunk count",
        gtf_filename,
    )
    return None


def load_gtf_region(index_path, gtf_filename, chrom, strand, lend, rend, margin):
    """``load_gtf`` restricted to a region, served from a tabix index.

    The fetch is widened by ``margin`` because that is the window the consumer
    asks about, not because a gene might be long. ``admissibility_offenders``
    tests each boundary against ``cut +/- margin``, so a gene lying wholly
    OUTSIDE the chunk but within the margin of its edge still blocks the cut. A
    fetch of ``[lend, rend]`` alone never returns that gene, and the check then
    reports zero offenders where the full scan reports one -- it does not fail,
    it stops being able to fire.

    Nothing wider is needed. Cut selection places boundaries only in the gaps
    between annotated islands, so a gene overlapping the chunk lies wholly
    inside it and all of its lines are inside the fetch; the derived span then
    equals the full scan's. The one case this cannot see -- a gene whose
    transcripts do not overlap each other, straddling a boundary with its far
    transcripts beyond the margin -- is a case cut selection cannot produce,
    because it blocks on the gene span it computed from the whole contig.
    """

    ingest = _GtfIngest(chrom, strand)
    fetch_lend = max(1, lend - margin)
    fetch_rend = rend + margin
    where = "{}@{}:{}-{}".format(gtf_filename, chrom, fetch_lend, fetch_rend)
    with pysam.TabixFile(index_path) as tabix:
        if chrom not in tabix.contigs:
            return ingest.finish()
        for line in tabix.fetch(chrom, fetch_lend - 1, fetch_rend):
            ingest.ingest(line, where)
    return ingest.finish()


def load_gtf_for_region(
    gtf_filename, chrom, strand, lend, rend, margin=DEFAULT_MARGIN, cache_dir=None
):
    """Annotation of one region, indexed where possible and scanned where not."""

    index_path = ensure_gtf_index(gtf_filename, cache_dir=cache_dir)
    if index_path is None:
        return load_gtf(gtf_filename, chrom, strand)
    return load_gtf_region(
        index_path, gtf_filename, chrom, strand, lend, rend, margin
    )


def load_gtf_for_contig(gtf_filename, chrom, strand="", cache_dir=None):
    """Annotation of one whole contig-strand, indexed where possible.

    No widening: the fetch already covers the contig, so nothing a gene could
    reach lies outside it.
    """

    index_path = ensure_gtf_index(gtf_filename, cache_dir=cache_dir)
    if index_path is None:
        return load_gtf(gtf_filename, chrom, strand)

    ingest = _GtfIngest(chrom, strand)
    where = "{}@{}".format(gtf_filename, chrom)
    with pysam.TabixFile(index_path) as tabix:
        if chrom not in tabix.contigs:
            return ingest.finish()
        for line in tabix.fetch(chrom):
            ingest.ingest(line, where)
    return ingest.finish()


def alignment_blocks(aln):
    """Aligned reference blocks, 1-based inclusive, split only at introns (N)."""

    blocks = []
    position = aln.reference_start  # 0-based
    block_start = position
    block_end = position
    for op, length in aln.cigartuples or ():
        if op in _CIGAR_REF_CONSUMING:
            block_end = position + length
            position = block_end
        elif op == _CIGAR_INTRON:
            if block_end > block_start:
                blocks.append((block_start + 1, block_end))
            position += length
            block_start = position
            block_end = position
    if block_end > block_start:
        blocks.append((block_start + 1, block_end))
    if not blocks:
        raise ExtractionError(
            "alignment {} consumes no reference bases".format(aln.query_name)
        )
    return blocks


def _is_nonprimary(aln):
    return bool(aln.flag & NONPRIMARY_FLAGS)


def _strand_matches(aln, strand):
    if not strand:
        return True
    return aln.is_forward if strand == "+" else aln.is_reverse


def retained_for_extraction(aln, strand, max_intron_length=None):
    """The alignments cut selection must cost, and emission must write.

    Primary alignments only, minus any carrying an intron the pipeline declines
    to model. This is a superset of what LRAA retains for quantification when
    secondary alignments are disallowed: LRAA further filters on mapping quality
    and per-base identity, so any alignment it could use is one this filter keeps.

    The intron rule comes from ``Util_funcs.has_disqualifying_long_intron``, the
    single implementation that also governs the strand split and read ingestion,
    so all three see one record set. It matters here beyond bookkeeping: an
    alignment carrying a 200 kb intron spans a fifth of a wiggle window, so
    counting one the pipeline has already discarded would push cut selection away
    from positions that are in fact free. ``max_intron_length`` defaults to the
    configured value; pass 0 to disable.
    """

    if aln.is_unmapped or _is_nonprimary(aln):
        return False
    if aln.is_duplicate or aln.is_qcfail:
        return False
    if max_intron_length is None:
        max_intron_length = LRAA_Globals.config["max_intron_length"]
    if Util_funcs.has_disqualifying_long_intron(aln, max_intron_length):
        return False
    return _strand_matches(aln, strand)


def blocked_cuts(lend, rend, margin):
    """The cut positions a record spanning [lend, rend] forbids."""

    return (lend - margin, rend + margin - 1)


class Island:
    """A maximal run of coordinate no boundary may divide.

    Annotated loci only. Alignments used to form islands too, back when a cut
    crossed by an alignment was refused; under the drop-overhang policy an
    alignment is a COST at a candidate position, not a prohibition, so it has no
    place in a structure whose only job is to say where cutting is forbidden.
    """

    __slots__ = (
        "lend",
        "rend",
        "content_lend",
        "content_rend",
        "gene_ids",
        "n_transcripts",
    )

    def __init__(self, lend, rend, content_lend, content_rend):
        self.lend = lend  # blocked-cut interval, margin included
        self.rend = rend
        self.content_lend = content_lend  # extent of the loci themselves
        self.content_rend = content_rend
        self.gene_ids = []
        self.n_transcripts = 0

    def absorb(self, lend, rend, content_lend, content_rend):
        self.rend = max(self.rend, rend)
        self.content_lend = min(self.content_lend, content_lend)
        self.content_rend = max(self.content_rend, content_rend)


def find_islands(annotation, chrom, strand="", margin=DEFAULT_MARGIN):
    """Indivisible annotated islands of one contig-strand, in coordinate order.

    Every same-strand gene locus falls inside exactly one island, so a boundary
    placed between islands leaves every locus wholly inside one chunk. Gene loci
    rather than individual transcripts: a gene is emitted whole or not at all
    (``genes_contained``), so splitting a gene between two isoforms would drop
    the gene from BOTH neighbours. Blocking on the gene span is therefore the
    strictly stronger form of "no cut inside an annotated isoform", and it is the
    form that matches what emission actually does.

    ``chrom`` is accepted for symmetry with the rest of the module; the
    annotation is already restricted to one contig by ``load_gtf``.
    """

    records = []  # (blocked_lend, blocked_rend, content_lend, content_rend, gene)
    for gene in annotation.genes.values():
        if strand and gene.strand != strand:
            continue
        lo, hi = blocked_cuts(gene.lend, gene.rend, margin)
        records.append((lo, hi, gene.lend, gene.rend, gene))

    records.sort(key=lambda r: (r[0], r[1]))

    islands = []
    for lo, hi, content_lo, content_hi, gene in records:
        if islands and lo <= islands[-1].rend + 1:
            # touching blocked intervals leave no admissible position between them
            island = islands[-1]
            island.absorb(lo, hi, content_lo, content_hi)
        else:
            island = Island(lo, hi, content_lo, content_hi)
            islands.append(island)
        island.gene_ids.append(gene.gene_id)
        island.n_transcripts += len(gene.transcript_ids)

    return islands


def cut_zones(islands, span_lend, span_rend):
    """Admissible cut positions inside ``[span_lend, span_rend]``, as ranges.

    A cut at ``b`` splits the span into ``[span_lend, b]`` and ``[b+1, span_rend]``,
    so ``b`` ranges over ``[span_lend, span_rend - 1]``.
    """

    zones = []
    cursor = span_lend
    for island in islands:
        if island.rend < cursor or island.lend > span_rend - 1:
            continue
        if island.lend > cursor:
            zones.append((cursor, min(island.lend - 1, span_rend - 1)))
        cursor = max(cursor, island.rend + 1)
    if cursor <= span_rend - 1:
        zones.append((cursor, span_rend - 1))
    return [(lo, hi) for lo, hi in zones if hi >= lo]


def overlap_offenders(bam_filename, annotation, chrom, strand, partitions):
    """Verify that no EMITTED read's OVERLAP set straddles a chunk boundary.

    Overlap, not splice compatibility: LRAA's majority-voting fallback assigns
    reads whose intron chain matches no annotated transcript, so overlap is the
    weaker and therefore safer notion to check. Returns a list of human-readable
    offences, empty when the partitioning is sound.

    Reads that straddle a boundary are skipped rather than reported: under the
    drop-overhang policy they are not emitted anywhere, so they have no chunk
    whose candidate set could be short. What must still hold is that a read that
    IS emitted sees every transcript it overlaps.
    """

    bounds = [p[0] for p in partitions]
    offenders = []
    with pysam.AlignmentFile(bam_filename, "rb") as bam:
        for aln in bam.fetch(chrom):
            if not retained_for_extraction(aln, strand):
                continue
            start = aln.reference_start + 1
            end = aln.reference_end
            index = bisect.bisect_right(bounds, start) - 1
            partition = partitions[index]
            if start < partition[0] or end > partition[1]:
                continue  # dropped for overhang, so it has no candidate set
            for transcript in annotation.transcripts_overlapping(start, end):
                if transcript.lend < partition[0] or transcript.rend > partition[1]:
                    offenders.append(
                        "read {} ({}:{}-{}) overlaps {} ({}-{}), which is not "
                        "inside its partition {}-{}".format(
                            aln.query_name,
                            chrom,
                            start,
                            end,
                            transcript.transcript_id,
                            transcript.lend,
                            transcript.rend,
                            partition[0],
                            partition[1],
                        )
                    )
    return offenders


def blocks_cut(lend, rend, cut, margin):
    """Does the record spanning [lend, rend] forbid a cut at ``cut``?

    Blocked positions are [lend - margin, rend + margin - 1], so the record
    blocks ``cut`` exactly when ``lend <= cut + margin`` and
    ``rend >= cut - margin + 1``. At margin 0 that is precisely "the record spans
    the junction between bases cut and cut+1", which is not an interval-overlap
    test -- hence the explicit predicate.
    """

    return lend <= cut + margin and rend >= cut - margin + 1


def admissibility_offenders(annotation, region, contig_length, margin=DEFAULT_MARGIN):
    """Annotated loci that forbid this region's boundaries, named.

    Loci only. An alignment crossing a boundary is dropped at emission, counted
    and named; a LOCUS crossing a boundary cannot be dropped, because
    ``genes_contained`` would omit it from both neighbours and its isoforms would
    disappear from the run. So this stays a hard refusal and takes no bam.

    A STRANDLESS region (``region.strand == ""``) is checked against the UNION:
    the strand test below falls through, and the annotation a strandless region
    loads carries both strands, so a boundary is refused if it cuts a locus on
    EITHER strand. That is the right rule and not merely a convenient
    fall-through -- a strandless chunk is the sole container for both strands'
    loci over its interval, so a locus it splits is lost exactly as a
    same-strand one would be.
    """

    edges = []
    if region.lend > 1:
        edges.append(("left", region.lend - 1))
    if region.rend < contig_length:
        edges.append(("right", region.rend))

    def _why(lend, rend, cut):
        if lend <= cut < rend:
            return " which straddles it"
        return " which is within the {} bp margin".format(margin)

    offenders = []
    for side, cut in edges:
        for gene in annotation.genes_overlapping(
            max(1, cut - margin), min(contig_length, cut + margin + 1)
        ):
            if region.strand and gene.strand != region.strand:
                continue
            if not blocks_cut(gene.lend, gene.rend, cut, margin):
                continue
            offenders.append(
                "{} edge cut at {} is blocked by gene locus {} ({}{}:{}-{}){}".format(
                    side,
                    cut,
                    gene.gene_id,
                    region.chrom,
                    gene.strand,
                    gene.lend,
                    gene.rend,
                    _why(gene.lend, gene.rend, cut),
                )
            )
    return offenders


def _write_mini_fasta(fasta, chrom, lend, rend, mini_contig_name, path, line_width=60):
    sequence = fasta.fetch(chrom, lend - 1, rend)
    with open(path, "wt") as ofh:
        print(">{}".format(mini_contig_name), file=ofh)
        for i in range(0, len(sequence), line_width):
            print(sequence[i : i + line_width], file=ofh)
    pysam.faidx(path)
    return len(sequence)


def _mini_header(source_header, mini_contig_name, mini_length):
    """Rebuild the SQ line. Inheriting the genome header via ``template=`` leaves
    the declared contig length longer than the contig actually written, which is
    inconsistent by construction."""

    header = source_header.to_dict()
    header["SQ"] = [{"SN": mini_contig_name, "LN": mini_length}]
    return pysam.AlignmentHeader.from_dict(header)


def _rebase_alignment(aln, offset, mini_header, mini_contig_name, mini_length):
    if aln.is_paired:
        raise ExtractionError(
            "alignment {} is flagged paired; this extractor supports single-end "
            "long-read alignments only and will not silently rewrite mate "
            "coordinates".format(aln.query_name)
        )

    new_start = aln.reference_start - offset
    new_end = aln.reference_end - offset
    if new_start < 0 or new_end > mini_length:
        raise ExtractionError(
            "alignment {} would land at {}-{} on a mini contig of 1-{}: it is not "
            "contained by the partition, so the boundary was not admissible".format(
                aln.query_name, new_start + 1, new_end, mini_length
            )
        )

    fields = aln.to_string().split("\t")
    fields[2] = mini_contig_name
    fields[3] = str(new_start + 1)
    fields[6] = "*"
    fields[7] = "0"
    fields[8] = "0"
    return pysam.AlignedSegment.fromstring("\t".join(fields), mini_header)


def _rebase_gtf_line(line, offset, mini_contig_name, mini_length, what):
    vals = line.rstrip("\n").split("\t")
    lend = int(vals[3]) - offset
    rend = int(vals[4]) - offset
    if lend < 1 or rend > mini_length:
        raise ExtractionError(
            "{} would land at {}-{} on a mini contig of 1-{}: it is not contained "
            "by the partition, so the boundary was not admissible".format(
                what, lend, rend, mini_length
            )
        )
    vals[0] = mini_contig_name
    vals[3] = str(lend)
    vals[4] = str(rend)
    return "\t".join(vals)


def _require_source_index(bam_filename):
    """Refuse to hand a consumer a source bam it cannot fetch a region out of.

    Under ``--reuse_source_bam`` the manifest names this file as the chunk's bam,
    and every downstream reader restricts it by contig, which needs the index.
    The region fetch above already needed one, so this cannot normally fail --
    but the manifest is a contract with a later stage in a later process, and an
    index deleted between them would surface as pysam's contextless
    ``fetch called on bamfile without index`` inside the strand split.
    """

    for suffix in (".bai", ".csi"):
        if os.path.exists(bam_filename + suffix):
            return
    raise ExtractionError(
        "--reuse_source_bam names {0} as the chunk's bam, but it has no .bai or "
        ".csi beside it. Every consumer of a reused source restricts it by "
        "contig, so the index is load-bearing rather than incidental. Index it "
        "(samtools index {0}) or drop --reuse_source_bam and let the chunk be "
        "extracted.".format(bam_filename)
    )


def extract_partition(
    genome_fa,
    bam,
    region,
    output_prefix,
    gtf=None,
    margin=DEFAULT_MARGIN,
    secondary_alignments="exclude",
    mini_contig_name=None,
    gtf_index_cache_dir=None,
    mixed_orientation_source="warn",
    max_intron_length=None,
    reuse_source_bam=False,
):
    """Extract one chunk. Returns the manifest, also written as JSON.

    Alignments that would extend past either boundary are DROPPED: counted,
    named into ``{output_prefix}.dropped_reads.txt``, and logged. Annotated loci
    are not droppable and a boundary that cuts one is refused outright.

    A region with NO strand suffix is a STRANDLESS chunk: both orientations are
    emitted into one bam over one mini contig, the manifest carries
    ``"strand": null`` and ``"strand_split_required": true``, and the split is a
    downstream step. ``mixed_orientation_source`` decides what a STRAND-SUFFIXED
    region does when the source bam still holds the other orientation -- see the
    refusal message for why that combination is a wrong answer rather than a
    filtered one. ``"warn"`` here, ``"refuse"`` at the CLI: filtering a mixed bam
    by raw orientation is a real capability with real callers, but nothing in the
    pipeline invokes it and inside the pipeline the combination can only be the
    ordering mistake.

    ``max_intron_length`` is the intron cap the RUN uses, not this module's
    default, and it has to be passed for the accounting to close. Under the
    strand-first pipeline the bam reaching this point had already been filtered
    at the run's value by the strand split, so reading the configured default
    here was invisible; a strandless chunk is cut from the RAW bam, so a value
    different from the default would make extraction retain a different record
    set than selection costed and than the in-chunk split will keep. ``None``
    reads ``LRAA_Globals.config``, 0 disables.

    ``reuse_source_bam`` is for the one chunk shape where the mini bam would be a
    RESTATEMENT of the source: a strandless chunk spanning its whole contig,
    ``[1, contig_length]``, offset 0, mini contig named and sized like the real
    one. There the emitted record stream is the source's own records on that
    contig, unmoved, minus what ``retained_for_extraction`` rejects -- and the
    tool that reads it next applies exactly that filter itself. So the copy buys
    nothing and costs the whole of the extraction: MEASURED on chrM of a 6.9 GB
    HG002 bam, 1,199,182 mapped records over 16,569 bp, the full extraction is
    67.5 s while the same fetch and the same predicates WITHOUT the rebase, the
    bgzf write and the index is 5.45 s -- 92 % of the cost is the copy.

    The region is still read and still counted, so the manifest states what a
    full extraction would have stated and remains this tool's own measurement
    rather than a number handed in from outside. What changes is that no bam is
    written and ``files.bam`` names the SOURCE, with
    ``bam_reused_from_source`` set so a consumer cannot mistake one for the
    other. Only whole-contig strandless regions qualify, and anything else is
    refused rather than approximated: a partial region needs the coordinate
    rebase, and a strand-suffixed one needs the orientation filter that only the
    written copy carries.

    ONE MEASURED DIFFERENCE, and it is in an intermediate rather than a result.
    ``_rebase_alignment`` normalises RNEXT/PNEXT/TLEN to ``*``/0/0 for an unpaired
    read; a reused source carries whatever the aligner wrote. On the bundled
    minigenome that is every one of 1,177 forward records differing in column 8
    alone -- PNEXT 0 against 1 -- with all ten other columns and every tag equal.
    Nothing downstream reads a mate coordinate for an unpaired read, and the
    merged ``quant.expr`` and ``quant.tracking`` come out identical across the two
    paths, which is the check that establishes it rather than the argument.
    """

    if secondary_alignments not in ("exclude", "reject"):
        raise ExtractionError(
            "secondary_alignments must be 'exclude' or 'reject', got {!r}".format(
                secondary_alignments
            )
        )
    if mixed_orientation_source not in ("warn", "refuse"):
        raise ExtractionError(
            "mixed_orientation_source must be 'warn' or 'refuse', got {!r}".format(
                mixed_orientation_source
            )
        )
    if margin < 0:
        raise ExtractionError("margin must be >= 0")
    if isinstance(region, str):
        region = parse_region(region)
    if mini_contig_name is None:
        mini_contig_name = region.chrom
    if reuse_source_bam and region.strand:
        raise ExtractionError(
            "--reuse_source_bam is defined only for a STRANDLESS whole-contig "
            "region, and {}{}:{}-{} names an orientation. A strand-suffixed chunk "
            "gets its orientation from the filter applied while the copy is "
            "written, so there is no copy to skip: reusing the source would hand "
            "the consumer both orientations.".format(
                region.chrom, region.strand, region.lend, region.rend
            )
        )
    if reuse_source_bam and mini_contig_name != region.chrom:
        raise ExtractionError(
            "--reuse_source_bam cannot rename the mini contig to {!r}: the reused "
            "records still carry {!r}, so the name has to be the source "
            "one.".format(mini_contig_name, region.chrom)
        )

    with pysam.FastaFile(genome_fa) as fasta:
        if region.chrom not in fasta.references:
            raise ExtractionError(
                "contig {} is absent from {}".format(region.chrom, genome_fa)
            )
        contig_length = fasta.get_reference_length(region.chrom)
        if region.rend > contig_length:
            raise ExtractionError(
                "requested partition {}:{}-{} runs past the end of {} ({} bp)".format(
                    region.chrom, region.lend, region.rend, region.chrom, contig_length
                )
            )
        if reuse_source_bam and (region.lend != 1 or region.rend != contig_length):
            raise ExtractionError(
                "--reuse_source_bam needs the region to span the whole contig, "
                "and {}:{}-{} covers {} bp of {}'s {}. A partial region rebases "
                "every alignment by {} and drops the ones that overhang, neither "
                "of which the source bam has had done to it.".format(
                    region.chrom,
                    region.lend,
                    region.rend,
                    region.rend - region.lend + 1,
                    region.chrom,
                    contig_length,
                    region.lend - 1,
                )
            )

        annotation = (
            load_gtf_for_region(
                gtf,
                region.chrom,
                region.strand,
                region.lend,
                region.rend,
                # the same margin admissibility_offenders applies below, so the
                # fetch covers exactly the coordinates the check will consult
                margin=margin,
                cache_dir=gtf_index_cache_dir,
            )
            if gtf
            else Annotation()
        )

        offenders = admissibility_offenders(annotation, region, contig_length, margin)
        if offenders:
            raise ExtractionError(
                "REJECTED: {}{}:{}-{} cuts an annotated locus at margin {} bp. {} "
                "blocking locus/loci; first {}: {}. Nothing was written. Unlike an "
                "alignment, a locus cannot be dropped and re-derived: a locus "
                "straddling a boundary is contained by neither neighbour, so BOTH "
                "would omit it and its isoforms would vanish from the run "
                "entirely. Choose a boundary that clears every annotated "
                "locus.".format(
                    region.chrom,
                    region.strand,
                    region.lend,
                    region.rend,
                    margin,
                    len(offenders),
                    min(3, len(offenders)),
                    "; ".join(offenders[:3]),
                )
            )

        offset = region.lend - 1
        mini_length = region.rend - region.lend + 1

        fasta_filename = "{}.fa".format(output_prefix)
        written = _write_mini_fasta(
            fasta,
            region.chrom,
            region.lend,
            region.rend,
            mini_contig_name,
            fasta_filename,
        )
        if written != mini_length:
            raise ExtractionError(
                "extracted {} bp of sequence for a {} bp partition".format(
                    written, mini_length
                )
            )

    counts = collections.Counter()
    dropped = []  # (name, start, end, side) in source coordinates
    # Set when a strand-suffixed region is being filtered out of a bam that still
    # carries the other orientation. Reported after the writer closes rather than
    # from inside it, so the partial chunk is removed rather than left to be
    # mistaken for output.
    mixed_orientation_offender = None
    # The region is READ and COUNTED either way -- only the writing is optional.
    # A manifest whose counts came from somewhere other than this loop would be
    # a different tool's claim about this chunk, and the downstream split
    # accounting checks itself against these numbers.
    bam_filename = (
        os.path.abspath(bam) if reuse_source_bam else "{}.bam".format(output_prefix)
    )
    with pysam.AlignmentFile(bam, "rb") as bamreader:
        mini_header = (
            None
            if reuse_source_bam
            else _mini_header(bamreader.header, mini_contig_name, mini_length)
        )
        nonprimary = []
        with (
            contextlib.nullcontext(None)
            if reuse_source_bam
            else pysam.AlignmentFile(bam_filename, "wb", header=mini_header)
        ) as bamwriter:
            for aln in bamreader.fetch(region.chrom, region.lend - 1, region.rend):
                if aln.is_unmapped:
                    continue
                if _is_nonprimary(aln) and _strand_matches(aln, region.strand):
                    counts["nonprimary_seen"] += 1
                    if len(nonprimary) < 5:
                        nonprimary.append(aln.query_name)
                # Retention is asked ORIENTATION-BLIND and the orientation tested
                # separately, because the two questions have different answers: a
                # record rejected for its flags or its introns says nothing about
                # which orientations the source bam holds, and that is exactly
                # what has to be known to tell a strand-split bam from a raw one.
                if not retained_for_extraction(aln, "", max_intron_length):
                    continue
                if not _strand_matches(aln, region.strand):
                    counts["opposite_orientation_seen"] += 1
                    if mixed_orientation_source == "refuse":
                        mixed_orientation_offender = aln.query_name
                        break
                    continue
                start = aln.reference_start + 1
                end = aln.reference_end
                if start < region.lend or end > region.rend:
                    # DROPPED, not truncated and not widened around. Truncating
                    # would invent an alignment the aligner never made; widening
                    # would break disjointness. The name is recorded so a pruned
                    # comparison baseline can be built from exactly this set --
                    # a count alone would not let anyone reproduce it.
                    sides = []
                    if start < region.lend:
                        sides.append("left")
                    if end > region.rend:
                        sides.append("right")
                    counts["alignments_dropped_overhang"] += 1
                    dropped.append((aln.query_name, start, end, "+".join(sides)))
                    continue
                counts["alignments_emitted"] += 1
                # Per orientation as well as in total, because for a strandless
                # chunk "both strands are in here" is the property the downstream
                # split depends on, and a single total cannot state it.
                if aln.is_forward:
                    counts["alignments_emitted_forward"] += 1
                else:
                    counts["alignments_emitted_reverse"] += 1
                if bamwriter is not None:
                    bamwriter.write(
                        _rebase_alignment(
                            aln, offset, mini_header, mini_contig_name, mini_length
                        )
                    )

    if mixed_orientation_offender is not None:
        # Nothing partial survives: a chunk bam holding the reads seen before the
        # offending one is not a smaller correct answer, it is a truncated one.
        # bam_filename is the SOURCE under --reuse_source_bam, and this branch is
        # unreachable there (it needs a strand-suffixed region, which reuse
        # refuses). Named explicitly anyway: the cost of getting that wrong once
        # is somebody's input bam.
        removable = [fasta_filename, fasta_filename + ".fai"]
        if not reuse_source_bam:
            removable.append(bam_filename)
        for path in removable:
            if os.path.exists(path):
                os.remove(path)
        raise ExtractionError(
            "REFUSED: region {}{}:{}-{} asks for one orientation but {} still "
            "holds the other (e.g. {}), so this bam has NOT been strand-separated. "
            "separate_bam_by_strand.py REWRITES is_reverse on every read whose "
            "inferred transcribed strand disagrees with the aligner, and the "
            "strand filter here reads the RAW flag -- filtering before that "
            "rewrite assigns every flipped read to the wrong chunk and still "
            "produces output that looks fine. Extract STRANDLESSLY (drop the {} "
            "suffix) and split the chunk afterwards, or point --bam at the "
            "strand-separated bam. Pass --mixed_orientation_source warn only if "
            "you mean to filter on aligner orientation. Nothing was kept.".format(
                region.chrom,
                region.strand,
                region.lend,
                region.rend,
                bam,
                mixed_orientation_offender,
                region.strand,
            )
        )

    if counts["opposite_orientation_seen"]:
        print(
            "NOTE: {} retained primary alignment(s) in {}{}:{}-{} carry the "
            "opposite orientation and were excluded by the strand filter. If this "
            "bam has not been strand-separated, that filter reads pre-flip flags "
            "and the exclusion is wrong; count recorded in the manifest.".format(
                counts["opposite_orientation_seen"],
                region.chrom,
                region.strand,
                region.lend,
                region.rend,
            ),
            file=sys.stderr,
        )

    if (
        not region.strand
        and counts["alignments_emitted"]
        and not (
            counts["alignments_emitted_forward"]
            and counts["alignments_emitted_reverse"]
        )
    ):
        # Reported, not refused: a sparse region can legitimately hold one
        # orientation. Over a real chunk it means the source was already split,
        # so the downstream per-chunk split will produce an empty partner.
        print(
            "NOTE: strandless chunk {}:{}-{} emitted {} forward and {} reverse "
            "alignment(s). A single orientation means the source bam was already "
            "strand-separated, so this chunk does not carry both strands.".format(
                region.chrom,
                region.lend,
                region.rend,
                counts["alignments_emitted_forward"],
                counts["alignments_emitted_reverse"],
            ),
            file=sys.stderr,
        )

    dropped_reads_filename = "{}.dropped_reads.txt".format(output_prefix)
    dropped_names = sorted({name for name, _, _, _ in dropped})
    with open(dropped_reads_filename, "wt") as ofh:
        for name in dropped_names:
            print(name, file=ofh)

    if dropped:
        print(
            "NOTE: dropped {} primary alignment(s) overhanging {}{}:{}-{}; names "
            "in {} (e.g. {}). Reason: the alignment is not contained by the "
            "chunk, and this extractor neither truncates nor widens.".format(
                len(dropped),
                region.chrom,
                region.strand,
                region.lend,
                region.rend,
                dropped_reads_filename,
                ", ".join(
                    "{} {}-{} ({})".format(*record) for record in dropped[:3]
                ),
            ),
            file=sys.stderr,
        )

    if counts["nonprimary_seen"]:
        message = (
            "{} secondary/supplementary alignment(s) overlap {}{}:{}-{} "
            "(e.g. {}).".format(
                counts["nonprimary_seen"],
                region.chrom,
                region.strand,
                region.lend,
                region.rend,
                ", ".join(nonprimary),
            )
        )
        if secondary_alignments == "reject":
            raise ExtractionError(
                message + " Primary alignments only are supported; rerun with "
                "--secondary_alignments exclude to drop them explicitly."
            )
        print(
            "NOTE: " + message + " Excluded from output; count recorded in the "
            "manifest.",
            file=sys.stderr,
        )

    if reuse_source_bam:
        # The source was fetched by region above, so it is already indexed -- and
        # indexing it again would write into the input directory, which is often
        # read-only and never this tool's to touch.
        _require_source_index(bam_filename)
    else:
        pysam.index(bam_filename)

    emitted_transcript_ids = []
    gtf_filename = None
    if gtf:
        gtf_filename = "{}.gtf".format(output_prefix)
        contained = sorted(
            annotation.genes_contained(region.lend, region.rend),
            key=lambda g: (g.lend, g.rend, g.gene_id),
        )
        with open(gtf_filename, "wt") as ofh:
            for gene in contained:
                for line in gene.lines:
                    print(
                        _rebase_gtf_line(
                            line,
                            offset,
                            mini_contig_name,
                            mini_length,
                            "gene {}".format(gene.gene_id),
                        ),
                        file=ofh,
                    )
                    counts["gtf_lines_emitted"] += 1
                for transcript_id in sorted(
                    gene.transcript_ids,
                    key=lambda t: (annotation.transcripts[t].lend, t),
                ):
                    transcript = annotation.transcripts[transcript_id]
                    # A transcript is emitted whole or not at all: every line it
                    # owns is rebased as one group, never filtered per-line, so a
                    # truncated model cannot be produced.
                    rebased = [
                        _rebase_gtf_line(
                            line,
                            offset,
                            mini_contig_name,
                            mini_length,
                            "transcript {}".format(transcript_id),
                        )
                        for line in transcript.lines
                    ]
                    for line in rebased:
                        print(line, file=ofh)
                    emitted_transcript_ids.append(transcript_id)
                    counts["gtf_transcripts_emitted"] += 1
                    counts["gtf_lines_emitted"] += len(rebased)

    manifest = {
        "chrom": region.chrom,
        # null, not "", when the chunk carries both orientations: a consumer
        # asking "which strand is this" must get an answer that cannot be
        # mistaken for one, and "" reads as a strand in a format string.
        "strand": region.strand or None,
        # The ordering constraint, stated in the artifact rather than left to the
        # caller's memory: this bam has NOT been strand-separated, so anything
        # that filters on orientation must run the split FIRST. is_reverse here
        # is still the aligner's, and the split rewrites it.
        "strand_split_required": not region.strand,
        "partition_lend": region.lend,
        "partition_rend": region.rend,
        "offset": offset,
        "margin": margin,
        # what to pass normalize_bam_by_strand.py --window_origin for this chunk:
        # the rebase offset in the GENOME frame, unreduced. Absolute 0-based ==
        # chunk-local 0-based + offset, so the chunk's depth-window grid is the
        # absolute grid. Every boundary is a multiple of depth_window, hence so
        # is this, hence no window straddles a cut.
        "window_origin": offset,
        "disjoint_hard_cut": True,
        # Whether files.bam is a mini bam this run WROTE or the SOURCE it read.
        # A consumer that splits, normalizes or quantifies the chunk needs to
        # know, because a reused source is not restricted to this contig: it
        # holds every other one too, and the reader has to say which it wants.
        "bam_reused_from_source": bool(reuse_source_bam),
        # True only when the chunk covers the whole contig, which is what makes
        # the reuse legal in the first place, and stated separately so the
        # condition can be checked rather than inferred from the flag.
        "spans_whole_contig": region.lend == 1 and region.rend == contig_length,
        "mini_contig_name": mini_contig_name,
        "mini_length": mini_length,
        "contig_length": contig_length,
        "emitted_gene_ids": sorted(
            g.gene_id for g in annotation.genes_contained(region.lend, region.rend)
        )
        if gtf
        else [],
        "emitted_transcript_ids": sorted(emitted_transcript_ids),
        "counts": {
            "alignments_emitted": counts["alignments_emitted"],
            "alignments_emitted_forward": counts["alignments_emitted_forward"],
            "alignments_emitted_reverse": counts["alignments_emitted_reverse"],
            "alignments_dropped_overhang": counts["alignments_dropped_overhang"],
            "opposite_orientation_excluded": counts["opposite_orientation_seen"],
            "nonprimary_excluded": counts["nonprimary_seen"],
            "gtf_transcripts_emitted": counts["gtf_transcripts_emitted"],
            "gtf_lines_emitted": counts["gtf_lines_emitted"],
        },
        "dropped_read_names": dropped_names,
        "files": {
            "fasta": fasta_filename,
            "bam": bam_filename,
            "gtf": gtf_filename,
            "dropped_reads": dropped_reads_filename,
        },
    }

    manifest_filename = "{}.partition.json".format(output_prefix)
    with open(manifest_filename, "wt") as ofh:
        json.dump(manifest, ofh, indent=2, sort_keys=True)
        print("", file=ofh)
    manifest["files"]["manifest"] = manifest_filename

    return manifest


def main(argv=None):

    parser = argparse.ArgumentParser(
        description="extract mini contig, bam and gtf for one chunk of a "
        "contig-strand. A boundary must lie outside every annotated locus, "
        "which cannot be dropped; an alignment that would extend past a "
        "boundary IS dropped, counted and named rather than truncated.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--genome_fa", type=str, required=True, help="genome fasta file"
    )

    parser.add_argument("--bam", type=str, required=True, help="bam aligned reads")

    parser.add_argument("--gtf", type=str, required=False, help="gtf annotation")

    parser.add_argument(
        "--region",
        type=str,
        required=True,
        help="the chunk, formatted chrom[+-]:lend-rend. Neither boundary may cut "
        "an annotated locus; choose them with select_contig_cut_points.py, which "
        "screens positions against the same rule this tool enforces.",
    )

    parser.add_argument(
        "--output_prefix",
        type=str,
        required=False,
        default=None,
        help="prefix for output files. Default: region string.",
    )

    parser.add_argument(
        "--margin",
        type=int,
        default=DEFAULT_MARGIN,
        help="clearance demanded either side of a cut, in bp. 0 means only that "
        "nothing crosses it; the default is 4x the largest boundary-snapping "
        "distance in LRAA_Globals (50 bp) so no snapping can reach across a cut.",
    )

    parser.add_argument(
        "--secondary_alignments",
        choices=("exclude", "reject"),
        default="exclude",
        help="secondary/supplementary alignments are out of scope: 'exclude' "
        "drops them from the output and records the count, "
        "'reject' fails if any are present.",
    )

    parser.add_argument(
        "--mini_contig_name",
        type=str,
        default=None,
        help="name of the extracted contig. Default: the source contig name.",
    )

    parser.add_argument(
        "--gtf_index_cache_dir",
        type=str,
        default=None,
        help="where to put the tabix index of --gtf when the GTF's own "
        "directory is not writable. The GTF itself is never modified.",
    )

    parser.add_argument(
        "--max_intron_length",
        type=int,
        default=LRAA_Globals.config["max_intron_length"],
        help="alignments carrying a longer intron are not extracted, matching "
        "the strand split and cut selection. Pass the value the RUN uses: a "
        "strandless chunk is cut from the raw bam, so a mismatch here retains a "
        "different record set than the selector costed. 0 disables.",
    )

    parser.add_argument(
        "--mixed_orientation_source",
        choices=("refuse", "warn"),
        default="refuse",
        help="what a STRAND-SUFFIXED --region does when the bam still holds the "
        "other orientation, i.e. has not been strand-separated. 'refuse' is the "
        "default because separate_bam_by_strand.py rewrites is_reverse and this "
        "tool's strand filter reads the raw flag, so filtering first assigns "
        "every flipped read to the wrong chunk while looking fine. Extract "
        "strandlessly and split afterwards. 'warn' only for deliberately "
        "filtering on aligner orientation.",
    )

    parser.add_argument(
        "--reuse_source_bam",
        action="store_true",
        default=False,
        help="for a STRANDLESS region spanning the whole contig, do not write a "
        "mini bam: count the region as usual and name --bam itself as the "
        "chunk's bam. At offset 0 over the full contig the mini bam is the "
        "source's own records unmoved, so the copy is pure cost -- measured 67.5 s "
        "against 5.45 s for the count alone on a 1.2 M-record contig. The "
        "manifest sets bam_reused_from_source, and the reader MUST restrict the "
        "source by contig: it holds every other contig too. Refused for a "
        "partial or strand-suffixed region.",
    )

    args = parser.parse_args(argv)
    region = parse_region(args.region)
    output_prefix = args.output_prefix
    if output_prefix is None:
        output_prefix = args.region.replace(",", "")

    manifest = extract_partition(
        genome_fa=args.genome_fa,
        bam=args.bam,
        region=region,
        output_prefix=output_prefix,
        gtf=args.gtf,
        margin=args.margin,
        secondary_alignments=args.secondary_alignments,
        mini_contig_name=args.mini_contig_name,
        gtf_index_cache_dir=args.gtf_index_cache_dir,
        mixed_orientation_source=args.mixed_orientation_source,
        max_intron_length=args.max_intron_length,
        reuse_source_bam=args.reuse_source_bam,
    )

    print(
        "chunk {}{}:{}-{} -> mini contig {} bp; {} primary alignments {}, "
        "{} dropped for overhang; {} transcripts emitted whole; margin {} bp".format(
            manifest["chrom"],
            manifest["strand"] or "",
            manifest["partition_lend"],
            manifest["partition_rend"],
            manifest["mini_length"],
            manifest["counts"]["alignments_emitted"],
            "counted, bam reused from the source"
            if manifest["bam_reused_from_source"]
            else "emitted",
            manifest["counts"]["alignments_dropped_overhang"],
            manifest["counts"]["gtf_transcripts_emitted"],
            manifest["margin"],
        ),
        file=sys.stderr,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
