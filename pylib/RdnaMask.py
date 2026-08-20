#!/usr/bin/env python3
"""rDNA-cassette read masking: exclude alignments landing on ribosomal-DNA-repeat
sequence before they ever reach splice-graph construction or quantification.

## Why this exists

A handful of reference loci (GRCh38's acrocentric chromosome short arms, plus
several unplaced/decoy scaffolds) carry collapsed, partial copies of the human
rDNA repeat unit. Reads genuinely derived from ribosomal RNA pile onto these loci
at extreme depth (hundreds of thousands of MAPQ-0 multi-mapping reads at a single
locus is observed in practice), driving a combinatorial explosion in candidate
multipaths and isoform assignment cost that is completely disproportionate to the
locus's real isoform content -- these loci consistently show the LOWEST discovered
yield of any comparably-costly locus in a run.

GTF annotation cannot identify these reads: GENCODE/Ensembl annotate only the
small mature-rRNA fragments (e.g. "5.8S") as genes at these positions, leaving the
much larger 18S/28S/external-and-internal-transcribed-spacer footprint -- where
the actual read pileup sits -- entirely unannotated. Sequence homology to the rDNA
repeat unit, not genome annotation, is what identifies these reads, hence a
dedicated masking step rather than a `gene_type` filter.

## What it does

At startup, LRAA aligns a bundled (or user-supplied) rDNA cassette FASTA against
`--genome` with `minimap2`. Every resulting hit span, padded by
``config["rdna_mask_pad"]`` bases, becomes an excluded region. Any alignment
whose reference blocks overlap an excluded region is treated exactly like any
other alignment `Util_funcs.quant_discard_reason` rejects (see there): it never
enters coverage normalization, splice-graph construction, or read-to-transcript
assignment.

## Where it runs

Once per LRAA invocation, against whatever `--genome` that invocation was given --
the whole genome for an unchunked run, or one chunk's own mini-genome under
`--chunk` (a separate `LRAA` process per chunk already runs this independently,
so the alignment is always local to what that process actually needs to search).
Building the mask is cheap relative to the run it protects: aligning the ~43 kb
cassette against a 10 Mb chunk genome is sub-second; against a full ~3.1 Gb human
genome it is on the order of a minute, still a rounding error in a run measured in
hours. See ``build_rdna_mask_bed`` for the on-disk cache that also avoids paying
the cost on every consumer of the same genome fasta.
"""

import hashlib
import logging
import os
import re
import subprocess
import tempfile

import intervaltree

import LRAA_Globals

# Util_funcs is NOT imported at module level: it imports THIS module to reach
# read_overlaps_mask from quant_discard_reason, so importing it back here would
# be circular. Deferred into the one function that needs it (_mask_cache_key)
# instead, so neither module has to care which one loads first.

logger = logging.getLogger(__name__)

# Repo-relative default: resources/human_rDNA_cassette.fa, alongside this file's
# own package (pylib/../resources), so it resolves the same way regardless of the
# caller's cwd -- exactly the convention LRAA itself uses for UTILDIR.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
DEFAULT_RDNA_CASSETTE_FASTA = os.path.join(
    _REPO_ROOT, "resources", "human_rDNA_cassette.fa"
)

# Bumped whenever the algorithm that turns a cassette-vs-genome alignment into a
# mask changes (the minimap2 preset, the padding rule, the merge rule). Named
# alongside the two file identities in every cache key below, for the same reason
# normalize_bam_by_strand.SPLICE_GRAPH_NORMALIZATION_METHOD is: nothing downstream
# can detect a stale mask on its own, so the key has to change instead.
MASK_METHOD_VERSION = "v2"


def resolve_rdna_fasta(rdna_mask_fasta_arg):
    """The cassette FASTA to align, or None if masking is disabled by argument.

    ``rdna_mask_fasta_arg`` is ``args.rdna_mask_fasta``: None means "use the
    bundled default", any other string is a caller-supplied override (a
    different organism's rDNA unit, for instance). Existence is checked here so
    a typo'd path fails at startup with a clear message rather than surfacing as
    an empty mask deep inside a chunk worker's log.
    """
    path = rdna_mask_fasta_arg or DEFAULT_RDNA_CASSETTE_FASTA
    path = os.path.realpath(path)
    if not os.path.exists(path):
        raise FileNotFoundError(
            "rDNA mask FASTA not found: {}. Pass an existing file to "
            "--rdna_mask_fasta, or omit it to use the bundled default at {}.".format(
                path, DEFAULT_RDNA_CASSETTE_FASTA
            )
        )
    return path


_CIGAR_OP_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_REF_CONSUMING_OPS = frozenset("MDN=X")


def _sam_hit_identity(fields):
    """Percent identity of one SAM record, or None if unmeasurable (no NM tag).

    Same computation as Util_funcs.alignment_per_id, reimplemented on raw SAM
    fields rather than a pysam.AlignedSegment: this parses minimap2's SAM
    output directly, before any bam/pysam object exists. NM is mismatches +
    inserted + deleted bases (edit distance), so this is the same "aligned
    bases minus edits, over aligned bases" identity every other consumer in
    this codebase reads off NM.
    """
    cigar = fields[5]
    aligned = sum(
        int(n) for n, op in _CIGAR_OP_RE.findall(cigar) if op in "MDN=X"
    )
    if aligned == 0:
        return None
    for tag in fields[11:]:
        if tag.startswith("NM:i:"):
            mismatches = int(tag[5:])
            return 100.0 * (aligned - mismatches) / aligned
    return None


def _sam_hit_spans(sam_path, min_length=0, min_per_id=0.0):
    """(contig, ref_start, ref_end) for every mapped record clearing the floors.

    0-based half-open, matching every other interval in this codebase.
    Unmapped records (minimap2 -a's placeholder for a query with no
    acceptable alignment) carry ``cigar == "*"`` and are skipped.

    Both floors guard the same failure mode: a short, low-identity hit that is
    coincidental rather than truly rDNA-homologous would otherwise become an
    excluded region indistinguishable from a real repeat-unit copy, at
    whatever scale the genome happens to produce chance matches. min_length
    alone would not catch a LONG low-identity alignment (asm20 tolerates
    substantial divergence by design, which is wanted for genuinely divergent
    paralogous copies but not for a long run of low-complexity sequence that
    merely resembles the cassette in the aggregate); min_per_id alone would
    not catch a SHORT, coincidentally-perfect match. Records with no NM tag
    (identity unmeasurable) are kept if long enough on length alone, matching
    Util_funcs.quant_discard_reason's own "absence of evidence is not
    evidence of a bad alignment" rule for identity floors.
    """
    spans = []
    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            rname, pos, cigar = fields[2], int(fields[3]), fields[5]
            if cigar == "*" or rname == "*":
                continue
            ref_len = sum(
                int(n) for n, op in _CIGAR_OP_RE.findall(cigar) if op in _REF_CONSUMING_OPS
            )
            if ref_len < min_length:
                continue
            if min_per_id > 0:
                per_id = _sam_hit_identity(fields)
                if per_id is not None and per_id < min_per_id:
                    continue
            spans.append((rname, pos - 1, pos - 1 + ref_len))
    return spans


def _merge_spans(spans, pad):
    """Pad each span and merge overlaps, per contig. Input order need not be sorted."""
    by_contig = {}
    for contig, s, e in spans:
        by_contig.setdefault(contig, []).append((max(0, s - pad), e + pad))

    merged = {}
    for contig, contig_spans in by_contig.items():
        contig_spans.sort()
        out = []
        for s, e in contig_spans:
            if out and s <= out[-1][1]:
                out[-1] = (out[-1][0], max(out[-1][1], e))
            else:
                out.append((s, e))
        merged[contig] = out
    return merged


def _mask_cache_key(genome_fasta, rdna_fasta, pad, min_length, min_per_id):
    import Util_funcs  # deferred: see the module docstring's note on the cycle

    genome_id = Util_funcs.file_identity_token(genome_fasta)
    rdna_id = Util_funcs.file_identity_token(rdna_fasta)
    digest = hashlib.blake2s(
        "{}|{}|{}|{}|{}|{}".format(
            MASK_METHOD_VERSION, genome_id, rdna_id, pad, min_length, min_per_id
        ).encode("utf-8"),
        digest_size=8,
    ).hexdigest()
    return digest


def build_rdna_mask_bed(
    genome_fasta,
    rdna_fasta,
    cache_dir="__rdna_mask_cache",
    pad=None,
    min_length=None,
    min_per_id=None,
    threads=1,
):
    """Align ``rdna_fasta`` against ``genome_fasta``; write/reuse a BED of hits.

    Checkpointed exactly like every other cached artifact in this codebase: the
    BED's name carries a digest of (mask method, genome identity, cassette
    identity, padding, length floor, identity floor), so a stale mask from a
    different genome, cassette, or quality floor can never be mistaken for a
    current one, and a second caller building the mask for the SAME genome
    (the per-chunk LRAA invocation and this chunk's own stage-4
    normalization, for instance) gets a cache hit instead of a second
    minimap2 run.

    ``min_length``/``min_per_id`` reject a hit before it can become part of
    the mask: without them, a single short, coincidentally-homologous
    alignment -- unavoidable at genome scale -- would exclude a region
    indistinguishable from a real rDNA-repeat-unit copy. See _sam_hit_spans
    for why both floors are needed together.

    Returns the BED path, or None if the cassette-vs-genome alignment found no
    hits at all (an ordinary genome with no rDNA-homologous sequence) -- an
    empty mask is not an error, and callers should treat None exactly like "no
    masking configured" rather than trying to load a BED that was never
    written.
    """
    if pad is None:
        pad = int(LRAA_Globals.config.get("rdna_mask_pad", 500))
    if min_length is None:
        min_length = int(LRAA_Globals.config.get("rdna_mask_min_hit_length", 200))
    if min_per_id is None:
        min_per_id = float(LRAA_Globals.config.get("rdna_mask_min_per_id", 80))

    # Not an error: a caller with no chunk genome on disk yet (a scheduling test
    # that never extracts one, for instance) gets "nothing to mask" rather than a
    # crash. Any REAL absence fails loudly and more informatively moments later,
    # at the pysam.FastaFile open every consumer of this genome also needs.
    if not os.path.exists(genome_fasta):
        logger.info(
            "rDNA mask: %s does not exist -- skipping (nothing to mask)",
            genome_fasta,
        )
        return None

    os.makedirs(cache_dir, exist_ok=True)
    key = _mask_cache_key(genome_fasta, rdna_fasta, pad, min_length, min_per_id)
    bed_path = os.path.join(cache_dir, "rdna_mask.{}.bed".format(key))
    checkpoint_path = bed_path + ".ok"
    empty_checkpoint_path = bed_path + ".empty.ok"

    if os.path.exists(empty_checkpoint_path):
        return None
    if os.path.exists(checkpoint_path):
        return bed_path

    # Two concurrent callers can both see "not done" above and both build (a
    # strandless chunk's two orientations do exactly this, against the SAME
    # chunk genome, on their own THREADS, not processes -- see chunk_worker).
    # tempfile.mkstemp, not a hand-rolled pid/id()-based name: two threads in
    # one process share a pid, and two objects created microseconds apart can
    # share an id() too, once CPython's allocator reuses the first one's
    # address the instant its refcount drops to zero. mkstemp's O_CREAT|O_EXCL
    # is the only uniqueness guarantee that holds across threads, processes,
    # and hosts.
    # Only the final bed_path and checkpoint are shared, and both are written
    # by a single os.replace/touch, so the loser of the race is a harmless,
    # byte-identical overwrite rather than a corrupted read.
    sam_fd, sam_path = tempfile.mkstemp(
        dir=cache_dir, prefix="rdna_mask.{}.".format(key), suffix=".sam"
    )
    os.close(sam_fd)
    tmp_bed_fd, tmp_bed_path = tempfile.mkstemp(
        dir=cache_dir, prefix="rdna_mask.{}.".format(key), suffix=".bed.tmp"
    )
    os.close(tmp_bed_fd)
    cmd = [
        "minimap2",
        "-a",
        "-x",
        "asm20",
        "--secondary=yes",
        "-p",
        "0.05",
        "-N",
        "50",
        "-t",
        str(max(1, int(threads))),
        genome_fasta,
        rdna_fasta,
    ]
    logger.info(
        "rDNA mask: aligning %s against %s (pad=%d bp)",
        os.path.basename(rdna_fasta),
        os.path.basename(genome_fasta),
        pad,
    )
    try:
        with open(sam_path, "w") as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)

        spans = _sam_hit_spans(sam_path, min_length=min_length, min_per_id=min_per_id)

        if not spans:
            logger.info(
                "rDNA mask: no hits against %s -- masking is a no-op for this genome",
                os.path.basename(genome_fasta),
            )
            open(empty_checkpoint_path, "w").close()
            return None

        merged = _merge_spans(spans, pad)
        total_regions = sum(len(v) for v in merged.values())
        total_bp = sum(e - s for v in merged.values() for s, e in v)
        with open(tmp_bed_path, "w") as out:
            for contig in sorted(merged):
                for s, e in merged[contig]:
                    out.write("{}\t{}\t{}\n".format(contig, s, e))
        os.replace(tmp_bed_path, bed_path)
        logger.info(
            "rDNA mask: %d region(s), %d bp total, across %d contig(s)",
            total_regions,
            total_bp,
            len(merged),
        )
        open(checkpoint_path, "w").close()
        return bed_path
    finally:
        for path in (sam_path, tmp_bed_path):
            try:
                os.remove(path)
            except OSError:
                pass


def load_mask_bed(bed_path):
    """A BED file (as ``build_rdna_mask_bed`` writes) into {contig: IntervalTree}.

    None in, None out: the caller-facing way to say "no mask configured", kept
    distinct from an IntervalTree that happens to be empty so a masking
    predicate can short-circuit on identity rather than iterate nothing.
    """
    if bed_path is None:
        return None
    mask = {}
    with open(bed_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            contig, s, e = line.split("\t")
            s, e = int(s), int(e)
            mask.setdefault(contig, intervaltree.IntervalTree())[s:e] = True
    return mask


def read_overlaps_mask(read, mask, min_overlap_bp=None):
    """True if ``read`` overlaps ``mask`` by at least ``min_overlap_bp``.

    ``mask`` is a {contig: IntervalTree} as ``load_mask_bed`` returns, or None
    (no mask configured -- always False). A read entirely on a contig absent
    from the mask is the common case for an ordinary genome and is rejected by
    the dict lookup alone, before any interval query runs.

    ``min_overlap_bp`` guards a case symmetric to the hit-quality floors in
    build_rdna_mask_bed: every excluded region is padded (config
    "rdna_mask_pad", default 500 bp) specifically to absorb alignment-boundary
    slop around a genuine rDNA-cassette hit, so a read whose alignment
    extends a handful of bases into that padding is far more likely an
    ordinary read from the adjacent unique sequence than one actually
    implicated in the locus's multi-mapping ambiguity -- discarding it
    outright would cost real read support for genes bordering a masked
    region for no corresponding gain in the mask's purpose. A read from
    genuinely inside a masked repeat copy overlaps by its whole aligned
    length and clears any threshold this small trivially, so the floor
    only ever spares boundary-adjacent reads, never the reads the mask
    exists to catch. None reads config["rdna_mask_min_overlap_bp"], the
    same "None defers to config" convention every other threshold here
    follows.
    """
    if not mask:
        return False
    tree = mask.get(read.reference_name)
    if not tree:
        return False
    if min_overlap_bp is None:
        min_overlap_bp = int(LRAA_Globals.config.get("rdna_mask_min_overlap_bp", 50))
    overlap = 0
    for block_start, block_end in read.get_blocks():
        for iv in tree.overlap(block_start, block_end):
            overlap += min(block_end, iv.end) - max(block_start, iv.begin)
            if overlap >= min_overlap_bp:
                return True
    return False
