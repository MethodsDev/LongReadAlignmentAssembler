#!/usr/bin/env python3

"""The equivalence gate for strandless chunking.

Strandless chunking cuts a contig BEFORE the strand split, so the split runs
inside each chunk in parallel instead of once over the whole BAM.  It is only
allowed if it computes the same answer as the strand-first arm, and "the same
answer" has two separable meanings.  They are checked by two subcommands, and the
cheap one comes first because it can kill the design on its own.

``strand-invariant`` -- the claim the whole design rests on, checkable without
quantifying anything.  Strand assignment is asserted to be PER-READ LOCAL:
``split_bam_by_strand`` decides a read's strand from that read alone, so the
strand it gets from a whole-contig split must equal the strand it gets from a
split run inside the chunk that contains it.  This subcommand extracts BOTH
assignments and compares them as a per-read mapping.  If they disagree, no
downstream number can be made to agree and the design is dead.  If they agree,
every remaining difference is partitioning, not strand.

    The assignments are READ BACK from the splitter's own two output BAMs rather
    than recomputed here.  A reimplementation would be a second parity notion and
    could agree with itself while disagreeing with the pipeline.  Which BAM the
    splitter wrote a record into IS its verdict, so that is what is compared.

``compare-arms`` -- two finished pipeline runs over the same substrate, one
strand-first and one strandless, compared at three levels: ``quant.expr`` byte
identity, ``quant.tracking`` row identity, and the per-read assignment sets.

    Byte identity is the bar, but the two arms select their cuts from different
    evidence -- one strand's coverage against both strands' union -- so their cut
    POSITIONS may differ, and different positions sever different reads.  A read
    severed in one arm and not the other is a legitimate input difference, not a
    defect.  So when the arms differ this subcommand asks the only question that
    distinguishes the two cases: is every differing row attributable to a read in
    the symmetric difference of the two severed sets?  If yes the divergence is
    the severing and nothing else.  If a transcript's numbers move with no
    severed-set read behind it, that is a real defect and it is named.

    Reads present in one arm and absent from the other are reported BY NAME with
    a reason resolved from the run's own artifacts, never as a count: the severed
    accounting differs between the two designs by construction, so a count
    difference cannot tell anyone whether it is correct.

THE ORDERING CONSTRAINT.  ``split_bam_by_strand`` REWRITES ``read.is_reverse``
when the strand it infers disagrees with the aligner, and the extractor's
``_strand_matches`` reads the raw flag.  So a chunk must be extracted
STRANDLESSLY and split AFTERWARDS; extracting ``chr1+:...`` from a raw BAM
assigns every flipped read to the wrong strand and still produces output that
looks fine.  ``strand-invariant`` is the check that catches it: run against a
strand-filtered extraction it reports exactly the flipped reads as disagreements.
``pylib/test_strandless_parity.py`` holds it to that.
"""

import argparse
import collections
import glob
import importlib.util
import json
import os
import subprocess
import sys

from importlib.machinery import SourceFileLoader

import pysam

_HERE = os.path.dirname(os.path.realpath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_HERE, ".."))

sys.path.insert(0, _HERE)
import ChunkedRun  # noqa: E402  (path insert must precede the import)
import LRAA_Globals  # noqa: E402


def _load_script(*relpath):
    """Import a repo script that is not importable by name."""

    path = os.path.join(REPO_ROOT, *relpath)
    loader = SourceFileLoader("strandless_parity_" + relpath[-1][:-3], path)
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


splitter = _load_script("util", "separate_bam_by_strand.py")
extractor = _load_script("util", "misc", "extract_contig_region_inputs.py")

SELECT_CUTS = os.path.join(REPO_ROOT, "util", "misc", "select_contig_cut_points.py")

# quant.tracking's mp_id is a per-process MultiPath counter, not data: two runs of
# the same input number their multipaths differently, so it is excluded from every
# content comparison here.  ChunkedRun.merge_and_translate says the same.
INCOMPARABLE_TRACKING_COLUMNS = ("mp_id",)


class ParityError(RuntimeError):
    pass


# ------------------------------------------------------------------ assignments


def split_assignments(
    bam,
    out_prefix,
    max_intron_length,
    infer_read_orient=False,
    genome_fa=None,
    gtf=None,
    offset=0,
):
    """Run the production strand splitter over ``bam`` and read back its verdict.

    Returns ``(assignments, counters, collisions)`` where ``assignments`` maps
    ``(read_name, reference_start)`` to ``'+'`` or ``'-'``.

    ``offset`` is added to every start so a chunk's rebased coordinates come back
    in the source frame; ``extract_contig_region_inputs`` rebases by subtracting
    exactly this value, so the sum is the original 0-based start.  The position is
    in the key because a read NAME is not a unique record identifier in general
    and silently overwriting one assignment with another would hide precisely the
    kind of disagreement this function exists to find.  Genuine collisions are
    returned rather than tolerated.
    """

    chrom_to_itree = None
    if infer_read_orient:
        if genome_fa is None:
            raise ParityError("inferred strand needs a genome FASTA")
        if gtf:
            chrom_to_itree = splitter.build_chrom_itrees(gtf)

    plus_bam = "{}.+.bam".format(out_prefix)
    minus_bam = "{}.-.bam".format(out_prefix)
    counters = splitter.split_bam_by_strand(
        bam,
        plus_bam,
        minus_bam,
        max_intron_length,
        infer_read_orient_flag=infer_read_orient,
        genome_fasta=genome_fa,
        chrom_to_itree=chrom_to_itree,
    )

    assignments = {}
    collisions = []
    for strand, path in (("+", plus_bam), ("-", minus_bam)):
        with pysam.AlignmentFile(path, "rb", check_sq=False) as fh:
            for aln in fh.fetch(until_eof=True):
                key = (aln.query_name, aln.reference_start + offset)
                if key in assignments:
                    collisions.append(
                        {
                            "read_name": key[0],
                            "reference_start": key[1],
                            "first_strand": assignments[key],
                            "second_strand": strand,
                        }
                    )
                    continue
                assignments[key] = strand
    return assignments, counters, collisions


def raw_flag_strands(bam, names):
    """``(name, start)`` -> aligner-flag strand, for the named reads only.

    Used to explain a disagreement: the splitter rewrites the flag when it infers
    a different strand, so the raw flag is only recoverable from the source BAM.
    """

    wanted = set(names)
    out = {}
    with pysam.AlignmentFile(bam, "rb") as fh:
        for aln in fh.fetch(until_eof=True):
            if aln.query_name in wanted:
                out[(aln.query_name, aln.reference_start)] = (
                    "+" if aln.is_forward else "-"
                )
    return out


# ---------------------------------------------------------------------- regions


Interval = collections.namedtuple("Interval", "chrom lend rend")


def intervals_from_cuts_json(path):
    """The selector's own ``<prefix>.cuts.json``, whatever strand produced it."""

    with open(path, "rt") as fh:
        selections = json.load(fh)
    intervals = []
    for selection in selections:
        for segment in selection["segments"]:
            intervals.append(
                Interval(selection["chrom"], segment["lend"], segment["rend"])
            )
    return intervals


def intervals_from_region_strings(region_strings):
    intervals = []
    for text in region_strings:
        region = extractor.parse_region(text)
        if region.strand:
            raise ParityError(
                "region {} names a strand. The invariant is about a STRANDLESS "
                "extraction: a strand-filtered extraction from a raw BAM filters "
                "on the aligner's flag before the splitter can rewrite it, which "
                "is the ordering mistake this check exists to catch. Drop the "
                "strand.".format(text)
            )
        intervals.append(Interval(region.chrom, region.lend, region.rend))
    return intervals


def chunk_dirs(chunks_root):
    """Chunk directories under ``chunks_root``, in name order."""

    return sorted(
        d
        for d in glob.glob(os.path.join(chunks_root, "*"))
        if os.path.exists(os.path.join(d, "chunk.partition.json"))
    )


def check_disjoint(intervals):
    """Overlapping intervals would double-count a read. Returns the overlaps."""

    overlaps = []
    by_chrom = collections.defaultdict(list)
    for interval in intervals:
        by_chrom[interval.chrom].append(interval)
    for chrom, items in by_chrom.items():
        items.sort(key=lambda i: (i.lend, i.rend))
        for left, right in zip(items, items[1:]):
            if right.lend <= left.rend:
                overlaps.append(
                    "{}:{}-{} overlaps {}:{}-{}".format(
                        chrom, left.lend, left.rend, chrom, right.lend, right.rend
                    )
                )
    return overlaps


# ----------------------------------------------------------- the invariant check


def check_strand_invariant(
    bam,
    genome_fa,
    work_dir,
    intervals=None,
    chunks_root=None,
    gtf=None,
    margin=None,
    max_intron_length=None,
    infer_read_orient=False,
    contig=None,
    keep_bams=False,
    gtf_index_cache_dir=None,
):
    """Whole-scope strand assignment against in-chunk strand assignment, per read.

    Exactly one of ``intervals`` (extracted here) or ``chunks_root`` (a real
    strandless run's ``chunks/`` directory, consumed as it stands) is used.
    Consuming a real run is the stronger check because it tests the artifact the
    pipeline actually produced; extracting here is self-contained and is what the
    unit tests use.

    Returns a report dict.  ``agreement`` is exact: ``agreed / compared``.
    """

    if (intervals is None) == (chunks_root is None):
        raise ParityError("give either intervals or chunks_root, not both")
    if max_intron_length is None:
        max_intron_length = LRAA_Globals.config["max_intron_length"]
    if margin is None:
        margin = extractor.DEFAULT_MARGIN
    os.makedirs(work_dir, exist_ok=True)

    # ---- the whole-scope side, computed with no knowledge of any chunk.
    whole, whole_counters, whole_collisions = split_assignments(
        bam,
        os.path.join(work_dir, "whole.strand"),
        max_intron_length,
        infer_read_orient=infer_read_orient,
        genome_fa=genome_fa,
        gtf=gtf,
    )

    # ---- the per-chunk side.
    chunks = []
    if chunks_root is not None:
        for cdir in chunk_dirs(chunks_root):
            with open(os.path.join(cdir, "chunk.partition.json"), "rt") as fh:
                manifest = json.load(fh)
            if manifest.get("strand"):
                raise ParityError(
                    "chunk {} carries strand {!r}. A strand-filtered chunk was "
                    "extracted on the raw flag before the splitter could rewrite "
                    "it; the invariant cannot be checked against it because the "
                    "chunk itself already committed the ordering mistake.".format(
                        cdir, manifest["strand"]
                    )
                )
            chunks.append(
                {
                    "chunk_id": os.path.basename(cdir),
                    "chrom": manifest["chrom"],
                    "lend": manifest["partition_lend"],
                    "rend": manifest["partition_rend"],
                    "offset": manifest["offset"],
                    "bam": os.path.join(cdir, "chunk.bam"),
                    "fasta": os.path.join(cdir, "chunk.fa"),
                    "gtf": os.path.join(cdir, "chunk.gtf"),
                    "dropped": manifest.get("dropped_read_names", []),
                }
            )
    else:
        overlaps = check_disjoint(intervals)
        if overlaps:
            raise ParityError(
                "the intervals are not disjoint, so a read could be compared "
                "twice: {}".format("; ".join(overlaps))
            )
        extract_root = os.path.join(work_dir, "extracted")
        os.makedirs(extract_root, exist_ok=True)
        for idx, interval in enumerate(intervals):
            chunk_id = "{}_{:02d}".format(interval.chrom, idx)
            cdir = os.path.join(extract_root, chunk_id)
            os.makedirs(cdir, exist_ok=True)
            manifest = extractor.extract_partition(
                genome_fa,
                bam,
                "{}:{}-{}".format(interval.chrom, interval.lend, interval.rend),
                os.path.join(cdir, "chunk"),
                gtf=gtf,
                margin=margin,
                secondary_alignments="exclude",
                gtf_index_cache_dir=gtf_index_cache_dir,
            )
            chunks.append(
                {
                    "chunk_id": chunk_id,
                    "chrom": interval.chrom,
                    "lend": interval.lend,
                    "rend": interval.rend,
                    "offset": manifest["offset"],
                    "bam": manifest["files"]["bam"],
                    "fasta": manifest["files"]["fasta"],
                    "gtf": manifest["files"]["gtf"],
                    "dropped": manifest["dropped_read_names"],
                }
            )

    overhang_dropped = {}  # read name -> chunk ids that dropped it
    for chunk in chunks:
        for name in chunk["dropped"]:
            overhang_dropped.setdefault(name, []).append(chunk["chunk_id"])

    agreed = 0
    disagreements = []
    chunk_only = []
    per_chunk = []
    flips_in_chunks = 0
    for chunk in chunks:
        chunk_assignments, counters, collisions = split_assignments(
            chunk["bam"],
            os.path.join(os.path.dirname(chunk["bam"]), "chunk.strand.parity"),
            max_intron_length,
            infer_read_orient=infer_read_orient,
            genome_fa=chunk["fasta"],
            gtf=chunk["gtf"] if (gtf and os.path.exists(chunk["gtf"])) else None,
            offset=chunk["offset"],
        )
        flips_in_chunks += counters["num_records_strand_flipped"]
        chunk_agreed = 0
        chunk_disagreed = 0
        for key, chunk_strand in chunk_assignments.items():
            whole_strand = whole.pop(key, None)
            if whole_strand is None:
                chunk_only.append(
                    {
                        "read_name": key[0],
                        "reference_start": key[1],
                        "chunk_id": chunk["chunk_id"],
                        "chunk_strand": chunk_strand,
                        "reason": "the whole-scope split has no verdict for this "
                        "record: it kept a record the whole-scope split "
                        "discarded, so the two filters do not agree",
                    }
                )
                continue
            if whole_strand == chunk_strand:
                agreed += 1
                chunk_agreed += 1
            else:
                chunk_disagreed += 1
                disagreements.append(
                    {
                        "read_name": key[0],
                        "reference_start": key[1],
                        "chunk_id": chunk["chunk_id"],
                        "whole_scope_strand": whole_strand,
                        "in_chunk_strand": chunk_strand,
                    }
                )
        per_chunk.append(
            {
                "chunk_id": chunk["chunk_id"],
                "region": "{}:{}-{}".format(
                    chunk["chrom"], chunk["lend"], chunk["rend"]
                ),
                "records_split_in_chunk": len(chunk_assignments),
                "agreed": chunk_agreed,
                "disagreed": chunk_disagreed,
                "records_strand_flipped_in_chunk": counters[
                    "num_records_strand_flipped"
                ],
                "collisions": collisions,
            }
        )
        if not keep_bams:
            for suffix in (".+.bam", ".-.bam"):
                path = os.path.join(
                    os.path.dirname(chunk["bam"]), "chunk.strand.parity" + suffix
                )
                if os.path.exists(path):
                    os.unlink(path)

    # whatever the chunks never claimed. Each one gets a reason, by name.
    covered = collections.defaultdict(list)
    for chunk in chunks:
        covered[chunk["chrom"]].append((chunk["lend"], chunk["rend"]))
    whole_only = []
    for (name, start), strand in whole.items():
        if name in overhang_dropped:
            reason = "overhang_dropped_at_" + ",".join(overhang_dropped[name])
        elif contig is not None and not _in_any(covered, contig, start):
            reason = "outside_every_chunk_interval"
        else:
            reason = "not_emitted_by_any_chunk"
        whole_only.append(
            {
                "read_name": name,
                "reference_start": start,
                "whole_scope_strand": strand,
                "reason": reason,
            }
        )
    whole_only.sort(key=lambda r: (r["read_name"], r["reference_start"]))

    if disagreements:
        flagged = raw_flag_strands(bam, [d["read_name"] for d in disagreements])
        for record in disagreements:
            record["aligner_flag_strand"] = flagged.get(
                (record["read_name"], record["reference_start"])
            )

    compared = agreed + len(disagreements)
    return {
        "mode": "inferred_orientation" if infer_read_orient else "aligner_flag",
        "bam": os.path.abspath(bam),
        "genome_fa": os.path.abspath(genome_fa),
        "gtf": os.path.abspath(gtf) if gtf else None,
        "chunks_root": os.path.abspath(chunks_root) if chunks_root else None,
        "max_intron_length": max_intron_length,
        "margin": margin,
        "chunks": len(chunks),
        "whole_scope_records_read": whole_counters["num_records"],
        "whole_scope_records_assigned": whole_counters["num_forward"]
        + whole_counters["num_reverse"],
        "whole_scope_records_strand_flipped": whole_counters[
            "num_records_strand_flipped"
        ],
        "whole_scope_collisions": whole_collisions,
        "in_chunk_records_strand_flipped": flips_in_chunks,
        "records_compared": compared,
        "records_agreed": agreed,
        "records_disagreed": len(disagreements),
        "agreement": (agreed / compared) if compared else None,
        "disagreements": disagreements,
        "records_only_in_chunks": chunk_only,
        "records_only_in_whole_scope": whole_only,
        "overhang_dropped_read_names": {
            name: ids for name, ids in sorted(overhang_dropped.items())
        },
        "per_chunk": per_chunk,
        "verdict": (
            "INVARIANT HOLDS"
            if compared and not disagreements and not chunk_only
            else "INVARIANT VIOLATED"
        ),
    }


def _in_any(covered, chrom, start0):
    for lend, rend in covered.get(chrom, ()):
        if lend - 1 <= start0 <= rend - 1:
            return True
    return False


def measure_ordering_violation(
    bam,
    genome_fa,
    stranded_region,
    work_dir,
    gtf=None,
    margin=None,
    max_intron_length=None,
):
    """What extracting a STRAND-FILTERED region from a raw BAM costs, by name.

    This is the mistake the design must not make, executed on purpose so its cost
    is a number rather than a warning.  ``retained_for_extraction`` filters on the
    RAW aligner flag, and ``split_bam_by_strand`` rewrites that flag whenever the
    strand it infers disagrees.  So extracting ``chr1+:...`` from a raw BAM keeps
    every read the aligner called forward and then treats all of them as forward
    -- including the ones whose splice dinucleotides say otherwise, which no later
    stage revisits.  The output looks fine.

    Returns a report naming every emitted read whose inferred strand contradicts
    the region's strand: those are the silently misassigned reads.
    """

    if max_intron_length is None:
        max_intron_length = LRAA_Globals.config["max_intron_length"]
    if margin is None:
        margin = extractor.DEFAULT_MARGIN
    region = extractor.parse_region(stranded_region)
    if not region.strand:
        raise ParityError(
            "{} names no strand, so there is no ordering mistake to "
            "measure".format(stranded_region)
        )
    os.makedirs(work_dir, exist_ok=True)
    manifest = extractor.extract_partition(
        genome_fa,
        bam,
        region,
        os.path.join(work_dir, "wrong_order"),
        gtf=gtf,
        margin=margin,
        secondary_alignments="exclude",
        # The extractor REFUSES this combination by default, and it is right to:
        # a strand-suffixed region against a bam still holding both orientations
        # is a wrong answer, not a filtered one. Overridden here because
        # committing the mistake is the whole point of this function.
        mixed_orientation_source="warn",
    )
    assignments, counters, _ = split_assignments(
        manifest["files"]["bam"],
        os.path.join(work_dir, "wrong_order.strand"),
        max_intron_length,
        infer_read_orient=True,
        genome_fa=manifest["files"]["fasta"],
        gtf=manifest["files"]["gtf"],
        offset=manifest["offset"],
    )
    misassigned = sorted(
        {name for (name, _), strand in assignments.items() if strand != region.strand}
    )
    return {
        "region": stranded_region,
        "records_emitted_by_the_strand_filter": manifest["counts"][
            "alignments_emitted"
        ],
        "records_the_split_would_move": counters["num_records_strand_flipped"],
        "misassigned_read_names": misassigned,
        "verdict": (
            "ORDERING MISTAKE IS DETECTABLE: {} of {} emitted records belong on "
            "the other strand".format(
                len(misassigned), manifest["counts"]["alignments_emitted"]
            )
            if misassigned
            else "no read in this region contradicts the aligner flag, so this "
            "region cannot demonstrate the mistake"
        ),
    }


# ------------------------------------------------------------------- arm loading


def load_arm(arm_dir):
    """The comparable artifacts of one pipeline output directory."""

    arm_dir = os.path.abspath(arm_dir)
    outputs_path = os.path.join(arm_dir, "outputs.json")
    if not os.path.exists(outputs_path):
        raise ParityError("{} has no outputs.json".format(arm_dir))
    with open(outputs_path, "rt") as fh:
        outputs = json.load(fh)
    if "chunked" not in outputs:
        raise ParityError(
            "{} holds no chunked arm; this compares two CHUNKED runs".format(arm_dir)
        )
    timing = {}
    timing_path = os.path.join(arm_dir, "timing.json")
    if os.path.exists(timing_path):
        with open(timing_path, "rt") as fh:
            timing = json.load(fh)
    return {
        "dir": arm_dir,
        "expr": outputs["chunked"]["quant_expr"],
        "tracking": outputs["chunked"]["quant_tracking"],
        "expr_rows": outputs["chunked"]["expr_rows"],
        "tracking_rows": outputs["chunked"]["tracking_rows"],
        "num_total_reads": outputs.get("num_total_reads"),
        "strandless": bool(timing.get("strandless_chunks")),
        "baseline_excluded_severed_reads": timing.get(
            "baseline_excluded_severed_reads"
        ),
        "timing": timing,
    }


def severed_sets(arm_dir):
    """What each arm lost, split by where it was lost.

    ``cut_selection`` is what the selector predicted it would cost; ``overhang``
    is what extraction actually dropped.  They are reported separately because
    they are produced by different code and a disagreement between them is itself
    a defect worth seeing.
    """

    predicted = ChunkedRun.severed_read_names(os.path.join(arm_dir, "cuts"))
    dropped = {}
    for path in sorted(
        glob.glob(os.path.join(arm_dir, "chunks", "*", "chunk.dropped_reads.txt"))
    ):
        chunk_id = os.path.basename(os.path.dirname(path))
        with open(path, "rt") as fh:
            for line in fh:
                name = line.strip()
                if name:
                    dropped.setdefault(name, []).append(chunk_id)
    return {"cut_selection": predicted, "overhang": dropped}


def chunk_input_names(arm_dir):
    """read name -> the chunk ids whose ``chunk.bam`` holds it.

    One streaming pass per chunk BAM.  Names, not positions: chunk BAMs carry
    rebased coordinates and the two arms cut in different places, so a position
    is not comparable between them.
    """

    where = {}
    for path in sorted(glob.glob(os.path.join(arm_dir, "chunks", "*", "chunk.bam"))):
        chunk_id = os.path.basename(os.path.dirname(path))
        with pysam.AlignmentFile(path, "rb") as fh:
            for aln in fh.fetch(until_eof=True):
                where.setdefault(aln.query_name, []).append(chunk_id)
    return where


def chunk_strand_placement(arm_dir):
    """read name -> ``chunk_id:strand`` for a strandless arm's post-split BAMs.

    Empty for a strand-first arm, whose chunks are strand-pure already and whose
    placement is therefore its chunk id.
    """

    placement = {}
    # The two names, exactly, rather than a wildcard: this module writes its own
    # chunk.strand.parity.{+,-}.bam beside them when asked to keep its
    # intermediates, and a wildcard would report "parity" as an orientation.
    for cdir in sorted(glob.glob(os.path.join(arm_dir, "chunks", "*"))):
        chunk_id = os.path.basename(cdir)
        for strand in ("+", "-"):
            path = os.path.join(cdir, "chunk.strand.{}.bam".format(strand))
            if not os.path.exists(path):
                continue
            with pysam.AlignmentFile(path, "rb", check_sq=False) as fh:
                for aln in fh.fetch(until_eof=True):
                    placement.setdefault(aln.query_name, []).append(
                        "{}:{}".format(chunk_id, strand)
                    )
    return placement


def explain_absence(name, evidence):
    """Why ``name`` is not in this arm's output. One reason, from its artifacts."""

    if name in evidence["overhang"]:
        return "overhang_dropped_at_" + ",".join(evidence["overhang"][name])
    if name in evidence["severed_predicted"]:
        return "severed_at_cut_selection_but_not_dropped_by_extraction"
    if name in evidence["chunk_inputs"]:
        return "in_chunk_input_{}_but_assigned_no_transcript".format(
            ",".join(evidence["chunk_inputs"][name])
        )
    return "absent_from_every_chunk_input"


# --------------------------------------------------------------- table handling


def load_expr(path):
    comments, header, rows = ChunkedRun.read_tsv(path)
    return header, rows


def load_tracking(path):
    comments, header, rows = ChunkedRun.read_tsv(path)
    return header, rows


def _column_index(header, name):
    try:
        return header.index(name)
    except ValueError:
        raise ParityError("no {} column in {}".format(name, header))


def _row_projector(header):
    """A row -> comparable-tuple function, with the column choice made once.

    Recomputing it per row costs a header scan per column per row, which on a
    million-row tracking table is the whole runtime.
    """

    keep = tuple(
        i
        for i, name in enumerate(header)
        if name not in INCOMPARABLE_TRACKING_COLUMNS
    )
    return lambda row: tuple(row[i] for i in keep)


def file_bytes_equal(a, b):
    """Byte equality of two files, decompressing ``.gz`` first.

    The gzip member carries an mtime, so two identical tables compress to
    different bytes; comparing the compressed form would report a difference that
    is not in the data.
    """

    opener = ChunkedRun.open_text
    with opener(a) as fa, opener(b) as fb:
        while True:
            ca = fa.read(1 << 20)
            cb = fb.read(1 << 20)
            if ca != cb:
                return False
            if not ca:
                return True


NUMERIC_EXPR_COLUMNS = (
    "uniq_reads",
    "all_reads",
    "isoform_fraction",
    "unique_gene_read_fraction",
    "TPM",
    "RPM_total_reads",
)
STRUCTURAL_EXPR_COLUMNS = ("gene_id", "exons", "introns", "splice_hash_code")


def _duplicate_ids(rows, column):
    """Ids appearing more than once, with their counts. Sorted, so it is diffable."""

    counts = collections.Counter(row[column] for row in rows)
    return {name: n for name, n in sorted(counts.items()) if n > 1}


def compare_expr(path_a, path_b):
    """Per-transcript comparison of two ``quant.expr`` tables."""

    header_a, rows_a = load_expr(path_a)
    header_b, rows_b = load_expr(path_b)
    if header_a != header_b:
        return {
            "byte_identical": False,
            "header_identical": False,
            "header_a": header_a,
            "header_b": header_b,
        }
    tid = _column_index(header_a, "transcript_id")

    # A transcript may appear once per table, and a table where it appears twice is
    # not a table this comparison can key on -- but it is also a FINDING, not a
    # reason to abort: a strandless chunk whose GTF carries both orientations,
    # handed unchanged to both of its per-orientation quant runs, emits every
    # transcript twice, each copy quantified against the other strand's reads. So
    # the duplicates are named and the per-transcript comparison below is skipped
    # rather than computed on an arbitrary one of each pair.
    duplicates_a = _duplicate_ids(rows_a, tid)
    duplicates_b = _duplicate_ids(rows_b, tid)
    if duplicates_a or duplicates_b:
        return {
            "byte_identical": file_bytes_equal(path_a, path_b),
            "header_identical": True,
            "rows_a": len(rows_a),
            "rows_b": len(rows_b),
            "duplicate_transcript_ids_a": duplicates_a,
            "duplicate_transcript_ids_b": duplicates_b,
            "content_identical": False,
            "not_comparable": (
                "transcript_id is not unique: {} duplicated in the strand-first "
                "table, {} in the strandless table. Each duplicate row was "
                "quantified in a different scope, so neither is a copy of the "
                "other and keying on transcript_id would silently pick one".format(
                    len(duplicates_a), len(duplicates_b)
                )
            ),
            "numeric": {},
            "structural": {},
            "transcripts_only_in_a": [],
            "transcripts_only_in_b": [],
            "row_order_differs": None,
        }

    index_a = {row[tid]: row for row in rows_a}
    index_b = {row[tid]: row for row in rows_b}
    only_a = sorted(set(index_a) - set(index_b))
    only_b = sorted(set(index_b) - set(index_a))
    shared = sorted(set(index_a) & set(index_b))

    numeric = {}
    for name in NUMERIC_EXPR_COLUMNS:
        if name not in header_a:
            continue
        col = header_a.index(name)
        worst = 0.0
        differing = []
        for transcript in shared:
            va = float(index_a[transcript][col] or 0)
            vb = float(index_b[transcript][col] or 0)
            if va != vb:
                differing.append(transcript)
                worst = max(worst, abs(va - vb))
        numeric[name] = {
            "transcripts_differing": len(differing),
            "max_abs_delta": worst,
            "differing_transcripts": differing,
        }
    structural = {}
    for name in STRUCTURAL_EXPR_COLUMNS:
        if name not in header_a:
            continue
        col = header_a.index(name)
        differing = [
            t for t in shared if index_a[t][col] != index_b[t][col]
        ]
        structural[name] = {
            "transcripts_differing": len(differing),
            "differing_transcripts": differing,
        }

    row_order_only = [row[tid] for row in rows_a] != [row[tid] for row in rows_b]
    any_value_differs = any(
        v["transcripts_differing"] for v in numeric.values()
    ) or any(v["transcripts_differing"] for v in structural.values())

    return {
        "byte_identical": file_bytes_equal(path_a, path_b),
        "header_identical": True,
        "rows_a": len(rows_a),
        "rows_b": len(rows_b),
        "transcripts_only_in_a": only_a,
        "transcripts_only_in_b": only_b,
        "row_order_differs": row_order_only,
        "numeric": numeric,
        "structural": structural,
        "content_identical": not (only_a or only_b or any_value_differs),
    }


def compare_tracking(path_a, path_b):
    """``quant.tracking`` compared as rows, then as a multiset, then by key.

    Order is checked SEPARATELY from content, and an ordering-only difference is
    proven rather than asserted: the two row multisets must be equal, which no
    content difference can satisfy.
    """

    header_a, rows_a = load_tracking(path_a)
    header_b, rows_b = load_tracking(path_b)
    if header_a != header_b:
        return {
            "byte_identical": False,
            "header_identical": False,
            "header_a": header_a,
            "header_b": header_b,
        }

    dropped = [c for c in INCOMPARABLE_TRACKING_COLUMNS if c in header_a]
    project = _row_projector(header_a)
    projected_a = [project(row) for row in rows_a]
    projected_b = [project(row) for row in rows_b]
    multiset_equal = collections.Counter(projected_a) == collections.Counter(
        projected_b
    )
    ordered_equal = projected_a == projected_b

    read_col = _column_index(header_a, "read_name")
    tx_col = _column_index(header_a, "transcript_id")
    gene_col = _column_index(header_a, "gene_id")
    frac_col = header_a.index("frac_assigned") if "frac_assigned" in header_a else None

    def keyed(rows):
        out = {}
        for row in rows:
            out[(row[read_col], row[tx_col], row[gene_col])] = row
        return out

    keys_a = keyed(rows_a)
    keys_b = keyed(rows_b)
    only_a = sorted(set(keys_a) - set(keys_b))
    only_b = sorted(set(keys_b) - set(keys_a))
    frac_differing = []
    frac_worst = 0.0
    if frac_col is not None:
        for key in set(keys_a) & set(keys_b):
            va = float(keys_a[key][frac_col] or 0)
            vb = float(keys_b[key][frac_col] or 0)
            if va != vb:
                frac_differing.append(key)
                frac_worst = max(frac_worst, abs(va - vb))

    return {
        "byte_identical": file_bytes_equal(path_a, path_b),
        "header_identical": True,
        "columns_excluded_as_incomparable": dropped,
        "rows_a": len(rows_a),
        "rows_b": len(rows_b),
        "rows_identical_in_order": ordered_equal,
        "row_multisets_equal": multiset_equal,
        "difference_is_ordering_only": multiset_equal and not ordered_equal,
        "keys_only_in_a": [list(k) for k in only_a],
        "keys_only_in_b": [list(k) for k in only_b],
        "frac_assigned_keys_differing": len(frac_differing),
        "frac_assigned_max_abs_delta": frac_worst,
        "frac_assigned_differing_keys": [list(k) for k in frac_differing],
        "reads_a": len({row[read_col] for row in rows_a}),
        "reads_b": len({row[read_col] for row in rows_b}),
        "read_names_a": {row[read_col] for row in rows_a},
        "read_names_b": {row[read_col] for row in rows_b},
    }


# ------------------------------------------------------------------ the verdict


def attribute_divergence(expr, tracking, severed_symmetric_difference):
    """Is every difference attributable to the reads only one arm ever saw?

    A transcript's numbers move when the reads assigned to it change, so a
    transcript that differs must have at least one read from the symmetric
    difference of the severed sets in its tracking rows.  One that differs
    WITHOUT one is a real defect: the arms saw the same reads for it and still
    disagreed.

    TPM is excluded and reported separately.  Its denominator is the total
    ``all_reads`` over the merged table (``ChunkedRun.merge_and_translate``, which
    reproduces ``LRAA:_merge_quant_expr_files``), so a single differing read
    rescales EVERY row.  Including it would mark the whole table unexplained on
    the strength of one severed read, which says nothing.
    """

    severed = set(severed_symmetric_difference)
    touched_transcripts = set()
    for key in tracking.get("keys_only_in_a", []) + tracking.get("keys_only_in_b", []):
        read_name, transcript_id = key[0], key[1]
        if read_name in severed:
            touched_transcripts.add(transcript_id)
    for key in tracking.get("frac_assigned_differing_keys", []):
        if key[0] in severed:
            touched_transcripts.add(key[1])

    unexplained = {}
    for name, stats in expr.get("numeric", {}).items():
        if name == "TPM":
            continue
        residue = [
            t for t in stats["differing_transcripts"] if t not in touched_transcripts
        ]
        if residue:
            unexplained[name] = residue
    for name, stats in expr.get("structural", {}).items():
        residue = [
            t for t in stats["differing_transcripts"] if t not in touched_transcripts
        ]
        if residue:
            unexplained[name] = residue

    unexplained_keys = [
        key
        for key in tracking.get("keys_only_in_a", [])
        + tracking.get("keys_only_in_b", [])
        if key[0] not in severed
    ]
    unexplained_fracs = [
        key
        for key in tracking.get("frac_assigned_differing_keys", [])
        if key[0] not in severed
    ]

    return {
        "severed_symmetric_difference_size": len(severed),
        "transcripts_touched_by_severed_reads": sorted(touched_transcripts),
        "unexplained_expr_columns": unexplained,
        "unexplained_tracking_keys": unexplained_keys,
        "unexplained_frac_assigned_keys": unexplained_fracs,
        "tpm_transcripts_differing": expr.get("numeric", {})
        .get("TPM", {})
        .get("transcripts_differing"),
        "tpm_max_abs_delta": expr.get("numeric", {}).get("TPM", {}).get("max_abs_delta"),
        "fully_explained_by_severing": not (
            unexplained or unexplained_keys or unexplained_fracs
        ),
    }


def compare_arms(strand_first_dir, strandless_dir, explain_reads=True):
    """The gate. Returns a report; ``verdict`` is the one-line answer."""

    arm_a = load_arm(strand_first_dir)
    arm_b = load_arm(strandless_dir)

    preconditions = []
    if arm_a["num_total_reads"] != arm_b["num_total_reads"]:
        preconditions.append(
            "num_total_reads differs: strand-first {} vs strandless {}. It is the "
            "-N passed identically to every chunk and it is RPM_total_reads' "
            "denominator, so the tables cannot be byte-identical while it "
            "differs, and no per-transcript comparison of RPM means anything. "
            "The denominator has to be the record set the arms both see.".format(
                arm_a["num_total_reads"], arm_b["num_total_reads"]
            )
        )

    expr = compare_expr(arm_a["expr"], arm_b["expr"])
    tracking = compare_tracking(arm_a["tracking"], arm_b["tracking"])

    severed_a = severed_sets(arm_a["dir"])
    severed_b = severed_sets(arm_b["dir"])
    lost_a = set(severed_a["overhang"]) | severed_a["cut_selection"]
    lost_b = set(severed_b["overhang"]) | severed_b["cut_selection"]
    symmetric = sorted(lost_a ^ lost_b)

    reads_a = tracking.pop("read_names_a", set())
    reads_b = tracking.pop("read_names_b", set())
    missing_from_b = sorted(reads_a - reads_b)
    missing_from_a = sorted(reads_b - reads_a)

    # The arms' INPUTS, compared by read name. This is the level-1 check of
    # docs/chunked_quant_parity_evaluation.md, in its stronger form: chunk BAMs
    # carry rebased coordinates and the two arms cut in different places, so a
    # position is not comparable between them but a name is. It costs one
    # streaming pass per chunk BAM, which is why it is behind ``explain_reads``.
    inputs_a = chunk_input_names(arm_a["dir"]) if explain_reads else {}
    inputs_b = chunk_input_names(arm_b["dir"]) if explain_reads else {}

    named_absences = {"absent_from_strandless": [], "absent_from_strand_first": []}
    if explain_reads and (missing_from_b or missing_from_a):
        evidence_b = {
            "overhang": severed_b["overhang"],
            "severed_predicted": severed_b["cut_selection"],
            "chunk_inputs": inputs_b,
        }
        evidence_a = {
            "overhang": severed_a["overhang"],
            "severed_predicted": severed_a["cut_selection"],
            "chunk_inputs": inputs_a,
        }
        placement_b = chunk_strand_placement(arm_b["dir"]) if missing_from_b else {}
        for name in missing_from_b:
            named_absences["absent_from_strandless"].append(
                {
                    "read_name": name,
                    "reason": explain_absence(name, evidence_b),
                    "strandless_strand_placement": placement_b.get(name),
                }
            )
        for name in missing_from_a:
            named_absences["absent_from_strand_first"].append(
                {"read_name": name, "reason": explain_absence(name, evidence_a)}
            )

    chunk_inputs = None
    if explain_reads:
        input_only_a = sorted(set(inputs_a) - set(inputs_b))
        input_only_b = sorted(set(inputs_b) - set(inputs_a))
        chunk_inputs = {
            "distinct_read_names_strand_first": len(inputs_a),
            "distinct_read_names_strandless": len(inputs_b),
            "names_identical": not input_only_a and not input_only_b,
            "only_in_strand_first": [
                {
                    "read_name": name,
                    "strand_first_chunks": inputs_a[name],
                    "reason": explain_absence(
                        name,
                        {
                            "overhang": severed_b["overhang"],
                            "severed_predicted": severed_b["cut_selection"],
                            "chunk_inputs": inputs_b,
                        },
                    ),
                }
                for name in input_only_a
            ],
            "only_in_strandless": [
                {
                    "read_name": name,
                    "strandless_chunks": inputs_b[name],
                    "reason": explain_absence(
                        name,
                        {
                            "overhang": severed_a["overhang"],
                            "severed_predicted": severed_a["cut_selection"],
                            "chunk_inputs": inputs_a,
                        },
                    ),
                }
                for name in input_only_b
            ],
        }

    # Every read only one arm severed, resolved individually. A count here would
    # be the least informative thing available: the two designs sever by
    # construction different sets, and what matters is whether the difference
    # reached quantification at all. A severed read that no arm assigned to any
    # transcript costs the comparison nothing; one that appears in a tracking
    # table is the reason a number moved.
    symmetric_detail = []
    for name in symmetric:
        symmetric_detail.append(
            {
                "read_name": name,
                "lost_by": "strand_first" if name in lost_a else "strandless",
                "in_strand_first_tracking": name in reads_a,
                "in_strandless_tracking": name in reads_b,
                "strand_first_overhang_chunks": severed_a["overhang"].get(name),
                "strandless_overhang_chunks": severed_b["overhang"].get(name),
                "reached_quantification_in_either_arm": (
                    name in reads_a or name in reads_b
                ),
            }
        )

    attribution = attribute_divergence(expr, tracking, symmetric)

    identical = expr["byte_identical"] and tracking["byte_identical"]
    # Non-uniqueness in either transcript table has to be its own verdict, because
    # the attribution below would otherwise pass it: with no per-transcript
    # comparison to draw from, nothing is "unexplained" and the arms would be
    # declared explained-by-severing while one table holds two rows per transcript.
    # An instrument that reports a clean result on an uninterpretable table is
    # worse than one that reports nothing.
    not_comparable = expr.get("not_comparable")
    if not_comparable and not identical:
        verdict = "NOT COMPARABLE: {}".format(not_comparable)
    elif identical:
        verdict = "EQUIVALENT: expr and tracking are byte-identical"
    elif (
        expr.get("content_identical")
        and tracking.get("rows_identical_in_order")
        and not missing_from_a
        and not missing_from_b
    ):
        # NOT an ordering difference: the rows are in the same order and agree
        # field for field once the incomparable columns are dropped. Calling this
        # "ordering" would misdescribe it, and the distinction is the whole reason
        # order and content are tested separately.
        verdict = (
            "EQUIVALENT UP TO {}: every expr value agrees and every tracking row "
            "agrees in the same order, differing only in {}, which is a "
            "per-process counter rather than data".format(
                ", ".join(tracking.get("columns_excluded_as_incomparable") or ["-"]),
                ", ".join(tracking.get("columns_excluded_as_incomparable") or ["-"]),
            )
        )
    elif (
        expr.get("content_identical")
        and tracking.get("row_multisets_equal")
        and not missing_from_a
        and not missing_from_b
    ):
        verdict = (
            "EQUIVALENT UP TO ROW ORDER: every expr value agrees, and the tracking "
            "row multisets are EQUAL while the row sequences are not -- which no "
            "content difference can satisfy, so the difference is order alone"
        )
    elif attribution["fully_explained_by_severing"] and symmetric:
        verdict = (
            "DIVERGENT, FULLY EXPLAINED BY SEVERING: every differing row involves "
            "a read in the symmetric difference of the two severed sets"
        )
    else:
        verdict = (
            "DIVERGENT AND NOT EXPLAINED BY SEVERING: at least one row differs "
            "with no severed-set read behind it"
        )
    if expr.get("duplicate_transcript_ids_a") or expr.get(
        "duplicate_transcript_ids_b"
    ):
        preconditions.append(
            "transcript_id is not unique in an expr table: {} duplicated in "
            "strand-first, {} in strandless. A transcript quantified twice in one "
            "run is a defect in that run whatever the arms agree on".format(
                len(expr.get("duplicate_transcript_ids_a") or {}),
                len(expr.get("duplicate_transcript_ids_b") or {}),
            )
        )

    return {
        "strand_first": {
            k: arm_a[k]
            for k in (
                "dir",
                "expr",
                "tracking",
                "expr_rows",
                "tracking_rows",
                "num_total_reads",
                "baseline_excluded_severed_reads",
            )
        },
        "strandless": {
            k: arm_b[k]
            for k in (
                "dir",
                "expr",
                "tracking",
                "expr_rows",
                "tracking_rows",
                "num_total_reads",
                "baseline_excluded_severed_reads",
            )
        },
        "preconditions_violated": preconditions,
        "expr": expr,
        "tracking": tracking,
        "severed": {
            "strand_first_cut_selection": sorted(severed_a["cut_selection"]),
            "strand_first_overhang": {
                k: v for k, v in sorted(severed_a["overhang"].items())
            },
            "strandless_cut_selection": sorted(severed_b["cut_selection"]),
            "strandless_overhang": {
                k: v for k, v in sorted(severed_b["overhang"].items())
            },
            "only_strand_first_lost": sorted(lost_a - lost_b),
            "only_strandless_lost": sorted(lost_b - lost_a),
            "symmetric_difference": symmetric,
            "symmetric_difference_detail": symmetric_detail,
        },
        "chunk_inputs": chunk_inputs,
        "reads_absent_by_name": named_absences,
        "attribution": attribution,
        "verdict": verdict,
        "passed": identical
        or (
            not not_comparable
            and bool(attribution["fully_explained_by_severing"])
        ),
    }


# ------------------------------------------------------------- end-to-end gate

PIPELINE = os.path.join(REPO_ROOT, "util", "misc", "run_chunked_quant_pipeline.py")


def _run_logged(cmd, log_path):
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "wt") as fh:
        print("$ " + " ".join(cmd), file=fh)
        fh.flush()
        return subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)


def run_gate(
    bam,
    genome_fa,
    gtf,
    out_dir,
    contig=None,
    HiFi=False,
    approx_MB_per_cut=None,
    approx_MB_per_cut_wiggle_window=None,
    cpu_budget=4,
    margin=None,
    max_intron_length=None,
    expect=None,
    skip_infer_invariant=False,
    skip_arms=False,
    ordering_cost_region=None,
    strand_first_arm="both",
):
    """The whole gate, cheapest and most falsifiable step first.

    Order is deliberate.  The strand-assignment invariant costs one extraction and
    two splits and can kill the design on its own, so it runs before either arm is
    quantified.  ``expect`` pins the strand-first arm's own figures, which is what
    makes a strandless number interpretable: without a control that reproduces
    known values, "strandless gives 553" says nothing about which arm moved.

    ``expect`` keys, all optional: ``expr_rows``, ``tracking_rows``,
    ``num_total_reads``, ``baseline_excluded_severed_reads``.

    ``strand_first_arm`` is ``"both"`` by default so the strand-first arm also
    produces the whole-contig control, which is what makes
    ``baseline_excluded_severed_reads`` available to pin.  Set it to ``"chunked"``
    where the control costs more than it tells you -- with nothing severed there is
    nothing for the control to subtract, and arm-to-arm equivalence never consults
    it.
    """

    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    report = {"out_dir": out_dir, "steps": [], "failures": []}

    common = []
    if contig:
        common += ["--contig", contig]
    if HiFi:
        common += ["--HiFi"]
    if approx_MB_per_cut is not None:
        common += ["--approx_MB_per_cut", str(approx_MB_per_cut)]
    if approx_MB_per_cut_wiggle_window is not None:
        common += [
            "--approx_MB_per_cut_wiggle_window",
            str(approx_MB_per_cut_wiggle_window),
        ]

    # ---- step 1: union cut selection, which is also the strandless geometry.
    cuts_prefix = os.path.join(out_dir, "union_cuts", "union")
    os.makedirs(os.path.dirname(cuts_prefix), exist_ok=True)
    selector = [
        sys.executable,
        SELECT_CUTS,
        "--bam",
        os.path.abspath(bam),
        "--genome_fa",
        os.path.abspath(genome_fa),
        "--gtf",
        os.path.abspath(gtf),
        "--gtf_index_cache_dir",
        os.path.join(out_dir, "gtf_index"),
        "--grid_origin",
        "0",
        "--min_per_id",
        "97" if HiFi else str(LRAA_Globals.config["min_per_id"]),
        "--min_mapping_quality",
        str(int(LRAA_Globals.config["min_mapping_quality_for_final_quant"])),
        "--output_prefix",
        cuts_prefix,
    ]
    # --strandless rather than merely omitting --strand: omitting it already
    # unions the two orientations' islands, but the declaration is what gets the
    # per-orientation denominator, the suffixless region strings, and a warning if
    # the bam turns out to hold one orientation -- which would mean an
    # already-split bam was passed and these "strandless" cuts serve one strand.
    selector += ["--strandless"]
    if contig:
        selector += ["--contig", contig]
    if approx_MB_per_cut is not None:
        selector += ["--approx_MB_per_cut", str(approx_MB_per_cut)]
    if approx_MB_per_cut_wiggle_window is not None:
        selector += [
            "--approx_MB_per_cut_wiggle_window",
            str(approx_MB_per_cut_wiggle_window),
        ]
    rc = _run_logged(selector, os.path.join(out_dir, "logs", "union_cuts.log"))
    if rc != 0:
        report["failures"].append("union cut selection failed; see logs/union_cuts.log")
        return report
    intervals = intervals_from_cuts_json(cuts_prefix + ".cuts.json")
    report["steps"].append({"step": "union_cuts", "intervals": len(intervals)})

    # ---- step 2: the invariant, on the aligner flag -- the mode that runs today.
    flag_report = check_strand_invariant(
        bam,
        genome_fa,
        os.path.join(out_dir, "invariant_flag"),
        intervals=intervals,
        gtf=gtf,
        contig=contig,
        margin=margin,
        max_intron_length=max_intron_length,
        gtf_index_cache_dir=os.path.join(out_dir, "gtf_index"),
    )
    report["invariant_aligner_flag"] = flag_report
    if flag_report["verdict"] != "INVARIANT HOLDS":
        report["failures"].append(
            "strand-assignment invariant violated on the aligner flag"
        )

    # ---- step 3: the invariant under inference, where a flip can happen at all.
    if not skip_infer_invariant:
        infer_report = check_strand_invariant(
            bam,
            genome_fa,
            os.path.join(out_dir, "invariant_infer"),
            chunks_root=os.path.join(out_dir, "invariant_flag", "extracted"),
            gtf=gtf,
            contig=contig,
            margin=margin,
            max_intron_length=max_intron_length,
            infer_read_orient=True,
        )
        report["invariant_inferred"] = infer_report
        if infer_report["verdict"] != "INVARIANT HOLDS":
            report["failures"].append(
                "strand-assignment invariant violated under inferred orientation"
            )

    # ---- step 4: what the wrong ordering would have cost, on this substrate.
    if ordering_cost_region:
        report["ordering_cost"] = measure_ordering_violation(
            bam,
            genome_fa,
            ordering_cost_region,
            os.path.join(out_dir, "ordering_cost"),
            gtf=gtf,
            margin=margin,
            max_intron_length=max_intron_length,
        )

    if skip_arms:
        return report

    # ---- step 5: the strand-first arm, with the control, so its severed
    # subtraction is exercised and its own figures are pinned.
    strand_first_dir = os.path.join(out_dir, "strand_first")
    rc = _run_logged(
        [sys.executable, PIPELINE, "--bam", os.path.abspath(bam),
         "--genome_fa", os.path.abspath(genome_fa), "--gtf", os.path.abspath(gtf),
         "--output_dir", strand_first_dir, "--arm", strand_first_arm,
         "--cpu_budget", str(cpu_budget)] + common,
        os.path.join(out_dir, "logs", "strand_first.log"),
    )
    if rc != 0:
        report["failures"].append("the strand-first arm failed; see logs/strand_first.log")
        return report
    arm_a = load_arm(strand_first_dir)
    report["strand_first_figures"] = {
        "expr_rows": arm_a["expr_rows"],
        "tracking_rows": arm_a["tracking_rows"],
        "num_total_reads": arm_a["num_total_reads"],
        "baseline_excluded_severed_reads": arm_a["baseline_excluded_severed_reads"],
    }
    for key, wanted in sorted((expect or {}).items()):
        got = report["strand_first_figures"].get(key)
        if got != wanted:
            report["failures"].append(
                "strand-first {} is {}, expected {}. The control has moved, so no "
                "strandless figure compared against it is interpretable".format(
                    key, got, wanted
                )
            )

    # ---- step 6: the strandless arm, on the SAME denominator.
    #
    # -N is taken from the strand-first arm rather than recomputed. Strandless
    # --arm chunked runs no stage 1 and so has no retained count to default to,
    # and RPM_total_reads is that number: letting each arm pick its own would
    # guarantee a difference in every row and prove nothing about chunking.
    strandless_dir = os.path.join(out_dir, "strandless")
    rc = _run_logged(
        [sys.executable, PIPELINE, "--bam", os.path.abspath(bam),
         "--genome_fa", os.path.abspath(genome_fa), "--gtf", os.path.abspath(gtf),
         "--output_dir", strandless_dir, "--arm", "chunked",
         "--strandless_chunks", "-N", str(arm_a["num_total_reads"]),
         "--cpu_budget", str(cpu_budget)] + common,
        os.path.join(out_dir, "logs", "strandless.log"),
    )
    if rc != 0:
        report["failures"].append("the strandless arm failed; see logs/strandless.log")
        return report
    arm_b = load_arm(strandless_dir)
    report["strandless_figures"] = {
        "expr_rows": arm_b["expr_rows"],
        "tracking_rows": arm_b["tracking_rows"],
        "num_total_reads": arm_b["num_total_reads"],
    }

    # ---- step 7: the invariant against the artifact the pipeline really wrote,
    # not one this module extracted for itself.
    real_chunks = os.path.join(strandless_dir, "chunks")
    if chunk_dirs(real_chunks):
        artifact_report = check_strand_invariant(
            bam,
            genome_fa,
            os.path.join(out_dir, "invariant_artifact"),
            chunks_root=real_chunks,
            gtf=gtf,
            contig=contig,
            margin=margin,
            max_intron_length=max_intron_length,
        )
        report["invariant_pipeline_artifact"] = artifact_report
        if artifact_report["verdict"] != "INVARIANT HOLDS":
            report["failures"].append(
                "the invariant fails against the chunks the pipeline itself wrote"
            )

    # ---- step 8: the equivalence comparison.
    arms = compare_arms(strand_first_dir, strandless_dir)
    report["arms"] = arms
    if arms["preconditions_violated"]:
        report["failures"].extend(arms["preconditions_violated"])
    if not arms["passed"]:
        report["failures"].append(arms["verdict"])
    report["verdict"] = arms["verdict"]
    report["passed"] = not report["failures"]
    return report


# -------------------------------------------------------------------------- CLI


def _print_invariant(report):
    print("STRAND-ASSIGNMENT INVARIANT [{}]".format(report["mode"]))
    print("  bam                        {}".format(report["bam"]))
    print("  chunks                     {}".format(report["chunks"]))
    print(
        "  whole-scope records        {} read, {} strand-assigned, {} flipped".format(
            report["whole_scope_records_read"],
            report["whole_scope_records_assigned"],
            report["whole_scope_records_strand_flipped"],
        )
    )
    print(
        "  in-chunk records flipped   {}".format(report["in_chunk_records_strand_flipped"])
    )
    print(
        "  compared                   {} records; agreed {}; disagreed {}".format(
            report["records_compared"],
            report["records_agreed"],
            report["records_disagreed"],
        )
    )
    if report["agreement"] is not None:
        print("  agreement                  {:.10f}".format(report["agreement"]))
    print(
        "  only in chunks             {}".format(len(report["records_only_in_chunks"]))
    )
    print(
        "  only in whole scope        {}".format(
            len(report["records_only_in_whole_scope"])
        )
    )
    reasons = collections.Counter(
        r["reason"].split("_at_")[0] for r in report["records_only_in_whole_scope"]
    )
    for reason, count in sorted(reasons.items()):
        print("      {:<42} {}".format(reason, count))
    for record in report["disagreements"][:50]:
        print(
            "  DISAGREE {} @{} chunk {}: whole {} vs chunk {} (aligner flag {})".format(
                record["read_name"],
                record["reference_start"],
                record["chunk_id"],
                record["whole_scope_strand"],
                record["in_chunk_strand"],
                record.get("aligner_flag_strand"),
            )
        )
    for record in report["records_only_in_chunks"][:50]:
        print(
            "  CHUNK-ONLY {} @{} chunk {} strand {}".format(
                record["read_name"],
                record["reference_start"],
                record["chunk_id"],
                record["chunk_strand"],
            )
        )
    print("  {}".format(report["verdict"]))


def _print_arms(report):
    print("ARM EQUIVALENCE")
    for label in ("strand_first", "strandless"):
        arm = report[label]
        print(
            "  {:<13} {} expr_rows {} tracking_rows {} N {}".format(
                label,
                arm["dir"],
                arm["expr_rows"],
                arm["tracking_rows"],
                arm["num_total_reads"],
            )
        )
    for note in report["preconditions_violated"]:
        print("  PRECONDITION  {}".format(note))
    expr = report["expr"]
    print(
        "  expr           byte_identical={} content_identical={} row_order_differs={}".format(
            expr["byte_identical"],
            expr.get("content_identical"),
            expr.get("row_order_differs"),
        )
    )
    for name, stats in sorted(expr.get("numeric", {}).items()):
        if stats["transcripts_differing"]:
            print(
                "      {:<26} {} transcripts differ, max abs delta {}".format(
                    name, stats["transcripts_differing"], stats["max_abs_delta"]
                )
            )
    for name, stats in sorted(expr.get("structural", {}).items()):
        if stats["transcripts_differing"]:
            print(
                "      {:<26} {} transcripts differ".format(
                    name, stats["transcripts_differing"]
                )
            )
    if expr.get("transcripts_only_in_a"):
        print(
            "      only in strand-first: {}".format(
                ", ".join(expr["transcripts_only_in_a"][:20])
            )
        )
    if expr.get("transcripts_only_in_b"):
        print(
            "      only in strandless:   {}".format(
                ", ".join(expr["transcripts_only_in_b"][:20])
            )
        )
    track = report["tracking"]
    print(
        "  tracking       byte_identical={} rows_in_order={} multisets_equal={} "
        "ordering_only={}".format(
            track["byte_identical"],
            track.get("rows_identical_in_order"),
            track.get("row_multisets_equal"),
            track.get("difference_is_ordering_only"),
        )
    )
    print(
        "      rows {} vs {}; keys only-A {} only-B {}; frac_assigned differing {}".format(
            track.get("rows_a"),
            track.get("rows_b"),
            len(track.get("keys_only_in_a", [])),
            len(track.get("keys_only_in_b", [])),
            track.get("frac_assigned_keys_differing"),
        )
    )
    inputs = report.get("chunk_inputs")
    if inputs:
        print(
            "  inputs         {} distinct read names strand-first, {} strandless; "
            "identical={}".format(
                inputs["distinct_read_names_strand_first"],
                inputs["distinct_read_names_strandless"],
                inputs["names_identical"],
            )
        )
        for record in inputs["only_in_strand_first"]:
            print(
                "      only strand-first fed it: {} (chunks {}) -> {}".format(
                    record["read_name"],
                    ",".join(record["strand_first_chunks"]),
                    record["reason"],
                )
            )
        for record in inputs["only_in_strandless"]:
            print(
                "      only strandless fed it:   {} (chunks {}) -> {}".format(
                    record["read_name"],
                    ",".join(record["strandless_chunks"]),
                    record["reason"],
                )
            )
    severed = report["severed"]
    print(
        "  severed        strand-first {} + {} overhang; strandless {} + {} "
        "overhang; symmetric difference {}".format(
            len(severed["strand_first_cut_selection"]),
            len(severed["strand_first_overhang"]),
            len(severed["strandless_cut_selection"]),
            len(severed["strandless_overhang"]),
            len(severed["symmetric_difference"]),
        )
    )
    for record in severed.get("symmetric_difference_detail", []):
        print(
            "      lost only by {:<12} {}  (in strand-first tracking: {}; in "
            "strandless tracking: {})".format(
                record["lost_by"],
                record["read_name"],
                record["in_strand_first_tracking"],
                record["in_strandless_tracking"],
            )
        )
    for label, records in sorted(report["reads_absent_by_name"].items()):
        for record in records:
            print(
                "  {} {} -> {}".format(
                    label.upper().replace("_", " "),
                    record["read_name"],
                    record["reason"],
                )
            )
    attribution = report["attribution"]
    for column, transcripts in sorted(
        attribution["unexplained_expr_columns"].items()
    ):
        print(
            "  UNEXPLAINED {} for {} transcripts: {}".format(
                column, len(transcripts), ", ".join(transcripts[:20])
            )
        )
    if attribution["tpm_transcripts_differing"]:
        print(
            "  TPM differs for {} transcripts, max abs delta {} -- globally coupled "
            "to the merged total, so one differing read moves every row".format(
                attribution["tpm_transcripts_differing"],
                attribution["tpm_max_abs_delta"],
            )
        )
    print("  {}".format(report["verdict"]))


def _print_gate(report):
    print("=" * 78)
    print("STRANDLESS CHUNKING GATE  {}".format(report["out_dir"]))
    print("=" * 78)
    for step in report["steps"]:
        print("  {}".format(step))
    for key in (
        "invariant_aligner_flag",
        "invariant_inferred",
        "invariant_pipeline_artifact",
    ):
        if key in report:
            print("")
            print("[{}]".format(key))
            _print_invariant(report[key])
    if "ordering_cost" in report:
        cost = report["ordering_cost"]
        print("")
        print("[ordering_cost] {}".format(cost["region"]))
        print(
            "  {} of {} emitted records would be misassigned by a strand-filtered "
            "extraction from the raw bam".format(
                len(cost["misassigned_read_names"]),
                cost["records_emitted_by_the_strand_filter"],
            )
        )
    if "strand_first_figures" in report:
        print("")
        print("[strand_first_figures] {}".format(report["strand_first_figures"]))
    if "strandless_figures" in report:
        print("[strandless_figures]   {}".format(report["strandless_figures"]))
    if "arms" in report:
        print("")
        _print_arms(report["arms"])
    print("")
    if report["failures"]:
        print("GATE FAILED:")
        for failure in report["failures"]:
            print("  - {}".format(failure))
    else:
        print("GATE PASSED")


def build_parser():
    parser = argparse.ArgumentParser(
        description="equivalence gate for strandless chunking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    inv = sub.add_parser(
        "strand-invariant",
        help="per-read agreement between the whole-scope and in-chunk strand split",
    )
    inv.add_argument("--bam", required=True)
    inv.add_argument("--genome_fa", required=True)
    inv.add_argument("--gtf", default=None)
    inv.add_argument("--work_dir", required=True)
    inv.add_argument("--contig", default=None)
    inv.add_argument(
        "--region",
        action="append",
        default=[],
        help="strandless region, repeatable, e.g. chr1:1-10000000",
    )
    inv.add_argument("--cuts_json", default=None, help="a selector <prefix>.cuts.json")
    inv.add_argument(
        "--chunks_root",
        default=None,
        help="a real strandless run's chunks/ directory, consumed as it stands",
    )
    inv.add_argument("--limit_regions", type=int, default=None)
    inv.add_argument("--margin", type=int, default=extractor.DEFAULT_MARGIN)
    inv.add_argument(
        "--max_intron_length",
        type=int,
        default=LRAA_Globals.config["max_intron_length"],
    )
    inv.add_argument(
        "--infer_read_orient",
        action="store_true",
        help="check the invariant for INFERRED orientation rather than the "
        "aligner flag. The pipeline does not pass this today, so it is the "
        "forward-looking case: inference is where a flip can happen at all",
    )
    inv.add_argument("--keep_bams", action="store_true")
    inv.add_argument("--gtf_index_cache_dir", default=None)
    inv.add_argument("--report", default=None, help="write the report as JSON here")

    arms = sub.add_parser(
        "compare-arms", help="diff a strand-first and a strandless pipeline run"
    )
    arms.add_argument("--strand_first_dir", required=True)
    arms.add_argument("--strandless_dir", required=True)
    arms.add_argument("--report", default=None)
    arms.add_argument(
        "--no_explain_reads",
        action="store_true",
        help="skip the per-read absence reasons, which stream every chunk BAM",
    )

    order = sub.add_parser(
        "ordering-cost",
        help="what a strand-filtered extraction from a raw BAM would misassign",
    )
    order.add_argument("--bam", required=True)
    order.add_argument("--genome_fa", required=True)
    order.add_argument("--gtf", default=None)
    order.add_argument("--work_dir", required=True)
    order.add_argument(
        "--region",
        required=True,
        help="a STRANDED region, e.g. chr1+:1-10000000 -- the mistake, on purpose",
    )
    order.add_argument("--margin", type=int, default=extractor.DEFAULT_MARGIN)
    order.add_argument(
        "--max_intron_length",
        type=int,
        default=LRAA_Globals.config["max_intron_length"],
    )
    order.add_argument("--report", default=None)

    whole = sub.add_parser(
        "gate",
        help="the whole gate: invariant, ordering cost, both arms, comparison",
    )
    whole.add_argument("--bam", required=True)
    whole.add_argument("--genome_fa", required=True)
    whole.add_argument("--gtf", required=True)
    whole.add_argument("--output_dir", required=True)
    whole.add_argument("--contig", default=None)
    whole.add_argument("--HiFi", action="store_true")
    whole.add_argument("--approx_MB_per_cut", type=float, default=None)
    whole.add_argument(
        "--approx_MB_per_cut_wiggle_window", type=float, default=None
    )
    whole.add_argument("--cpu_budget", type=int, default=4)
    whole.add_argument("--margin", type=int, default=None)
    whole.add_argument("--max_intron_length", type=int, default=None)
    whole.add_argument(
        "--ordering_cost_region",
        default=None,
        help="a STRANDED region to measure the ordering mistake against, e.g. "
        "chr1+:1-10000000",
    )
    whole.add_argument("--skip_infer_invariant", action="store_true")
    whole.add_argument(
        "--strand_first_arm",
        choices=("chunked", "both"),
        default="both",
        help="whether the strand-first run also produces the whole-contig "
        "control. Only 'both' can pin baseline_excluded_severed_reads; "
        "arm-to-arm equivalence never reads it",
    )
    whole.add_argument(
        "--skip_arms",
        action="store_true",
        help="stop after the invariant; the arms are the expensive half",
    )
    whole.add_argument(
        "--expect",
        default=None,
        help="comma-separated key=value pins on the STRAND-FIRST arm, e.g. "
        "expr_rows=555,tracking_rows=541,num_total_reads=2266,"
        "baseline_excluded_severed_reads=5",
    )
    whole.add_argument("--report", default=None)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.command == "strand-invariant":
        intervals = None
        chunks_root = args.chunks_root
        if chunks_root is None:
            if args.cuts_json:
                intervals = intervals_from_cuts_json(args.cuts_json)
            elif args.region:
                intervals = intervals_from_region_strings(args.region)
            else:
                raise ParityError(
                    "need one of --region, --cuts_json or --chunks_root"
                )
            if args.limit_regions:
                intervals = intervals[: args.limit_regions]
        report = check_strand_invariant(
            args.bam,
            args.genome_fa,
            args.work_dir,
            intervals=intervals,
            chunks_root=chunks_root,
            gtf=args.gtf,
            margin=args.margin,
            max_intron_length=args.max_intron_length,
            infer_read_orient=args.infer_read_orient,
            contig=args.contig,
            keep_bams=args.keep_bams,
            gtf_index_cache_dir=args.gtf_index_cache_dir,
        )
        _print_invariant(report)
        ok = report["verdict"] == "INVARIANT HOLDS"
    elif args.command == "ordering-cost":
        report = measure_ordering_violation(
            args.bam,
            args.genome_fa,
            args.region,
            args.work_dir,
            gtf=args.gtf,
            margin=args.margin,
            max_intron_length=args.max_intron_length,
        )
        print("EXTRACTION ORDERING COST [{}]".format(report["region"]))
        print(
            "  emitted by the raw-flag strand filter   {}".format(
                report["records_emitted_by_the_strand_filter"]
            )
        )
        print(
            "  records the split would move            {}".format(
                report["records_the_split_would_move"]
            )
        )
        print(
            "  distinct reads silently misassigned     {}".format(
                len(report["misassigned_read_names"])
            )
        )
        for name in report["misassigned_read_names"][:20]:
            print("      {}".format(name))
        print("  {}".format(report["verdict"]))
        # A demonstration, not a gate: finding misassignments is the POINT.
        ok = True
    elif args.command == "gate":
        expect = {}
        for item in (args.expect or "").split(","):
            if not item.strip():
                continue
            key, _, value = item.partition("=")
            expect[key.strip()] = int(value)
        report = run_gate(
            args.bam,
            args.genome_fa,
            args.gtf,
            args.output_dir,
            contig=args.contig,
            HiFi=args.HiFi,
            approx_MB_per_cut=args.approx_MB_per_cut,
            approx_MB_per_cut_wiggle_window=args.approx_MB_per_cut_wiggle_window,
            cpu_budget=args.cpu_budget,
            margin=args.margin,
            max_intron_length=args.max_intron_length,
            expect=expect,
            skip_infer_invariant=args.skip_infer_invariant,
            skip_arms=args.skip_arms,
            ordering_cost_region=args.ordering_cost_region,
            strand_first_arm=args.strand_first_arm,
        )
        _print_gate(report)
        ok = not report["failures"]
    else:
        report = compare_arms(
            args.strand_first_dir,
            args.strandless_dir,
            explain_reads=not args.no_explain_reads,
        )
        _print_arms(report)
        ok = report["passed"] and not report["preconditions_violated"]

    if args.report:
        with open(args.report, "wt") as fh:
            json.dump(report, fh, indent=2, sort_keys=True, default=sorted)
            print("", file=fh)
        print("report: {}".format(args.report))
    return 0 if ok else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except ParityError as err:
        print("\nGATE FAILED\n{}".format(err), file=sys.stderr)
        sys.exit(2)
