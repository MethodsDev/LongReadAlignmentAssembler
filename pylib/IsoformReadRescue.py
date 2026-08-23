#!/usr/bin/env python3

import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from contextlib import contextmanager
import hashlib

import pysam
from intervaltree import IntervalTree

import LRAA_Globals
import Util_funcs
from MultiPath import MultiPath
from Pretty_alignment import Pretty_alignment


logger = logging.getLogger(__name__)

# Why a read is offered to transcriptome rescue. Shared vocabulary rather than four
# string literals per call site, because the batch path and the streaming path must
# report the same populations under the same names for their candidate sets to be
# comparable at all -- a typo in one of them reads as a category the other never
# produced. low_perID deliberately reuses Util_funcs.quant_discard_reason's key, since
# that predicate is what puts a read in this category.
RESCUE_CANDIDATE_LOW_PER_ID = "low_perID"
RESCUE_CANDIDATE_SPACER_PATH = "spacer_path"
RESCUE_CANDIDATE_NO_GRAPH_PATH = "no_graph_path"
RESCUE_CANDIDATE_UNASSIGNED_TO_TARGETS = "unassigned_to_targets"

RESCUE_CANDIDATE_CATEGORIES = (
    RESCUE_CANDIDATE_LOW_PER_ID,
    RESCUE_CANDIDATE_SPACER_PATH,
    RESCUE_CANDIDATE_NO_GRAPH_PATH,
    RESCUE_CANDIDATE_UNASSIGNED_TO_TARGETS,
)


def _no_rescue_read_details():
    """The shape _parse_rescue_alignments returns, with nothing in it.

    One definition rather than one per early return: the two early returns below used to
    disagree about which keys they carried, so which keys a caller could rely on depended
    on which way rescue happened to bail out.
    """
    return {
        "read_name_to_multipaths": {},
        "alignment_rejections": {},
        "reads_locality_declined": 0,
        "reads_locality_displaced": 0,
    }


def rescue_unassigned_reads_to_transcriptome(
    splice_graph,
    transcripts,
    contig_seq_str,
    bam_file,
    contig_acc,
    region_lend,
    region_rend,
    read_names,
    read_path_mapper=None,
    return_read_details=False,
):
    """
    Rescue previously unassigned reads by aligning them to local transcript sequences and
    projecting accepted transcript hits back into splice-graph node paths.

    This is the only transcriptome-alignment path there is. Rescue is on or off; there is
    no mode that assigns every read by transcript alignment, and none that adjudicates a
    whole-genome assignment against a whole-transcriptome one. Rescue is always a retry
    for reads the genomic graph could not place, always confined to the locus the read's
    own genome alignment sits in, and always judged against what that alignment already
    explains.

    Acceptance rules:
    - mapped alignment only
    - the target's exons must overlap the read's own genome alignment blocks. The genome
      alignment is trusted for LOCALITY but not for optimality: rescue runs because the
      alignment fits badly, so it may refine placement within the locus the read sits in
      and may not move the read to another one. Measured on chr21 HG002 quant-only
      ref-guided, 22 of the 33 (read, transcript) credits rescue made went to models the
      read does not overlap, the farthest by 26.6 Mb. Coordinate overlap only -- no
      minimum overlap, no fraction of the read, no intron-chain or path-compatibility
      reasoning, since a suboptimal alignment's structure is not evidence about anything.
    - no reference-skip (N) cigar operations
    - after Pretty_alignment-style small-gap merging, exactly one merged transcript block remains
    - if multiple top-scoring transcript hits survive for a read, they must all project to the same node path
    """

    minimap2_exe = shutil.which("minimap2")
    if minimap2_exe is None:
        # Skipping here would drop every read this call was meant to recover and
        # let the run report the smaller number as if it were the answer.
        raise RuntimeError(
            "minimap2 not found in PATH, so transcriptome rescue cannot run. "
            "Install minimap2, or disable it with "
            "--no_rescue_unassigned_reads_via_transcriptome_alignment."
        )

    transcript_models = _build_transcript_models(
        splice_graph,
        transcripts,
        contig_seq_str,
    )
    if not transcript_models:
        if return_read_details:
            return [], _no_rescue_read_details()
        return []

    (
        read_name_to_seq,
        read_name_to_genome_explained,
        read_name_to_genome_gap_id,
        read_name_to_allowed_target_ids,
    ) = _collect_read_sequences(
        bam_file,
        contig_acc,
        region_lend,
        region_rend,
        read_names,
        splice_graph.get_contig_strand(),
        # Locality is unconditional and is not a mode: every rescue is confined to
        # targets the read's own genome alignment overlaps. It is one additional
        # rejection layered onto the rules already here -- same offered reads, same
        # thresholds, same genome baseline, same order.
        exon_overlap_index=_build_exon_overlap_index(transcript_models),
    )
    # Reads that lack a graph path get no exemption from the genome baselines above.
    # Lacking a graph path says the graph cannot represent the read's structure; it says
    # nothing about whether the read's genome alignment is measurable. Measured on chr20
    # HG002, the rescued reads that lacked a path had genome alignments explaining a
    # median 99.8% of the read at MAPQ 60, so the baseline is both available and
    # meaningful for them. Declining a rescue that explains less than that baseline costs
    # nothing either: the read simply stays unassigned, which is where it already was.

    if not read_name_to_seq:
        logger.info(
            "[%s%s] transcriptome rescue skipped: no reads with retrievable sequences",
            splice_graph.get_contig_acc(),
            splice_graph.get_contig_strand(),
        )
        if return_read_details:
            return [], _no_rescue_read_details()
        return []

    # Keep the scratch tree inside the run's per-(contig,strand) temp dir so that it is
    # swept by the normal contigtmp cleanup (and by `rm -rf __*` in the test Makefiles)
    # rather than accumulating in the working directory. Falls back to cwd when the
    # rescue is invoked outside a contig worker (e.g. unit tests).
    tmp_parent = os.environ.get("LRAA_TMP_DIR") or os.getcwd()
    if not os.path.isdir(tmp_parent):
        tmp_parent = os.getcwd()
    tmp_dir = tempfile.mkdtemp(
        prefix=f"tx_rescue.{contig_acc}.{splice_graph.get_contig_strand()}.",
        dir=tmp_parent,
    )
    keep_tmp = bool(LRAA_Globals.DEBUG) or bool(
        LRAA_Globals.config.get("no_cleanup", False)
    )
    try:
        transcript_fa = os.path.join(tmp_dir, "transcripts.fa")
        reads_fa = os.path.join(tmp_dir, "reads.fa")
        rescue_sam = os.path.join(tmp_dir, "rescue.sam")

        _write_transcript_fasta(transcript_fa, transcript_models)
        _write_reads_fasta(reads_fa, read_name_to_seq)
        _run_minimap2_transcriptome_alignment(
            transcript_fa, reads_fa, rescue_sam, minimap2_exe
        )

        rescued_mps, read_details = _parse_rescue_alignments(
            rescue_sam,
            splice_graph,
            transcript_models,
            read_path_mapper=read_path_mapper,
            read_name_to_allowed_target_ids=read_name_to_allowed_target_ids,
            read_name_to_genome_explained=read_name_to_genome_explained,
            read_name_to_genome_gap_id=read_name_to_genome_gap_id,
        )

        # Rescue outcomes: how many reads were offered, how many came back, and for those
        # that did not, which rule declined them. locality_displaced_reads is the one to
        # watch -- reads whose best transcriptome hit sat at a locus the read is not at,
        # which is exactly what the locality criterion refuses.
        logger.info(
            "[%s%s] transcriptome rescue: offered=%d rescued=%d locality_declined_reads=%d locality_displaced_reads=%d rejections=%s",
            splice_graph.get_contig_acc(),
            splice_graph.get_contig_strand(),
            len(read_name_to_seq),
            len(rescued_mps),
            read_details.get("reads_locality_declined", 0),
            read_details.get("reads_locality_displaced", 0),
            ",".join(
                f"{reason}={count}"
                for reason, count in sorted(
                    read_details.get("alignment_rejections", {}).items()
                )
            )
            or "none",
        )

        if return_read_details:
            return rescued_mps, read_details

        return rescued_mps
    finally:
        if not keep_tmp:
            shutil.rmtree(tmp_dir, ignore_errors=True)


def _run_minimap2_transcriptome_alignment(transcript_fa, reads_fa, rescue_sam, minimap2_exe):
    preset = _resolve_rescue_minimap2_preset()
    try:
        minimap_threads = max(1, int(LRAA_Globals.config.get("tool_threads", 1)))
    except Exception:
        minimap_threads = 1
    cmd = [
        minimap2_exe,
        "-a",
        "-t",
        str(minimap_threads),
        "--secondary=yes",
        "-N",
        "50",
        "-f",
        str(_resolve_rescue_minimap2_filter_fraction()),
    ]
    if preset:
        cmd.extend(["-x", str(preset)])
    cmd.extend([transcript_fa, reads_fa])
    with open(rescue_sam, "wt") as sam_fh:
        subprocess.run(cmd, check=True, stdout=sam_fh)


def _collect_read_sequences(
    bam_file,
    contig_acc,
    region_lend,
    region_rend,
    target_read_names=None,
    target_strand=None,
    exon_overlap_index=None,
):
    """Read sequences plus what each read's own genome alignment establishes about it.

    Three things come off the one primary record in a single pass, because all three are
    properties of that record and a second pass over the region would cost the same scan
    again: the sequence to realign, the baseline the rescue has to beat, and -- when an
    exon overlap index is supplied -- which targets the read's aligned blocks touch.
    """
    remaining = None if target_read_names is None else set(target_read_names)
    read_name_to_seq = {}
    read_name_to_genome_explained = {}
    read_name_to_genome_gap_id = {}
    read_name_to_allowed_target_ids = None if exon_overlap_index is None else {}
    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:
        if region_lend is not None and region_rend is not None:
            fetch_iter = bam_reader.fetch(contig_acc, max(int(region_lend) - 1, 0), int(region_rend))
        else:
            fetch_iter = bam_reader.fetch(contig_acc)
        for read in fetch_iter:
            if read.is_unmapped or read.is_supplementary or read.is_secondary:
                # Only a primary record is guaranteed to carry the full read sequence;
                # hard-clipped supplementary records would supply a truncated one.
                continue
            if target_strand is not None:
                if read.is_forward and target_strand != "+":
                    continue
                if read.is_reverse and target_strand != "-":
                    continue
            read_name = Util_funcs.get_read_name_include_sc_encoding(read)
            if remaining is not None and read_name not in remaining:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            if read.is_reverse:
                seq = _reverse_complement(seq)
            if read_name not in read_name_to_seq:
                read_name_to_seq[read_name] = seq
                # What this read's genome alignment already explains, so a rescue that
                # would explain less of it can be declined. The gap-aware identity is
                # carried alongside because explained bases cannot see a rescue that
                # agrees over a shorter span, having skipped part of its target.
                genome_explained = _explained_read_bases(read)
                if genome_explained is not None:
                    read_name_to_genome_explained[read_name] = genome_explained
                genome_gap_id = _gap_aware_identity(read)
                if genome_gap_id is not None:
                    read_name_to_genome_gap_id[read_name] = genome_gap_id
                if read_name_to_allowed_target_ids is not None:
                    # The genome alignment is trusted for LOCALITY but not assumed
                    # optimal. Rescue exists because that alignment fits badly, so it
                    # refines placement WITHIN the locus the read already occupies and
                    # never relocates the read to a different one. That is the whole of
                    # the criterion, which is why plain coordinate overlap of the
                    # aligned blocks against the target's exons is sufficient, and why
                    # nothing structural is wanted: a suboptimal alignment's intron
                    # chain is not evidence about anything.
                    #
                    # An empty set is recorded rather than skipped. The read stays in
                    # the offered population and is still realigned exactly as before;
                    # it simply has no target it may be credited to. Dropping it here
                    # instead would change which reads rescue is asked about, and this
                    # is meant to add one rejection, not reshape the candidate set.
                    read_name_to_allowed_target_ids[read_name] = set(
                        _get_alignment_overlapping_targets(read, exon_overlap_index)
                    )
                if remaining is not None:
                    remaining.discard(read_name)
            if remaining is not None and not remaining:
                break
    return (
        read_name_to_seq,
        read_name_to_genome_explained,
        read_name_to_genome_gap_id,
        read_name_to_allowed_target_ids,
    )


def _build_transcript_models(
    splice_graph, transcripts, contig_seq_str, include_monoexonic=False
):
    transcript_models = {}
    for transcript in transcripts:
        if not include_monoexonic and not transcript.has_introns():
            continue
        simple_path = transcript.get_simple_path()
        path_no_boundaries = [
            node_id
            for node_id in simple_path
            if not re.match("TSS:|POLYA:", node_id)
        ]
        exon_nodes_genomic = [node_id for node_id in path_no_boundaries if node_id.startswith("E:")]
        if not exon_nodes_genomic:
            continue

        exon_nodes_transcript_order = list(exon_nodes_genomic)
        if transcript.get_strand() == "-":
            exon_nodes_transcript_order = list(reversed(exon_nodes_transcript_order))

        tx_coord_map = []
        tx_pos = 1
        for node_id in exon_nodes_transcript_order:
            node_obj = splice_graph.get_node_obj_via_id(node_id)
            lend, rend = node_obj.get_coords()
            seg_len = rend - lend + 1
            tx_coord_map.append(
                {
                    "tx_lend": tx_pos,
                    "tx_rend": tx_pos + seg_len - 1,
                    "node_id": node_id,
                    "genomic_lend": lend,
                    "genomic_rend": rend,
                }
            )
            tx_pos += seg_len

        exon_index_in_path = {
            node_id: idx for idx, node_id in enumerate(path_no_boundaries)
        }
        gene_id = transcript.get_gene_id()
        transcript_id = transcript.get_transcript_id()
        target_id = f"{gene_id}^{transcript_id}"
        transcript_models[target_id] = {
            "transcript": transcript,
            "gene_id": gene_id,
            "transcript_id": transcript.get_transcript_id(),
            "target_id": target_id,
            "sequence": _build_transcript_sequence(transcript, contig_seq_str),
            "genomic_exon_segments": transcript.get_exon_segments(),
            "path_no_boundaries": path_no_boundaries,
            "tx_coord_map": tx_coord_map,
            "path_index": exon_index_in_path,
        }
    return transcript_models


def _build_transcript_sequence(transcript, contig_seq_str):
    exon_segments = transcript.get_exon_segments()
    exon_seqs = [contig_seq_str[lend - 1 : rend] for lend, rend in exon_segments]
    if transcript.get_strand() == "-":
        exon_seqs = [_reverse_complement(seq) for seq in reversed(exon_seqs)]
    return "".join(exon_seqs)


def _reverse_complement(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _write_transcript_fasta(transcript_fa, transcript_models):
    with open(transcript_fa, "wt") as ofh:
        for target_id, model in transcript_models.items():
            print(f">{target_id}", file=ofh)
            print(model["sequence"], file=ofh)


def _write_reads_fasta(reads_fa, read_name_to_seq):
    with open(reads_fa, "wt") as ofh:
        for read_name, seq in read_name_to_seq.items():
            print(f">{read_name}", file=ofh)
            print(seq, file=ofh)


@contextmanager
def _rescue_alignment_records(source):
    """Alignment records from a SAM path, or an already-open iterable of them.

    Streaming rescue builds pysam records from mappy hits instead of reading a SAM
    file, and has to run the very same acceptance loop over them: a second
    implementation of those rules is exactly how the two rescue paths would drift
    apart without anything reporting it.
    """
    if isinstance(source, (str, bytes, os.PathLike)):
        with pysam.AlignmentFile(source, "r") as sam_reader:
            yield sam_reader.fetch(until_eof=True)
    else:
        yield source


def _parse_rescue_alignments(
    rescue_alignments,
    splice_graph,
    transcript_models,
    read_path_mapper=None,
    read_name_to_allowed_target_ids=None,
    read_name_to_genome_explained=None,
    read_name_to_genome_gap_id=None,
):
    read_to_hits = defaultdict(list)
    # Why alignments were declined, in rescue's own terms. Rescue previously reported
    # only how many reads came back, which leaves an unrescued read indistinguishable
    # from an unoffered one and gives no way to see which rule is doing the work --
    # locality above all, since it is the one that can silently do nothing.
    rejections = defaultdict(int)
    # Best score among a read's locality-declined alignments, so the report can name the
    # population that matters: reads whose BEST transcriptome hit was somewhere the read
    # is not. A read with only a lower-scoring non-local hit lost nothing.
    locality_declined_best_score = {}
    min_per_id = float(_resolve_rescue_min_per_id())
    min_aligned_frac = float(
        LRAA_Globals.config.get("rescue_unassigned_min_aligned_read_frac", 0) or 0
    )
    max_indel_length = int(
        LRAA_Globals.config.get("rescue_unassigned_max_indel_length", 0) or 0
    )
    with _rescue_alignment_records(rescue_alignments) as alignment_records:
        for read in alignment_records:
            if read.is_unmapped or read.is_supplementary:
                rejections["unmapped_or_supplementary"] += 1
                continue
            if read.reference_name not in transcript_models:
                rejections["target_not_modelled"] += 1
                continue
            if any(code == 3 for code, _ in (read.cigartuples or [])):
                rejections["reference_skip"] += 1
                continue
            # A long indel is structural disagreement, not sequencing error: the read
            # skips (D) or carries (I) a stretch the target does not share, which is
            # what an exon-level difference looks like once the alignment is against a
            # cDNA rather than the genome. Nothing else here sees it -- the aligned
            # fraction counts query bases, so a deletion costs it nothing; the
            # explained-bases baseline excludes deleted bases by construction; and the
            # single-merged-segment test cannot fire because a deletion emits its own
            # genome segment and bridges its own gap.
            #
            # Calibrated against reads realigned to the transcript their genome intron
            # chain exactly matches, where every indel is platform error: the largest
            # deletion seen was 32 on PacBio HiFi and 45 on ONT cDNA.
            if max_indel_length > 0 and any(
                code in (1, 2) and length >= max_indel_length
                for code, length in (read.cigartuples or [])
            ):
                rejections["long_indel"] += 1
                continue
            if not _passes_percent_identity(read, min_per_id):
                rejections["low_per_id"] += 1
                continue
            # The genome-baseline test below only constrains reads that already had a
            # genome alignment. Reads with no graph path -- the bulk of what rescue
            # retries -- reach here with nothing bounding how much of the read the
            # target accounts for: percent identity is measured over the aligned
            # portion alone and is blind to clipping, so a 200bp end of a 3kb read at
            # 85% identity would count as full support for the isoform. Require the
            # alignment to span essentially the whole read. Measured as aligned length
            # over read length rather than as matched bases, so ONT/PacBio error rates
            # do not make the requirement unsatisfiable; mismatches are already
            # constrained by rescue_unassigned_min_per_id.
            if min_aligned_frac > 0:
                read_length = read.infer_read_length() or read.query_length
                aligned_length = read.query_alignment_length
                if read_length and aligned_length is not None:
                    if (aligned_length / read_length) < min_aligned_frac:
                        rejections["low_aligned_read_frac"] += 1
                        continue
            if read_name_to_genome_explained:
                # An alignment explaining less of the read than the genome already does
                # cannot correct anything; it only adds a competing path shaped like the
                # target. It is not merely additive either: an accepted rescue detaches
                # the read from the path it already had, in
                # LRAA._detach_rescued_reads_from_original_paths, and a path left with no
                # reads is dropped -- so admitting a worse alignment actively withdraws
                # support from the better structure.
                #
                # A read with no genome alignment has no baseline here and is
                # unconstrained. Reads lacking a graph path are NOT exempt: they keep
                # their baseline like any other read.
                genome_explained = read_name_to_genome_explained.get(read.query_name)
                if genome_explained is not None:
                    rescue_explained = _explained_read_bases(read)
                    if (
                        rescue_explained is not None
                        and rescue_explained < genome_explained
                    ):
                        rejections["explains_less_than_genome"] += 1
                        continue
            if read_name_to_genome_gap_id:
                # Explained bases cannot separate a clean spliced genome alignment from
                # a transcript alignment that skips part of its target: deleted bases
                # are not read bases, so both explain the read equally and the test
                # above ties. Gap-aware identity charges the skipped span, so the
                # rescue must agree with its target at least as well as the genome
                # alignment agrees with the genome, over the extent each covers.
                genome_gap_id = read_name_to_genome_gap_id.get(read.query_name)
                if genome_gap_id is not None:
                    rescue_gap_id = _gap_aware_identity(read)
                    if rescue_gap_id is not None and rescue_gap_id < genome_gap_id:
                        rejections["agrees_worse_than_genome"] += 1
                        continue
            # Locality, last of the content rules and the only one about WHERE the
            # target is rather than how well it fits. Everything above has already
            # passed, so an alignment declined here was a viable placement refused for
            # sitting somewhere the read does not: that is what makes the counters below
            # mean something. Still above the merge and projection, which are the
            # expensive steps.
            #
            # The genome alignment is trusted for locality and not for optimality --
            # rescue runs precisely because it fits badly -- so the test is plain
            # coordinate overlap of the read's aligned blocks against the target's exons,
            # computed once per read in _collect_read_sequences. No minimum overlap, no
            # fraction of the read, no intron-chain or path-compatibility reasoning; a
            # suboptimal alignment's structure is not evidence about anything.
            if read_name_to_allowed_target_ids is not None and (
                read.query_name not in read_name_to_allowed_target_ids
                or read.reference_name
                not in read_name_to_allowed_target_ids[read.query_name]
            ):
                rejections["locality"] += 1
                score = _alignment_score(read)
                prior = locality_declined_best_score.get(read.query_name)
                if prior is None or score > prior:
                    locality_declined_best_score[read.query_name] = score
                continue
            merged_segments = Pretty_alignment.read_to_pretty_alignment_segments(read)
            if len(merged_segments) != 1:
                rejections["not_one_merged_block"] += 1
                continue
            merged_lend, merged_rend = merged_segments[0]
            left_soft_clipping, right_soft_clipping = _get_soft_clipping_lengths(read)
            projected_path = _project_alignment_to_graph_path(
                transcript_models[read.reference_name],
                merged_lend,
                merged_rend,
                read_path_mapper=read_path_mapper,
                left_soft_clipping=left_soft_clipping,
                right_soft_clipping=right_soft_clipping,
            )
            if not projected_path:
                rejections["no_graph_projection"] += 1
                continue
            read_to_hits[read.query_name].append(
                {
                    "score": _alignment_score(read),
                    "path": tuple(projected_path),
                }
            )

    rescued_mps = []
    read_name_to_multipaths = defaultdict(list)
    for read_name, hits in read_to_hits.items():
        read_key = _normalize_read_identifier(read_name)
        if read_key is None:
            continue
        best_score = max(hit["score"] for hit in hits)
        best_hits = [hit for hit in hits if hit["score"] == best_score]
        # Every top-scoring hit has to imply the same graph path. A read whose best
        # alignments disagree about its structure is not evidence for either of them.
        projected_paths = {hit["path"] for hit in best_hits}
        if len(projected_paths) != 1:
            continue
        mp_paths = [list(path_tuple) for path_tuple in sorted(projected_paths)]
        multipath = MultiPath(
            splice_graph,
            mp_paths,
            read_types={"PacBio"},
            read_names={read_name},
            read_count=1,
        )
        rescued_mps.append(multipath)
        read_name_to_multipaths[read_key].append(multipath)

    # Reads whose placement locality actually changed: some alignment that cleared every
    # content rule was declined for locality, and it outscored whatever survived -- or
    # nothing survived, so the read would have been credited non-locally or not at all.
    # This is the size of the defect. Counted after the loop because "best" is not known
    # until every alignment has been seen.
    locality_displaced = 0
    for declined_read_name, declined_score in locality_declined_best_score.items():
        kept = read_to_hits.get(declined_read_name)
        if not kept or declined_score > max(hit["score"] for hit in kept):
            locality_displaced += 1

    return rescued_mps, {
        "read_name_to_multipaths": dict(read_name_to_multipaths),
        "alignment_rejections": dict(rejections),
        "reads_locality_declined": len(locality_declined_best_score),
        "reads_locality_displaced": locality_displaced,
    }


def _build_exon_overlap_index(transcript_models):
    exon_overlap_index = IntervalTree()
    for target_id, model in transcript_models.items():
        for lend, rend in model["genomic_exon_segments"]:
            exon_overlap_index[lend - 1 : rend] = target_id
    return exon_overlap_index


def _get_alignment_overlapping_targets(read, exon_overlap_index):
    target_id_to_overlap_bp = defaultdict(int)
    for block_lend, block_rend in read.get_blocks():
        for interval in exon_overlap_index.overlap(block_lend, block_rend):
            overlap_bp = min(block_rend, interval.end) - max(block_lend, interval.begin)
            if overlap_bp > 0:
                target_id_to_overlap_bp[interval.data] += overlap_bp
    return target_id_to_overlap_bp


def _passes_percent_identity(read, min_per_id):
    """Percent identity for acceptance, measured gap-aware.

    The previous formula divided matched-minus-NM by the matched count alone. NM
    counts inserted and deleted bases while the denominator excludes them, so the
    result was not an identity: it could go negative, and it scored an alignment with
    zero substitutions and a 592-base insertion at 0%. It was also the only gate that
    saw indels at all, which made a real defect load-bearing.

    _gap_aware_identity() charges gap bases to the denominator instead, which is a
    standard definition, cannot go negative, and preserves the behaviour that mattered:
    across 313 rescue alignments the two disagreed on 2 acceptance decisions at the
    HiFi threshold of 97, both within 0.001 of it.
    """
    per_id_fraction = _gap_aware_identity(read)
    if per_id_fraction is None:
        return True
    return (per_id_fraction * 100.0) >= min_per_id


def _resolve_rescue_minimap2_preset():
    preset = LRAA_Globals.config.get("rescue_unassigned_minimap2_preset", "auto")
    if preset in (None, "", "auto"):
        return None
    return str(preset)


def _resolve_rescue_minimap2_filter_fraction():
    filter_fraction = LRAA_Globals.config.get(
        "rescue_unassigned_minimap2_filter_fraction", 0
    )
    if filter_fraction is None:
        return 0
    return filter_fraction


def _resolve_rescue_min_per_id():
    min_per_id = LRAA_Globals.config.get("rescue_unassigned_min_per_id", None)
    if min_per_id is None:
        return float(LRAA_Globals.config["min_per_id"])
    return float(min_per_id)


def _alignment_score(read):
    """The aligner's AS, or a stand-in when it is absent.

    Deliberately uses the RAW edit count, not Util_funcs.substitution_count():
    this stands in for an alignment score, and a score should be charged for
    indels as well as substitutions. That is the opposite of what
    _explained_read_bases() and _gap_aware_identity() need, which is why they
    delegate and this does not.

    Charged via Util_funcs.alignment_edit_count(), which converts nM to an edit
    distance by adding I and D. Reading the raw tag instead made this score blind
    to indels whenever the producer wrote nM, so two candidates differing only in
    indel content tied -- the exact comparison this exists to decide. Same
    alignment, `10M5I10M`: 18 under nM:i:2 against 13 under the equivalent
    NM:i:7. Both now score 13.
    """
    if read.has_tag("AS"):
        return int(read.get_tag("AS"))

    edit_count = Util_funcs.alignment_edit_count(read)
    if edit_count is None:
        edit_count = 0
    # Shared primitive: M + = + X. Computing it inline took M alone unless M was
    # exactly zero, which under-counts the 0.32% of records mixing M with =/X.
    aligned = Util_funcs.aligned_base_count(read)
    return int(aligned) - int(edit_count)


def _explained_read_bases(read):
    """Read bases this alignment accounts for as matches.

    Clipped bases are excluded, so a target that disagrees with part of a read cannot
    score well by clipping the disagreement away; percent identity, measured over the
    aligned portion only, is blind to exactly that. Computed the same way for a genome
    alignment and for a transcriptome alignment of the same read, so the two compare.
    """
    cigar_stats = read.get_cigar_stats()
    # Shared primitive: M + = + X. Computing it inline took M alone unless M was
    # exactly zero, which under-counts the 0.32% of records mixing M with =/X.
    aligned = Util_funcs.aligned_base_count(read)
    if aligned <= 0:
        return None

    # Only substitutions reduce the read bases an alignment explains: inserted
    # bases are already outside the aligned count, and deleted bases are not read
    # bases at all. Delegated because the indel subtraction is correct for NM and
    # WRONG for nM, which already excludes indels -- doing it inline here treated
    # `10M5I10M` with nM:i:2 as having zero substitutions.
    substitutions = Util_funcs.substitution_count(read)
    if substitutions is None:
        return int(aligned)
    return max(0, int(aligned) - substitutions)


def _gap_aware_identity(read):
    """Matched read bases over the alignment's full extent, gaps included.

    _explained_read_bases() answers "how much of the read did this account for" and
    deliberately ignores deleted reference bases, since those are not read bases. That
    makes it blind to exon-level disagreement: a read skipping 74 bases the target
    contains explains exactly as much of itself as a clean spliced genome alignment
    does, and ties the baseline.

    This counts the skipped bases in the denominator, so an alignment that agrees with
    its target over a shorter span scores lower. Reference skips (N) are excluded --
    an intron is not a gap in the read's agreement with the genome, which is what makes
    a genome alignment comparable to a transcriptome one here.

    Returns None when NM is unavailable, matching _explained_read_bases().
    """
    substitutions = Util_funcs.substitution_count(read)
    if substitutions is None:
        return None

    cigar_stats = read.get_cigar_stats()
    # Shared primitive: M + = + X. Computing it inline took M alone unless M was
    # exactly zero, which under-counts the 0.32% of records mixing M with =/X.
    aligned = Util_funcs.aligned_base_count(read)
    insertions = int(cigar_stats[0][1])
    deletions = int(cigar_stats[0][2])
    span = int(aligned) + insertions + deletions
    if span <= 0:
        return None

    return max(0.0, float(int(aligned) - substitutions) / float(span))



def _normalize_read_identifier(value):
    try:
        if isinstance(value, int):
            return int(value)
        if isinstance(value, str):
            try:
                return int(value)
            except Exception:
                pass
            name_store = getattr(LRAA_Globals, "READ_NAME_STORE", None)
            if name_store is not None:
                rid = name_store.get_or_add(value)
                if rid is not None:
                    return int(rid)
            digest = hashlib.sha1(value.encode("utf-8", "ignore")).hexdigest()
            return int(digest[:16], 16)
    except Exception:
        return None
    return None


def _get_soft_clipping_lengths(read):
    cigar_tuples = read.cigartuples or []
    left_soft_clipping = cigar_tuples[0][1] if cigar_tuples and cigar_tuples[0][0] == 4 else 0
    right_soft_clipping = (
        cigar_tuples[-1][1] if cigar_tuples and cigar_tuples[-1][0] == 4 else 0
    )
    return left_soft_clipping, right_soft_clipping


def _project_alignment_to_graph_path(
    model,
    tx_lend,
    tx_rend,
    read_path_mapper=None,
    left_soft_clipping=None,
    right_soft_clipping=None,
):
    if read_path_mapper is not None:
        genomic_segments = _project_interval_to_genomic_segments(
            model, tx_lend, tx_rend
        )
        if not genomic_segments:
            return None
        genomic_left_soft_clipping = left_soft_clipping
        genomic_right_soft_clipping = right_soft_clipping
        if model["transcript"].get_strand() == "-":
            genomic_left_soft_clipping, genomic_right_soft_clipping = (
                right_soft_clipping,
                left_soft_clipping,
            )
        return read_path_mapper(
            genomic_segments,
            refine_TSS_simple_path=True,
            refine_PolyA_simple_path=True,
            snap_nearby_boundary_features=True,
            left_soft_clipping=genomic_left_soft_clipping,
            right_soft_clipping=genomic_right_soft_clipping,
        )

    return _project_interval_to_path(model, tx_lend, tx_rend)


def _project_interval_to_genomic_segments(model, tx_lend, tx_rend):
    genomic_segments = []
    strand = model["transcript"].get_strand()
    for seg_info in model["tx_coord_map"]:
        seg_lend = seg_info["tx_lend"]
        seg_rend = seg_info["tx_rend"]
        if seg_rend < tx_lend or seg_lend > tx_rend:
            continue

        overlap_tx_lend = max(seg_lend, tx_lend)
        overlap_tx_rend = min(seg_rend, tx_rend)

        offset_lend = overlap_tx_lend - seg_lend
        offset_rend = overlap_tx_rend - seg_lend

        if strand == "+":
            genomic_lend = seg_info["genomic_lend"] + offset_lend
            genomic_rend = seg_info["genomic_lend"] + offset_rend
        else:
            genomic_rend = seg_info["genomic_rend"] - offset_lend
            genomic_lend = seg_info["genomic_rend"] - offset_rend

        genomic_segments.append((genomic_lend, genomic_rend))

    if not genomic_segments:
        return None

    genomic_segments.sort(key=lambda x: (x[0], x[1]))
    return _merge_contiguous_genomic_segments(genomic_segments)


def _merge_contiguous_genomic_segments(genomic_segments):
    if not genomic_segments:
        return genomic_segments

    merged_segments = []
    curr_lend, curr_rend = genomic_segments[0]
    for lend, rend in genomic_segments[1:]:
        if lend <= curr_rend + 1:
            curr_rend = max(curr_rend, rend)
        else:
            merged_segments.append((curr_lend, curr_rend))
            curr_lend, curr_rend = lend, rend
    merged_segments.append((curr_lend, curr_rend))
    return merged_segments


def _project_interval_to_path(model, tx_lend, tx_rend):
    overlapped_exon_nodes = []
    for seg_info in model["tx_coord_map"]:
        seg_lend = seg_info["tx_lend"]
        seg_rend = seg_info["tx_rend"]
        node_id = seg_info["node_id"]
        if not (seg_rend < tx_lend or seg_lend > tx_rend):
            overlapped_exon_nodes.append(node_id)
    if not overlapped_exon_nodes:
        return None

    node_positions = sorted(model["path_index"][node_id] for node_id in overlapped_exon_nodes)
    start_idx = node_positions[0]
    end_idx = node_positions[-1]
    projected = model["path_no_boundaries"][start_idx : end_idx + 1]
    if not projected:
        return None
    return projected


# -----------------------------------------------------------------------------
# Transcriptome Read Rescue Implementation Notes
# -----------------------------------------------------------------------------
#
# High-level purpose
# ------------------
# Transcriptome read rescue is used to recover read evidence that is not usable
# from the genome-aligned BAM path alone. The core idea is:
#
# 1. Build local transcript sequences from the input transcript models already
#    represented in the current splice graph.
# 2. Extract candidate read sequences from the genome BAM.
# 3. Align those reads to the local transcript FASTA with minimap2.
# 4. Project accepted transcriptome alignments back into genome coordinates.
# 5. Convert the projected genome intervals into LRAA MultiPath objects so they
#    can be quantified by the same downstream graph/EM machinery as ordinary
#    genome-derived read paths.
#
# Entry point
# -----------
# rescue_unassigned_reads_to_transcriptome() is the only one. Rescue is on or off:
# it always operates on an explicit read-name set -- the reads the genome read-path
# process could not place -- against the transcript models of the current splice
# graph, and always confined to the targets the read's own genome alignment
# overlaps. There is no mode that assigns every read by transcript alignment, and
# none that adjudicates a whole-genome assignment against a whole-transcriptome
# one.
#
# Transcript model construction
# -----------------------------
# _build_transcript_models() derives a transcript FASTA and a transcript-to-
# genome coordinate map from the current Transcript objects and splice graph.
# Boundary nodes are intentionally removed at this stage:
#
#     path_no_boundaries = simple_path without TSS:/POLYA: nodes
#
# This is because the transcriptome alignment itself is projected through exon
# coordinates first. TSS/PolyA handling is added later when the projected genome
# segments are converted back into a splice-graph read path.
#
# Read sequence extraction and orientation
# ----------------------------------------
# _collect_read_sequences() pulls read sequences from the genome BAM. If the BAM
# record is on the reverse strand, the read sequence is reverse-complemented before
# transcriptome alignment. This makes the extracted sequence correspond to the
# transcript/cDNA orientation. The same pass records the two genome baselines a
# rescue has to beat and the targets the read's aligned blocks touch, because all
# three are properties of the one primary record.
#
# Transcriptome alignment
# -----------------------
# _run_minimap2_transcriptome_alignment() runs minimap2 against the local
# transcript FASTA with secondary alignments enabled. The rescue alignment uses
# the NATIVE TOOL thread count from LRAA_Globals.config["tool_threads"],
# retains up to 50 secondary alignments (-N 50), and uses the configured minimap2
# filter fraction, defaulting to -f 0 for permissive rescue sensitivity.
#
# Alignment acceptance
# --------------------
# _parse_rescue_alignments() filters transcriptome alignments before producing
# read paths. Accepted rescue alignments must:
#
# - be mapped and non-supplementary;
# - map to one of the local transcript models;
# - overlap the read's own genome alignment blocks (locality, unconditional);
# - contain no reference-skip (N) cigar operation;
# - pass the configured percent identity threshold;
# - collapse to exactly one merged transcript block after Pretty_alignment-style
#   small-gap merging.
#
# For each read, only best-scoring transcriptome alignments are retained, and the
# read is rescued only when all of them project to the same graph path.
#
# Projection back to genome/read paths
# ------------------------------------
# _project_alignment_to_graph_path() is where transcript-aligned intervals are
# turned into LRAA graph paths.
#
# Normal LRAA path:
#     If read_path_mapper is provided, the transcript interval is first projected
#     to genomic exon segments with _project_interval_to_genomic_segments().
#     Those genomic segments are then passed to read_path_mapper() with:
#
#         refine_TSS_simple_path=True
#         refine_PolyA_simple_path=True
#         snap_nearby_boundary_features=True
#         left_soft_clipping=<from transcriptome alignment>
#         right_soft_clipping=<from transcriptome alignment>
#
#     This is the genome-equivalent rescue mode. It means rescued transcriptome
#     alignments use the same splice-graph read-path construction logic used for
#     genome alignments, including TSS and PolyA boundary-node refinement.
#
#     Both LRAA execution paths that rescue provide read_path_mapper: final
#     quantification, and the early transcriptome rescue during isoform
#     reconstruction.
#
#     Therefore, normal LRAA quant/rescue operation uses this genome-equivalent
#     path and can incorporate TSS/PolyA nodes through the shared mapper.
#
# Fallback path:
#     If read_path_mapper is not provided, which is not expected in normal LRAA
#     execution and should only occur in direct helper calls, legacy code, or
#     narrowly scoped tests, _project_interval_to_path() maps the transcript
#     interval directly onto the transcript model's exon/intron nodes with
#     path_no_boundaries. In that fallback mode, TSS and PolyA nodes are not
#     added because the shared genome read-path mapper is not being used.
#
# Practical consequence
# ---------------------
# In current LRAA quant/rescue usage, read_path_mapper is supplied. Rescue paths
# should therefore include TSS and PolyA nodes identically to genome-derived read
# paths whenever the projected genomic segments and soft-clipping evidence meet
# the existing splice-graph boundary-refinement criteria. The exon/intron-only
# fallback exists for direct helper usage, not as the intended production mode.
