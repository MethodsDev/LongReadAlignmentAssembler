#!/usr/bin/env python3

import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict

import pysam

import LRAA_Globals
import Util_funcs
from MultiPath import MultiPath
from Pretty_alignment import Pretty_alignment


logger = logging.getLogger(__name__)


def rescue_unassigned_reads_to_transcriptome(
    splice_graph,
    transcripts,
    contig_seq_str,
    bam_file,
    contig_acc,
    region_lend,
    region_rend,
    read_names,
):
    """
    Rescue previously unassigned reads by aligning them to local transcript sequences and
    projecting accepted transcript hits back into splice-graph node paths.

    Acceptance rules:
    - mapped alignment only
    - no reference-skip (N) cigar operations
    - after Pretty_alignment-style small-gap merging, exactly one merged transcript block remains
    - if multiple top-scoring transcript hits survive for a read, they must all project to the same node path
    """

    if not read_names:
        return []

    minimap2_exe = shutil.which("minimap2")
    if minimap2_exe is None:
        logger.warning("minimap2 not found in PATH; skipping transcriptome rescue")
        return []

    transcript_models = _build_transcript_models(splice_graph, transcripts, contig_seq_str)
    if not transcript_models:
        return []

    read_name_to_seq = _collect_read_sequences(
        bam_file, contig_acc, region_lend, region_rend, read_names
    )
    if not read_name_to_seq:
        logger.info(
            "[%s%s] transcriptome rescue skipped: no failed reads with retrievable sequences",
            splice_graph.get_contig_acc(),
            splice_graph.get_contig_strand(),
        )
        return []

    tmp_dir = tempfile.mkdtemp(
        prefix=f"tx_rescue.{contig_acc}.{splice_graph.get_contig_strand()}.",
        dir=os.getcwd(),
    )
    keep_tmp = bool(LRAA_Globals.DEBUG)
    try:
        transcript_fa = os.path.join(tmp_dir, "transcripts.fa")
        reads_fa = os.path.join(tmp_dir, "reads.fa")
        rescue_sam = os.path.join(tmp_dir, "rescue.sam")

        _write_transcript_fasta(transcript_fa, transcript_models)
        _write_reads_fasta(reads_fa, read_name_to_seq)
        _run_minimap2_transcriptome_alignment(transcript_fa, reads_fa, rescue_sam, minimap2_exe)

        rescued_mps = _parse_rescue_alignments(
            rescue_sam, splice_graph, transcript_models
        )

        logger.info(
            "[%s%s] transcriptome rescue: requested=%d rescued=%d",
            splice_graph.get_contig_acc(),
            splice_graph.get_contig_strand(),
            len(read_name_to_seq),
            len(rescued_mps),
        )

        return rescued_mps
    finally:
        if not keep_tmp:
            shutil.rmtree(tmp_dir, ignore_errors=True)


def _run_minimap2_transcriptome_alignment(transcript_fa, reads_fa, rescue_sam, minimap2_exe):
    preset = LRAA_Globals.config.get("rescue_unassigned_minimap2_preset", "map-hifi")
    cmd = [
        minimap2_exe,
        "-a",
        "-x",
        str(preset),
        "--secondary=yes",
        "-N",
        "50",
        transcript_fa,
        reads_fa,
    ]
    with open(rescue_sam, "wt") as sam_fh:
        subprocess.run(cmd, check=True, stdout=sam_fh)


def _collect_read_sequences(bam_file, contig_acc, region_lend, region_rend, target_read_names):
    remaining = set(target_read_names)
    read_name_to_seq = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:
        if region_lend is not None and region_rend is not None:
            fetch_iter = bam_reader.fetch(contig_acc, max(int(region_lend) - 1, 0), int(region_rend))
        else:
            fetch_iter = bam_reader.fetch(contig_acc)
        for read in fetch_iter:
            read_name = Util_funcs.get_read_name_include_sc_encoding(read)
            if read_name not in remaining:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            if read_name not in read_name_to_seq:
                read_name_to_seq[read_name] = seq
                remaining.discard(read_name)
            if not remaining:
                break
    return read_name_to_seq


def _build_transcript_models(splice_graph, transcripts, contig_seq_str):
    transcript_models = {}
    for transcript in transcripts:
        if not transcript.has_introns():
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
            tx_coord_map.append((tx_pos, tx_pos + seg_len - 1, node_id))
            tx_pos += seg_len

        exon_index_in_path = {
            node_id: idx for idx, node_id in enumerate(path_no_boundaries)
        }
        transcript_models[transcript.get_transcript_id()] = {
            "transcript": transcript,
            "sequence": _build_transcript_sequence(transcript, contig_seq_str),
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
        for transcript_id, model in transcript_models.items():
            print(f">{transcript_id}", file=ofh)
            print(model["sequence"], file=ofh)


def _write_reads_fasta(reads_fa, read_name_to_seq):
    with open(reads_fa, "wt") as ofh:
        for read_name, seq in read_name_to_seq.items():
            print(f">{read_name}", file=ofh)
            print(seq, file=ofh)


def _parse_rescue_alignments(rescue_sam, splice_graph, transcript_models):
    read_to_hits = defaultdict(list)
    min_per_id = float(LRAA_Globals.config.get("rescue_unassigned_min_per_id", 80.0))
    with pysam.AlignmentFile(rescue_sam, "r") as sam_reader:
        for read in sam_reader.fetch(until_eof=True):
            if read.is_unmapped or read.is_supplementary:
                continue
            if read.reference_name not in transcript_models:
                continue
            if any(code == 3 for code, _ in (read.cigartuples or [])):
                continue
            if not _passes_percent_identity(read, min_per_id):
                continue
            merged_segments = Pretty_alignment.read_to_pretty_alignment_segments(read)
            if len(merged_segments) != 1:
                continue
            merged_lend, merged_rend = merged_segments[0]
            projected_path = _project_interval_to_path(
                transcript_models[read.reference_name], merged_lend, merged_rend
            )
            if not projected_path:
                continue
            read_to_hits[read.query_name].append(
                {
                    "score": _alignment_score(read),
                    "path": tuple(projected_path),
                    "transcript_id": read.reference_name,
                }
            )

    rescued_mps = []
    for read_name, hits in read_to_hits.items():
        best_score = max(hit["score"] for hit in hits)
        best_hits = [hit for hit in hits if hit["score"] == best_score]
        projected_paths = {hit["path"] for hit in best_hits}
        if len(projected_paths) != 1:
            continue
        rescued_mps.append(
            MultiPath(
                splice_graph,
                [list(projected_paths.pop())],
                read_types={"PacBio"},
                read_names={read_name},
                read_count=1,
            )
        )

    return rescued_mps


def _passes_percent_identity(read, min_per_id):
    mismatch_count = None
    if read.has_tag("NM"):
        mismatch_count = int(read.get_tag("NM"))
    elif read.has_tag("nM"):
        mismatch_count = int(read.get_tag("nM"))
    if mismatch_count is None:
        return True

    cigar_stats = read.get_cigar_stats()
    aligned_base_count = cigar_stats[0][0]
    if aligned_base_count == 0:
        aligned_base_count = cigar_stats[0][7] + cigar_stats[0][8]
    if aligned_base_count <= 0:
        return False
    per_id = 100.0 - (mismatch_count / aligned_base_count) * 100.0
    return per_id >= min_per_id


def _alignment_score(read):
    if read.has_tag("AS"):
        return int(read.get_tag("AS"))

    mismatch_count = 0
    if read.has_tag("NM"):
        mismatch_count = int(read.get_tag("NM"))
    elif read.has_tag("nM"):
        mismatch_count = int(read.get_tag("nM"))
    cigar_stats = read.get_cigar_stats()
    aligned_base_count = cigar_stats[0][0]
    if aligned_base_count == 0:
        aligned_base_count = cigar_stats[0][7] + cigar_stats[0][8]
    return int(aligned_base_count) - int(mismatch_count)


def _project_interval_to_path(model, tx_lend, tx_rend):
    overlapped_exon_nodes = []
    for seg_lend, seg_rend, node_id in model["tx_coord_map"]:
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
