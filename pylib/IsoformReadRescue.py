#!/usr/bin/env python3

import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
import hashlib

import pysam
from intervaltree import IntervalTree

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
    read_path_mapper=None,
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

    return build_transcriptome_alignment_multipaths(
        splice_graph,
        transcripts,
        contig_seq_str,
        bam_file,
        contig_acc,
        region_lend,
        region_rend,
        read_names=read_names,
        read_path_mapper=read_path_mapper,
        target_strand=splice_graph.get_contig_strand(),
        include_monoexonic=False,
        require_unique_path_across_best_hits=True,
        log_label="transcriptome rescue",
    )


def build_transcriptome_alignment_multipaths(
    splice_graph,
    transcripts,
    contig_seq_str,
    bam_file,
    contig_acc,
    region_lend,
    region_rend,
    read_names=None,
    read_path_mapper=None,
    target_strand=None,
    include_monoexonic=False,
    require_unique_path_across_best_hits=True,
    split_multipaths_by_gene=False,
    genome_target_gating=False,
    log_label="transcriptome alignment",
    return_read_details=False,
):
    minimap2_exe = shutil.which("minimap2")
    if minimap2_exe is None:
        logger.warning("minimap2 not found in PATH; skipping %s", log_label)
        return []

    transcript_models = _build_transcript_models(
        splice_graph,
        transcripts,
        contig_seq_str,
        include_monoexonic=include_monoexonic,
    )
    if not transcript_models:
        if return_read_details:
            return [], {
                "read_name_to_multipaths": {},
                "read_name_to_best_score": {},
                "read_name_to_primary_score": {},
                "requested_read_names": set(),
            }
        return []

    candidate_rows = None
    read_name_to_allowed_target_ids = None
    read_name_to_primary_score = {}
    read_name_to_primary_per_id = {}
    if genome_target_gating:
        candidate_tsv = None
        (
            read_name_to_seq,
            read_name_to_allowed_target_ids,
            gating_stats,
            candidate_rows,
            read_name_to_primary_score,
            read_name_to_primary_per_id,
        ) = (
            _collect_genome_gated_read_targets(
                bam_file,
                contig_acc,
                region_lend,
                region_rend,
                transcript_models,
                target_read_names=read_names,
                target_strand=target_strand,
                include_secondary_candidates=LRAA_Globals.config.get(
                    "allow_secondary_alignments", False
                ),
            )
        )
        primary_considered = gating_stats["primary_considered"]
        primary_retained = gating_stats["primary_retained_for_fastq"]
        retain_frac = (
            primary_retained / primary_considered if primary_considered > 0 else 0.0
        )
        avg_targets_per_read = (
            gating_stats["candidate_target_links"] / gating_stats["candidate_reads"]
            if gating_stats["candidate_reads"] > 0
            else 0.0
        )
        logger.info(
            "[%s%s] %s genome gating: primary_retained=%d/%d (%.3f), candidate_reads=%d, candidate_alignments(primary=%d secondary=%d), avg_targets_per_read=%.2f",
            splice_graph.get_contig_acc(),
            splice_graph.get_contig_strand(),
            log_label,
            primary_retained,
            primary_considered,
            retain_frac,
            gating_stats["candidate_reads"],
            gating_stats["candidate_primary_alignments"],
            gating_stats["candidate_secondary_alignments"],
            avg_targets_per_read,
        )
    else:
        read_name_to_seq = _collect_read_sequences(
            bam_file, contig_acc, region_lend, region_rend, read_names, target_strand
        )

    if not read_name_to_seq:
        logger.info(
            "[%s%s] %s skipped: no reads with retrievable sequences",
            splice_graph.get_contig_acc(),
            splice_graph.get_contig_strand(),
            log_label,
        )
        if return_read_details:
            return [], {
                "read_name_to_multipaths": {},
                "read_name_to_best_score": {},
                "read_name_to_primary_score": read_name_to_primary_score,
                "read_name_to_best_per_id": {},
                "read_name_to_primary_per_id": read_name_to_primary_per_id,
                "requested_read_names": set(),
            }
        return []

    tmp_dir = tempfile.mkdtemp(
        prefix=f"tx_rescue.{contig_acc}.{splice_graph.get_contig_strand()}.",
        dir=os.getcwd(),
    )
    keep_tmp = bool(LRAA_Globals.DEBUG) or bool(
        LRAA_Globals.config.get("no_cleanup", False)
    )
    try:
        transcript_fa = os.path.join(tmp_dir, "transcripts.fa")
        reads_fa = os.path.join(tmp_dir, "reads.fa")
        rescue_sam = os.path.join(tmp_dir, "rescue.sam")
        candidate_tsv = os.path.join(tmp_dir, "genome_target_candidates.tsv")

        _write_transcript_fasta(transcript_fa, transcript_models)
        _write_reads_fasta(reads_fa, read_name_to_seq)
        if candidate_rows is not None:
            _write_candidate_tsv(
                candidate_tsv,
                candidate_rows,
                retained_read_names=set(read_name_to_seq.keys()),
            )
        _run_minimap2_transcriptome_alignment(
            transcript_fa, reads_fa, rescue_sam, minimap2_exe
        )

        rescued_mps, read_details = _parse_rescue_alignments(
            rescue_sam,
            splice_graph,
            transcript_models,
            read_path_mapper=read_path_mapper,
            require_unique_path_across_best_hits=require_unique_path_across_best_hits,
            split_multipaths_by_gene=split_multipaths_by_gene,
            read_name_to_allowed_target_ids=read_name_to_allowed_target_ids,
        )

        if return_read_details:
            read_details["read_name_to_primary_score"] = {
                _normalize_read_identifier(read_name): score
                for read_name, score in read_name_to_primary_score.items()
                if _normalize_read_identifier(read_name) is not None
            }
            read_details["read_name_to_primary_per_id"] = {
                _normalize_read_identifier(read_name): per_id
                for read_name, per_id in read_name_to_primary_per_id.items()
                if _normalize_read_identifier(read_name) is not None
            }
            read_details["requested_read_names"] = {
                _normalize_read_identifier(read_name)
                for read_name in read_name_to_seq.keys()
                if _normalize_read_identifier(read_name) is not None
            }

        logger.info(
            "[%s%s] %s: requested=%d rescued=%d",
            splice_graph.get_contig_acc(),
            splice_graph.get_contig_strand(),
            log_label,
            len(read_name_to_seq),
            len(rescued_mps),
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
        minimap_threads = max(
            1, int(LRAA_Globals.config.get("num_threads_per_worker", 1))
        )
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
):
    remaining = None if target_read_names is None else set(target_read_names)
    read_name_to_seq = {}
    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:
        if region_lend is not None and region_rend is not None:
            fetch_iter = bam_reader.fetch(contig_acc, max(int(region_lend) - 1, 0), int(region_rend))
        else:
            fetch_iter = bam_reader.fetch(contig_acc)
        for read in fetch_iter:
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
                if remaining is not None:
                    remaining.discard(read_name)
            if remaining is not None and not remaining:
                break
    return read_name_to_seq


def _collect_genome_gated_read_targets(
    bam_file,
    contig_acc,
    region_lend,
    region_rend,
    transcript_models,
    target_read_names=None,
    target_strand=None,
    include_secondary_candidates=False,
):
    exon_overlap_index = _build_exon_overlap_index(transcript_models)
    remaining = None if target_read_names is None else set(target_read_names)
    read_name_to_seq = {}
    read_name_to_allowed_target_ids = defaultdict(set)
    read_name_to_primary_score = {}
    read_name_to_primary_per_id = {}
    candidate_rows = []
    min_per_id = float(_resolve_rescue_min_per_id())

    stats = {
        "primary_considered": 0,
        "primary_retained_for_fastq": 0,
        "candidate_reads": 0,
        "candidate_target_links": 0,
        "candidate_primary_alignments": 0,
        "candidate_secondary_alignments": 0,
    }

    with pysam.AlignmentFile(bam_file, "rb") as bam_reader:
        if region_lend is not None and region_rend is not None:
            fetch_iter = bam_reader.fetch(
                contig_acc, max(int(region_lend) - 1, 0), int(region_rend)
            )
        else:
            fetch_iter = bam_reader.fetch(contig_acc)

        for read in fetch_iter:
            if target_strand is not None:
                if read.is_forward and target_strand != "+":
                    continue
                if read.is_reverse and target_strand != "-":
                    continue

            if read.is_unmapped or read.is_supplementary:
                continue
            if read.is_secondary and not include_secondary_candidates:
                continue
            if read.is_paired and not read.is_proper_pair:
                continue
            if read.is_duplicate or read.is_qcfail:
                continue
            if read.mapping_quality < int(LRAA_Globals.config["min_mapping_quality"]):
                continue
            passes_per_id = _passes_percent_identity(read, min_per_id)
            if read.is_secondary and not passes_per_id:
                continue
            if not read.is_secondary and not passes_per_id:
                # Keep low-perID primary alignments available for transcriptome
                # arbitration/rescue when the genome-derived path is unusable.
                # Secondary candidate generation remains stricter to limit noise.
                pass
            elif not passes_per_id:
                continue

            read_name = Util_funcs.get_read_name_include_sc_encoding(read)
            if remaining is not None and read_name not in remaining:
                continue

            if not read.is_secondary:
                stats["primary_considered"] += 1

            target_id_to_overlap_bp = _get_alignment_overlapping_targets(
                read, exon_overlap_index
            )
            if not target_id_to_overlap_bp:
                continue

            if read.is_secondary:
                stats["candidate_secondary_alignments"] += 1
            else:
                stats["candidate_primary_alignments"] += 1

            if not read.is_secondary:
                seq = read.query_sequence
                if seq and read_name not in read_name_to_seq:
                    if read.is_reverse:
                        seq = _reverse_complement(seq)
                    read_name_to_seq[read_name] = seq
                    stats["primary_retained_for_fastq"] += 1
                if remaining is not None:
                    remaining.discard(read_name)
                read_name_to_primary_score[read_name] = _alignment_score(read)
                read_name_to_primary_per_id[read_name] = _alignment_per_id_fraction(read)

            for target_id, overlap_bp in target_id_to_overlap_bp.items():
                read_name_to_allowed_target_ids[read_name].add(target_id)
                model = transcript_models[target_id]
                candidate_rows.append(
                    {
                        "read_name": read_name,
                        "target_id": target_id,
                        "gene_id": model["gene_id"],
                        "transcript_id": model["transcript_id"],
                        "is_primary": int(not read.is_secondary),
                        "is_secondary": int(read.is_secondary),
                        "mapq": int(read.mapping_quality),
                        "genomic_lend": int(read.reference_start) + 1,
                        "genomic_rend": int(read.reference_end),
                        "overlap_bp": int(overlap_bp),
                    }
                )

            if remaining is not None and not remaining:
                # Do not break here: secondary alignments for previously seen primary reads
                # may still appear later in coordinate-sorted BAMs and should contribute candidates.
                pass

    retained_read_names = set(read_name_to_seq.keys())
    read_name_to_allowed_target_ids = {
        read_name: target_ids
        for read_name, target_ids in read_name_to_allowed_target_ids.items()
        if read_name in retained_read_names
    }
    stats["candidate_reads"] = len(read_name_to_allowed_target_ids)
    stats["candidate_target_links"] = sum(
        len(target_ids) for target_ids in read_name_to_allowed_target_ids.values()
    )

    return (
        read_name_to_seq,
        read_name_to_allowed_target_ids,
        stats,
        candidate_rows,
        read_name_to_primary_score,
        read_name_to_primary_per_id,
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


def _write_candidate_tsv(candidate_tsv, candidate_rows, retained_read_names=None):
    retained = None if retained_read_names is None else set(retained_read_names)
    with open(candidate_tsv, "wt") as ofh:
        print(
            "\t".join(
                [
                    "read_name",
                    "target_id",
                    "gene_id",
                    "transcript_id",
                    "is_primary",
                    "is_secondary",
                    "mapq",
                    "genomic_lend",
                    "genomic_rend",
                    "overlap_bp",
                ]
            ),
            file=ofh,
        )
        for row in candidate_rows:
            if retained is not None and row["read_name"] not in retained:
                continue
            print(
                "\t".join(
                    [
                        str(row["read_name"]),
                        str(row["target_id"]),
                        str(row["gene_id"]),
                        str(row["transcript_id"]),
                        str(row["is_primary"]),
                        str(row["is_secondary"]),
                        str(row["mapq"]),
                        str(row["genomic_lend"]),
                        str(row["genomic_rend"]),
                        str(row["overlap_bp"]),
                    ]
                ),
                file=ofh,
            )


def _parse_rescue_alignments(
    rescue_sam,
    splice_graph,
    transcript_models,
    read_path_mapper=None,
    require_unique_path_across_best_hits=True,
    split_multipaths_by_gene=False,
    read_name_to_allowed_target_ids=None,
):
    read_to_hits = defaultdict(list)
    min_per_id = float(_resolve_rescue_min_per_id())
    with pysam.AlignmentFile(rescue_sam, "r") as sam_reader:
        for read in sam_reader.fetch(until_eof=True):
            if read.is_unmapped or read.is_supplementary:
                continue
            if read.reference_name not in transcript_models:
                continue
            if (
                read_name_to_allowed_target_ids is not None
                and (
                    read.query_name not in read_name_to_allowed_target_ids
                    or read.reference_name
                    not in read_name_to_allowed_target_ids[read.query_name]
                )
            ):
                continue
            if any(code == 3 for code, _ in (read.cigartuples or [])):
                continue
            if not _passes_percent_identity(read, min_per_id):
                continue
            merged_segments = Pretty_alignment.read_to_pretty_alignment_segments(read)
            if len(merged_segments) != 1:
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
                continue
            read_to_hits[read.query_name].append(
                {
                    "score": _alignment_score(read),
                    "per_id": _alignment_per_id_fraction(read),
                    "path": tuple(projected_path),
                    "target_id": read.reference_name,
                    "transcript_id": transcript_models[read.reference_name][
                        "transcript_id"
                    ],
                    "gene_id": transcript_models[read.reference_name]["gene_id"],
                }
            )

    rescued_mps = []
    read_name_to_multipaths = defaultdict(list)
    read_name_to_best_score = {}
    read_name_to_best_per_id = {}
    for read_name, hits in read_to_hits.items():
        read_key = _normalize_read_identifier(read_name)
        if read_key is None:
            continue
        best_score = max(hit["score"] for hit in hits)
        best_hits = [hit for hit in hits if hit["score"] == best_score]
        read_name_to_best_score[read_key] = best_score
        read_name_to_best_per_id[read_key] = max(
            hit.get("per_id", 0.0) for hit in best_hits
        )
        if split_multipaths_by_gene:
            gene_to_hits = defaultdict(list)
            for hit in best_hits:
                gene_to_hits[hit["gene_id"]].append(hit)
            for gene_hits in gene_to_hits.values():
                projected_paths = {hit["path"] for hit in gene_hits}
                if require_unique_path_across_best_hits and len(projected_paths) != 1:
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
        else:
            projected_paths = {hit["path"] for hit in best_hits}
            if require_unique_path_across_best_hits and len(projected_paths) != 1:
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

    return rescued_mps, {
        "read_name_to_multipaths": dict(read_name_to_multipaths),
        "read_name_to_best_score": read_name_to_best_score,
        "read_name_to_best_per_id": read_name_to_best_per_id,
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
    per_id_fraction = _alignment_per_id_fraction(read)
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


def _alignment_per_id_fraction(read):
    mismatch_count = None
    if read.has_tag("NM"):
        mismatch_count = int(read.get_tag("NM"))
    elif read.has_tag("nM"):
        mismatch_count = int(read.get_tag("nM"))
    if mismatch_count is None:
        return None

    cigar_stats = read.get_cigar_stats()
    aligned_base_count = cigar_stats[0][0]
    if aligned_base_count == 0:
        aligned_base_count = cigar_stats[0][7] + cigar_stats[0][8]
    if aligned_base_count <= 0:
        return 0.0

    return max(0.0, float(aligned_base_count - mismatch_count) / float(aligned_base_count))


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
        genomic_segments = _project_interval_to_genomic_segments(model, tx_lend, tx_rend)
        if not genomic_segments:
            return None
        return read_path_mapper(
            genomic_segments,
            refine_TSS_simple_path=True,
            refine_PolyA_simple_path=True,
            snap_nearby_boundary_features=True,
            left_soft_clipping=left_soft_clipping,
            right_soft_clipping=right_soft_clipping,
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
