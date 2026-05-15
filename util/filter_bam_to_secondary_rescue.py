#!/usr/bin/env python3

import argparse
import csv
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass

import pysam


logging.basicConfig(
    format="%(asctime)-15s %(levelname)s %(module)s.%(funcName)s:\n\t%(message)s\n",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class AlignmentMetrics:
    as_score: float
    ms: float
    nm_rate: float
    aligned_query_len: int
    max_intron: int


def _run(cmd):
    logger.info("Running: %s", " ".join(cmd))
    subprocess.check_call(cmd)


def _numeric_tag(read, tag_name):
    if not read.has_tag(tag_name):
        return math.nan
    try:
        value = read.get_tag(tag_name)
    except Exception:
        return math.nan
    return float(value) if isinstance(value, (int, float)) else math.nan


def _alignment_metrics(read):
    aligned_query_len = read.query_alignment_length or 0
    max_intron = 0
    for op, length in read.cigartuples or []:
        if op == 3:
            max_intron = max(max_intron, length)

    as_score = _numeric_tag(read, "AS")
    ms = _numeric_tag(read, "ms")
    nm = _numeric_tag(read, "NM")
    nm_rate = math.nan
    if aligned_query_len > 0 and not math.isnan(nm):
        nm_rate = nm / aligned_query_len

    return AlignmentMetrics(
        as_score=as_score,
        ms=ms,
        nm_rate=nm_rate,
        aligned_query_len=aligned_query_len,
        max_intron=max_intron,
    )


def _metrics_are_complete(metrics):
    if metrics is None:
        return False
    if metrics.aligned_query_len <= 0:
        return False
    return not (
        math.isnan(metrics.as_score)
        or math.isnan(metrics.ms)
        or math.isnan(metrics.nm_rate)
    )


def should_rescue_secondary(primary_metrics, secondary_metrics, secondary_rank_by_as):
    """
    Return (keep_boolean, branch_name) for the secondary-rescue heuristic.

    Branches follow the SIM_BAMs-derived rule:
      delta_nm_rate <= 0
      AND (
        (delta_ms > -0.5 AND rank_by_AS == 1)
        OR
        (delta_ms <= -0.5 AND delta_nm_rate <= -0.005
         AND log2((secondary_max_intron + 1) / (primary_max_intron + 1)) > 0.03)
      )
    """
    if (
        not _metrics_are_complete(primary_metrics)
        or not _metrics_are_complete(secondary_metrics)
        or secondary_rank_by_as is None
    ):
        return False, "missing_metrics"

    delta_ms = secondary_metrics.ms - primary_metrics.ms
    delta_nm_rate = secondary_metrics.nm_rate - primary_metrics.nm_rate
    log2_max_intron_ratio = math.log2(
        (secondary_metrics.max_intron + 1) / (primary_metrics.max_intron + 1)
    )

    if delta_nm_rate > 0:
        return False, "failed_delta_nm_rate"

    if delta_ms > -0.5 and secondary_rank_by_as == 1:
        return True, "branch1"

    if (
        delta_ms <= -0.5
        and delta_nm_rate <= -0.005
        and log2_max_intron_ratio > 0.03
    ):
        return True, "branch2"

    return False, "failed_heuristic"


def _iter_name_groups(bam_reader):
    current_name = None
    current_group = []
    for read in bam_reader.fetch(until_eof=True):
        if current_name is None:
            current_name = read.query_name
            current_group = [read]
            continue
        if read.query_name == current_name:
            current_group.append(read)
            continue
        yield current_name, current_group
        current_name = read.query_name
        current_group = [read]
    if current_name is not None:
        yield current_name, current_group


def _sort_key_for_rank(read):
    metrics = _alignment_metrics(read)
    as_score = metrics.as_score
    if math.isnan(as_score):
        as_score = -math.inf
    return (
        -as_score,
        read.reference_name or "",
        read.reference_start if read.reference_start is not None else -1,
        read.cigarstring or "",
        read.flag,
    )


def _choose_primary_record(non_secondary_records):
    complete = []
    incomplete = []
    for read in non_secondary_records:
        metrics = _alignment_metrics(read)
        row = (read, metrics)
        if _metrics_are_complete(metrics):
            complete.append(row)
        else:
            incomplete.append(row)

    candidates = complete or incomplete
    if not candidates:
        return None, None

    def key(item):
        read, metrics = item
        as_score = metrics.as_score
        if math.isnan(as_score):
            as_score = -math.inf
        ms = metrics.ms
        if math.isnan(ms):
            ms = -math.inf
        return (
            -as_score,
            -ms,
            read.reference_name or "",
            read.reference_start if read.reference_start is not None else -1,
            read.cigarstring or "",
            read.flag,
        )

    read, metrics = sorted(candidates, key=key)[0]
    return read, metrics


def _filter_group(read_name, group, bam_writer, stats):
    stats["read_groups"] += 1
    stats["input_records"] += len(group)

    mapped_non_supp_records = [
        r for r in group if (not r.is_unmapped) and (not r.is_supplementary)
    ]
    non_secondary_records = [r for r in mapped_non_supp_records if not r.is_secondary]
    secondary_records = [r for r in mapped_non_supp_records if r.is_secondary]

    stats["non_secondary_records"] += len(non_secondary_records)
    stats["secondary_candidates"] += len(secondary_records)

    for read in group:
        if read.is_unmapped:
            stats["dropped_unmapped"] += 1
        elif read.is_supplementary:
            stats["dropped_supplementary"] += 1

    for read in non_secondary_records:
        bam_writer.write(read)
        stats["kept_non_secondary"] += 1

    _primary_read, primary_metrics = _choose_primary_record(non_secondary_records)
    if secondary_records and primary_metrics is None:
        stats["groups_with_secondaries_without_primary"] += 1

    ranked_secondaries = sorted(secondary_records, key=_sort_key_for_rank)
    secondary_id_to_rank = {id(read): idx for idx, read in enumerate(ranked_secondaries, 1)}

    for read in secondary_records:
        secondary_metrics = _alignment_metrics(read)
        rank = secondary_id_to_rank.get(id(read))
        keep, branch = should_rescue_secondary(primary_metrics, secondary_metrics, rank)
        if keep:
            bam_writer.write(read)
            stats["kept_secondary_heuristic"] += 1
            stats[f"kept_secondary_{branch}"] += 1
        else:
            stats["dropped_secondary"] += 1
            stats[f"dropped_secondary_{branch}"] += 1


def _write_summary(summary_path, stats):
    secondary_candidates = stats["secondary_candidates"]
    kept_secondary = stats["kept_secondary_heuristic"]
    kept_secondary_frac = (
        kept_secondary / secondary_candidates if secondary_candidates else 0.0
    )
    fieldnames = [
        "read_groups",
        "input_records",
        "non_secondary_records",
        "kept_non_secondary",
        "secondary_candidates",
        "kept_secondary_heuristic",
        "kept_secondary_branch1",
        "kept_secondary_branch2",
        "dropped_secondary",
        "dropped_secondary_missing_metrics",
        "dropped_secondary_failed_delta_nm_rate",
        "dropped_secondary_failed_heuristic",
        "dropped_unmapped",
        "dropped_supplementary",
        "groups_with_secondaries_without_primary",
        "secondary_kept_frac",
    ]
    row = {field: stats.get(field, 0) for field in fieldnames}
    row["secondary_kept_frac"] = f"{kept_secondary_frac:.6f}"
    with open(summary_path, "wt", newline="") as ofh:
        writer = csv.DictWriter(ofh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter BAM by read name to retain all mapped non-supplementary primary "
            "alignments and only secondary alignments passing the LRAA secondary-rescue heuristic."
        )
    )
    parser.add_argument("--input_bam", required=True, help="coordinate-sorted input BAM")
    parser.add_argument(
        "--output_bam", required=True, help="coordinate-sorted filtered output BAM"
    )
    parser.add_argument(
        "--summary_tsv",
        default=None,
        help="optional summary TSV path (default: <output_bam>.summary.tsv)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="threads for samtools name-sort / coordinate-sort / index",
    )
    parser.add_argument(
        "--workdir",
        default=None,
        help="working directory for temporary BAMs (default: temporary directory)",
    )
    args = parser.parse_args()

    threads = max(1, int(args.threads))
    input_bam = os.path.abspath(args.input_bam)
    output_bam = os.path.abspath(args.output_bam)
    summary_tsv = (
        os.path.abspath(args.summary_tsv)
        if args.summary_tsv
        else os.path.abspath(output_bam + ".summary.tsv")
    )
    out_dir = os.path.dirname(output_bam) or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)

    if args.workdir:
        workdir = os.path.abspath(args.workdir)
        os.makedirs(workdir, exist_ok=True)
        cleanup_workdir = False
    else:
        workdir = tempfile.mkdtemp(prefix="lraa_secondary_rescue.")
        cleanup_workdir = True

    try:
        name_sorted_bam = os.path.join(workdir, "name_sorted.bam")
        filtered_unsorted_bam = os.path.join(workdir, "filtered.unsorted.bam")
        sort_prefix = os.path.join(workdir, "samtools_tmp")

        _run(
            [
                "samtools",
                "sort",
                "-n",
                "-@",
                str(threads),
                "-T",
                sort_prefix + ".name",
                "-o",
                name_sorted_bam,
                input_bam,
            ]
        )

        stats = {
            "read_groups": 0,
            "input_records": 0,
            "non_secondary_records": 0,
            "kept_non_secondary": 0,
            "secondary_candidates": 0,
            "kept_secondary_heuristic": 0,
            "kept_secondary_branch1": 0,
            "kept_secondary_branch2": 0,
            "dropped_secondary": 0,
            "dropped_secondary_missing_metrics": 0,
            "dropped_secondary_failed_delta_nm_rate": 0,
            "dropped_secondary_failed_heuristic": 0,
            "dropped_unmapped": 0,
            "dropped_supplementary": 0,
            "groups_with_secondaries_without_primary": 0,
        }

        with pysam.AlignmentFile(name_sorted_bam, "rb") as bam_in, pysam.AlignmentFile(
            filtered_unsorted_bam, "wb", template=bam_in
        ) as bam_out:
            for read_name, group in _iter_name_groups(bam_in):
                _filter_group(read_name, group, bam_out, stats)

        _run(
            [
                "samtools",
                "sort",
                "-@",
                str(threads),
                "-T",
                sort_prefix + ".coord",
                "-o",
                output_bam,
                filtered_unsorted_bam,
            ]
        )
        _run(["samtools", "index", "-@", str(threads), output_bam])

        _write_summary(summary_tsv, stats)
        secondary_kept_frac = (
            stats["kept_secondary_heuristic"] / stats["secondary_candidates"]
            if stats["secondary_candidates"] > 0
            else 0.0
        )
        logger.info(
            "Filtered BAM complete: groups=%d input_records=%d kept_non_secondary=%d secondary_candidates=%d kept_secondary_heuristic=%d kept_secondary_branch1=%d kept_secondary_branch2=%d dropped_secondary=%d dropped_secondary_missing_metrics=%d dropped_secondary_failed_delta_nm_rate=%d dropped_secondary_failed_heuristic=%d dropped_supplementary=%d secondary_kept_frac=%.4f summary=%s",
            stats["read_groups"],
            stats["input_records"],
            stats["kept_non_secondary"],
            stats["secondary_candidates"],
            stats["kept_secondary_heuristic"],
            stats["kept_secondary_branch1"],
            stats["kept_secondary_branch2"],
            stats["dropped_secondary"],
            stats["dropped_secondary_missing_metrics"],
            stats["dropped_secondary_failed_delta_nm_rate"],
            stats["dropped_secondary_failed_heuristic"],
            stats["dropped_supplementary"],
            secondary_kept_frac,
            summary_tsv,
        )

    finally:
        if cleanup_workdir:
            shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as exc:
        logger.error("Command failed with exit code %s", exc.returncode)
        sys.exit(exc.returncode)
