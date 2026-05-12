#!/usr/bin/env python3

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile

import pysam


logging.basicConfig(
    format="%(asctime)-15s %(levelname)s %(module)s.%(funcName)s:\n\t%(message)s\n",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


def _run(cmd):
    logger.info("Running: %s", " ".join(cmd))
    subprocess.check_call(cmd)


def _get_alignment_score(read):
    if read.has_tag("AS"):
        try:
            return int(read.get_tag("AS"))
        except Exception:
            return None
    return None


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


def _filter_group(read_name, group, bam_writer, stats):
    stats["read_groups"] += 1
    stats["input_records"] += len(group)

    mapped_non_supp_records = [
        r for r in group if (not r.is_unmapped) and (not r.is_supplementary)
    ]
    scored_records = [
        (r, _get_alignment_score(r))
        for r in mapped_non_supp_records
    ]
    scored_records = [(r, score) for (r, score) in scored_records if score is not None]

    best_as = None
    if scored_records:
        best_as = max(score for (_, score) in scored_records)
    elif mapped_non_supp_records:
        stats["groups_without_scored_alignment"] += 1

    for read in group:
        if read.is_unmapped:
            stats["dropped_unmapped"] += 1
            continue

        if read.is_supplementary:
            stats["dropped_supplementary"] += 1
            continue

        if read.is_secondary:
            stats["secondary_records"] += 1
        else:
            stats["non_secondary_records"] += 1

        aln_as = _get_alignment_score(read)
        if best_as is None:
            stats["dropped_no_best_as"] += 1
            continue
        if aln_as is None:
            stats["dropped_no_as"] += 1
            continue

        if aln_as == best_as:
            bam_writer.write(read)
            stats["kept_best_as_tied"] += 1
            if read.is_secondary:
                stats["kept_secondary_best_as_tied"] += 1
            else:
                stats["kept_non_secondary_best_as_tied"] += 1
        else:
            stats["dropped_lower_as"] += 1


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter BAM by read name to retain only mapped non-supplementary alignments "
            "whose AS score ties the best AS observed for that read. Primary/secondary "
            "flags are ignored for selection."
        )
    )
    parser.add_argument("--input_bam", required=True, help="coordinate-sorted input BAM")
    parser.add_argument("--output_bam", required=True, help="coordinate-sorted filtered output BAM")
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
    out_dir = os.path.dirname(output_bam) or os.getcwd()
    os.makedirs(out_dir, exist_ok=True)

    if args.workdir:
        workdir = os.path.abspath(args.workdir)
        os.makedirs(workdir, exist_ok=True)
        cleanup_workdir = False
    else:
        workdir = tempfile.mkdtemp(prefix="lraa_tied_secondaries.")
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
            "secondary_records": 0,
            "kept_best_as_tied": 0,
            "kept_non_secondary_best_as_tied": 0,
            "kept_secondary_best_as_tied": 0,
            "dropped_lower_as": 0,
            "dropped_no_best_as": 0,
            "dropped_no_as": 0,
            "dropped_unmapped": 0,
            "dropped_supplementary": 0,
            "groups_without_scored_alignment": 0,
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

        kept_total = stats["kept_best_as_tied"]
        secondary_kept_frac = (
            stats["kept_secondary_best_as_tied"] / stats["secondary_records"]
            if stats["secondary_records"] > 0
            else 0.0
        )
        logger.info(
            "Filtered BAM complete: groups=%d input_records=%d non_secondary_records=%d secondary_records=%d kept_total=%d kept_non_secondary_best_as_tied=%d kept_secondary_best_as_tied=%d dropped_lower_as=%d dropped_no_best_as=%d dropped_no_as=%d dropped_supplementary=%d secondary_kept_frac=%.4f",
            stats["read_groups"],
            stats["input_records"],
            stats["non_secondary_records"],
            stats["secondary_records"],
            kept_total,
            stats["kept_non_secondary_best_as_tied"],
            stats["kept_secondary_best_as_tied"],
            stats["dropped_lower_as"],
            stats["dropped_no_best_as"],
            stats["dropped_no_as"],
            stats["dropped_supplementary"],
            secondary_kept_frac,
        )
        if stats["groups_without_scored_alignment"]:
            logger.warning(
                "Pathological read groups observed: groups_without_scored_alignment=%d",
                stats["groups_without_scored_alignment"],
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
