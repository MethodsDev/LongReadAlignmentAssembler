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

    primary_records = [
        r for r in group if (not r.is_secondary) and (not r.is_unmapped) and (not r.is_supplementary)
    ]
    primary_as = None
    if primary_records:
        primary_scores = [_get_alignment_score(r) for r in primary_records]
        primary_scores = [x for x in primary_scores if x is not None]
        if primary_scores:
            primary_as = max(primary_scores)
        if len(primary_records) > 1:
            stats["multi_primary_groups"] += 1
    else:
        stats["groups_without_primary"] += 1

    for read in group:
        if read.is_unmapped:
            stats["dropped_unmapped"] += 1
            continue

        if not read.is_secondary:
            bam_writer.write(read)
            stats["kept_non_secondary"] += 1
            continue

        stats["secondary_records"] += 1
        sec_as = _get_alignment_score(read)
        if primary_as is None:
            stats["dropped_secondary_no_primary"] += 1
            continue
        if sec_as is None:
            stats["dropped_secondary_no_as"] += 1
            continue
        if sec_as == primary_as:
            bam_writer.write(read)
            stats["kept_secondary_tied"] += 1
        else:
            stats["dropped_secondary_untied"] += 1


def main():
    parser = argparse.ArgumentParser(
        description="Filter BAM to retain all non-secondary records plus only those secondary records whose AS score ties the primary alignment score for the same read."
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
            "secondary_records": 0,
            "kept_non_secondary": 0,
            "kept_secondary_tied": 0,
            "dropped_secondary_untied": 0,
            "dropped_secondary_no_primary": 0,
            "dropped_secondary_no_as": 0,
            "dropped_unmapped": 0,
            "groups_without_primary": 0,
            "multi_primary_groups": 0,
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

        kept_total = stats["kept_non_secondary"] + stats["kept_secondary_tied"]
        secondary_kept_frac = (
            stats["kept_secondary_tied"] / stats["secondary_records"]
            if stats["secondary_records"] > 0
            else 0.0
        )
        logger.info(
            "Filtered BAM complete: groups=%d input_records=%d kept_total=%d kept_non_secondary=%d kept_secondary_tied=%d dropped_secondary_untied=%d dropped_secondary_no_primary=%d dropped_secondary_no_as=%d secondary_kept_frac=%.4f",
            stats["read_groups"],
            stats["input_records"],
            kept_total,
            stats["kept_non_secondary"],
            stats["kept_secondary_tied"],
            stats["dropped_secondary_untied"],
            stats["dropped_secondary_no_primary"],
            stats["dropped_secondary_no_as"],
            secondary_kept_frac,
        )
        if stats["groups_without_primary"] or stats["multi_primary_groups"]:
            logger.warning(
                "Pathological read groups observed: groups_without_primary=%d multi_primary_groups=%d",
                stats["groups_without_primary"],
                stats["multi_primary_groups"],
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
