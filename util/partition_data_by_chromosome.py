#!/usr/bin/env python3
"""Split BAM/FASTA/GTF inputs into per-chromosome files in a single pass."""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, Iterable, List, Optional

import pysam


if not logging.getLogger().handlers:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(processName)s] %(levelname)s: %(message)s")

LOGGER = logging.getLogger(__name__)


def _samtools_threads(cli_value: Optional[int] = None) -> int:
    if cli_value is not None:
        return max(1, cli_value)

    env_value = os.environ.get("PARTITION_SAMTOOLS_THREADS")
    if env_value:
        try:
            return max(1, int(env_value))
        except ValueError:
            LOGGER.warning(
                "Invalid PARTITION_SAMTOOLS_THREADS value %r; falling back to 1 thread",
                env_value,
            )
    return 1


def _write_empty_bam(
    path: str, chromosomes: List[str], lengths: Optional[Dict[str, int]] = None
) -> None:
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {
                "SN": chrom,
                "LN": max(1, int(lengths.get(chrom, 1))) if lengths else 1,
            }
            for chrom in chromosomes
        ],
    }
    with pysam.AlignmentFile(path, "wb", header=header):
        pass


def _normalize_path(path: Optional[str]) -> Optional[str]:
    if path is None:
        return None
    trimmed = path.strip()
    if not trimmed or trimmed.lower() in {"none", "null"}:
        return None
    return trimmed


def _unique_ordered(items: Iterable[str]) -> List[str]:
    seen = set()
    ordered: List[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            ordered.append(item)
    return ordered


def _clean_output_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)
    for entry in os.scandir(path):
        if entry.is_file():
            os.remove(entry.path)


def _maybe_index_bam(bam_path: str) -> None:
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM not found: {bam_path}")
    bai_candidates = [bam_path + ".bai", os.path.splitext(bam_path)[0] + ".bai"]
    if any(os.path.exists(p) for p in bai_candidates):
        return
    pysam.index(bam_path)


def _ensure_fasta_index(fasta_path: str) -> None:
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA not found: {fasta_path}")
    fai_path = fasta_path + ".fai"
    if os.path.exists(fai_path):
        return
    pysam.faidx(fasta_path)


def _collect_mapped_counts(bam_path: str) -> Dict[str, int]:
    stats_output = pysam.idxstats(bam_path)
    counts: Dict[str, int] = {}
    for line in stats_output.strip().splitlines():
        fields = line.split("\t")
        if len(fields) != 4:
            continue
        chrom, _length, mapped, _unmapped = fields
        if chrom == "*":
            continue
        try:
            counts[chrom] = int(mapped)
        except ValueError:
            counts[chrom] = 0
    return counts


def _collect_chromosomes_from_bam(bam_path: str) -> List[str]:
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        return list(bam.references)


def _collect_chromosomes_from_fasta(fasta_path: str) -> List[str]:
    _ensure_fasta_index(fasta_path)
    with pysam.FastaFile(fasta_path) as fasta:
        return list(fasta.references)


def _collect_chromosomes_from_gtf(gtf_path: str) -> List[str]:
    chroms: List[str] = []
    open_fn = gzip.open if gtf_path.endswith(".gz") else open
    with open_fn(gtf_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0].strip()
            if chrom:
                chroms.append(chrom)
    return chroms


def _partition_bam(
    bam_path: Optional[str],
    chromosomes: List[str],
    out_dir: str,
    samtools_threads: Optional[int] = None,
) -> None:
    _clean_output_dir(out_dir)
    if bam_path is None or not os.path.exists(bam_path):
        for chrom in chromosomes:
            _write_empty_bam(os.path.join(out_dir, f"{chrom}.bam"), [chrom])
        return

    _maybe_index_bam(bam_path)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        references = set(bam.references)
        chrom_lengths = dict(zip(bam.references, bam.lengths))

    mapped_counts = _collect_mapped_counts(bam_path)
    threads = _samtools_threads(samtools_threads)

    for chrom in chromosomes:
        out_path = os.path.join(out_dir, f"{chrom}.bam")
        if chrom not in references:
            LOGGER.info("BAM partition: chromosome %s missing from BAM header; writing empty stub", chrom)
            _write_empty_bam(out_path, [chrom], chrom_lengths)
            continue

        mapped = mapped_counts.get(chrom, 0)
        if mapped == 0:
            LOGGER.info("BAM partition: chromosome %s has no mapped reads; writing empty stub", chrom)
            _write_empty_bam(out_path, [chrom], chrom_lengths)
            continue

        LOGGER.info(
            "BAM partition: extracting chromosome %s with %d mapped alignments using %d thread(s)",
            chrom,
            mapped,
            threads,
        )
        start_time = time.time()
        try:
            pysam.view(
                "-@",
                str(threads),
                "-h",
                "-b",
                "-o",
                out_path,
                bam_path,
                chrom,
                catch_stdout=False,
            )
        except pysam.SamtoolsError as exc:  # pragma: no cover - htslib surface error
            LOGGER.warning(
                "BAM partition: samtools view failed for %s (%s); writing empty stub",
                chrom,
                exc,
            )
            _write_empty_bam(out_path, [chrom], chrom_lengths)
            continue

        elapsed = time.time() - start_time
        rate = mapped / (elapsed / 60.0) if elapsed > 0 else 0.0
        LOGGER.info(
            "BAM partition: wrote %s with %d alignments in %.2f minutes (%.1f alignments/min)",
            chrom,
            mapped,
            elapsed / 60.0,
            rate,
        )


def _write_fasta_record(handle, chrom: str, sequence: str) -> None:
    handle.write(f">{chrom}\n")
    for idx in range(0, len(sequence), 60):
        handle.write(sequence[idx : idx + 60] + "\n")


def _partition_fasta(
    fasta_path: Optional[str],
    chromosomes: List[str],
    out_dir: str,
) -> None:
    _clean_output_dir(out_dir)
    if fasta_path is None or not os.path.exists(fasta_path):
        for chrom in chromosomes:
            with open(os.path.join(out_dir, f"{chrom}.genome.fasta"), "wt") as handle:
                handle.write(f">{chrom}\n")
        return

    _ensure_fasta_index(fasta_path)

    outputs = {
        chrom: open(os.path.join(out_dir, f"{chrom}.genome.fasta"), "wt")
        for chrom in chromosomes
    }
    try:
        seen: Dict[str, bool] = {chrom: False for chrom in chromosomes}
        logged_chroms: set[str] = set()

        # Stream through the FASTA once to avoid repeated random-access fetches.
        with pysam.FastxFile(fasta_path) as fasta_handle:
            for entry in fasta_handle:
                chrom_name = entry.name.split()[0]
                if chrom_name not in outputs or seen.get(chrom_name):
                    continue
                if chrom_name not in logged_chroms:
                    LOGGER.info("FASTA partition: encountered chromosome %s", chrom_name)
                    logged_chroms.add(chrom_name)
                _write_fasta_record(outputs[chrom_name], chrom_name, entry.sequence)
                seen[chrom_name] = True

        for chrom, handle in outputs.items():
            if not seen[chrom]:
                handle.write(f">{chrom}\n")
    finally:
        for handle in outputs.values():
            handle.close()


def _partition_gtf(
    gtf_path: Optional[str],
    chromosomes: List[str],
    out_dir: str,
) -> None:
    _clean_output_dir(out_dir)

    outputs = {
        chrom: open(os.path.join(out_dir, f"{chrom}.annot.gtf"), "wt")
        for chrom in chromosomes
    }
    has_records: Dict[str, bool] = {chrom: False for chrom in chromosomes}
    logged_chroms: set[str] = set()

    if gtf_path is not None:
        if not os.path.exists(gtf_path):
            raise FileNotFoundError(f"GTF not found: {gtf_path}")
        open_fn = gzip.open if gtf_path.endswith(".gz") else open
        with open_fn(gtf_path, "rt") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                chrom = line.split("\t", 1)[0].strip()
                sink = outputs.get(chrom)
                if sink is not None:
                    if chrom not in logged_chroms:
                        LOGGER.info("GTF partition: encountered chromosome %s", chrom)
                        logged_chroms.add(chrom)
                    sink.write(line)
                    has_records[chrom] = True

    for chrom, handle in outputs.items():
        if not has_records[chrom]:
            handle.write("# no gtf records\n")
        handle.close()


def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Partition BAM/FASTA/GTF by chromosome")
    parser.add_argument("--input-bam", dest="input_bam", type=str, default=None)
    parser.add_argument("--genome-fasta", dest="genome_fasta", type=str, default=None)
    parser.add_argument("--annot-gtf", dest="annot_gtf", type=str, default=None)
    parser.add_argument("--chromosomes", nargs="*", help="Ordered chromosome list", default=None)
    parser.add_argument("--bam-out-dir", default="split_bams")
    parser.add_argument("--fasta-out-dir", default="split_fastas")
    parser.add_argument("--gtf-out-dir", default="split_gtfs")
    parser.add_argument(
        "--samtools-threads",
        type=int,
        default=None,
        help="Threads to use for samtools view (default pulls from PARTITION_SAMTOOLS_THREADS or 1)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = _parse_args(argv)

    input_bam = _normalize_path(args.input_bam)
    genome_fasta = _normalize_path(args.genome_fasta)
    annot_gtf = _normalize_path(args.annot_gtf)

    chromosome_inputs: List[str] = []
    if args.chromosomes:
        for entry in args.chromosomes:
            chromosome_inputs.extend(part for part in entry.replace(",", " ").split() if part)

    chromosomes: List[str]
    if chromosome_inputs:
        chromosomes = _unique_ordered(chromosome_inputs)
    else:
        collected: List[str] = []
        if input_bam:
            collected.extend(_collect_chromosomes_from_bam(input_bam))
        if genome_fasta:
            collected.extend(_collect_chromosomes_from_fasta(genome_fasta))
        if annot_gtf:
            collected.extend(_collect_chromosomes_from_gtf(annot_gtf))
        chromosomes = _unique_ordered(collected)

    if not chromosomes:
        raise ValueError("No chromosomes supplied or detected.")

    partition_jobs = (
        ("BAM", _partition_bam, (input_bam, chromosomes, args.bam_out_dir, args.samtools_threads)),
        ("FASTA", _partition_fasta, (genome_fasta, chromosomes, args.fasta_out_dir)),
        ("GTF", _partition_gtf, (annot_gtf, chromosomes, args.gtf_out_dir)),
    )

    with ProcessPoolExecutor(max_workers=len(partition_jobs)) as executor:
        future_to_name = {
            executor.submit(func, *func_args): job_name for job_name, func, func_args in partition_jobs
        }
        for future in as_completed(future_to_name):
            job_name = future_to_name[future]
            try:
                future.result()
            except Exception as exc:
                raise RuntimeError(f"{job_name} partition step failed") from exc

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - surface friendly message
        print(f"partition_data_by_chromosome: {exc}", file=sys.stderr)
        raise
