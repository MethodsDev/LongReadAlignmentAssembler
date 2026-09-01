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


def _extract_one_contig(bam_path, chrom, out_path, threads, chrom_lengths, mapped):
    """One contig's extraction: its own output file, nothing shared with any other.

    Split out of the loop so it can run in a worker. Everything it needs is passed in,
    because the header and the mapped counts are read ONCE by the caller -- re-reading
    the header of a 188 GB bam per worker would spend more than the fan-out saves.
    """

    if mapped == 0:
        LOGGER.info("BAM partition: chromosome %s has no mapped reads; writing empty stub", chrom)
        _write_empty_bam(out_path, [chrom], chrom_lengths)
        return

    LOGGER.info(
        "BAM partition: extracting chromosome %s with %d mapped alignments using %d additional thread(s)",
        chrom,
        mapped,
        threads,
    )
    start_time = time.time()
    try:
        # --no-PG: this writes one output BAM per chromosome, and samtools appends one
        # @PG record per existing chain tip on every write. With the shared
        # bam_for_sg's 2,727,296-record header that added ~302,848 records here.
        # `samtools view -b` also copies the source header into every output, so each
        # per-chromosome BAM inherited that whole chain -- ~1.1 GiB as uncompressed SAM
        # header TEXT, far smaller on disk once BGZF-compressed, but replicated 25
        # times per cluster: a chrY.bam holding only 26,725 alignments measured 87.9 MB
        # on disk, and across the 325 per-chromosome SG outputs roughly 27 GB of the
        # 44 GB total was duplicated header rather than alignment data.
        pysam.view(
            "--no-PG",
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
        return

    elapsed = time.time() - start_time
    rate = mapped / (elapsed / 60.0) if elapsed > 0 else 0.0
    LOGGER.info(
        "BAM partition: wrote %s with %d alignments in %.2f minutes (%.1f alignments/min)",
        chrom,
        mapped,
        elapsed / 60.0,
        rate,
    )


def _plan_bam_partition(
    bam_path: Optional[str],
    chromosomes: List[str],
    out_dir: str,
    label: str,
    samtools_threads: Optional[int] = None,
):
    """Prepare `out_dir` and return the contig extractions that still need running.

    Stubs -- absent bam, contig missing from the header, contig with no mapped reads --
    are written HERE, because they cost nothing and leaving them in the work list would
    make a worker slot wait on a file write. What comes back is only real extraction.

    A planner rather than an executor so both bam kinds can be flattened into ONE pool.
    Giving each kind its own pool of N cannot honour a budget of N: floor division
    strands capacity at odd budgets, and a floor of 1 per kind runs 2 concurrent
    samtools when the budget says 1.
    """

    _clean_output_dir(out_dir)
    if bam_path is None or not os.path.exists(bam_path):
        for chrom in chromosomes:
            _write_empty_bam(os.path.join(out_dir, f"{chrom}.bam"), [chrom])
        return []

    _maybe_index_bam(bam_path)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        references = set(bam.references)
        chrom_lengths = dict(zip(bam.references, bam.lengths))

    mapped_counts = _collect_mapped_counts(bam_path)
    threads = _samtools_threads(samtools_threads)

    work = []
    for chrom in chromosomes:
        out_path = os.path.join(out_dir, f"{chrom}.bam")
        if chrom not in references:
            LOGGER.info(
                "%s partition: chromosome %s missing from BAM header; writing empty stub",
                label,
                chrom,
            )
            _write_empty_bam(out_path, [chrom], chrom_lengths)
            continue
        mapped = mapped_counts.get(chrom, 0)
        if mapped == 0:
            LOGGER.info(
                "%s partition: chromosome %s has no mapped reads; writing empty stub",
                label,
                chrom,
            )
            _write_empty_bam(out_path, [chrom], chrom_lengths)
            continue
        work.append((bam_path, chrom, out_path, threads, chrom_lengths, mapped))

    LOGGER.info(
        "%s partition: %d contig(s) to extract at %d additional samtools thread(s) each",
        label,
        len(work),
        threads,
    )
    return work


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
    parser.add_argument("--bam-for-sg", dest="bam_for_sg", type=str, default=None)
    parser.add_argument("--genome-fasta", dest="genome_fasta", type=str, default=None)
    parser.add_argument("--annot-gtf", dest="annot_gtf", type=str, default=None)
    parser.add_argument("--chromosomes", nargs="*", help="Ordered chromosome list", default=None)
    parser.add_argument("--bam-out-dir", default="split_bams")
    parser.add_argument("--bam-for-sg-out-dir", dest="bam_for_sg_out_dir", default="split_bams_for_sg")
    parser.add_argument("--fasta-out-dir", default="split_fastas")
    parser.add_argument("--gtf-out-dir", default="split_gtfs")
    parser.add_argument(
        "--samtools-threads",
        type=int,
        default=None,
        help="Threads to use for samtools view (default pulls from PARTITION_SAMTOOLS_THREADS or 1)",
    )
    # The contig loop in _partition_bam is the whole cost of this script: FASTA and
    # GTF finish in seconds and bam_for_sg is often absent, so the four-way job pool
    # is one serial pass over the largest bam in the pipeline. MEASURED on a 188 GB
    # library (1.51 B mapped reads, 25 contigs) on a 28-core host: 27 minutes and
    # still running with ONE core busy, on track to exceed an hour before a single
    # shard starts.
    #
    # Contigs are independent -- separate output files, no shared state, and
    # --samtools-threads already sizes each invocation -- so they fan out. Default 1
    # keeps today's behaviour, because this task also runs once per cluster (14 to 32
    # times in a single-cell run) where a wide pool would drain the box, which is why
    # its cpu reservation was reduced to 5 in the first place. The caller that owns
    # the big top-level partition raises this and its reservation together.
    parser.add_argument(
        "--num-workers",
        dest="num_workers",
        type=int,
        default=1,
        help="contigs to extract concurrently (default 1, the historical serial "
        "pass). --samtools-threads is samtools' -@, which is ADDITIONAL threads, so "
        "each worker needs samtools_threads + 1 runnable threads and a caller's cpu "
        "reservation must cover num_workers * (samtools_threads + 1)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = _parse_args(argv)

    input_bam = _normalize_path(args.input_bam)
    bam_for_sg = _normalize_path(args.bam_for_sg)
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

    # FASTA and GTF first, in their own 2-wide pool, because they are single-threaded
    # and short next to the bam work, and running them alongside the bam pool keeps
    # today's overlap. They are the "+ 2" the caller's reservation carries; they are
    # NOT drawn from the bam worker budget, so that budget stays exactly the number of
    # concurrent samtools processes.
    light_jobs = (
        ("FASTA", _partition_fasta, (genome_fasta, chromosomes, args.fasta_out_dir)),
        ("GTF", _partition_gtf, (annot_gtf, chromosomes, args.gtf_out_dir)),
    )

    # Both bam kinds flattened into ONE work list and ONE pool, so `num_workers` is
    # the exact count of concurrent samtools processes -- not a per-kind figure that
    # doubles when bam_for_sg exists, and not a floor division that strands capacity
    # at odd budgets or runs 2 when the budget says 1.
    bam_work = []
    bam_work += _plan_bam_partition(
        input_bam, chromosomes, args.bam_out_dir, "BAM", args.samtools_threads
    )
    bam_work += _plan_bam_partition(
        bam_for_sg, chromosomes, args.bam_for_sg_out_dir, "BAM_FOR_SG", args.samtools_threads
    )

    # Longest first, ACROSS both kinds, by the mapped count already read from each
    # index. Callers pass main_chromosomes read-count descending, which orders each
    # kind correctly on its own -- but concatenating them puts every splice-graph
    # contig behind every primary one, so the largest SG contig would start only after
    # the whole primary bam finished and would then run alone at the tail. Sorting the
    # flattened list is what actually makes the pool's last job a small one.
    bam_work.sort(key=lambda unit: unit[5], reverse=True)

    workers = max(1, min(args.num_workers, len(bam_work))) if bam_work else 1
    threads_each = _samtools_threads(args.samtools_threads)
    LOGGER.info(
        "partition: %d contig extraction(s) across both bam kind(s), %d at a time, "
        "%d additional samtools thread(s) each -> up to %d runnable thread(s) for the "
        "bam work plus 2 for FASTA/GTF",
        len(bam_work),
        workers,
        threads_each,
        workers * (threads_each + 1),
    )

    failures = []
    with ProcessPoolExecutor(max_workers=len(light_jobs)) as light_pool:
        light_futures = {
            light_pool.submit(func, *func_args): name for name, func, func_args in light_jobs
        }

        with ProcessPoolExecutor(max_workers=workers) as bam_pool:
            bam_futures = {}
            for unit in bam_work:
                bam_futures[bam_pool.submit(_extract_one_contig, *unit)] = unit[1]
            for future in as_completed(bam_futures):
                chrom = bam_futures[future]
                try:
                    future.result()
                except Exception as exc:
                    failures.append(("contig {}".format(chrom), exc))

        for future in as_completed(light_futures):
            name = light_futures[future]
            try:
                future.result()
            except Exception as exc:
                failures.append((name, exc))

    if failures:
        what, exc = failures[0]
        raise RuntimeError(
            "{} partition step failed ({} failure(s) total)".format(what, len(failures))
        ) from exc

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # pragma: no cover - surface friendly message
        print(f"partition_data_by_chromosome: {exc}", file=sys.stderr)
        raise
