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
    """samtools' -@, which is ADDITIONAL threads on top of the main one.

    An explicit 0 is honoured rather than raised to 1: it is samtools' own default and
    the only value that fits one extraction plus the two light jobs inside a 3-core
    grant. Absent any value the default is 1, unchanged.
    """

    if cli_value is not None:
        return max(0, cli_value)

    env_value = os.environ.get("PARTITION_SAMTOOLS_THREADS")
    if env_value:
        try:
            return max(0, int(env_value))
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


def _cgroup_cpu_quota(root: str = "/sys/fs/cgroup") -> Optional[int]:
    """Cores this process was GRANTED by its cgroup, or None if no quota is set.

    This is the number a container runtime actually enforces, which is what makes it
    worth reading: a WDL `cpu` declaration is a request, and the backend may grant
    less. miniwdl's swarm backend applies it as `--limit-cpu`, Cromwell sets a per-VM
    quota, and both surface here.

    ``root`` is a parameter so this is testable without a container.

    v2 puts "<quota> <period>" in cpu.max ("max" meaning unlimited); v1 splits it
    across cpu.cfs_quota_us (-1 unlimited) and cpu.cfs_period_us. Rounded DOWN to whole
    cores, then floored at 1, because a fractional grant still has to run something.
    """

    v2 = os.path.join(root, "cpu.max")
    try:
        with open(v2, "rt") as fh:
            quota_s, period_s = fh.read().split()[:2]
        if quota_s != "max":
            period = int(period_s)
            if period > 0:
                return max(1, int(quota_s) // period)
        return None
    except (OSError, ValueError):
        pass

    try:
        with open(os.path.join(root, "cpu", "cpu.cfs_quota_us"), "rt") as fh:
            quota = int(fh.read().strip())
        with open(os.path.join(root, "cpu", "cpu.cfs_period_us"), "rt") as fh:
            period = int(fh.read().strip())
        if quota > 0 and period > 0:
            return max(1, quota // period)
    except (OSError, ValueError):
        pass

    return None


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
    # The cores the CALLER reserved, so the pool can hold itself inside them. This is
    # the contract; affinity is only a secondary net, because a backend may constrain
    # scheduling shares without constraining a cpuset: miniwdl reports "cpu adjusted to
    # host limit" and hands the task fewer cores, yet the task can still SEE every
    # core, so sched_getaffinity would not bind there. Passing the reservation makes
    # the cap independent of whether the backend enforces one.
    parser.add_argument(
        "--reserved-cpu",
        dest="reserved_cpu",
        type=int,
        default=None,
        help="cores the caller reserved for this task. The pool is held so that "
        "num_workers * (samtools_threads + 1) + 2 stays within it, so the script "
        "cannot oversubscribe its own reservation even where nothing enforces it",
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

    # The cores this process actually HAS, resolved before any work is planned, because
    # the effective -@ is derived from it and every planned unit carries that value into
    # its samtools invocation. Deciding it afterwards would log one number and run
    # another.
    #
    # The CGROUP QUOTA is the only number the runtime both grants and enforces. It is
    # how a cpu reservation is actually applied under docker: miniwdl's swarm backend
    # passes --limit-cpu, Cromwell/Terra sets a quota per VM, and both land as cpu.max
    # (v2) or cpu.cfs_quota_us (v1) inside the task. That makes it the right primary
    # signal, because it reflects what was GRANTED rather than what was asked for --
    # the case where miniwdl reports "cpu adjusted to host limit" and hands the task
    # fewer cores than the WDL requested. MEASURED inside `docker run --cpus=8`:
    # sched_getaffinity still reports all 16 host cores, the quota reports 8.
    #
    # AFFINITY catches the cpuset form of the same grant (taskset, Apptainer with one).
    try:
        visible = len(os.sched_getaffinity(0))
    except AttributeError:  # pragma: no cover - non-Linux
        visible = os.cpu_count() or 1
    granted = _cgroup_cpu_quota()
    binding_grant = min([c for c in (granted, visible) if c])

    requested_threads = _samtools_threads(args.samtools_threads)
    threads_each = requested_threads

    # A GRANT below the floor is a fact about the environment, not a caller mistake, so
    # it is adapted to rather than refused: failing here would abort a run that can
    # still complete, just slower. -@ is the one lever with no correctness meaning, so
    # it gives way first. Capping the POOL alone does not cover this -- one worker at
    # -@ 4 is 5 threads, plus the 2 light jobs, under an enforced 3-core quota.
    if binding_grant < requested_threads + 3:
        threads_each = max(0, binding_grant - 3)
        LOGGER.warning(
            "partition: --samtools-threads reduced from %d to %d -- only %d core(s) "
            "granted and one worker plus the 2 FASTA/GTF jobs needs %d at the "
            "requested value. Throughput only; -@ is ADDITIONAL threads",
            requested_threads,
            threads_each,
            binding_grant,
            requested_threads + 3,
        )
        if binding_grant < 3:
            LOGGER.warning(
                "partition: %d core(s) cannot hold one extraction plus the 2 FASTA/GTF "
                "jobs even at -@ 0; running 3 process(es) anyway, which is the minimum "
                "this script has",
                binding_grant,
            )

    # Both bam kinds flattened into ONE work list and ONE pool, so `num_workers` is
    # the exact count of concurrent samtools processes -- not a per-kind figure that
    # doubles when bam_for_sg exists, and not a floor division that strands capacity
    # at odd budgets or runs 2 when the budget says 1.
    bam_work = []
    # threads_each, not args.samtools_threads: the value resolved against the grant
    # above is the one each unit must carry into its samtools invocation.
    bam_work += _plan_bam_partition(
        input_bam, chromosomes, args.bam_out_dir, "BAM", threads_each
    )
    bam_work += _plan_bam_partition(
        bam_for_sg, chromosomes, args.bam_for_sg_out_dir, "BAM_FOR_SG", threads_each
    )

    # Longest first, ACROSS both kinds, by the mapped count already read from each
    # index. Callers pass main_chromosomes read-count descending, which orders each
    # kind correctly on its own -- but concatenating them puts every splice-graph
    # contig behind every primary one, so the largest SG contig would start only after
    # the whole primary bam finished and would then run alone at the tail. Sorting the
    # flattened list is what actually makes the pool's last job a small one.
    bam_work.sort(key=lambda unit: unit[5], reverse=True)

    # THREE caps on the POOL, because no one of them is sufficient. The grant numbers
    # come from above; the RESERVATION is self-consistency only -- it holds the pool
    # inside what the caller PROMISED the scheduler even where nothing enforces it, and
    # cannot detect a grant smaller than the request, which is what the other two are
    # for.

    def _affordable(cores):
        # minus the 2 light FASTA/GTF jobs, divided by the runnable threads a worker
        # needs (-@ is ADDITIONAL, so + 1 for its main thread)
        return max(1, (cores - 2) // (threads_each + 1))

    caps = {}
    if granted:
        caps["cgroup quota({} core)".format(granted)] = _affordable(granted)
    caps["affinity({} core)".format(visible)] = _affordable(visible)
    if args.reserved_cpu is not None and args.reserved_cpu > 0:
        caps["reservation({} core)".format(args.reserved_cpu)] = _affordable(args.reserved_cpu)

    # A reservation smaller than the floor cannot be honoured by capping: one worker
    # is 5 runnable threads at -@ 4 and the two light jobs are 2 more, so 7 is the
    # minimum this script can run inside. Below that, max(1, ...) would floor the pool
    # at one worker and STILL exceed the reservation -- the same violation the cap
    # exists to prevent, just quieter. Refused rather than silently overrun, because
    # the caller has two real fixes (raise the reservation, or lower
    # --samtools-threads) and the script cannot choose between them.
    minimum_cpu = threads_each + 3
    if args.reserved_cpu is not None and 0 < args.reserved_cpu < minimum_cpu:
        raise ValueError(
            "--reserved-cpu {} cannot run this task: one worker needs {} runnable "
            "thread(s) at --samtools-threads {} and the FASTA/GTF jobs need 2 more, "
            "so the floor is {}. Raise the reservation or lower "
            "--samtools-threads.".format(
                args.reserved_cpu,
                threads_each + 1,
                threads_each,
                minimum_cpu,
            )
        )

    requested = max(1, args.num_workers)
    ceiling = min(caps.values())
    workers = max(1, min(requested, len(bam_work), ceiling)) if bam_work else 1

    if bam_work and workers < min(requested, len(bam_work)):
        binding = [k for k, v in caps.items() if v == ceiling]
        LOGGER.warning(
            "partition: --num-workers %d reduced to %d by %s -- each worker needs %d "
            "runnable thread(s), so the requested pool would oversubscribe",
            requested,
            workers,
            " and ".join(sorted(binding)),
            threads_each + 1,
        )

    LOGGER.info(
        "partition: %d contig extraction(s) across both bam kind(s), %d at a time, "
        "%d additional samtools thread(s) each -> up to %d runnable thread(s) for the "
        "bam work plus 2 for FASTA/GTF, against %d visible core(s)",
        len(bam_work),
        workers,
        threads_each,
        workers * (threads_each + 1),
        visible,
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
