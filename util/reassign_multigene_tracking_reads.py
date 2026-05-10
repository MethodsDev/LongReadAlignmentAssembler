#!/usr/bin/env python3

import argparse
import csv
import gzip
import logging
import math
import shutil
from collections import defaultdict, deque
from statistics import median


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


def open_maybe_gzip(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Reallocate cross-gene multimapping reads using a constrained EM while "
            "keeping non-cross-gene tracking fractions fixed."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--quant_expr", required=True, help="input merged LRAA quant.expr file"
    )
    parser.add_argument(
        "--tracking",
        required=True,
        help="input merged LRAA quant.tracking[.gz] file",
    )
    parser.add_argument(
        "--output_expr", required=True, help="output corrected quant.expr file"
    )
    parser.add_argument(
        "--output_tracking",
        required=True,
        help="output corrected quant.tracking[.gz] file",
    )
    parser.add_argument(
        "--max_iter",
        type=int,
        default=250,
        help="maximum EM iterations for each affected gene component",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=1e-9,
        help="L1 convergence tolerance on transcript counts within each component",
    )
    parser.add_argument(
        "--min_expr",
        type=float,
        default=1e-8,
        help="minimum abundance floor used during EM updates",
    )
    return parser.parse_args()


def parse_float(value, default=0.0):
    if value is None or value == "":
        return default
    return float(value)


def format_count(value):
    return f"{value:.1f}"


def format_prob(value):
    return f"{value:.3f}"


def tx_key_from_row(row):
    return (row["gene_id"], row["transcript_id"])


def read_expr_rows(path):
    with open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"Error, no header found in quant expr file: {path}")
        rows = list(reader)
        fieldnames = list(reader.fieldnames)
    return rows, fieldnames


def read_tracking_rows(path):
    with open_maybe_gzip(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"Error, no header found in tracking file: {path}")
        rows = list(reader)
        fieldnames = list(reader.fieldnames)
    return rows, fieldnames


def derive_rpm_scale(expr_rows):
    ratios = []
    for row in expr_rows:
        all_reads = parse_float(row.get("all_reads"), 0.0)
        rpm_total_reads = parse_float(row.get("RPM_total_reads"), 0.0)
        if all_reads > 0 and math.isfinite(rpm_total_reads):
            ratios.append(rpm_total_reads / all_reads)
    if not ratios:
        return 0.0
    return median(ratios)


def build_components(affected_reads, read_to_genes):
    gene_to_reads = defaultdict(set)
    for read_name in affected_reads:
        for gene_id in read_to_genes[read_name]:
            gene_to_reads[gene_id].add(read_name)

    components = []
    remaining_genes = set(gene_to_reads)
    while remaining_genes:
        seed_gene = next(iter(remaining_genes))
        genes = set()
        reads = set()
        queue = deque([("gene", seed_gene)])
        while queue:
            node_type, node_id = queue.popleft()
            if node_type == "gene":
                if node_id in genes:
                    continue
                genes.add(node_id)
                remaining_genes.discard(node_id)
                for read_name in gene_to_reads[node_id]:
                    if read_name not in reads:
                        queue.append(("read", read_name))
            else:
                if node_id in reads:
                    continue
                reads.add(node_id)
                for gene_id in read_to_genes[node_id]:
                    if gene_id not in genes:
                        queue.append(("gene", gene_id))
        components.append((genes, reads))

    return components


def run_component_em(
    component_reads,
    read_candidates,
    component_transcripts,
    fixed_counts,
    init_counts,
    max_iter,
    tol,
    min_expr,
):
    counts = {
        tx_key: max(init_counts.get(tx_key, 0.0), fixed_counts.get(tx_key, 0.0), min_expr)
        for tx_key in component_transcripts
    }
    read_to_fracs = {}

    for _ in range(max_iter):
        ambiguous_counts = defaultdict(float)
        new_read_to_fracs = {}

        for read_name in component_reads:
            candidates = read_candidates[read_name]
            denom = 0.0
            numerators = {}
            for tx_key, weight in candidates.items():
                abundance = max(counts.get(tx_key, 0.0), min_expr)
                numer = max(weight, 0.0) * abundance
                numerators[tx_key] = numer
                denom += numer

            if denom <= 0:
                fallback_total = sum(max(weight, 0.0) for weight in candidates.values())
                if fallback_total <= 0:
                    fallback_total = float(len(candidates))
                    numerators = {tx_key: 1.0 for tx_key in candidates}
                else:
                    numerators = {
                        tx_key: max(weight, 0.0) for tx_key, weight in candidates.items()
                    }
                denom = fallback_total

            read_fracs = {}
            for tx_key, numer in numerators.items():
                frac = numer / denom if denom > 0 else 0.0
                read_fracs[tx_key] = frac
                ambiguous_counts[tx_key] += frac

            new_read_to_fracs[read_name] = read_fracs

        new_counts = {}
        for tx_key in component_transcripts:
            new_counts[tx_key] = fixed_counts.get(tx_key, 0.0) + ambiguous_counts.get(
                tx_key, 0.0
            )

        delta = sum(abs(new_counts[tx_key] - counts[tx_key]) for tx_key in component_transcripts)
        counts = new_counts
        read_to_fracs = new_read_to_fracs

        if delta <= tol:
            break

    return counts, read_to_fracs


def write_tracking(path, fieldnames, rows):
    with open_maybe_gzip(path, "wt") as ofh:
        writer = csv.DictWriter(
            ofh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_expr(path, fieldnames, rows):
    with open(path, "wt", newline="") as ofh:
        writer = csv.DictWriter(
            ofh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()

    expr_rows, expr_fieldnames = read_expr_rows(args.quant_expr)
    tracking_rows, tracking_fieldnames = read_tracking_rows(args.tracking)

    expr_by_tx = {}
    gene_to_transcripts = defaultdict(list)
    for row in expr_rows:
        tx_key = tx_key_from_row(row)
        if tx_key in expr_by_tx:
            raise RuntimeError(f"Duplicate transcript row found in quant expr: {tx_key}")
        expr_by_tx[tx_key] = row
        gene_to_transcripts[row["gene_id"]].append(tx_key)

    read_to_indices = defaultdict(list)
    read_to_genes = defaultdict(set)
    for idx, row in enumerate(tracking_rows):
        read_name = row["read_name"]
        read_to_indices[read_name].append(idx)
        read_to_genes[read_name].add(row["gene_id"])

    affected_reads = {
        read_name for read_name, gene_ids in read_to_genes.items() if len(gene_ids) > 1
    }
    logger.info(
        "Loaded %d expr rows, %d tracking rows, and %d cross-gene reads",
        len(expr_rows),
        len(tracking_rows),
        len(affected_reads),
    )

    if not affected_reads:
        shutil.copyfile(args.quant_expr, args.output_expr)
        if args.tracking != args.output_tracking:
            shutil.copyfile(args.tracking, args.output_tracking)
        logger.info("No cross-gene reads found; copied inputs to outputs unchanged.")
        return

    fixed_counts = defaultdict(float)
    for row in tracking_rows:
        if row["read_name"] in affected_reads:
            continue
        fixed_counts[tx_key_from_row(row)] += parse_float(row.get("frac_assigned"), 0.0)

    expr_init_counts = {
        tx_key_from_row(row): parse_float(row.get("all_reads"), 0.0) for row in expr_rows
    }

    components = build_components(affected_reads, read_to_genes)
    logger.info("Found %d affected gene/read components", len(components))

    for component_index, (component_genes, component_reads) in enumerate(components, start=1):
        component_transcripts = []
        for gene_id in sorted(component_genes):
            component_transcripts.extend(gene_to_transcripts.get(gene_id, []))

        read_candidates = {}
        for read_name in component_reads:
            candidates = {}
            for idx in read_to_indices[read_name]:
                row = tracking_rows[idx]
                tx_key = tx_key_from_row(row)
                if tx_key in candidates:
                    raise RuntimeError(
                        f"Duplicate transcript candidate encountered for read {read_name}: {tx_key}"
                    )
                candidates[tx_key] = parse_float(row.get("read_weight"), 1.0)
            read_candidates[read_name] = candidates

        logger.info(
            "Component %d: %d genes, %d transcripts, %d cross-gene reads",
            component_index,
            len(component_genes),
            len(component_transcripts),
            len(component_reads),
        )

        _, read_to_fracs = run_component_em(
            component_reads,
            read_candidates,
            component_transcripts,
            fixed_counts,
            expr_init_counts,
            args.max_iter,
            args.tol,
            args.min_expr,
        )

        for read_name in component_reads:
            read_fracs = read_to_fracs[read_name]
            for idx in read_to_indices[read_name]:
                tx_key = tx_key_from_row(tracking_rows[idx])
                tracking_rows[idx]["frac_assigned"] = format_prob(read_fracs[tx_key])

    updated_counts = defaultdict(float)
    updated_uniq_reads = defaultdict(int)
    for row in tracking_rows:
        tx_key = tx_key_from_row(row)
        frac = parse_float(row.get("frac_assigned"), 0.0)
        updated_counts[tx_key] += frac
        if frac >= 0.9995:
            updated_uniq_reads[tx_key] += 1

    gene_counts = defaultdict(float)
    total_reported_read_count = 0.0
    for tx_key, count in updated_counts.items():
        gene_counts[tx_key[0]] += count
        total_reported_read_count += count

    rpm_scale = derive_rpm_scale(expr_rows)
    logger.info(
        "Updated total reported read count = %.3f; RPM_total_reads scale factor = %.6f",
        total_reported_read_count,
        rpm_scale,
    )

    for row in expr_rows:
        tx_key = tx_key_from_row(row)
        gene_id = row["gene_id"]
        all_reads = updated_counts.get(tx_key, 0.0)
        gene_total = gene_counts.get(gene_id, 0.0)
        uniq_reads = updated_uniq_reads.get(tx_key, 0)
        isoform_fraction = all_reads / gene_total if gene_total > 0 else 0.0
        unique_gene_read_fraction = uniq_reads / gene_total if gene_total > 0 else 0.0
        tpm = all_reads / total_reported_read_count * 1e6 if total_reported_read_count > 0 else 0.0
        rpm_total_reads = all_reads * rpm_scale if rpm_scale > 0 else 0.0

        if "uniq_reads" in row:
            row["uniq_reads"] = str(int(uniq_reads))
        if "all_reads" in row:
            row["all_reads"] = format_count(all_reads)
        if "isoform_fraction" in row:
            row["isoform_fraction"] = format_prob(isoform_fraction)
        if "unique_gene_read_fraction" in row:
            row["unique_gene_read_fraction"] = format_prob(unique_gene_read_fraction)
        if "TPM" in row:
            row["TPM"] = format_prob(tpm)
        if "RPM_total_reads" in row:
            row["RPM_total_reads"] = format_prob(rpm_total_reads)

    write_expr(args.output_expr, expr_fieldnames, expr_rows)
    write_tracking(args.output_tracking, tracking_fieldnames, tracking_rows)
    logger.info("Wrote corrected expr to %s", args.output_expr)
    logger.info("Wrote corrected tracking to %s", args.output_tracking)


if __name__ == "__main__":
    main()
