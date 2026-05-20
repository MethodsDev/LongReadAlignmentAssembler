#!/usr/bin/env python3

import argparse
import csv
import gzip
import logging
import math
import os
import shutil
import subprocess
import tempfile
from collections import defaultdict
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


def iter_non_comment_lines(fh):
    for line in fh:
        if line.startswith("#"):
            continue
        yield line


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Requantify connected cross-gene multimapping components while keeping "
            "unaffected tracking fractions fixed."
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
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.01,
        help="base alpha regularization applied to ambiguous component reads",
    )
    parser.add_argument(
        "--tmp_dir",
        default=None,
        help="temporary directory for external sort intermediates",
    )
    parser.add_argument(
        "--sort_buffer",
        default="4G",
        help="buffer size passed to GNU sort -S for tracking-file sorting",
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
        reader = csv.DictReader(iter_non_comment_lines(fh), delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError(f"Error, no header found in quant expr file: {path}")
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


def run_component_em(
    component_reads,
    read_candidates,
    component_transcripts,
    init_counts,
    max_iter,
    tol,
    min_expr,
    base_alpha,
):
    counts = {
        tx_key: max(init_counts.get(tx_key, 0.0), min_expr)
        for tx_key in component_transcripts
    }
    read_to_fracs = {}
    transcript_alphas = defaultdict(float)
    for read_name in component_reads:
        candidates = read_candidates[read_name]
        if len(candidates) > 1:
            for tx_key in candidates:
                transcript_alphas[tx_key] += base_alpha

    for _ in range(max_iter):
        assigned_counts = defaultdict(float)
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
                assigned_counts[tx_key] += frac

            new_read_to_fracs[read_name] = read_fracs

        new_expression_counts = {}
        for tx_key in component_transcripts:
            new_expression_counts[tx_key] = (
                assigned_counts.get(tx_key, 0.0)
                + transcript_alphas.get(tx_key, 0.0)
                + min_expr
            )

        delta = sum(
            abs(new_expression_counts[tx_key] - counts[tx_key])
            for tx_key in component_transcripts
        )
        counts = new_expression_counts
        read_to_fracs = new_read_to_fracs

        if delta <= tol:
            break

    final_counts = defaultdict(float)
    for read_fracs in read_to_fracs.values():
        for tx_key, frac in read_fracs.items():
            final_counts[tx_key] += frac

    return final_counts, read_to_fracs


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


class DisjointSet:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def add(self, item):
        if item not in self.parent:
            self.parent[item] = item
            self.rank[item] = 0

    def find(self, item):
        self.add(item)
        if self.parent[item] != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, left, right):
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.rank[left_root] < self.rank[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        if self.rank[left_root] == self.rank[right_root]:
            self.rank[left_root] += 1


def read_tracking_fieldnames(path):
    with open_maybe_gzip(path, "rt") as fh:
        header = ""
        for line in fh:
            if line.startswith("#"):
                continue
            header = line
            break
    if not header:
        raise RuntimeError(f"Error, no header found in tracking file: {path}")
    return header.rstrip("\r\n").split("\t")


def _sort_env():
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    return env


def sort_tracking_body(input_path, output_path, fieldnames, sort_field, tmp_dir, sort_buffer):
    sort_col = fieldnames.index(sort_field) + 1
    cmd = [
        "sort",
        "-t",
        "\t",
        "-S",
        sort_buffer,
        "-T",
        tmp_dir,
        "-k",
        f"{sort_col},{sort_col}",
    ]
    logger.info("Sorting tracking rows by %s: %s", sort_field, " ".join(cmd))
    with open(output_path, "wt", newline="") as ofh:
        proc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=ofh,
            text=True,
            env=_sort_env(),
        )
        assert proc.stdin is not None
        with open_maybe_gzip(input_path, "rt") as ifh:
            header = ""
            for line in ifh:
                if line.startswith("#"):
                    continue
                header = line
                break
            if not header:
                raise RuntimeError(f"Error, no header found in tracking file: {input_path}")
            for line in ifh:
                if line.startswith("#"):
                    continue
                proc.stdin.write(line)
        proc.stdin.close()
        returncode = proc.wait()
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, cmd)


def sort_tsv_body(input_path, output_path, key_specs, tmp_dir, sort_buffer):
    cmd = ["sort", "-t", "\t", "-S", sort_buffer, "-T", tmp_dir]
    for key_spec in key_specs:
        cmd.extend(["-k", key_spec])
    logger.info("Sorting intermediate TSV: %s", " ".join(cmd))
    with open(input_path, "rt", newline="") as ifh, open(output_path, "wt", newline="") as ofh:
        subprocess.check_call(cmd, stdin=ifh, stdout=ofh, env=_sort_env())


def iter_tracking_groups(sorted_tracking_path, fieldnames):
    with open(sorted_tracking_path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, fieldnames=fieldnames, delimiter="\t")
        current_read = None
        group = []
        for row in reader:
            read_name = row["read_name"]
            if current_read is None:
                current_read = read_name
            if read_name != current_read:
                yield current_read, group
                current_read = read_name
                group = []
            group.append(row)
        if current_read is not None:
            yield current_read, group


def dedupe_tracking_group(rows):
    key_to_row = {}
    duplicate_rows = 0
    for row in rows:
        key = (row["read_name"], row["gene_id"], row["transcript_id"])
        if key not in key_to_row:
            key_to_row[key] = dict(row)
            continue

        duplicate_rows += 1
        kept_row = key_to_row[key]
        if "read_weight" in row and "read_weight" in kept_row:
            kept_row["read_weight"] = format_prob(
                max(
                    parse_float(kept_row.get("read_weight"), 1.0),
                    parse_float(row.get("read_weight"), 1.0),
                )
            )
    return list(key_to_row.values()), duplicate_rows


def scan_read_sorted_tracking(read_sorted_path, fieldnames, cross_read_genes_path):
    dsu = DisjointSet()
    total_rows = 0
    kept_rows = 0
    duplicate_rows = 0
    affected_read_count = 0

    with open(cross_read_genes_path, "wt", newline="") as cross_fh:
        writer = csv.writer(cross_fh, delimiter="\t", lineterminator="\n")
        for read_name, group in iter_tracking_groups(read_sorted_path, fieldnames):
            total_rows += len(group)
            deduped_rows, group_duplicate_rows = dedupe_tracking_group(group)
            duplicate_rows += group_duplicate_rows
            kept_rows += len(deduped_rows)

            genes = sorted({row["gene_id"] for row in deduped_rows})
            if len(genes) > 1:
                affected_read_count += 1
                first_gene = genes[0]
                dsu.add(first_gene)
                for gene_id in genes[1:]:
                    dsu.union(first_gene, gene_id)
                for gene_id in genes:
                    writer.writerow([read_name, gene_id])
                continue

    if duplicate_rows:
        logger.info(
            "Collapsed %d duplicate tracking rows while scanning %d rows",
            duplicate_rows,
            total_rows,
        )

    return {
        "fieldnames": fieldnames,
        "total_rows": total_rows,
        "kept_rows": kept_rows,
        "affected_read_count": affected_read_count,
        "dsu": dsu,
    }


def write_read_components(cross_read_genes_path, read_components_path, dsu):
    root_to_component_id = {}
    component_genes = defaultdict(set)
    component_count = 0
    affected_reads = 0

    with open(cross_read_genes_path, "rt", newline="") as ifh, open(
        read_components_path, "wt", newline=""
    ) as ofh:
        reader = csv.reader(ifh, delimiter="\t")
        writer = csv.writer(ofh, delimiter="\t", lineterminator="\n")
        current_read = None
        genes = []
        for row in reader:
            read_name, gene_id = row
            if current_read is None:
                current_read = read_name
            if read_name != current_read:
                root = dsu.find(genes[0])
                if root not in root_to_component_id:
                    component_count += 1
                    root_to_component_id[root] = component_count
                component_id = root_to_component_id[root]
                writer.writerow([current_read, component_id])
                component_genes[component_id].update(genes)
                affected_reads += 1
                current_read = read_name
                genes = []
            genes.append(gene_id)

        if current_read is not None:
            root = dsu.find(genes[0])
            if root not in root_to_component_id:
                component_count += 1
                root_to_component_id[root] = component_count
            component_id = root_to_component_id[root]
            writer.writerow([current_read, component_id])
            component_genes[component_id].update(genes)
            affected_reads += 1

    logger.info(
        "Found %d affected gene/read components across %d cross-gene reads",
        component_count,
        affected_reads,
    )
    return component_genes


def write_component_candidate_rows(
    read_sorted_path,
    fieldnames,
    gene_to_component,
    candidates_path,
):
    candidate_rows = 0
    component_read_count = 0
    fixed_counts = defaultdict(float)
    fixed_uniq_reads = defaultdict(int)

    with open(candidates_path, "wt", newline="") as ofh:
        writer = csv.writer(ofh, delimiter="\t", lineterminator="\n")
        for read_name, group in iter_tracking_groups(read_sorted_path, fieldnames):
            deduped_rows, _ = dedupe_tracking_group(group)
            component_ids = sorted(
                {
                    gene_to_component[row["gene_id"]]
                    for row in deduped_rows
                    if row["gene_id"] in gene_to_component
                }
            )

            if not component_ids:
                for row in deduped_rows:
                    tx_key = tx_key_from_row(row)
                    frac = parse_float(row.get("frac_assigned"), 0.0)
                    fixed_counts[tx_key] += frac
                    if frac >= 0.9995:
                        fixed_uniq_reads[tx_key] += 1
                continue

            if len(component_ids) > 1:
                raise RuntimeError(
                    f"Read {read_name} spans multiple cross-gene EM components: {component_ids}"
                )

            component_id = component_ids[0]
            component_read_count += 1
            for row in deduped_rows:
                if row["gene_id"] not in gene_to_component:
                    raise RuntimeError(
                        f"Read {read_name} has candidate gene {row['gene_id']} outside affected component {component_id}"
                    )
                writer.writerow(
                    [
                        component_id,
                        read_name,
                        row["gene_id"],
                        row["transcript_id"],
                        row.get("read_weight", "1.0"),
                    ]
                )
                candidate_rows += 1

    logger.info(
        "Wrote %d component read/transcript candidate rows across %d component reads",
        candidate_rows,
        component_read_count,
    )
    return fixed_counts, fixed_uniq_reads


def iter_candidate_components(sorted_candidates_path):
    with open(sorted_candidates_path, "rt", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        current_component_id = None
        read_candidates = {}
        for row in reader:
            component_id = int(row[0])
            read_name = row[1]
            tx_key = (row[2], row[3])
            read_weight = parse_float(row[4], 1.0)
            if current_component_id is None:
                current_component_id = component_id
            if component_id != current_component_id:
                yield current_component_id, read_candidates
                current_component_id = component_id
                read_candidates = {}
            read_candidates.setdefault(read_name, {})[tx_key] = read_weight
        if current_component_id is not None:
            yield current_component_id, read_candidates


def run_component_ems_to_fraction_file(
    sorted_candidates_path,
    fractions_path,
    component_genes,
    gene_to_transcripts,
    fixed_counts,
    fixed_uniq_reads,
    expr_init_counts,
    max_iter,
    tol,
    min_expr,
    base_alpha,
):
    updated_counts = defaultdict(float)
    updated_counts.update(fixed_counts)
    updated_uniq_reads = defaultdict(int)
    updated_uniq_reads.update(fixed_uniq_reads)

    with open(fractions_path, "wt", newline="") as ofh:
        writer = csv.writer(ofh, delimiter="\t", lineterminator="\n")
        for component_index, (component_id, read_candidates) in enumerate(
            iter_candidate_components(sorted_candidates_path), start=1
        ):
            component_reads = list(read_candidates)
            component_transcripts = []
            for gene_id in sorted(component_genes.get(component_id, [])):
                component_transcripts.extend(gene_to_transcripts.get(gene_id, []))

            logger.info(
                "Component %d: %d genes, %d transcripts, %d component reads",
                component_index,
                len(component_genes.get(component_id, [])),
                len(component_transcripts),
                len(component_reads),
            )

            _, read_to_fracs = run_component_em(
                component_reads,
                read_candidates,
                component_transcripts,
                expr_init_counts,
                max_iter,
                tol,
                min_expr,
                base_alpha,
            )

            for read_name, read_fracs in read_to_fracs.items():
                for tx_key, frac in read_fracs.items():
                    weight = read_candidates[read_name][tx_key]
                    writer.writerow(
                        [
                            read_name,
                            tx_key[0],
                            tx_key[1],
                            format_prob(frac),
                            format_prob(weight),
                        ]
                    )
                    updated_counts[tx_key] += frac
                    if frac >= 0.9995:
                        updated_uniq_reads[tx_key] += 1

    return updated_counts, updated_uniq_reads


def iter_fraction_groups(sorted_fractions_path):
    with open(sorted_fractions_path, "rt", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        current_read = None
        frac_map = {}
        for row in reader:
            read_name = row[0]
            tx_key = (row[1], row[2])
            if current_read is None:
                current_read = read_name
            if read_name != current_read:
                yield current_read, frac_map
                current_read = read_name
                frac_map = {}
            frac_map[tx_key] = (row[3], row[4])
        if current_read is not None:
            yield current_read, frac_map


def write_final_tracking_from_read_sorted(
    read_sorted_path,
    output_path,
    fieldnames,
    sorted_fractions_path,
):
    fraction_iter = iter_fraction_groups(sorted_fractions_path)
    try:
        fraction_read_name, fraction_map = next(fraction_iter)
    except StopIteration:
        fraction_read_name = None
        fraction_map = {}

    with open_maybe_gzip(output_path, "wt") as ofh:
        writer = csv.DictWriter(
            ofh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()

        for read_name, group in iter_tracking_groups(read_sorted_path, fieldnames):
            while fraction_read_name is not None and fraction_read_name < read_name:
                try:
                    fraction_read_name, fraction_map = next(fraction_iter)
                except StopIteration:
                    fraction_read_name = None
                    fraction_map = {}
                    break

            deduped_rows, _ = dedupe_tracking_group(group)
            active_fraction_map = fraction_map if fraction_read_name == read_name else {}

            for row in deduped_rows:
                tx_key = tx_key_from_row(row)
                if tx_key in active_fraction_map:
                    row["frac_assigned"], row_weight = active_fraction_map[tx_key]
                    if "read_weight" in row:
                        row["read_weight"] = row_weight
                writer.writerow(row)

            if fraction_read_name == read_name:
                try:
                    fraction_read_name, fraction_map = next(fraction_iter)
                except StopIteration:
                    fraction_read_name = None
                    fraction_map = {}


def main():
    args = parse_args()

    expr_rows, expr_fieldnames = read_expr_rows(args.quant_expr)
    tracking_fieldnames = read_tracking_fieldnames(args.tracking)

    expr_by_tx = {}
    gene_to_transcripts = defaultdict(list)
    for row in expr_rows:
        tx_key = tx_key_from_row(row)
        if tx_key in expr_by_tx:
            raise RuntimeError(f"Duplicate transcript row found in quant expr: {tx_key}")
        expr_by_tx[tx_key] = row
        gene_to_transcripts[row["gene_id"]].append(tx_key)

    with tempfile.TemporaryDirectory(prefix="reassign_multigene_tracking.", dir=args.tmp_dir) as tmp_dir:
        read_sorted_path = os.path.join(tmp_dir, "tracking.by_read.tsv")
        cross_read_genes_path = os.path.join(tmp_dir, "cross_read_genes.tsv")
        read_components_path = os.path.join(tmp_dir, "cross_read_components.tsv")
        candidates_unsorted_path = os.path.join(tmp_dir, "component_candidates.unsorted.tsv")
        candidates_sorted_path = os.path.join(tmp_dir, "component_candidates.sorted.tsv")
        fractions_unsorted_path = os.path.join(tmp_dir, "read_fractions.unsorted.tsv")
        fractions_sorted_path = os.path.join(tmp_dir, "read_fractions.sorted.tsv")

        sort_tracking_body(
            args.tracking,
            read_sorted_path,
            tracking_fieldnames,
            "read_name",
            tmp_dir,
            args.sort_buffer,
        )

        tracking_state = scan_read_sorted_tracking(
            read_sorted_path,
            tracking_fieldnames,
            cross_read_genes_path,
        )

        logger.info(
            "Loaded %d expr rows, streamed %d tracking rows (%d unique candidates), and found %d cross-gene reads",
            len(expr_rows),
            tracking_state["total_rows"],
            tracking_state["kept_rows"],
            tracking_state["affected_read_count"],
        )

        if tracking_state["affected_read_count"] == 0:
            open(fractions_sorted_path, "wt").close()
            shutil.copyfile(args.quant_expr, args.output_expr)
            write_final_tracking_from_read_sorted(
                read_sorted_path,
                args.output_tracking,
                tracking_fieldnames,
                fractions_sorted_path,
            )
            logger.info(
                "No cross-gene reads found; copied expr and wrote deduplicated tracking output."
            )
            return

        component_genes = write_read_components(
            cross_read_genes_path,
            read_components_path,
            tracking_state["dsu"],
        )
        gene_to_component = {
            gene_id: component_id
            for component_id, genes in component_genes.items()
            for gene_id in genes
        }

        fixed_counts, fixed_uniq_reads = write_component_candidate_rows(
            read_sorted_path,
            tracking_fieldnames,
            gene_to_component,
            candidates_unsorted_path,
        )
        sort_tsv_body(
            candidates_unsorted_path,
            candidates_sorted_path,
            ["1,1n", "2,2"],
            tmp_dir,
            args.sort_buffer,
        )

        expr_init_counts = {
            tx_key_from_row(row): parse_float(row.get("all_reads"), 0.0) for row in expr_rows
        }

        updated_counts, updated_uniq_reads = run_component_ems_to_fraction_file(
            candidates_sorted_path,
            fractions_unsorted_path,
            component_genes,
            gene_to_transcripts,
            fixed_counts,
            fixed_uniq_reads,
            expr_init_counts,
            args.max_iter,
            args.tol,
            args.min_expr,
            args.alpha,
        )
        sort_tsv_body(
            fractions_unsorted_path,
            fractions_sorted_path,
            ["1,1", "2,2", "3,3"],
            tmp_dir,
            args.sort_buffer,
        )

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
        write_final_tracking_from_read_sorted(
            read_sorted_path,
            args.output_tracking,
            tracking_fieldnames,
            fractions_sorted_path,
        )
        logger.info("Wrote corrected expr to %s", args.output_expr)
        logger.info("Wrote corrected tracking to %s", args.output_tracking)


if __name__ == "__main__":
    main()
