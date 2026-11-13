#!/usr/bin/env python3
"""Per-cluster pseudobulk quant.expr files from ONE pooled single-cell quantification.

The cluster-guided single-cell workflow quantifies every cell against the final catalog
in a single LRAA quant-only pass (the pooled, all-cells EM). This script splits that one
assignment by cell cluster: each read row of the pooled tracking file is credited to the
cluster of its cell barcode, and the rows are summed per (cluster, isoform). No EM is run
here, so every cluster's ambiguous reads keep the fractions the pooled EM gave them.

Columns are computed exactly as LRAA's streaming quant.expr writer computes them
(pylib/StreamingQuant.py write_expr), restricted to the cluster's reads:

    uniq_reads                 read rows with is_unique == 1
    all_reads                  sum of frac_assigned
    isoform_fraction           all_reads / the gene's all_reads in the cluster
    unique_gene_read_fraction  uniq_reads / the gene's all_reads in the cluster
    TPM                        all_reads / the cluster's all_reads, x 1e6
    RPM_total_reads            see below
    uniq_FSM_reads             read rows with is_unique == 1 and is_FSM == 1
    has_FSM_read               1 if any read row has is_FSM == 1

The remaining columns (exons, introns, splice_hash_code and any others) are copied from
the pooled quant.expr, which also fixes the isoform list: every cluster file lists every
isoform of the pooled file, zeros included, sorted by transcript_id as LRAA writes them.

RPM_total_reads: LRAA divides by the mapped reads in the bam it quantified. A cluster's
own mapped-read total is not known here, so the denominator is the cluster's assigned
reads (its all_reads total), which makes this column equal TPM. A comment line in each
file says so.

Reads from barcodes that are in no cluster are counted in the pooled quantification (and
in the single-cell matrices built from it) but in no cluster file. Their totals are
reported in the --partition_summary table.

Before writing anything the script checks that the tracking reproduces the pooled
quant.expr: per isoform, the per-cluster sums plus the unclustered remainder must equal
the pooled uniq_reads, uniq_FSM_reads and has_FSM_read exactly; a mismatch means the two
files are not from the same run and stops the script. all_reads is compared to the pooled
file's one-decimal rounding and only warned about: LRAA's default streaming writer reports
the tracking's sum, while the EM report path (stream_reads off) reports EM's own counts.

Output files are named <output_prefix>.<cluster>.LRAA.quant-only.quant.expr, the names
the former per-cluster quant runs wrote, so readers that take the cluster from the file
name (util/sc/site_usage/build_site_event_read_tracks.py) keep working.
"""

import argparse
import collections
import gzip
import logging
import os
import shutil
import subprocess
import sys

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s : %(levelname)s : %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

UNCLUSTERED = "__unclustered__"
# all_reads is written with one decimal, so a correct sum can sit up to 0.05 away
ALL_READS_TOLERANCE = 0.051


def open_text(path):
    """Text handle on a plain or gzipped file; gzip goes through pigz when available."""
    if not path.endswith(".gz"):
        return open(path, "rt")
    if shutil.which("pigz"):
        proc = subprocess.Popen(["pigz", "-dc", path], stdout=subprocess.PIPE, text=True)
        return _ProcReader(proc, path)
    return gzip.open(path, "rt")


class _ProcReader:
    """File-like wrapper that checks the decompressor's exit status on close."""

    def __init__(self, proc, path):
        self.proc = proc
        self.path = path

    def __iter__(self):
        return iter(self.proc.stdout)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.proc.stdout.close()
        rc = self.proc.wait()
        if rc != 0 and exc[0] is None:
            raise RuntimeError(f"pigz -dc {self.path} exited {rc}")
        return False


def read_cell_clusters(path):
    """cell_barcode -> cluster, parsed as util/sc/partition_bam_by_cell_cluster.py does:
    two tab-separated fields per line; the exact header row 'cell_barcode<TAB>cluster'
    on line 1 is skipped; other lines without two fields are skipped."""
    cell_to_cluster = dict()
    with open(path, "rt") as fh:
        for lineno, line in enumerate(fh, start=1):
            vals = line.rstrip().split("\t")
            if len(vals) != 2:
                continue
            barcode, cluster = vals
            if lineno == 1 and (barcode, cluster) == ("cell_barcode", "cluster"):
                continue
            cell_to_cluster[barcode] = cluster
    if not cell_to_cluster:
        raise SystemExit(f"no cell -> cluster assignments read from {path}")
    return cell_to_cluster


def read_pooled_expr(path):
    """(comment lines, header fields, rows as lists) of the pooled quant.expr."""
    comments, header, rows = [], None, []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
                continue
            vals = line.rstrip("\n").split("\t")
            if header is None:
                header = vals
                continue
            rows.append(vals)
    if header is None:
        raise SystemExit(f"no header line in {path}")
    for col in ("gene_id", "transcript_id", "uniq_reads", "all_reads", "isoform_fraction",
                "unique_gene_read_fraction", "TPM", "RPM_total_reads", "uniq_FSM_reads",
                "has_FSM_read"):
        if col not in header:
            raise SystemExit(f"{path}: column {col} missing")
    return comments, header, rows


class Totals:
    __slots__ = ("frac", "uniq", "uniq_fsm", "any_fsm")

    def __init__(self):
        self.frac = 0.0
        self.uniq = 0
        self.uniq_fsm = 0
        self.any_fsm = False


def tally_tracking(path, cell_to_cluster, known_transcripts):
    """(cluster, transcript_id) -> Totals over the tracking rows, plus per-group counts of
    read rows and distinct barcodes. Barcodes in no cluster go to UNCLUSTERED."""
    totals = collections.defaultdict(Totals)
    barcodes_seen = collections.defaultdict(set)
    col = None
    n_rows = 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            vals = line.rstrip("\n").split("\t")
            if vals[0] == "gene_id":
                if col is None:
                    col = {name: i for i, name in enumerate(vals)}
                    for name in ("transcript_id", "read_name", "frac_assigned", "is_unique", "is_FSM"):
                        if name not in col:
                            raise SystemExit(f"{path}: tracking column {name} missing")
                    i_tx, i_read = col["transcript_id"], col["read_name"]
                    i_frac, i_uniq, i_fsm = col["frac_assigned"], col["is_unique"], col["is_FSM"]
                continue
            if col is None:
                raise SystemExit(f"{path}: data row before the column header")
            n_rows += 1
            transcript_id = vals[i_tx]
            if transcript_id not in known_transcripts:
                raise SystemExit(
                    f"{path}: transcript {transcript_id} is not in the pooled quant.expr; "
                    "the tracking and quant.expr must come from the same quantification")
            barcode = vals[i_read].split("^", 1)[0]
            cluster = cell_to_cluster.get(barcode, UNCLUSTERED)
            barcodes_seen[cluster].add(barcode)
            t = totals[(cluster, transcript_id)]
            t.frac += float(vals[i_frac])
            if vals[i_fsm] == "1":
                t.any_fsm = True
            if vals[i_uniq] == "1":
                t.uniq += 1
                if vals[i_fsm] == "1":
                    t.uniq_fsm += 1
    if col is None:
        raise SystemExit(f"{path}: no tracking column header found")
    logger.info("read %d tracking rows", n_rows)
    return totals, barcodes_seen


def check_against_pooled(header, rows, totals):
    """Per isoform, clusters + unclustered must reproduce the pooled quant.expr."""
    ix = {name: i for i, name in enumerate(header)}
    by_tx = collections.defaultdict(Totals)
    for (cluster, tx), t in totals.items():
        s = by_tx[tx]
        s.frac += t.frac
        s.uniq += t.uniq
        s.uniq_fsm += t.uniq_fsm
        s.any_fsm = s.any_fsm or t.any_fsm
    bad, off = [], []
    for r in rows:
        tx = r[ix["transcript_id"]]
        s = by_tx.get(tx, Totals())
        if (int(r[ix["uniq_reads"]]) != s.uniq
                or int(r[ix["uniq_FSM_reads"]]) != s.uniq_fsm
                or (r[ix["has_FSM_read"]] == "1") != s.any_fsm):
            bad.append((tx, r[ix["uniq_reads"]], s.uniq, r[ix["uniq_FSM_reads"]], s.uniq_fsm,
                        r[ix["has_FSM_read"]], int(s.any_fsm)))
        pooled_all = float(r[ix["all_reads"]])
        if abs(pooled_all - s.frac) > ALL_READS_TOLERANCE + 1e-4 * pooled_all:
            off.append((tx, r[ix["all_reads"]], f"{s.frac:.3f}"))
    if off:
        # The streaming writer (LRAA's default, stream_reads) reports all_reads as the sum
        # of the tracking's frac_assigned, so this should not happen there. The EM report
        # path (stream_reads off) reports EM's own counts, which can differ slightly.
        for o in off[:10]:
            logger.warning("all_reads differs, pooled quant.expr vs tracking sum: %s %s / %s", *o)
        logger.warning("%d of %d isoforms: pooled all_reads differs from the tracking's sum of "
                       "frac_assigned; cluster values are tracking sums", len(off), len(rows))
    if bad:
        for b in bad[:10]:
            logger.error("pooled vs tracking mismatch: %s uniq %s/%s uniqFSM %s/%s hasFSM %s/%s", *b)
        raise SystemExit(
            f"{len(bad)} isoforms: the tracking does not reproduce the pooled quant.expr's "
            "uniq_reads, uniq_FSM_reads or has_FSM_read; the two files must come from the "
            "same quantification")


def write_cluster_expr(path, comments, header, rows, cluster, totals):
    ix = {name: i for i, name in enumerate(header)}
    gene_of = {r[ix["transcript_id"]]: r[ix["gene_id"]] for r in rows}
    gene_frac = collections.defaultdict(float)
    cluster_total = 0.0
    for r in rows:
        t = totals.get((cluster, r[ix["transcript_id"]]))
        if t is not None:
            gene_frac[gene_of[r[ix["transcript_id"]]]] += t.frac
            cluster_total += t.frac
    zero = Totals()
    with open(path, "wt") as ofh:
        for c in comments:
            print(c, file=ofh)
        print(f"# cluster {cluster}: derived from the pooled (all cells) quantification by "
              "summing the cluster's cells; RPM_total_reads is per million of the cluster's "
              "assigned reads (equal to TPM)", file=ofh)
        print("\t".join(header), file=ofh)
        for r in sorted(rows, key=lambda r: r[ix["transcript_id"]]):
            tx = r[ix["transcript_id"]]
            t = totals.get((cluster, tx), zero)
            gene_total = gene_frac.get(gene_of[tx], 0.0)
            out = list(r)
            out[ix["uniq_reads"]] = f"{t.uniq}"
            out[ix["all_reads"]] = f"{t.frac:.1f}"
            out[ix["isoform_fraction"]] = f"{(t.frac / gene_total if gene_total > 0 else 0.0):.3f}"
            out[ix["unique_gene_read_fraction"]] = f"{(t.uniq / gene_total if gene_total > 0 else 0.0):.3f}"
            tpm = t.frac / cluster_total * 1e6 if cluster_total > 0 else 0.0
            out[ix["TPM"]] = f"{tpm:.3f}"
            out[ix["RPM_total_reads"]] = f"{tpm:.3f}"
            out[ix["uniq_FSM_reads"]] = f"{t.uniq_fsm}"
            out[ix["has_FSM_read"]] = "1" if t.any_fsm else "0"
            print("\t".join(out), file=ofh)
    return cluster_total


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tracking", required=True,
                        help="pooled LRAA quant.tracking (.gz ok)")
    parser.add_argument("--quant_expr", required=True,
                        help="pooled LRAA quant.expr from the same quantification")
    parser.add_argument("--cell_clusters", required=True,
                        help="cell_barcode <tab> cluster (as given to partition_bam_by_cell_cluster.py)")
    parser.add_argument("--output_prefix", required=True,
                        help="files are written as <output_prefix>.<cluster>.LRAA.quant-only.quant.expr")
    parser.add_argument("--output_dir", default=".", help="directory for the per-cluster files")
    parser.add_argument("--partition_summary", default=None,
                        help="optional TSV: per cluster (and unclustered) barcodes with reads, "
                             "read rows' assigned reads, unique reads")
    args = parser.parse_args()

    cell_to_cluster = read_cell_clusters(args.cell_clusters)
    clusters = sorted(set(cell_to_cluster.values()))
    if UNCLUSTERED in clusters:
        raise SystemExit(f"a cluster is named {UNCLUSTERED}, which this script reserves")
    logger.info("%d barcodes in %d clusters", len(cell_to_cluster), len(clusters))

    comments, header, rows = read_pooled_expr(args.quant_expr)
    known = {r[header.index("transcript_id")] for r in rows}
    logger.info("%d isoforms in the pooled quant.expr", len(rows))

    totals, barcodes_seen = tally_tracking(args.tracking, cell_to_cluster, known)
    check_against_pooled(header, rows, totals)

    os.makedirs(args.output_dir, exist_ok=True)
    assigned = dict()
    for cluster in clusters:
        path = os.path.join(args.output_dir, f"{args.output_prefix}.{cluster}.LRAA.quant-only.quant.expr")
        assigned[cluster] = write_cluster_expr(path, comments, header, rows, cluster, totals)
        logger.info("cluster %s: %.1f assigned reads -> %s", cluster, assigned[cluster], path)

    group_frac = collections.defaultdict(float)
    group_uniq = collections.defaultdict(int)
    for (c, _), t in totals.items():
        group_frac[c] += t.frac
        group_uniq[c] += t.uniq
    unclustered = group_frac[UNCLUSTERED]
    pooled = sum(group_frac.values())
    logger.info("unclustered barcodes: %d, %.1f assigned reads (%.2f%% of %.1f)",
                len(barcodes_seen.get(UNCLUSTERED, ())), unclustered,
                100.0 * unclustered / pooled if pooled else 0.0, pooled)

    if args.partition_summary:
        with open(args.partition_summary, "wt") as ofh:
            print("# Reads of the pooled quantification by cell cluster. assigned_reads: sum of "
                  "frac_assigned; unique_reads: read rows assigned to one isoform only. The "
                  "__unclustered__ row holds barcodes in no cluster: counted in the single-cell "
                  "matrices, in no cluster's pseudobulk.", file=ofh)
            print("\t".join(["cluster", "barcodes_with_reads", "assigned_reads", "unique_reads",
                             "fraction_of_assigned_reads"]), file=ofh)
            for cluster in clusters + [UNCLUSTERED]:
                a, u = group_frac[cluster], group_uniq[cluster]
                print("\t".join([cluster, str(len(barcodes_seen.get(cluster, ()))), f"{a:.1f}", str(u),
                                 f"{(a / pooled if pooled else 0.0):.6f}"]), file=ofh)

    return 0


if __name__ == "__main__":
    sys.exit(main())
