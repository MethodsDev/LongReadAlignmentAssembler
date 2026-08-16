#!/usr/bin/env python3
import sys
import os
import argparse

sys.path.insert(
    0,
    os.path.sep.join(
        [os.path.dirname(os.path.realpath(__file__)), "..", "..", "pylib"]
    ),
)

from ResourceMonitor import summarize_resource_log  # noqa: E402


def summarize_file(path):
    """One implementation, shared with LRAA's per-worker summary.

    This used to be a second copy. The two agreed until one was corrected, and
    the surviving copy went on reporting every work unit shorter than the write
    interval as lasting 0.0 seconds.
    """

    summary = summarize_resource_log(path)
    if summary is None:
        print(f"ERROR reading {path}", file=sys.stderr)
    return summary


def _load_timeseries(path):
    """Return (elapsed, total_rss, total_cpu, note) arrays for plotting."""
    elapsed = []
    total_rss = []
    total_cpu = []
    note_val = None
    try:
        with open(path, "rt") as fh:
            header = next(fh).rstrip("\n").split("\t")
            idx = {k: i for i, k in enumerate(header)}
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                try:
                    el = float(parts[idx.get("elapsed_sec", 1)])
                    rss = float(parts[idx.get("rss_mb", 2)])
                    cpu = float(parts[idx.get("cpu_percent", 3)])
                    rss_ch = float(parts[idx.get("rss_mb_children", 4)])
                    cpu_ch = float(parts[idx.get("cpu_percent_children", 5)])
                    note_val = parts[idx.get("note", -1)] if idx.get("note", -1) != -1 else note_val
                except Exception:
                    continue
                elapsed.append(el)
                total_rss.append(rss + rss_ch)
                total_cpu.append(cpu + cpu_ch)
    except Exception as e:
        print(f"ERROR reading {path}: {e}", file=sys.stderr)
        return None
    return elapsed, total_rss, total_cpu, (note_val or "")


def _plot_timeseries(path, outdir, fmt="png", dpi=120, title=None):
    """Plot RSS (MB) and CPU (%) over elapsed time for a single resources.tsv file."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"Plotting disabled: matplotlib not available ({e})", file=sys.stderr)
        return None

    ts = _load_timeseries(path)
    if not ts:
        return None
    elapsed, total_rss, total_cpu, note = ts
    if not elapsed:
        return None
    base = os.path.splitext(os.path.basename(path))[0]
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{base}.{fmt}")

    fig, ax1 = plt.subplots(figsize=(8, 4))
    ax1.plot(elapsed, total_rss, color="tab:blue", label="RSS (MB)")
    ax1.set_xlabel("Elapsed (s)")
    ax1.set_ylabel("RSS (MB)", color="tab:blue")
    ax1.tick_params(axis='y', labelcolor='tab:blue')

    ax2 = ax1.twinx()
    ax2.plot(elapsed, total_cpu, color="tab:orange", alpha=0.7, label="CPU (%)")
    ax2.set_ylabel("CPU (%)", color="tab:orange")
    ax2.tick_params(axis='y', labelcolor='tab:orange')

    ttl = title
    if not ttl:
        # try to derive contig/strand from note (e.g., worker:chr1:+)
        ttl = note if note else base
    fig.suptitle(ttl)

    # build a combined legend
    lines, labels = [], []
    for ax in (ax1, ax2):
        l, lab = ax.get_legend_handles_labels()
        lines.extend(l)
        labels.extend(lab)
    fig.legend(lines, labels, loc='upper right')

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outfile, dpi=dpi)
    plt.close(fig)
    return outfile


def main():
    parser = argparse.ArgumentParser(
        description="Summarize LRAA resource monitor logs and optionally plot per-(contig,strand) usage over time"
    )
    parser.add_argument("files", nargs='+', help="resources.tsv files to summarize/plot")
    parser.add_argument("--plot", action="store_true", help="generate per-file plots (RSS+CPU over time)")
    parser.add_argument("--outdir", default=None, help="directory to write plots (default: alongside each input file)")
    parser.add_argument("--format", default="png", choices=["png", "pdf", "svg"], help="plot image format")
    parser.add_argument("--dpi", type=int, default=120, help="plot DPI for raster formats")
    args = parser.parse_args()

    print("\t".join(["file", "samples", "duration_sec", "peak_rss_mb_total", "avg_rss_mb_total", "peak_cpu_percent_total", "avg_cpu_percent_total"]))
    for p in args.files:
        res = summarize_file(p)
        if not res:
            continue
        print("\t".join([
            res["file"],
            str(res["samples"]),
            f"{res['duration_sec']:.1f}",
            f"{res['peak_rss_mb_total']:.1f}",
            f"{res['avg_rss_mb_total']:.1f}",
            f"{res['peak_cpu_percent_total']:.1f}",
            f"{res['avg_cpu_percent_total']:.1f}",
        ]))

        if args.plot:
            outdir = args.outdir if args.outdir else os.path.dirname(os.path.abspath(p))
            _plot_timeseries(p, outdir=outdir, fmt=args.format, dpi=args.dpi)


if __name__ == "__main__":
    main()
