#!/usr/bin/env python3
import sys


def summarize_file(path):
    peak_rss = 0.0
    peak_cpu = 0.0
    sum_rss = 0.0
    sum_cpu = 0.0
    n = 0
    first_ts = None
    last_ts = None
    try:
        with open(path, "rt") as fh:
            header = next(fh).rstrip("\n").split("\t")
            idx = {k: i for i, k in enumerate(header)}
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                try:
                    ts = float(parts[idx.get("epoch_ts", 0)])
                    rss = float(parts[idx.get("rss_mb", 2)])
                    cpu = float(parts[idx.get("cpu_percent", 3)])
                    rss_ch = float(parts[idx.get("rss_mb_children", 4)])
                    cpu_ch = float(parts[idx.get("cpu_percent_children", 5)])
                except Exception:
                    continue
                total_rss = rss + rss_ch
                total_cpu = cpu + cpu_ch
                peak_rss = max(peak_rss, total_rss)
                peak_cpu = max(peak_cpu, total_cpu)
                sum_rss += total_rss
                sum_cpu += total_cpu
                n += 1
                if first_ts is None:
                    first_ts = ts
                last_ts = ts
    except Exception as e:
        print(f"ERROR reading {path}: {e}", file=sys.stderr)
        return None
    avg_rss = (sum_rss / n) if n else 0.0
    avg_cpu = (sum_cpu / n) if n else 0.0
    duration = (last_ts - first_ts) if (first_ts is not None and last_ts is not None) else 0.0
    return {
        "file": path,
        "samples": n,
        "duration_sec": duration,
        "peak_rss_mb_total": peak_rss,
        "avg_rss_mb_total": avg_rss,
        "peak_cpu_percent_total": peak_cpu,
        "avg_cpu_percent_total": avg_cpu,
    }


def main():
    if len(sys.argv) < 2:
        print("usage: summarize_resource_logs.py <resources.tsv> [more.tsv ...]", file=sys.stderr)
        sys.exit(1)
    print("\t".join(["file", "samples", "duration_sec", "peak_rss_mb_total", "avg_rss_mb_total", "peak_cpu_percent_total", "avg_cpu_percent_total"]))
    for p in sys.argv[1:]:
        res = summarize_file(p)
        if not res:
            continue
        print(
            "\t".join(
                [
                    res["file"],
                    str(res["samples"]),
                    f"{res['duration_sec']:.1f}",
                    f"{res['peak_rss_mb_total']:.1f}",
                    f"{res['avg_rss_mb_total']:.1f}",
                    f"{res['peak_cpu_percent_total']:.1f}",
                    f"{res['avg_cpu_percent_total']:.1f}",
                ]
            )
        )


if __name__ == "__main__":
    main()
