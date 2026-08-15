#!/usr/bin/env python3
"""Response surface figure: MARD against alpha, per sample, both 3p settings.

One panel per stratum, samples drawn individually rather than averaged. The
whole finding is that averaging across strata cancels two real and opposing
effects, so a mean curve here would hide the result.

alpha=0 has no place on a log axis; it is drawn one decade below the smallest
positive grid point and the axis label says so. Writing it as 1e-4 in the data
would misstate the spacing.
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_config as C  # noqa: E402

ZERO_AT = 1e-4  # where alpha=0 is drawn; one decade below the 3e-4 grid point


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics", default=os.path.join(C.BASE, "results", "grid_metrics.tsv"))
    ap.add_argument("--out", default=os.path.join(C.BASE, "results", "alpha_response_curves.png"))
    ap.add_argument("--metric", default="mard")
    args = ap.parse_args()

    df = pd.read_csv(args.metrics, sep="\t")
    df["alpha"] = pd.to_numeric(df["alpha"], errors="coerce")
    df[args.metric] = pd.to_numeric(df[args.metric], errors="coerce")
    df["x"] = df["alpha"].where(df["alpha"] > 0, ZERO_AT)
    df["stratum"] = df["sample"].map(
        lambda s: "SIRV E1" if "_E1_" in s else "SIRV E2" if "_E2_" in s else "MORF"
    )

    strata = [s for s in ("SIRV E1", "SIRV E2", "MORF") if (df.stratum == s).any()]
    fig, axes = plt.subplots(1, len(strata), figsize=(5.2 * len(strata), 4.4), squeeze=False)
    colors = plt.cm.tab10.colors
    for ax, st in zip(axes[0], strata):
        sub = df[df.stratum == st]
        for i, (smp, g) in enumerate(sub.groupby("sample")):
            for w, style, mk in ((1, "-", "o"), (0, "--", None)):
                gg = g[g.threep == w].sort_values("x")
                if gg.empty:
                    continue
                ax.plot(gg["x"], gg[args.metric], style, color=colors[i % 10],
                        marker=mk, ms=3, lw=1.4,
                        label=smp.replace("CL_", "").replace("_sirv", "") if w else None)
        ax.axvline(C.DEFAULT_ALPHA, color="0.6", lw=0.8, ls=":")
        ax.set_xscale("log")
        ax.set_xlabel(f"EM_alpha   (alpha=0 drawn at {ZERO_AT:g}; dotted line = shipped default)")
        ax.set_ylabel(args.metric)
        ax.set_title(f"{st}   solid = 3' weighting on, dashed = off")
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
