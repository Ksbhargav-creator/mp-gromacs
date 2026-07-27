#!/usr/bin/env python3
"""
plot_instrumented.py  --  plot the runtime force-loop product-magnitude
distribution recorded by the instrumented GROMACS build (nga_range_stats.hpp).

Reads nga_gromacs_stats.csv (operand, binade_exp, count) and renders one
log2-binade histogram panel per operand -- the runtime companion to the static
operand histograms from plot_histograms.py, on the same magnitude axis. Overlays
the exponent reach of float32 / float16 / posit<16,1> / posit<32,2> so under/
overflow risk for each candidate format is read straight off the plot.

Usage:
    python3 plot_instrumented.py nga_gromacs_stats.csv --label benchMEM --out figs/
"""
import argparse, os, csv
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# min-normal binade (log2) for each format's smallest representable magnitude
FORMAT_MINEXP = {
    "float16":     (-14, "#d1495b"),
    "posit<16,1>": (-28, "#2e8b57"),
    "float32":     (-126, "#3a7ca5"),
    "posit<32,2>": (-120, "#6a4c93"),
}


def load(path):
    data = defaultdict(lambda: defaultdict(int))
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            data[row["operand"]][int(row["binade_exp"])] += int(row["count"])
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", help="nga_gromacs_stats.csv from the instrumented run")
    ap.add_argument("--label", default="instrumented")
    ap.add_argument("--out", default="figs")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    data = load(args.csv)
    ops = list(data.keys())
    ncol = 2
    nrow = int(np.ceil(len(ops) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.2 * nrow))
    axes = np.atleast_1d(axes).ravel()

    for ax, op in zip(axes, ops):
        exps = sorted(data[op])
        counts = [data[op][e] for e in exps]
        ax.bar(exps, counts, width=0.9, color="#404040")
        span = (max(exps) - min(exps)) if exps else 0
        ax.set_title(f"{op}\nbinade span = {span} (~{span*0.301:.1f} decades)",
                     fontsize=9)
        ax.set_xlabel("binade  (exponent of |value|, base 2)")
        ax.set_ylabel("count"); ax.grid(True, alpha=0.2)
        ymax = ax.get_ylim()[1]
        lo = min(exps) if exps else -1
        for i, (name, (mexp, col)) in enumerate(FORMAT_MINEXP.items()):
            if mexp > lo - 8:
                ax.axvline(mexp, color=col, ls="--", lw=1, alpha=0.7)
                ax.text(mexp, ymax * (0.95 - 0.1 * i), f" {name} min",
                        color=col, fontsize=6.5, rotation=90, va="top")

    for ax in axes[len(ops):]:
        ax.axis("off")
    fig.suptitle(f"GROMACS instrumented force-loop product magnitudes  --  {args.label}",
                 fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = os.path.join(args.out, f"instrumented_histograms_{args.label}.png")
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
