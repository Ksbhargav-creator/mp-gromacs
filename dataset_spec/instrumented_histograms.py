#!/usr/bin/env python3
"""
instrumented_histograms.py -- value-distribution histograms from the RUNTIME instrumented
GROMACS/CHARMM stats (Method A), in fp32 and fp64.

Input
-----
The aggregate stats CSVs written by the instrumented binary:
    runs_instr/<APP>/<fp32|fp64>/<system>/nga_gromacs_stats.csv
each a per-operand log2 (binade) histogram:  operand,binade_exp,count
where binade_exp = floor(log2|x|). Counts are summed across all systems of an app so the
result is the app-wide runtime distribution of each quantity.

Output
------
One figure per app: a small-multiple grid, one panel per operand, plotting
    number of operands  (y)   vs   log10(magnitude)  (x)
with fp32 and fp64 overlaid. x is log10 (binade_exp is converted via log10(2)); the
dashed line marks |x| = 1 (posit's peak-accuracy point).

Note: the binade histogram is coarse (one bin per power of two), so fp32 and fp64 nearly
coincide here -- that is expected and is itself the point (both precisions see the same
value distribution). The true per-bit fp32-vs-fp64 error needs the per-interaction trace
(the paired error harness), not these aggregate bins.

Usage
-----
    python3 dataset_spec/instrumented_histograms.py                 # both apps, defaults
    python3 dataset_spec/instrumented_histograms.py --runs runs_instr --out docs
"""
import argparse
import csv
import glob
import math
import os
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOG10_2 = math.log10(2.0)

# a readable operand order (operands not listed are appended alphabetically)
ORDER = ["coulomb_qq", "coulomb_term", "rinv",
         "lj_c6", "lj_c12", "lj_rinv6",
         "lj_term6", "lj_term12", "lj_force",
         "nb_fscal", "nb_force_acc"]

PREC_STYLE = {  # precision -> (color, how to draw)
    "fp64": ("#1f77b4", "bars"),
    "fp32": ("#d62728", "step"),
}


def load_app(runs_dir, app):
    """-> {precision: {operand: {binade_exp: summed_count}}}"""
    out = {}
    for prec in ("fp32", "fp64"):
        per_op = defaultdict(lambda: defaultdict(int))
        files = glob.glob(os.path.join(runs_dir, app, prec, "*", "nga_gromacs_stats.csv"))
        for f in files:
            with open(f) as fh:
                for row in csv.DictReader(fh):
                    try:
                        e = int(row["binade_exp"]); c = int(row["count"])
                    except (KeyError, ValueError):
                        continue
                    per_op[row["operand"]][e] += c
        if files:
            out[prec] = per_op
    return out


def operand_list(app_data):
    ops = set()
    for prec in app_data.values():
        ops |= set(prec.keys())
    ordered = [o for o in ORDER if o in ops] + sorted(o for o in ops if o not in ORDER)
    return ordered


def plot_app(app, app_data, out_dir):
    ops = operand_list(app_data)
    if not ops:
        print(f"[{app}] no operands -- skipping")
        return
    ncol = 3
    nrow = math.ceil(len(ops) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.7 * nrow))
    axes = axes.ravel()

    for ax, op in zip(axes, ops):
        for prec, (color, how) in PREC_STYLE.items():
            hist = app_data.get(prec, {}).get(op)
            if not hist:
                continue
            exps = sorted(hist)
            xs = [e * LOG10_2 for e in exps]          # log10(magnitude)
            ys = [hist[e] for e in exps]
            if how == "bars":
                ax.bar(xs, ys, width=LOG10_2, align="edge", color=color,
                       alpha=0.55, label=prec)
            else:
                ax.step([x + LOG10_2 for x in xs], ys, where="pre",
                        color=color, lw=1.3, label=prec)
        ax.axvline(0.0, color="k", ls="--", lw=0.7)
        ax.set_title(op, fontsize=9)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.15)
    # hide any unused panels
    for ax in axes[len(ops):]:
        ax.set_visible(False)

    # shared labels + one legend
    fig.supxlabel(r"$\log_{10}$(magnitude)", fontsize=10)
    fig.supylabel("Number of operands", fontsize=10)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", fontsize=9)
    fig.suptitle(f"{app}: runtime operand/product distributions (fp32 vs fp64)", fontsize=12)
    fig.tight_layout(rect=(0.02, 0.02, 1, 0.96))

    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"hist_runtime_{app}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    # quick numeric summary
    total = {p: sum(sum(h.values()) for h in d.values()) for p, d in app_data.items()}
    print(f"[{app}] {len(ops)} operands; pooled counts " +
          ", ".join(f"{p}={n:,}" for p, n in total.items()))
    print(f"   wrote {out}")


def plot_app_pooled(app, app_data, out_dir):
    """ALL instrumented quantities aggregated into one distribution for the app:
    number of instrumented values (y) vs log10(magnitude) (x), fp32 and fp64 overlaid."""
    if not app_data:
        print(f"[{app}] no data -- skipping pooled figure")
        return
    # sum counts over every operand, per precision -> {prec: {binade_exp: total_count}}
    pooled = {}
    for prec, per_op in app_data.items():
        agg = defaultdict(int)
        for hist in per_op.values():
            for e, c in hist.items():
                agg[e] += c
        pooled[prec] = agg

    fig, ax = plt.subplots(figsize=(10, 5))
    for prec, (color, how) in PREC_STYLE.items():
        agg = pooled.get(prec)
        if not agg:
            continue
        exps = sorted(agg)
        xs = [e * LOG10_2 for e in exps]         # log10(magnitude)
        ys = [agg[e] for e in exps]
        if how == "bars":
            ax.bar(xs, ys, width=LOG10_2, align="edge", color=color, alpha=0.55, label=prec)
        else:
            ax.step([x + LOG10_2 for x in xs], ys, where="pre", color=color, lw=1.4, label=prec)
    ax.axvline(0.0, color="k", ls="--", lw=0.9, label="|x| = 1")
    ax.set_xlabel(r"$\log_{10}$(magnitude) of instrumented value", fontsize=11)
    ax.set_ylabel("Number of instrumented values", fontsize=11)
    ax.set_title(f"{app}: pooled runtime value distribution (all instrumented quantities)",
                 fontsize=12)
    ax.ticklabel_format(axis="y", style="plain")
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _pos: f"{int(v):,}"))
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(alpha=0.15)
    fig.tight_layout()

    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"hist_runtime_pooled_{app}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    totals = {p: sum(a.values()) for p, a in pooled.items()}
    print(f"[{app}] pooled ALL operands; total instrumented values " +
          ", ".join(f"{p}={n:,}" for p, n in totals.items()))
    print(f"   wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default="runs_instr", help="root of instrumented run outputs")
    ap.add_argument("--out", default="docs", help="output directory for figures")
    ap.add_argument("--apps", nargs="+", default=["GROMACS", "CHARMM"])
    ap.add_argument("--pooled", action="store_true",
                    help="one figure per app pooling ALL instrumented quantities together")
    args = ap.parse_args()

    for app in args.apps:
        data = load_app(args.runs, app)
        if args.pooled:
            plot_app_pooled(app, data, args.out)
        else:
            plot_app(app, data, args.out)


if __name__ == "__main__":
    main()
