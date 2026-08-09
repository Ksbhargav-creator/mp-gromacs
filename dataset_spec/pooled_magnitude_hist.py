#!/usr/bin/env python3
"""
pooled_magnitude_hist.py -- pooled magnitude histogram: COUNT vs log10(|operand|).

This is a different view from histograms.py. histograms.py is density-normalized and
per-kernel; this one POOLS every operand value (across all systems and all quantities)
for an app into a single distribution and plots the raw *count* against log10(|value|).
It answers one blunt question: across the whole app, how many operand values sit in each
decade of magnitude -- and how many fall outside what fp32 can even represent?

fp32 focus
----------
GROMACS stores force-field parameters in single precision, so we cast each pooled value
to float32 before taking |x| and log10. That cast is what makes these literally "fp32
operands"; it also lets us count the values that FLUSH TO ZERO (underflow) or OVERFLOW in
fp32 -- operands the format cannot hold, which is part of the posit motivation.

Source
------
The raw operand arrays in  <run-dir>/*/operands.npz  written by scripts/nga_gromacs.py.
Defaults: GROMACS -> runs/ , CHARMM -> runs_charmm/ .

Usage
-----
    python3 dataset_spec/pooled_magnitude_hist.py                 # both apps, defaults
    python3 dataset_spec/pooled_magnitude_hist.py --params-only   # drop sampled products
    python3 dataset_spec/pooled_magnitude_hist.py --binwidth 0.25 --out docs
"""
import argparse
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Arrays that are sampled products / coordinates rather than raw force-field operands.
# --params-only excludes these so the pool is just the distinct FF parameters.
DERIVED = {"coulomb_qq_product", "position"}

# fp32 magnitude landmarks (IEEE binary32)
FP32_MIN_NORMAL = np.float32(np.finfo(np.float32).tiny)   # ~1.18e-38
FP32_MAX        = np.float32(np.finfo(np.float32).max)    # ~3.40e38


def pool_app(run_dir, params_only=False):
    """Pool every operand value under run_dir/*/operands.npz into one fp32 magnitude array.

    Returns (log10_abs_values, composition, n_underflow, n_overflow, n_files)."""
    files = sorted(glob.glob(os.path.join(run_dir, "*", "operands.npz")))
    chunks, comp = [], {}
    n_uf = n_of = 0
    for f in files:
        d = np.load(f)
        for k in d.files:
            if params_only and k in DERIVED:
                continue
            a = np.asarray(d[k], dtype=np.float64).ravel()
            a = a[np.isfinite(a) & (a != 0.0)]
            if a.size == 0:
                continue
            a32 = np.abs(a).astype(np.float32)            # <-- the fp32 cast
            n_uf += int(np.sum(a32 == np.float32(0.0)))   # underflowed to zero
            n_of += int(np.sum(~np.isfinite(a32)))        # overflowed to inf
            good = a32[(a32 > 0) & np.isfinite(a32)]
            if good.size:
                chunks.append(good.astype(np.float64))
                comp[k] = comp.get(k, 0) + int(good.size)
    pooled = np.concatenate(chunks) if chunks else np.empty(0)
    return pooled, comp, n_uf, n_of, len(files)          # raw fp32 magnitudes


def make_log_bins(mag, binwidth):
    """Log-spaced magnitude bin EDGES (so a log x-axis shows real values, e.g. 0.2)."""
    lo = np.floor(np.log10(mag.min()) / binwidth) * binwidth
    hi = np.ceil(np.log10(mag.max()) / binwidth) * binwidth
    return 10.0 ** np.arange(lo, hi + binwidth, binwidth)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gromacs-dir", default="runs", help="dir of GROMACS operands.npz runs")
    ap.add_argument("--charmm-dir", default="runs_charmm", help="dir of CHARMM operands.npz runs")
    ap.add_argument("--binwidth", type=float, default=0.5, help="log10 bin width (decades)")
    ap.add_argument("--params-only", action="store_true",
                    help="pool only distinct FF parameters (drop sampled coulomb_qq_product/position)")
    ap.add_argument("--out", default="docs", help="output directory for the figure")
    args = ap.parse_args()

    apps = [
        ("GROMACS", args.gromacs_dir, "#1f77b4"),
        ("CHARMM",  args.charmm_dir,  "#d62728"),
    ]
    pooled = {}
    for name, d, _ in apps:
        mag, comp, n_uf, n_of, nf = pool_app(d, args.params_only)
        pooled[name] = (mag, comp, n_uf, n_of, nf)
        print(f"[{name}] {nf} run(s), {mag.size:,} pooled fp32 operand values"
              f"  ({'params only' if args.params_only else 'all operands'})")
        if mag.size:
            lo10, hi10 = np.log10(mag.min()), np.log10(mag.max())
            print(f"   magnitude span: {mag.min():.2e} .. {mag.max():.2e}  "
                  f"({hi10-lo10:.1f} decades)")
        if n_uf or n_of:
            print(f"   fp32 cannot represent: {n_uf:,} underflow->0, {n_of:,} overflow->inf")
        top = sorted(comp.items(), key=lambda kv: -kv[1])[:6]
        print("   top contributors: " + ", ".join(f"{k}={v:,}" for k, v in top))

    fp32_lo, fp32_hi = np.log10(FP32_MIN_NORMAL), np.log10(FP32_MAX)
    os.makedirs(args.out, exist_ok=True)
    tag = "_paramsonly" if args.params_only else ""

    # One SEPARATE histogram per app -- GROMACS and CHARMM never share a figure.
    for name, _, color in apps:
        mag = pooled[name][0]
        if mag.size == 0:
            print(f"[{name}] no operands pooled -- skipping figure")
            continue
        log10 = np.log10(mag)                       # x = log10 of the operand magnitude
        lo = np.floor(log10.min() / args.binwidth) * args.binwidth
        hi = np.ceil(log10.max() / args.binwidth) * args.binwidth
        bins = np.arange(lo, hi + args.binwidth, args.binwidth)

        fig, ax = plt.subplots(figsize=(9, 5))
        # y = number of operands ; x = log10(operand magnitude)
        ax.hist(log10, bins=bins, color=color, alpha=0.85, edgecolor="white", linewidth=0.3)
        ax.axvline(0.0, color="k", ls="--", lw=0.9, label="|x| = 1")
        ax.set_xlim(lo - 1, hi + 1)
        ax.set_xlabel(r"$\log_{10}$(operand magnitude), fp32")
        ax.set_ylabel("Number of operands")
        # show full integer counts with thousands separators (no "1e6" offset)
        ax.ticklabel_format(axis="y", style="plain")
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _pos: f"{int(v):,}"))
        subtitle = "FF parameters only" if args.params_only else "all operands"
        ax.set_title(f"{name}: pooled fp32 operand magnitudes  ({subtitle})")
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(alpha=0.15, which="both")
        fig.text(0.5, 0.005,
                 f"fp32 normal range is 1e{fp32_lo:.0f}..1e{fp32_hi:.0f}; "
                 f"all operands fall well inside it",
                 ha="center", fontsize=8, style="italic", color="gray")
        fig.tight_layout(rect=(0, 0.02, 1, 1))

        out = os.path.join(args.out, f"hist_pooled_fp32_{name}{tag}.png")
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
