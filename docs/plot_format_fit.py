#!/usr/bin/env python3
"""Fraction bits (precision) vs value magnitude for several posit formats, with the
measured GROMACS operand ranges overlaid -- read straight from the dataset records so
the figure regenerates from committed data (no hardcoded numbers).

Usage:
    python3 docs/plot_format_fit.py dataset/records/gromacs.nonbonded_coulomb.*.json
"""
import glob, json, math, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def regime_bits(k):
    return k + 2 if k >= 0 else -k + 1

def frac(n, es, s):
    return max(0, n - 1 - regime_bits(math.floor(s / (2**es))) - es)

def band_from_record(path):
    """Return list of (label, s_min, s_max) for each quantity in a record."""
    r = json.load(open(path))
    sysname = r["input_case"]["name"]
    out = []
    for q in r["quantities"]:
        s = q["stats"]
        if "abs_min" in s and s["abs_min"] > 0:
            out.append((f"{sysname} {q['name'].replace('coulomb_qq_product','|q_i q_j|')}",
                        math.floor(math.log2(s["abs_min"])), math.ceil(math.log2(s["abs_max"]))))
    return out

def main(paths):
    scales = np.arange(-26, 27)
    formats = [("posit16, es=1", 16, 1, "#d1495b"),
               ("posit16, es=2", 16, 2, "#e6a817"),
               ("posit32, es=2", 32, 2, "#3a7ca5")]
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    for label, n, es, c in formats:
        ax.plot(scales, [frac(n, es, s) for s in scales], color=c, lw=2, label=label)

    bands = []
    for p in paths:
        bands += band_from_record(p)
    # keep it readable: dedup and cap
    seen, uniq = set(), []
    for b in bands:
        key = (b[1], b[2], b[0])
        if key not in seen:
            seen.add(key); uniq.append(b)
    colors = ["#2e8b57", "#6a4c93", "#404040", "#b5651d", "#1b6ca8", "#8a2be2"]
    for i, (lab, lo, hi) in enumerate(uniq):
        c = colors[i % len(colors)]
        y = 30 - 2.6 * i
        ax.axvspan(lo, hi, color=c, alpha=0.05)
        ax.annotate("", xy=(lo, y), xytext=(hi, y),
                    arrowprops=dict(arrowstyle="<->", color=c, lw=1.5))
        ax.text((lo + hi) / 2, y + 0.5, lab, color=c, ha="center", fontsize=7.5)

    ax.axhline(23, color="#3a7ca5", ls=":", lw=1); ax.text(-25.5, 23.4, "fp32 acc (~7 dec)", color="#3a7ca5", fontsize=8)
    ax.axhline(10, color="#d1495b", ls=":", lw=1); ax.text(-25.5, 10.4, "fp16 acc (~3.3 dec)", color="#d1495b", fontsize=8)
    ax.axvline(0, color="0.7", lw=0.8); ax.text(0.2, 1, "|x|=1", color="0.5", fontsize=8)
    ax.set_xlabel("value magnitude  (scale = log2 |x|)")
    ax.set_ylabel("fraction bits kept  (precision)")
    ax.set_title("Posit tapered accuracy vs. measured GROMACS operand ranges")
    ax.set_ylim(0, 33); ax.set_xlim(-26, 26); ax.grid(alpha=0.25); ax.legend(loc="upper right")
    fig.tight_layout()
    out = "docs/format_fit_nonbonded.png"
    fig.savefig(out, dpi=150)
    print("wrote", out, "from", len(uniq), "operand bands across", len(paths), "records")

if __name__ == "__main__":
    args = sys.argv[1:]
    paths = []
    for a in args:
        paths += sorted(glob.glob(a)) if any(c in a for c in "*?[]") else [a]
    if not paths:
        # sensible default: the two extremes
        paths = ["dataset/records/gromacs.nonbonded_coulomb.benchmem_82k.json",
                 "dataset/records/gromacs.nonbonded_coulomb.benchrib.json"]
    main(paths)
