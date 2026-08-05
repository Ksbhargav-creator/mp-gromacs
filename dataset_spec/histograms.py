#!/usr/bin/env python3
"""
histograms.py -- value-distribution histograms per kernel and per app, built from the
committed records only (no raw arrays needed).

Method (agreed):
  * axis = log2|x| (binary scale); we histogram magnitudes, on a common integer-binade grid
  * each quantity's stored `log2_hist` is re-binned onto that common grid
  * each quantity is normalized to a DENSITY (area = 1) so sample-count differences
    (per-atom charges vs per-type c6/c12) don't bias the picture
  * per kernel  = mean of its quantities' densities (each quantity an equal voice)
  * per app     = mean of its kernels' densities   (each kernel an equal voice)
  * posit representable ranges + the |x|=1 accuracy peak are overlaid, so "the app's
    values lie in posit's dynamic range" is shown directly.

Usage:
    python3 dataset_spec/histograms.py dataset/records/*.json --out docs
"""
import argparse, glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

KERNEL_ORDER = ["nonbonded_coulomb", "nonbonded_lj", "bonded",
                "lincs_constraint", "integration", "pme_reciprocal"]

def rebin(src_edges, src_counts, target_edges):
    """Redistribute src_counts onto target bins by fractional overlap (area-conserving)."""
    out = np.zeros(len(target_edges) - 1)
    te_lo, te_hi = target_edges[:-1], target_edges[1:]
    for c, a, b in zip(src_counts, src_edges[:-1], src_edges[1:]):
        if c == 0 or b <= a:
            continue
        ov = np.clip(np.minimum(b, te_hi) - np.maximum(a, te_lo), 0, None)
        out += c * ov / (b - a)
    return out

def density(counts):
    s = counts.sum()
    return counts / s if s > 0 else counts

def load(paths):
    """-> {(app, kernel): [ (qname, edges, counts) ]}, and global scale bounds."""
    data, lo, hi = {}, np.inf, -np.inf
    for p in paths:
        r = json.load(open(p))
        app, kern = r["app"]["name"], r["kernel"]["name"]
        for q in r["quantities"]:
            h = q["stats"].get("log2_hist")
            if not h or not h["counts"]:
                continue
            e = np.array(h["edges"], float); c = np.array(h["counts"], float)
            data.setdefault((app, kern), []).append((q["name"], e, c))
            lo, hi = min(lo, e[0]), max(hi, e[-1])
    return data, lo, hi

def grid(lo, hi, per_binade):
    lo, hi = np.floor(lo), np.ceil(hi)
    return np.arange(lo, hi + 1.0 / per_binade, 1.0 / per_binade)

def kernel_density(quantities, edges):
    """mean of each quantity's re-binned, normalized density."""
    ds = [density(rebin(e, c, edges)) for (_, e, c) in quantities]
    return np.mean(ds, axis=0) if ds else np.zeros(len(edges) - 1)

def posit_overlay(ax, ymax):
    """shade representable ranges; mark accuracy peak at |x|=1."""
    # max representable scale = (2^es)*(n-2)
    for (n, es, col, lab) in [(16, 1, "#d1495b", "posit16,es1"),
                              (16, 2, "#e6a817", "posit16,es2"),
                              (32, 2, "#3a7ca5", "posit32,es2")]:
        s = (2 ** es) * (n - 2)
        ax.axvline(-s, color=col, ls="--", lw=1, alpha=0.7)
        ax.axvline(s,  color=col, ls="--", lw=1, alpha=0.7)
        ax.text(-s, ymax * 0.9, f" {lab} min", color=col, rotation=90, va="top", fontsize=6.5)
    ax.axvline(0, color="0.5", lw=1)
    ax.text(0.3, ymax * 0.97, "|x|=1  (accuracy peak)", color="0.4", fontsize=7, va="top")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("records", nargs="+")
    ap.add_argument("--out", default="docs")
    ap.add_argument("--per-binade", type=int, default=2, help="bins per binade")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    paths = []
    for a in args.records:
        paths += sorted(glob.glob(a)) if any(c in a for c in "*?[]") else [a]
    data, lo, hi = load(paths)
    edges = grid(lo, hi, args.per_binade)
    centers = 0.5 * (edges[:-1] + edges[1:])

    apps = sorted({app for (app, _) in data})
    # per-kernel density for each (app,kernel)
    kdens = {ak: kernel_density(q, edges) for ak, q in data.items()}

    # ---- Figure 1: per-kernel densities (one app), kernels overlaid ----------
    for app in apps:
        ks = [k for k in KERNEL_ORDER if (app, k) in kdens] + \
             [k for (a, k) in kdens if a == app and k not in KERNEL_ORDER]
        fig, ax = plt.subplots(figsize=(10, 5.2))
        for k in ks:
            ax.plot(centers, kdens[(app, k)], lw=1.8, label=k)
        ymax = ax.get_ylim()[1]
        posit_overlay(ax, ymax)
        ax.set_xlabel("value magnitude  (scale = log2 |x|)")
        ax.set_ylabel("density (per-kernel, area = 1)")
        ax.set_title(f"{app}: value distribution per computational kernel")
        ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.25)
        fig.tight_layout(); fig.savefig(f"{args.out}/hist_kernels_{app}.png", dpi=150)
        print(f"wrote {args.out}/hist_kernels_{app}.png")

    # ---- Figure 2: per-app densities (mean of kernels), apps overlaid --------
    fig, ax = plt.subplots(figsize=(10, 5.2))
    for app in apps:
        ker = [kdens[(app, k)] for (a, k) in kdens if a == app]
        appd = np.mean(ker, axis=0)
        ax.fill_between(centers, appd, alpha=0.25)
        ax.plot(centers, appd, lw=2, label=f"{app} ({len(ker)} kernels)")
    ymax = ax.get_ylim()[1]
    posit_overlay(ax, ymax)
    ax.set_xlabel("value magnitude  (scale = log2 |x|)")
    ax.set_ylabel("density (per-app, area = 1)")
    ax.set_title("Application value distributions vs. posit dynamic range")
    ax.legend(fontsize=9); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(f"{args.out}/hist_apps.png", dpi=150)
    print(f"wrote {args.out}/hist_apps.png  (apps: {', '.join(apps)})")

if __name__ == "__main__":
    main()
