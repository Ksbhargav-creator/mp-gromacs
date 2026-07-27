#!/usr/bin/env python3
"""
plot_histograms.py  --  render GROMACS operand value-distribution histograms.

Produces the HPC analogue of the LeNet/MNIST value-distribution slide: for each
operand (partial charge, mass, Coulomb q_i*q_j product, LJ c6/c12, and -- from
the instrumented build -- runtime force-loop product magnitudes) it plots the
distribution of log10(|value|). The x-axis is decades of magnitude; annotated
guide bands show the exponent reach of the number formats under consideration
(IEEE float32, float16/bfloat16, and a posit<16,1>/<32,2> tapered range) so the
plot directly makes the "custom NGA" case: where the mass of the distribution
sits vs. what each format can represent without over/underflow.

Usage:
    python3 plot_histograms.py operands.npz --summary summary.json --out figs/
"""
import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Approximate representable |value| reach (decades of normalized magnitude).
FORMAT_RANGES = {
    "float32":       (-37.8, 38.5,  "#3a7ca5"),
    "float16":       (-4.9,  4.8,   "#d1495b"),
    "bfloat16":      (-37.8, 38.5,  "#8a5a44"),
    "posit<16,1>":   (-8.4,  8.4,   "#2e8b57"),
    "posit<32,2>":   (-36.7, 36.7,  "#6a4c93"),
}
PRETTY = {
    "charge": "Partial charge q_i  [e]",
    "mass": "Atomic mass m_i  [amu]",
    "position": "Coordinate  [nm]",
    "coulomb_qq_product": "Coulomb operand |q_i q_j|  [e^2]",
    "lj_c6": "Lennard-Jones c6",
    "lj_c12": "Lennard-Jones c12",
    "product_magnitude": "Force-loop product magnitude (instrumented)",
    "lj_term": "LJ term magnitude (instrumented)",
    "coulomb_term": "Coulomb term magnitude (instrumented)",
}


def _log10_abs(v):
    v = np.asarray(v, dtype=np.float64).ravel()
    v = v[np.isfinite(v)]
    v = np.abs(v[v != 0.0])
    return np.log10(v) if v.size else np.empty(0)


def _panel(ax, values, title, show_formats=True):
    lg = _log10_abs(values)
    if lg.size == 0:
        ax.text(0.5, 0.5, "no nonzero data", ha="center", va="center",
                transform=ax.transAxes, color="0.5")
        ax.set_title(title, fontsize=10)
        return
    lo, hi = np.floor(lg.min()), np.ceil(lg.max())
    bins = np.linspace(lo, hi, max(30, int((hi - lo) * 6)))
    ax.hist(lg, bins=bins, color="#404040", alpha=0.85, edgecolor="none")
    dr = lg.max() - lg.min()
    ax.set_title(f"{title}\ndynamic range = {dr:.1f} decades  (n={lg.size:,})",
                 fontsize=9)
    ax.set_xlabel("log10 |value|"); ax.set_ylabel("count")
    ax.grid(True, alpha=0.2)
    if show_formats:
        ymax = ax.get_ylim()[1]
        for i, (name, (flo, fhi, col)) in enumerate(FORMAT_RANGES.items()):
            # only draw the underflow edge that is near the data
            if flo > lo - 3:
                ax.axvline(flo, color=col, ls="--", lw=1, alpha=0.7)
                ax.text(flo, ymax * (0.95 - 0.09 * i), f" {name} min",
                        color=col, fontsize=6.5, rotation=90, va="top")


def ordered_keys(data):
    keys = [k for k in ("charge", "mass", "position", "coulomb_qq_product",
                        "lj_c6", "lj_c12", "product_magnitude",
                        "lj_term", "coulomb_term") if k in data.files]
    keys += [k for k in data.files if k not in keys]
    return keys


def plot_single(values, title, label, out_path):
    """One operand, one figure -- sized for a single slide."""
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    _panel(ax, values, title)
    fig.suptitle(f"{label}", fontsize=12, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print("wrote", out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz", help="operands.npz from extract_operands.py")
    ap.add_argument("--summary", default=None, help="summary.json (for the title label)")
    ap.add_argument("--out", default="figs", help="output directory for PNGs")
    ap.add_argument("--label", default=None)
    ap.add_argument("--separate", action="store_true",
                    help="write one PNG per operand (for one-per-slide) instead of a grid")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    data = np.load(args.npz)
    label = args.label
    if args.summary and os.path.exists(args.summary):
        with open(args.summary) as f:
            label = label or json.load(f).get("label")
    label = label or os.path.basename(os.path.dirname(os.path.abspath(args.npz)))

    keys = ordered_keys(data)

    if args.separate:
        for k in keys:
            out = os.path.join(args.out, f"operand_histogram_{label}_{k}.png")
            plot_single(data[k], PRETTY.get(k, k), label, out)
        return

    ncol = 2
    nrow = int(np.ceil(len(keys) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.2 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, k in zip(axes, keys):
        _panel(ax, data[k], PRETTY.get(k, k))
    for ax in axes[len(keys):]:
        ax.axis("off")
    fig.suptitle(f"GROMACS value-distribution / dynamic range  --  {label}",
                 fontsize=13, y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    out = os.path.join(args.out, f"operand_histograms_{label}.png")
    fig.savefig(out, dpi=140)
    print("wrote", out)


if __name__ == "__main__":
    main()
