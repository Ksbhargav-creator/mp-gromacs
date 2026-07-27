#!/usr/bin/env python3
"""
extract_operands.py  --  mp-spice / NGA side-task (GROMACS operand distributions)

Purpose
-------
Read the *numeric operands* that feed GROMACS's force computation directly out
of a compiled run input (.tpr), WITHOUT rebuilding or instrumenting GROMACS.
This is the cheap "first histogram" path (analogue of the LeNet weight/activation
value-distribution slide): the distribution of partial charges q_i, masses m_i,
and -- if a `gmx dump` text file is supplied -- the Lennard-Jones c6/c12
coefficients that a custom posit/tapered-precision format would be sized against.

Two sources, by design:
  * charges, masses, (positions if a structure file is given): pulled from the
    .tpr binary in pure Python via MDAnalysis' TPRParser. No gmx binary needed.
    (Validated against MDAnalysisTests real .tpr files: adk_oplsaa, cobrotoxin.)
  * LJ c6/c12 nonbonded parameter matrix (nbfp): parsed from the text output of
    `gmx dump -s system.tpr` (see parse_gmx_dump.py). This is the authoritative
    route for LJ operands; MDAnalysis does not expose the full nbfp matrix.
    On the instrumented build, the same nbfp matrix is also dumped at forcerec
    init (see instrumentation/), which is the preferred cross-check.

Output
------
<out>/operands.npz     raw arrays (charge, mass, |q_i*q_j| sample, c6, c12, ...)
<out>/summary.json     per-operand DynamicRangeStats (mirrors the SPICE side)

The DynamicRangeStats block intentionally mirrors the fields recorded by
include/sw/mp_spice/klu_study.hpp on the SPICE side (min_abs, max_abs, decades,
log2 regime histogram) so the two apps can be compared on the same axes.
"""
import argparse, json, os, sys
import numpy as np


def dynamic_range_stats(values, name, nbins_log2=64):
    """Mirror of the SPICE-side DynamicRangeStats / product_magnitude_stats.

    Records the magnitude spread of a set of operands the way a posit regime
    selection cares about: absolute min/max of nonzero entries, span in decades
    and in binades (log2), plus a log2(|x|) histogram (the "regime histogram").
    """
    v = np.asarray(values, dtype=np.float64).ravel()
    finite = v[np.isfinite(v)]
    nonzero = finite[finite != 0.0]
    a = np.abs(nonzero)
    stats = {
        "name": name,
        "count": int(finite.size),
        "count_zero": int(np.sum(finite == 0.0)),
        "count_nonzero": int(nonzero.size),
        "signed_min": float(finite.min()) if finite.size else None,
        "signed_max": float(finite.max()) if finite.size else None,
    }
    if a.size:
        amin, amax = float(a.min()), float(a.max())
        stats.update({
            "abs_min": amin,
            "abs_max": amax,
            "decades": float(np.log10(amax / amin)) if amin > 0 else None,
            "binades_log2": float(np.log2(amax / amin)) if amin > 0 else None,
            "pct": {p: float(np.percentile(a, p)) for p in (1, 5, 50, 95, 99)},
        })
        # log2 regime histogram (posit regime field lives here)
        lo, hi = np.floor(np.log2(amin)), np.ceil(np.log2(amax))
        edges = np.arange(lo, hi + 1.0, max(1.0, (hi - lo) / nbins_log2))
        if edges.size < 2:
            edges = np.array([lo, lo + 1.0])
        counts, _ = np.histogram(np.log2(a), bins=edges)
        stats["log2_hist"] = {
            "edges": edges.tolist(),
            "counts": counts.astype(int).tolist(),
        }
    return stats


def sample_pair_products(charges, max_pairs=2_000_000, seed=0):
    """Distribution of |q_i * q_j| -- the operand product the Coulomb kernel forms
    (numerator of q_i q_j / r before the 1/r factor). This previews, from static
    data alone, the product-magnitude distribution that the instrumented kernel
    logs at runtime with the real 1/r included."""
    rng = np.random.default_rng(seed)
    q = np.asarray(charges, dtype=np.float64)
    q = q[q != 0.0]
    n = q.size
    if n < 2:
        return np.empty(0)
    total = n * (n - 1) // 2
    if total <= max_pairs:
        i, j = np.triu_indices(n, k=1)
    else:  # random sample of unordered pairs
        i = rng.integers(0, n, size=max_pairs)
        j = rng.integers(0, n, size=max_pairs)
        keep = i != j
        i, j = i[keep], j[keep]
    return np.abs(q[i] * q[j])


def extract_arrays(tpr, struct=None, gmx_dump=None, max_pairs=2_000_000, do_pairs=True):
    """Pull the operand value arrays out of a .tpr. Returns (arrays, meta).

    Reusable core shared by this script's CLI and the batch driver nga_gromacs.py.
    arrays maps operand name -> np.ndarray of values; meta carries n_atoms + notes.
    """
    try:
        import MDAnalysis as mda
    except ImportError:
        raise SystemExit("MDAnalysis is required: pip install 'MDAnalysis>=2.7'")
    import warnings
    warnings.filterwarnings("ignore")

    u = mda.Universe(tpr, struct) if struct else mda.Universe(tpr)
    ag = u.atoms
    charges = np.asarray(ag.charges, dtype=np.float64)
    masses = np.asarray(ag.masses, dtype=np.float64)
    arrays = {"charge": charges, "mass": masses}
    notes = []

    try:
        arrays["position"] = np.asarray(ag.positions, dtype=np.float64).ravel()
    except Exception:
        notes.append("no coordinates in input (pass a structure file for positions)")

    if do_pairs:
        arrays["coulomb_qq_product"] = sample_pair_products(charges, max_pairs=max_pairs)

    if gmx_dump:
        try:
            from parse_gmx_dump import parse_nbfp
        except ImportError:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from parse_gmx_dump import parse_nbfp
        c6, c12 = parse_nbfp(gmx_dump)
        if c6.size:
            arrays["lj_c6"], arrays["lj_c12"] = c6, c12
        else:
            notes.append("gmx dump supplied but no c6/c12 parsed (check format/version)")

    return arrays, {"n_atoms": int(ag.n_atoms), "notes": notes}


def build_summary(label, tpr, arrays, meta):
    return {
        "label": label,
        "n_atoms": meta["n_atoms"],
        "tpr": os.path.abspath(tpr),
        "notes": meta.get("notes", []),
        "operands": {k: dynamic_range_stats(v, k) for k, v in arrays.items()},
    }


def print_range_table(summary, indent="  "):
    print(f"{indent}operand              count      |min|        |max|      decades")
    print(indent + "-" * 62)
    for k, s in summary["operands"].items():
        if "abs_min" in s:
            print(f"{indent}{k:<18} {s['count_nonzero']:>8}  {s['abs_min']:.3e}  "
                  f"{s['abs_max']:.3e}   {s['decades']:>6.2f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tpr", required=True, help="GROMACS .tpr run input")
    ap.add_argument("--struct", default=None,
                    help="optional coordinate file (.gro/.pdb/.trr) for positions")
    ap.add_argument("--gmx-dump", default=None,
                    help="optional text output of `gmx dump -s system.tpr` for c6/c12")
    ap.add_argument("--out", default="operands_out", help="output directory")
    ap.add_argument("--label", default=None, help="system label (e.g. benchMEM)")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    label = args.label or os.path.splitext(os.path.basename(args.tpr))[0]

    arrays, meta = extract_arrays(args.tpr, args.struct, args.gmx_dump)
    print(f"[{label}] n_atoms={meta['n_atoms']}  "
          f"nonzero_charges={np.count_nonzero(arrays['charge'])}")
    for n in meta["notes"]:
        print(f"  ({n})")

    np.savez_compressed(os.path.join(args.out, "operands.npz"), **arrays)
    summary = build_summary(label, args.tpr, arrays, meta)
    with open(os.path.join(args.out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    print()
    print_range_table(summary)
    print(f"\n  wrote {args.out}/operands.npz and summary.json")


if __name__ == "__main__":
    main()
