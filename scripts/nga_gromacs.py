#!/usr/bin/env python3
"""
nga_gromacs.py -- operand-distribution study driver (GROMACS side of the NGA work).

The GROMACS counterpart of the mp-spice study harness: point it at one or many
.tpr files and it pulls the numeric operands straight out of each (charges,
masses, coordinates, |q_i q_j| Coulomb products, and LJ c6/c12 when a gmx dump is
provided), records DynamicRangeStats per operand, writes per-system results, and
emits a cross-system comparison table -- all with no GROMACS build required
(pure-Python .tpr read via MDAnalysis).

Examples
--------
    # single system
    ./nga_gromacs.py benchMEM.tpr

    # batch: several systems + a whole directory, with plots + comparison
    ./nga_gromacs.py benchMEM.tpr benchRIB.tpr systems/ --plot --out study1

    # include LJ c6/c12 by pointing at matching `gmx dump` text files
    #   (files matched to each .tpr by basename, e.g. benchMEM.dump.txt)
    ./nga_gromacs.py *.tpr --gmx-dump-dir dumps/ --plot

Outputs (under --out, default: nga_study/)
    <label>/operands.npz     raw operand value arrays
    <label>/summary.json     DynamicRangeStats per operand
    <label>/*.png            histograms (if --plot)
    comparison.csv           dynamic range (decades) per operand, per system
"""
import argparse, glob, json, os, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from extract_operands import extract_arrays, build_summary, print_range_table

HERE = os.path.dirname(os.path.abspath(__file__))
OPERAND_ORDER = ["charge", "mass", "position", "coulomb_qq_product", "lj_c6", "lj_c12"]


def resolve_tprs(patterns):
    """Expand file paths, globs, and directories into a flat list of .tpr files."""
    out = []
    for p in patterns:
        if os.path.isdir(p):
            out += sorted(glob.glob(os.path.join(p, "*.tpr")))
        elif any(c in p for c in "*?[]"):
            out += sorted(glob.glob(p))
        else:
            out.append(p)
    seen, uniq = set(), []
    for t in out:
        if t not in seen and os.path.exists(t):
            seen.add(t); uniq.append(t)
        elif not os.path.exists(t):
            print(f"  ! skipping missing file: {t}")
    return uniq


def match_sidecar(tpr, directory, suffixes):
    """Find a matching struct/dump file for a .tpr by shared basename."""
    if not directory:
        return None
    stem = os.path.splitext(os.path.basename(tpr))[0]
    for suf in suffixes:
        cand = os.path.join(directory, stem + suf)
        if os.path.exists(cand):
            return cand
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tprs", nargs="+", help=".tpr files, globs, or directories")
    ap.add_argument("--out", default="nga_study", help="output directory")
    ap.add_argument("--struct-dir", default=None,
                    help="dir of coordinate files (matched to each .tpr by basename)")
    ap.add_argument("--gmx-dump-dir", default=None,
                    help="dir of `gmx dump` text files (basename-matched) for c6/c12")
    ap.add_argument("--plot", action="store_true", help="also render histograms")
    ap.add_argument("--no-pairs", action="store_true",
                    help="skip the |q_i q_j| pairwise-product sampling (faster)")
    ap.add_argument("--max-pairs", type=int, default=2_000_000)
    args = ap.parse_args()

    tprs = resolve_tprs(args.tprs)
    if not tprs:
        sys.exit("no .tpr files found")
    os.makedirs(args.out, exist_ok=True)

    print("=" * 68)
    print(f" mp-nga-gromacs :: operand-distribution study   ({len(tprs)} system(s))")
    print("=" * 68)

    comparison = []
    for tpr in tprs:
        label = os.path.splitext(os.path.basename(tpr))[0]
        struct = match_sidecar(tpr, args.struct_dir, (".gro", ".pdb", ".trr"))
        dump = match_sidecar(tpr, args.gmx_dump_dir, (".dump.txt", ".dump", ".txt"))
        sysdir = os.path.join(args.out, label)
        os.makedirs(sysdir, exist_ok=True)

        print(f"\n--- {label} ---")
        print(f"  tpr:    {tpr}")
        if struct: print(f"  struct: {struct}")
        if dump:   print(f"  dump:   {dump}")

        arrays, meta = extract_arrays(tpr, struct, dump,
                                      max_pairs=args.max_pairs,
                                      do_pairs=not args.no_pairs)
        for n in meta["notes"]:
            print(f"  ({n})")
        print(f"  n_atoms={meta['n_atoms']}  nonzero_charges={np.count_nonzero(arrays['charge'])}\n")

        np.savez_compressed(os.path.join(sysdir, "operands.npz"), **arrays)
        summary = build_summary(label, tpr, arrays, meta)
        with open(os.path.join(sysdir, "summary.json"), "w") as f:
            json.dump(summary, f, indent=2)
        print_range_table(summary)

        if args.plot:
            subprocess.run([sys.executable, os.path.join(HERE, "plot_histograms.py"),
                            os.path.join(sysdir, "operands.npz"),
                            "--summary", os.path.join(sysdir, "summary.json"),
                            "--out", sysdir, "--label", label], check=False)

        row = {"system": label, "n_atoms": meta["n_atoms"]}
        for k, s in summary["operands"].items():
            row[k] = s.get("decades")
        comparison.append(row)

    # ---- cross-system comparison table (decades of dynamic range) ----
    cols = [c for c in OPERAND_ORDER if any(c in r for r in comparison)]
    csv_path = os.path.join(args.out, "comparison.csv")
    with open(csv_path, "w") as f:
        f.write("system,n_atoms," + ",".join(cols) + "\n")
        for r in comparison:
            f.write(f"{r['system']},{r['n_atoms']}," +
                    ",".join(f"{r.get(c):.3f}" if r.get(c) is not None else ""
                             for c in cols) + "\n")

    print("\n" + "=" * 68)
    print(" CROSS-SYSTEM DYNAMIC RANGE  (decades of |value|)")
    print("=" * 68)
    hdr = f"{'system':<18}{'atoms':>10}  " + "".join(f"{c[:11]:>12}" for c in cols)
    print(hdr); print("-" * len(hdr))
    for r in comparison:
        line = f"{r['system']:<18}{r['n_atoms']:>10}  " + "".join(
            (f"{r.get(c):>12.2f}" if r.get(c) is not None else f"{'-':>12}") for c in cols)
        print(line)
    print(f"\nwrote {csv_path}")
    print(f"per-system results under {args.out}/<system>/")


if __name__ == "__main__":
    main()
