#!/usr/bin/env python3
"""
parse_gmx_dump.py -- extract bonded/nonbonded/constraint operands from `gmx dump`.

Generate the input first:
    gmx dump -s system.tpr > system.dump.txt 2>/dev/null

The ffparams section lists one `functype[N]=NAME, key= val, key= val, ...` line per
distinct interaction type. We parse those generically (any NAME, any keys) so it works
across GROMACS versions and force fields, then map the known families to quantities:

    LJ_SR                    -> lj_c6, lj_c12          (nonbonded_lj)
    BONDS/G96BONDS/HARMONIC  -> bond_b0, bond_k        (bonded)
    ANGLES/G96ANGLES/UREY_.. -> angle_theta0, angle_k  (bonded)
    CONSTR/CONSTRNC/SETTLE   -> constraint_length      (lincs_constraint)

If a family isn't picked up, run  `parse_gmx_dump.py --diagnose file.dump.txt`  to see the
exact NAMEs and field keys present, and extend the *_NAMES / field lists below.
"""
import os
import re
import sys
import numpy as np
from collections import defaultdict

FLOAT = r"[-+]?(?:\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|nan|inf)"
FUNC_RE = re.compile(r"functype\[\d+\]\s*=\s*([A-Za-z0-9_]+)\s*,?(.*)")
KV_RE = re.compile(r"([A-Za-z0-9_]+)\s*=\s*(" + FLOAT + r")")


def parse_functypes(path):
    """name -> list of {field: value} dicts, one per functype entry of that name."""
    out = defaultdict(list)
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            m = FUNC_RE.search(line)
            if not m:
                continue
            kv = {k: float(v) for k, v in KV_RE.findall(m.group(2))}
            out[m.group(1)].append(kv)
    return out


def _collect(entries, keys):
    """Pull the first matching key from each entry -> array."""
    vals = []
    for e in entries:
        for k in keys:
            if k in e:
                vals.append(e[k]); break
    return np.array(vals, dtype=np.float64)


BOND_NAMES   = ("BONDS", "G96BONDS", "HARMONIC", "MORSE", "CUBICBONDS", "FENEBONDS")
ANGLE_NAMES  = ("ANGLES", "G96ANGLES", "UREY_BRADLEY", "LINEAR_ANGLES", "RESTRANGLES",
                "CROSS_BOND_BONDS", "QUARTIC_ANGLES")
PDIH_NAMES   = ("PDIHS", "PIDIHS", "FOURDIHS")
IDIH_NAMES   = ("IDIHS",)
CONSTR_NAMES = ("CONSTR", "CONSTRNC")


def parse_nbfp(path):
    """LJ c6/c12 per pair type (nonbonded_lj)."""
    ft = parse_functypes(path)
    c6 = _collect(ft.get("LJ_SR", []), ("c6",))
    c12 = _collect(ft.get("LJ_SR", []), ("c12",))
    return c6, c12


def parse_bonded(path):
    """bonded operands: bond, angle, proper- and improper-dihedral params.
    (Empty arrays are fine -- e.g. no bond_* when constraints=all-bonds converts
    every bond to a CONSTR; that stiffness is then in the lincs_constraint kernel.)"""
    ft = parse_functypes(path)
    out = {k: [] for k in ("bond_b0", "bond_k", "angle_theta0", "angle_k",
                            "dih_k", "dih_phase", "idih_k", "idih_xi0")}

    def take(entry, dst, keys):
        for k in keys:
            if k in entry:
                out[dst].append(entry[k]); return

    for nm in BOND_NAMES:
        for e in ft.get(nm, []):
            take(e, "bond_b0", ("b0A", "b0", "bA"))
            take(e, "bond_k",  ("cbA", "cb", "kb", "kA"))
    for nm in ANGLE_NAMES:
        for e in ft.get(nm, []):
            take(e, "angle_theta0", ("thetaA", "thA", "theta0", "tA"))
            take(e, "angle_k",      ("ctA", "ct", "cthetaA", "kthetaA", "ktheta"))
    for nm in PDIH_NAMES:
        for e in ft.get(nm, []):
            take(e, "dih_k",     ("cpA", "cp", "kphi"))
            take(e, "dih_phase", ("phiA", "phi0", "phi"))
    for nm in IDIH_NAMES:
        for e in ft.get(nm, []):
            take(e, "idih_k",   ("cxA", "cx", "kxi"))
            take(e, "idih_xi0", ("xiA", "xi0", "xi"))
    return {k: np.array(v, dtype=np.float64) for k, v in out.items()}


def parse_constraints(path):
    """constraint_length array (CONSTR + SETTLE d_OH / d_HH)."""
    ft = parse_functypes(path)
    d = []
    for nm in CONSTR_NAMES:
        for e in ft.get(nm, []):
            for k in ("dA", "d", "distA"):
                if k in e: d.append(e[k]); break
    for e in ft.get("SETTLE", []):
        for k in ("doh", "dOH", "d_OH", "dhh", "dHH", "d_HH"):
            if k in e: d.append(e[k])
    return {"constraint_length": np.array(d)}


def diagnose(path):
    ft = parse_functypes(path)
    print(f"functype families in {os.path.basename(path)}  ({sum(len(v) for v in ft.values())} entries):")
    for name, entries in sorted(ft.items()):
        fields = sorted({k for e in entries for k in e})
        print(f"  {name:16} {len(entries):7} entries   fields: {fields}")


def _report(name, a):
    nz = np.abs(a[np.isfinite(a) & (a != 0)])
    if nz.size:
        print(f"  {name:16} {a.size:5}  |min| {nz.min():.3e}  |max| {nz.max():.3e}  "
              f"decades {np.log10(nz.max()/nz.min()):.2f}")
    else:
        print(f"  {name:16} {a.size:5}  (empty)")


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--diagnose" in args:
        for p in [a for a in args if a != "--diagnose"]:
            diagnose(p)
        sys.exit()
    if len(args) != 1:
        sys.exit("usage: parse_gmx_dump.py [--diagnose] system.dump.txt")
    path = args[0]
    c6, c12 = parse_nbfp(path)
    _report("lj_c6", c6); _report("lj_c12", c12)
    for k, v in parse_bonded(path).items():
        _report(k, v)
    for k, v in parse_constraints(path).items():
        _report(k, v)
    if c6.size == 0:
        print("  (no LJ_SR parsed -- run with --diagnose to inspect the functype format)")
