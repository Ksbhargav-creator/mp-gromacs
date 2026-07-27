#!/usr/bin/env python3
"""
parse_gmx_dump.py  --  extract Lennard-Jones c6/c12 operands from `gmx dump`.

Usage on your research machine (needs the gmx binary from the build):
    gmx dump -s system.tpr > system.dump.txt 2>/dev/null
    python3 parse_gmx_dump.py system.dump.txt

Why this route: MDAnalysis reads charges/masses straight from the .tpr, but it
does not expose the full nonbonded parameter matrix (nbfp). `gmx dump` prints the
ffparams block, whose leading LJ_SR functype entries ARE the atnr x atnr nbfp
matrix of per-type-pair (c6, c12) coefficients -- the LJ operands a posit taper
would be sized against.

Format note (verify against YOUR GROMACS version once):
GROMACS prints these as lines resembling
    functype[0]=LJ_SR, c6= 2.61587e-03, c12= 2.63795e-06
The parser below is deliberately format-tolerant: it captures any line that
carries both a `c6=` and a `c12=` float. If your version differs, print a few
lines of the dump (`grep -m3 c6 system.dump.txt`) and adjust FLOAT/LINE_RE.
The instrumented build also dumps nbfp at forcerec init (see instrumentation/),
which is the version-proof cross-check for this parser.
"""
import re
import sys
import numpy as np

FLOAT = r"[-+]?(?:\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?|nan|inf)"
LINE_RE = re.compile(r"c6\s*=\s*(" + FLOAT + r").*?c12\s*=\s*(" + FLOAT + r")")


def parse_nbfp(path):
    """Return (c6_array, c12_array) for every LJ pair entry in a gmx dump file."""
    c6, c12 = [], []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if "c6" not in line or "c12" not in line:
                continue
            m = LINE_RE.search(line)
            if m:
                try:
                    c6.append(float(m.group(1)))
                    c12.append(float(m.group(2)))
                except ValueError:
                    pass
    return np.array(c6, dtype=np.float64), np.array(c12, dtype=np.float64)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: parse_gmx_dump.py system.dump.txt")
    c6, c12 = parse_nbfp(sys.argv[1])
    print(f"parsed {c6.size} LJ pair entries")
    for name, a in (("c6", c6), ("c12", c12)):
        nz = np.abs(a[a != 0])
        if nz.size:
            print(f"  {name}: |min| {nz.min():.3e}  |max| {nz.max():.3e}  "
                  f"decades {np.log10(nz.max()/nz.min()):.2f}")
