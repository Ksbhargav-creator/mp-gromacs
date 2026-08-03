#!/usr/bin/env python3
"""
format_fit.py -- from a record's measured value range, compute the posit bit budget
(regime / exponent / fraction) needed, and recommend a format.

Posit<n,es> layout:  sign(1) | regime(variable) | exponent(<=es) | fraction(rest)
  value ~ useed^k * 2^e * (1+frac),  useed = 2^(2^es),  scale = 2^es * k + e
  - the regime+exponent fields encode the binary SCALE (exponent) of the value
  - fraction bits give precision; they are whatever remains after regime+exponent
  - near |x|=1 the regime is shortest (2 bits) so precision peaks (tapered accuracy)

Given a measured operand spanning [abs_min, abs_max]:
  s_min = floor(log2 abs_min),  s_max = ceil(log2 abs_max)     # binary scales
  span  = s_max - s_min   (== binades_log2)                    # RANGE requirement
Regime length at scale s:  k = floor(s / 2^es);  R = k+2 (k>=0) or -k+1 (k<0)
Fraction bits at scale s:  f(s) = n - 1 - R - es
Guarantee f_target bits across the whole range -> use the WORST (longest-regime) end:
  n_min = f_target + 1 + es + max(R(s_min), R(s_max))
Decimals of accuracy (Gustafson):  d ~ (f+1) * log10(2)
"""
import json, math, sys, glob, os
LOG10_2 = 0.30102999566

def regime_bits(k):
    return k + 2 if k >= 0 else -k + 1

def scales(abs_min, abs_max):
    return math.floor(math.log2(abs_min)), math.ceil(math.log2(abs_max))

def frac_at(n, es, s):
    k = math.floor(s / (2 ** es))
    return n - 1 - regime_bits(k) - es

def min_n(f_target, es, s_min, s_max):
    R = max(regime_bits(math.floor(s_min/(2**es))),
            regime_bits(math.floor(s_max/(2**es))))
    return f_target + 1 + es + R

def decimals(f):
    return (f + 1) * LOG10_2

def analyse(name, abs_min, abs_max):
    s_min, s_max = scales(abs_min, abs_max)
    span = s_max - s_min
    print(f"\n  {name}: |x| in [{abs_min:.3e}, {abs_max:.3e}]  "
          f"scales [{s_min}, {s_max}]  span = {span} binades")
    # posit needs its exponent reach to cover the range:  (2^es)*(n-2) >= max|s|
    for es in (1, 2, 3):
        row = []
        for f_target, tag in ((10, "fp16-acc ~3.3dec"), (23, "fp32-acc ~7.2dec")):
            n = min_n(f_target, es, s_min, s_max)
            row.append(f"{tag}: n>={n}")
        # precision a standard posit16/32 would actually deliver at the worst scale
        f16 = min(frac_at(16, es, s_min), frac_at(16, es, s_max))
        f32 = min(frac_at(32, es, s_min), frac_at(32, es, s_max))
        print(f"    es={es}:  " + " | ".join(row)
              + f"   [posit16 worst {max(f16,0)}f/{decimals(max(f16,0)):.1f}dec,"
                f" posit32 worst {max(f32,0)}f/{decimals(max(f32,0)):.1f}dec]")

def best_posit(abs_min, abs_max, f_target):
    """Smallest posit<n,es> that keeps >= f_target fraction bits across the range
    AND whose exponent reach covers the range. Returns (n, es)."""
    s_min, s_max = scales(abs_min, abs_max)
    best = None
    for es in (0, 1, 2, 3):
        n = min_n(f_target, es, s_min, s_max)
        reach = (2 ** es) * (n - 2)                       # max |scale| the format spans
        if reach < max(abs(s_min), abs(s_max)):
            continue
        if best is None or n < best[0]:
            best = (n, es)
    return best

def from_record(path):
    r = json.load(open(path))
    print("="*88); print(f" {r['id']}   (kernel={r['kernel']['name']}, motif={r['motif']})")
    print("="*88)
    for q in r["quantities"]:
        s = q["stats"]
        if "abs_min" in s and s["abs_min"] > 0:
            analyse(f"{q['name']} [{q['role']}]", s["abs_min"], s["abs_max"])

def recommend_table(paths):
    """Markdown table: the format the dataset recommends per (system, quantity)."""
    print("| system | quantity | role | span (binades) | posit (fp16-class, ~3 dec) | posit (fp32-class, ~7 dec) |")
    print("|---|---|---|---|---|---|")
    for p in paths:
        r = json.load(open(p))
        sysname = r["input_case"]["name"]
        for q in r["quantities"]:
            s = q["stats"]
            if "abs_min" not in s or s["abs_min"] <= 0:
                continue
            span = round(s["binades_log2"])
            n16, es16 = best_posit(s["abs_min"], s["abs_max"], 10)
            n32, es32 = best_posit(s["abs_min"], s["abs_max"], 23)
            print(f"| {sysname} | {q['name']} | {q['role']} | {span} "
                  f"| posit&lt;{n16},{es16}&gt; | posit&lt;{n32},{es32}&gt; |")

KERNEL_ORDER = ["nonbonded_coulomb", "nonbonded_lj", "bonded",
                "lincs_constraint", "integration", "pme_reciprocal"]

def recommend_by_kernel(paths):
    """Per-kernel governing format: the widest requirement across all its quantities
    and all input cases (a kernel is ported as one unit, so the hardest quantity wins)."""
    kern = {}
    for p in paths:
        r = json.load(open(p))
        for q in r["quantities"]:
            s = q["stats"]
            if "abs_min" in s and s["abs_min"] > 0:
                kern.setdefault(r["kernel"]["name"], []).append(
                    (s["abs_min"], s["abs_max"], s["binades_log2"]))
    ks = [k for k in KERNEL_ORDER if k in kern] + [k for k in kern if k not in KERNEL_ORDER]
    print("| kernel | n quantities | max span (binades) | fp16-class (~3 dec) | fp32-class (~7 dec) |")
    print("|---|---|---|---|---|")
    for k in ks:
        e = kern[k]
        def gov(f):
            bn, be = 0, 0
            for amn, amx, _ in e:
                n, es = best_posit(amn, amx, f)
                if n > bn: bn, be = n, es
            return bn, be
        n16, es16 = gov(10); n32, es32 = gov(23)
        print(f"| {k} | {len(e)} | {max(x[2] for x in e):.1f} "
              f"| posit&lt;{n16},{es16}&gt; | posit&lt;{n32},{es32}&gt; |")

if __name__ == "__main__":
    args = sys.argv[1:]
    mode = ("--kernels" if "--kernels" in args else
            "--table" if "--table" in args else None)
    args = [a for a in args if a not in ("--table", "--kernels")]
    paths = []
    for a in args:
        paths += sorted(glob.glob(a)) if any(c in a for c in "*?[]") else [a]
    if not paths:
        sys.exit("usage: format_fit.py [--table|--kernels] <record.json> ...")
    if mode == "--kernels":
        recommend_by_kernel(paths)
    elif mode == "--table":
        recommend_table(paths)
    else:
        for p in paths:
            from_record(p)
