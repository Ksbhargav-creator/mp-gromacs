#!/usr/bin/env python3
"""Fraction bits (precision) vs value magnitude for several posit formats, with the
measured GROMACS nonbonded-Coulomb operand ranges overlaid. Shows the data sits where
posit precision peaks, and how much precision each format keeps across the range."""
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def regime_bits(k):
    return k + 2 if k >= 0 else -k + 1

def frac(n, es, s):
    return max(0, n - 1 - regime_bits(math.floor(s / (2**es))) - es)

scales = np.arange(-26, 27)
formats = [("posit16, es=1", 16, 1, "#d1495b"),
           ("posit16, es=2", 16, 2, "#e6a817"),
           ("posit32, es=2", 32, 2, "#3a7ca5")]

fig, ax = plt.subplots(figsize=(9.5, 5.2))
for label, n, es, c in formats:
    ax.plot(scales, [frac(n, es, s) for s in scales], color=c, lw=2, label=label)

# measured operand ranges (scale = log2|x|), from the records
bands = [("benchMEM charge",        -5,  1, "#2e8b57"),
         ("benchMEM |q_i q_j|",     -9,  2, "#6a4c93"),
         ("benchRIB |q_i q_j|",    -21,  2, "#404040")]
for i, (lab, lo, hi, c) in enumerate(bands):
    y = 30 - 3*i
    ax.axvspan(lo, hi, ymin=0, ymax=1, color=c, alpha=0.06)
    ax.annotate("", xy=(lo, y), xytext=(hi, y),
                arrowprops=dict(arrowstyle="<->", color=c, lw=1.6))
    ax.text((lo+hi)/2, y+0.6, lab, color=c, ha="center", fontsize=8)

ax.axhline(23, color="#3a7ca5", ls=":", lw=1); ax.text(-25.5, 23.4, "fp32 acc (~7 dec)", color="#3a7ca5", fontsize=8)
ax.axhline(10, color="#d1495b", ls=":", lw=1); ax.text(-25.5, 10.4, "fp16 acc (~3.3 dec)", color="#d1495b", fontsize=8)
ax.axvline(0, color="0.7", lw=0.8); ax.text(0.2, 1, "|x|=1", color="0.5", fontsize=8)

ax.set_xlabel("value magnitude  (scale = log2 |x|)")
ax.set_ylabel("fraction bits kept  (precision)")
ax.set_title("Posit tapered accuracy vs. GROMACS nonbonded-Coulomb operand ranges")
ax.set_ylim(0, 33); ax.set_xlim(-26, 26); ax.grid(alpha=0.25); ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig("docs/format_fit_nonbonded.png", dpi=150)
print("wrote docs/format_fit_nonbonded.png")
