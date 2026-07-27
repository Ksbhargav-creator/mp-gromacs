# mp-gromacs — GROMACS value-distribution / dynamic-range histograms

Build histograms for molecular dynamics — the distribution of the real numbers a
GROMACS force computation actually forms — to make the case for a **custom NGA**
(a posit/tapered-precision format sized to an app's real value distribution rather
than a generic one). 

The pipeline mirrors the SPICE side (`include/sw/mp_spice/klu_study.hpp`,
`DynamicRangeStats` / `product_magnitude_stats`).

## Instrumentation of GROMACS

A `.tpr` is a compiled GROMACS run, not a readable value file like a SPICE `.mtx`.
But that only blocks the *intermediate products* of the force loop — **the numeric
operands are fully recoverable without any rebuild.** So the work splits into two
paths of very different cost:

| Path | What it captures |
|------|------------------|
| **Operands** | charges, masses, LJ c6/c12, coords, `q_iq_j` products |
| **Force-loop products** | runtime `q_iq_j/r`, `c6/r^6`, `c12/r^12` magnitudes |
