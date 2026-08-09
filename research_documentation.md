# ABSTRACT DRAFT

    In this paper, we curated real number distributions needed for the entire range of some HPC apps. Posits 
    are termed as next generation arithmetic because of its tapered accuracy. Unlike IEEE floats, posits can be 
    customized according to the dynamic range of a HPC app. But when we actually try to port an app, there will be a 
    need for a lot of experiments to customize posits to the app. This is especially a problem without any reference 
    dataset because we will have to start the experiments from scratch each time we port a HPC app. We addressed this 
    problem by creating a motif-indexed dataset. These datasets will characterize the numerical behaviour of the different 
    compute phases of the HPC app. Every phase is indexed by a computational motif(n\_body\_methods, spectral_methods etc...)
    and it's operands, intermediate products, accumulation length and condition number of summation are recorded by instrumentation.
    A arithmetic format is recommended according to the instrumented data. In this paper, we worked on three HPC applications - GROMACS,
    CHARMM GUI and SPICE. With this dataset we will have more of a straightforward process compared to the prior trial and error approach.



# What GROMACS is

GROMACS is a molecular-dynamics engine. It simulates how a collection of atoms
physically moves over time: each step it computes the force on every atom, advances
the atoms forward by a tiny timestep, then recomputes and repeats. After millions of
these force loops we obtain a trajectory of the molecules moving. The benchmark
`.tpr` files (https://www.mpinat.mpg.de/grubmueller/bench) contain each atom's
position, mass, partial charge (q_i), and the force-field parameters that define how
the atoms interact.

**Operands.** The numeric quantities that feed the force computation are the per-atom
partial charges q_i (units of e), masses m_i (amu), coordinates (nm), and the
Lennard-Jones coefficients c6/c12 per atom-type pair. GROMACS's unit system (nm, amu,
e, kJ/mol) is molecular-scale, so these operands are naturally O(1) in magnitude.

**The math (force field).** The total potential energy is a sum of *bonded* terms
(bonds, angles, dihedrals) and *nonbonded* terms between atom pairs: Lennard-Jones
`V = c12/r^12 − c6/r^6` (van der Waals) and Coulomb `V = f·q_i·q_j / r`
(electrostatics). The force on an atom is the negative gradient of this energy.

**Numerical kernels.** The expensive kernels are 
(1) the short-range nonbonded kernels
— the nbnxm cluster/Verlet-list kernels that evaluate LJ + short-range Coulomb over
neighbour pairs — and (2) PME (Particle Mesh Ewald), an FFT-based Poisson solver for
long-range electrostatics (the one stage resembling a linear solver). Time integration
uses the leap-frog Verlet integrator (a = F/m, then update velocity and position), and
constraint solvers (LINCS, SETTLE for water) hold the fastest bonds rigid so a larger
timestep stays stable. Note GROMACS already runs in *mixed precision* by default:
single-precision positions/velocities/forces, double-precision energy/virial
accumulation — so the field already accepts that precision matters here.


# What CHARMM-GUI is

CHARMM-GUI (https://www.charmm-gui.org) is a web-based **input generator** — not a
simulation engine and not a dataset. Given a structure (e.g. a PDB), it builds a
complete simulation system (topology, coordinates, and force-field files) for a chosen
force field, commonly CHARMM36m, and exports it in the input formats of several MD
engines (GROMACS, NAMD, AMBER, OpenMM, CHARMM). Because it performs no runtime
arithmetic, there is nothing to instrument *inside* CHARMM-GUI. To obtain CHARMM-side
value distributions we use CHARMM-GUI to generate a CHARMM36m system and then run that
system through the same pipeline — either by exporting GROMACS format (`grompp` →
`.tpr`) or by reading the CHARMM-native PSF topology directly (MDAnalysis reads PSF
charges and masses). "GROMACS vs CHARMM" is therefore a comparison of two force fields
under identical tooling.


# What we instrumented in GROMACS

We haven't instrumented it yet in the traditional sense, but extracted the operands instead.
We extract the operands directly from each compiled `.tpr` run input using a
pure-Python reader (MDAnalysis) — no GROMACS rebuild required. For each system we
record the partial charges q_i, atomic masses m_i, and a sampled distribution of the
pairwise Coulomb operand |q_i·q_j| (the product the Coulomb kernel forms), compute
per-operand dynamic-range statistics, and plot a log-magnitude histogram for each —
the HPC analogue of the LeNet/MNIST value-distribution plot. (LJ c6/c12 panels are
available when a `gmx dump` of the tpr is supplied; kernel-level instrumentation to
record runtime force-loop products is built but not yet run.)

The distributions show the operands lie inside posit's high-precision zone. On
benchMEM, ~100% of charges and masses and ~98% of Coulomb products fall within
posit16's golden zone (|x| ∈ [1/16, 16]), with medians near unity (charge 0.41,
mass 1.008, |q_i·q_j| 0.29).

**Experimentation results** (dynamic range in decades = log10(max/min) of nonzero |x|):

| System | Class | N atoms | charge | mass | Coulomb product |
|---|---|---|---|---|---|
| benchMEM | Kutzner std MD — membrane protein | 81,743 | 1.53 | 1.55 | 2.96 |
| benchRIB | Kutzner std MD — ribosome in water | 2,136,412 | 3.46 | 1.81 | 6.48 |
| benchPEP | Kutzner std MD — peptides in water | 12,495,503 | 2.51 | 1.55 | 2.89 |
| benchPEP-h | Kutzner std MD — peptides (H-mass repart.) | 12,495,503 | 2.51 | 1.55 | 2.89 |
| cmet_eq | Kutzner binding-affinity (FEP, equilibration) | 67,291 | 3.15 | 1.55 | 5.71 |
| cmet_ti | Kutzner binding-affinity (FEP, thermodynamic integ.) | 67,291 | 3.15 | 1.55 | 5.71 |

**Findings.** 

(1) Across the whole size ladder (82K → 12.5M atoms) operand ranges stay
bounded (mass ~1.5–1.8 decades, charge ~1.5–3.5) and centred near unity — a good posit
fit at every scale. 

(2) The range is a property of the system / force field, not of run size or variant:
benchPEP and benchPEP-h are identical, cmet_eq and cmet_ti are identical, and the 12.5M 
benchPEP (charge 2.51, Coulomb 2.89) is actually *narrower* than the 2.1M benchRIB (3.46, 6.48).

(3) The widest ranges appear in the ribosome (benchRIB, Coulomb 6.48 decades — highly charged 
nucleic-acid backbone) and the binding-affinity systems (cmet, 5.71) — the wide-product regime
where a single narrow posit cannot suffice and a quire / mixed scheme matters most. 

(4) The same molecule spans 2.3 decades of charge under OPLS vs 1.6 under CHARMM — 
direct evidence for a per-force-field tailored format. Caveat: these are *operand* 
ranges; the runtime *products* (LJ r^-12 tower, Coulomb q_i·q_j/r) span even wider 
and are where mixed precision will actually be tested.

To actually instrument GROMACS, we need to have it on our local machine. I downloaded
the tar ball and tried to build GROMACS. But I got an error since Mac compiler doesn't 
have OpenMP, so I downloaded it without OpenMP(`-DMGX_OPENMP=OFF`).(This turned out to
be a mistake) 

Now we need to decide where we need to hook the instrumentation in the kernels. I copied
instrumentation/nga_range_stats.hpp to gromacs src code at it's root as we would need it
to instrument the kernels. A brief overview of nga_range_stats.hpp - It's similar to dynamic_
range.hpp in SPICE and records instrumented values. 

GROMACS's real type is float or double depending on the GMX_DOUBLE CMake option, so the same
instrumented source compiled twice will give us fp32 and fp64 arithmetic as the engine performs
it. From these values, we can make our value-distribution histograms.

As for the Instrumentation method, it starts with two header files - nga_range_stats.hpp and 
nga_trace.hpp. A brief overview on these header files

| Header              | Flag           | Macro                                  | Purpose                                                                                                                               |
|---------------------|----------------|----------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| nga_range_stats.hpp | NGA_INSTRUMENT | NGA_RECORD(name, value),NGA_FLUSH(stem)| aggregate per-quantity dynamic-range stats (min/max/decades + log2 regime histogram) |
| nga_trace.hpp       | NGA_TRACE      |NGA_TRACE_ROW(i,j,qq,rinv,c6,c12), NGA_TRACE_CLOSE()|one CSV row per interaction; input to the paired fp32/fp64 error harness|

# What we instrumented in CHARMM GUI

Not yet run through a CHARMM-GUI-generated system. Preliminary result: Claude verified the
same pipeline ingests CHARMM-native PSF topology directly — a CHARMM adenylate-kinase
PSF gave 1.60 / 1.50 / 3.20 decades for charge / mass / Coulomb product (the CHARMM row
above), confirming CHARMM operands reach the pipeline unchanged. Planned next step:
generate a CHARMM36m system in CHARMM-GUI (Membrane or Solution Builder), export it in
GROMACS format (`grompp` → `.tpr`) or read the PSF directly, and run the identical
driver; the only missing adapter is a CHARMM parameter-file (ε/Rmin → c6/c12) parser
for the LJ panels.


# Mixed-precision plan

Future direction: Build a trace and replay offline system. Use the instrumentation
(nga_range_stats.hpp) we already have and record the operands and products of the force
loops in double precision(the ground truth). We will then replay these operations using 
the recorded inputs offline in a posit+quire library(universal) and we can compare this
with the ground truth. This is easier and a much more approachable step than rewriting
GROMACS kernels for posit+quire.
