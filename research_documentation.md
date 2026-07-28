# ABSTRACT DRAFT

    Tapered-precision formats such as posits promise accuracy at lower bit-width than
    IEEE floating point, but only when their non-uniform precision matches an
    application's actual operand magnitudes. Many HPC and AI applications have operand
    magnitudes that fall inside posit's high-precision zone, motivating us to test these
    applications under mixed precision. We work with the molecular-dynamics engine
    GROMACS and the input generator CHARMM-GUI.
    
    We built a pipeline that extracts the operand distributions feeding the GROMACS and
    CHARMM force computations — partial charges, atomic masses, and pairwise Coulomb
    products |q_i·q_j| — directly from compiled run inputs, and records their dynamic
    range. On the Kutzner benchMEM benchmark (81,743-atom membrane protein) these
    operands occupy a bounded range (~1.5 decades for charges, ~3 for Coulomb products),
    with nearly 100% of values inside posit16's high-precision zone. Across the full
    Kutzner size ladder (82K–12.5M atoms) the operands stay bounded and centred near
    unity, but the width is set by the system and force field rather than size: the same
    molecule spans 2.3 decades of charge under OPLS vs 1.6 under CHARMM, and highly
    charged systems (ribosome, protein–ligand binding) push the Coulomb products wider
    (to ~6.5 decades) — evidence the distribution is a force-field and chemistry
    property, and thus a candidate for a tailored format. These results indicate the MD force loop is a
    promising target for posit-based mixed precision, our immediate next step.


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
| adk (OPLS) * | validation — adenylate kinase | 47,681 | 2.32 | 1.50 | 3.56 |
| cobrotoxin * | validation | 19,385 | 2.32 | 1.55 | 3.56 |
| adk (CHARMM PSF) * | validation | 3,341 | 1.60 | 1.50 | 3.20 |

\* Validation systems from the MDAnalysis test set, used to validate the pipeline and
illustrate force-field dependence. All six Kutzner benchmarks are now run: the full
standard-MD size ladder (benchMEM 82K → benchRIB 2.1M → benchPEP 12.5M) plus the cmet
protein-ligand binding-affinity set.

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


# What we instrumented in CHARMM GUI

Not yet run through a CHARMM-GUI-generated system. Preliminary result: we verified the
same pipeline ingests CHARMM-native PSF topology directly — a CHARMM adenylate-kinase
PSF gave 1.60 / 1.50 / 3.20 decades for charge / mass / Coulomb product (the CHARMM row
above), confirming CHARMM operands reach the pipeline unchanged. Planned next step:
generate a CHARMM36m system in CHARMM-GUI (Membrane or Solution Builder), export it in
GROMACS format (`grompp` → `.tpr`) or read the PSF directly, and run the identical
driver; the only missing adapter is a CHARMM parameter-file (ε/Rmin → c6/c12) parser
for the LJ panels.


# Mixed-precision plan

*To be written.*
