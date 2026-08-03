#!/usr/bin/env python3
"""
kernel_math_walkthrough.py -- the GROMACS force math, worked by hand, and how it
maps to a dataset record. Deliberately tiny (4 atoms) so every number is checkable.

The recipe this teaches (works for ANY kernel / app):
  1. Write the kernel's math.
  2. Tag every symbol as an OPERAND (an input the kernel reads), a PRODUCT (an
     intermediate it forms), or an ACCUMULATOR (a running sum).
  3. Measure the magnitude distribution of each -> that's what goes in the record.
"""
import numpy as np

# GROMACS unit system: length nm, charge e, mass amu, energy kJ/mol.
ONE_4PI_EPS0 = 138.935458   # kJ mol^-1 nm e^-2  (GROMACS's Coulomb prefactor f)

# ---- a tiny system: 2 "water-like" molecules, 4 atoms ---------------------
# OPERANDS read straight from the topology (.tpr):
names  = ["O1", "H1", "O2", "H2"]
charge = np.array([-0.834, 0.417, -0.834, 0.417])      # q_i, units e   (OPERAND)
mass   = np.array([15.999, 1.008, 15.999, 1.008])      # m_i, units amu (OPERAND)
pos    = np.array([[0.00, 0.00, 0.00],                 # coordinates, nm (OPERAND)
                   [0.09, 0.00, 0.00],
                   [0.30, 0.00, 0.00],
                   [0.39, 0.00, 0.00]])
# LJ params per atom type (OPERANDS). Only O has LJ here (H has none), water-ish:
c6  = {"O": 2.617e-3, "H": 0.0}     # kJ mol^-1 nm^6   (OPERAND)
c12 = {"O": 2.634e-6, "H": 0.0}     # kJ mol^-1 nm^12  (OPERAND)
atype = ["O", "H", "O", "H"]

print("="*74)
print(" ONE PAIR, FULLY WORKED  (check this by hand)")
print("="*74)
i, j = 0, 2                        # O1 -- O2
d   = pos[i] - pos[j]             # displacement vector (nm)
r2  = d @ d                       # squared distance   (PRODUCT)
r   = np.sqrt(r2)                 # distance           (PRODUCT)
rinv = 1.0 / r                    # 1/r                (PRODUCT)
qq  = ONE_4PI_EPS0 * charge[i] * charge[j]     # f*q_i*q_j  (PRODUCT: the Coulomb operand product)
Vcoul = qq * rinv                              # Coulomb energy term (PRODUCT)
fscal_coul = qq * rinv * rinv * rinv           # Coulomb force scalar = f qq / r^3 (multiplies d)
c6ij  = c6[atype[i]]  if atype[i]==atype[j] else np.sqrt(c6[atype[i]]*c6[atype[j]])
c12ij = c12[atype[i]] if atype[i]==atype[j] else np.sqrt(c12[atype[i]]*c12[atype[j]])
rinv2 = rinv*rinv
rinv6 = rinv2**3                               # 1/r^6  (PRODUCT)
Vlj = c12ij*rinv6*rinv6 - c6ij*rinv6           # LJ 12-6 energy term (PRODUCT)
FrLJ = 12.0*c12ij*rinv6*rinv6 - 6.0*c6ij*rinv6 # LJ force prefactor  (PRODUCT)
fscal_lj = FrLJ * rinv2                         # LJ force scalar (multiplies d)

print(f" pair {names[i]}-{names[j]}:")
print(f"   q_i={charge[i]:+.3f} e   q_j={charge[j]:+.3f} e            (operands)")
print(f"   r   = |pos_i - pos_j| = {r:.4f} nm                        (product)")
print(f"   1/r = {rinv:.4f} nm^-1                                    (product)")
print(f"   qq  = f*q_i*q_j = {ONE_4PI_EPS0:.3f} * {charge[i]:+.3f} * {charge[j]:+.3f} = {qq:.4f}  (product)")
print(f"   Coulomb energy  = qq/r          = {Vcoul:.4f} kJ/mol      (product)")
print(f"   Coulomb fscal   = qq/r^3        = {fscal_coul:.4f}         (product)")
print(f"   c6={c6ij:.3e}  c12={c12ij:.3e}                            (operands)")
print(f"   1/r^6 = {rinv6:.4e}                                       (product)")
print(f"   LJ energy = c12/r^12 - c6/r^6   = {Vlj:.4e} kJ/mol        (product)")

# ---- now ALL pairs, and the per-atom force ACCUMULATION --------------------
print("\n" + "="*74)
print(" ALL PAIRS  +  PER-ATOM FORCE ACCUMULATION")
print("="*74)
N = len(names)
force = np.zeros((N, 3))                       # ACCUMULATOR: sum of pair forces per atom
qq_list, coul_list, lj_list = [], [], []
print(f" {'pair':<8}{'r(nm)':>8}{'qq':>10}{'coul E':>10}{'LJ E':>12}")
for i in range(N):
    for j in range(i+1, N):
        d = pos[i]-pos[j]; r = np.sqrt(d@d); rinv=1/r
        qq = ONE_4PI_EPS0*charge[i]*charge[j]
        Vcoul = qq*rinv
        fcoul = qq*rinv**3
        c6ij  = np.sqrt(c6[atype[i]]*c6[atype[j]])
        c12ij = np.sqrt(c12[atype[i]]*c12[atype[j]])
        rinv6 = rinv**6
        Vlj = c12ij*rinv6*rinv6 - c6ij*rinv6
        FrLJ = 12*c12ij*rinv6*rinv6 - 6*c6ij*rinv6
        fpair = (fcoul + FrLJ*rinv**2) * d      # total pair force vector on i (and -d on j)
        force[i] += fpair                        # <-- accumulation onto atom i
        force[j] -= fpair                        # <-- and onto atom j
        qq_list.append(abs(qq)); coul_list.append(abs(Vcoul)); lj_list.append(abs(Vlj))
        print(f" {names[i]}-{names[j]:<5}{r:>8.3f}{qq:>10.3f}{Vcoul:>10.3f}{Vlj:>12.3e}")
print("\n per-atom force magnitude |F_i| (the ACCUMULATOR result):")
for i in range(N):
    print(f"   {names[i]}: |F| = {np.linalg.norm(force[i]):.3f} kJ mol^-1 nm^-1  "
          f"(sum over {N-1} neighbours)")

# ---- how these become a dataset record ------------------------------------
def decades(vals):
    a = np.abs(np.asarray(vals, float)); a = a[a>0]
    return float(np.log10(a.max()/a.min())) if a.size else 0.0

print("\n" + "="*74)
print(" WHAT THE RECORD MEASURES  (magnitude spread of each quantity)")
print("="*74)
rows = [
    ("charge",             "operand",     "nonbonded_coulomb", np.abs(charge)),
    ("mass",               "operand",     "integration",       mass),
    ("coulomb_qq_product", "product",     "nonbonded_coulomb", qq_list),
    ("coulomb_term",       "product",     "nonbonded_coulomb", coul_list),
    ("lj_term",            "product",     "nonbonded_lj",      lj_list),
    ("force_accumulation", "accumulator", "integration",       [np.linalg.norm(f) for f in force]),
]
print(f" {'quantity':<20}{'role':<12}{'kernel':<20}{'decades':>8}")
for name, role, kern, vals in rows:
    print(f" {name:<20}{role:<12}{kern:<20}{decades(vals):>8.2f}")
print("\n Each row above -> one 'quantity' in a record; the record's `kernel` field")
print(" groups the rows that belong to the same kernel; `motif` tags that kernel.")
