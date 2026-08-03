#!/usr/bin/env python3
"""
summary_to_records.py -- convert operand summary.json files into dataset records.

Each summary.json (produced by nga_gromacs.py) mixes operands that belong to
*different* computational kernels: partial charge and the q_i q_j product belong
to the nonbonded-Coulomb kernel, c6/c12 to the nonbonded-LJ kernel, mass to the
integration (leap-frog update) kernel. The dataset unit is (kernel, input), so
this tool splits one summary.json into ONE record per computational kernel,
attaching the schema metadata (app, kernel, motif, input_case, provenance) and
carrying each operand's DynamicRangeStats across unchanged.

Only the layers that static operand extraction can honestly fill are populated:
metadata + quantities. The accumulation and format_fit layers are left absent
(a note records why) because they require the instrumented-trace and replay steps.

Usage:
    python3 summary_to_records.py <summary.json | dir | glob> ... [options]
    python3 summary_to_records.py ../study_benchMEM --out ../dataset/records

Options:
    --app-version STR   GROMACS version string (default 2025.3)
    --source-url URL    input-case provenance URL (default: Kutzner benchmark page)
    --license STR       input-case license (default CC-BY-4.0)
    --out DIR           output directory (default: dataset/records)
    --no-validate       skip JSON-Schema validation of emitted records
"""
import argparse, datetime, glob, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMA_PATH = os.path.join(HERE, "kernel_characterization.schema.json")

# quantity name -> (kernel, role, units)
OPERAND_MAP = {
    # nonbonded-Coulomb (static + runtime)
    "charge":             ("nonbonded_coulomb", "operand", "e"),
    "coulomb_qq_product": ("nonbonded_coulomb", "product", "e^2"),
    "rinv":               ("nonbonded_coulomb", "product", "nm^-1"),
    "coulomb_term":       ("nonbonded_coulomb", "product", "kJ mol^-1"),
    # nonbonded-LJ
    "lj_c6":              ("nonbonded_lj",       "operand", "kJ mol^-1 nm^6"),
    "lj_c12":             ("nonbonded_lj",       "operand", "kJ mol^-1 nm^12"),
    "lj_rinv6":           ("nonbonded_lj",       "product", "nm^-6"),
    "lj_term6":           ("nonbonded_lj",       "product", "kJ mol^-1"),
    "lj_term12":          ("nonbonded_lj",       "product", "kJ mol^-1"),
    "lj_force":           ("nonbonded_lj",       "product", "kJ mol^-1 nm^-1"),
    # PME reciprocal
    "pme_qgrid":          ("pme_reciprocal",     "operand", "e"),
    "pme_influence":      ("pme_reciprocal",     "product", "1"),
    "pme_struct_factor":  ("pme_reciprocal",     "product", "e^2"),
    "pme_recip_term":     ("pme_reciprocal",     "product", "kJ mol^-1"),
    "pme_recip_energy":   ("pme_reciprocal",     "accumulator", "kJ mol^-1"),
    # bonded  (static operands from gmx dump: k, b0, theta0; runtime: dr, energy, force)
    "bond_k":             ("bonded",             "operand", "kJ mol^-1 nm^-2"),
    "bond_b0":            ("bonded",             "operand", "nm"),
    "angle_k":            ("bonded",             "operand", "kJ mol^-1 rad^-2"),
    "angle_theta0":       ("bonded",             "operand", "deg"),
    "dih_k":              ("bonded",             "operand", "kJ mol^-1"),
    "dih_phase":          ("bonded",             "operand", "deg"),
    "idih_k":             ("bonded",             "operand", "kJ mol^-1 rad^-2"),
    "idih_xi0":           ("bonded",             "operand", "deg"),
    "ub_k":               ("bonded",             "operand", "kJ mol^-1 nm^-2"),  # CHARMM Urey-Bradley 1-3
    "ub_r0":              ("bonded",             "operand", "nm"),
    # CMAP -- CHARMM phi/psi backbone correction map (grid interpolation)
    "cmap_grid":          ("cmap",               "operand", "kJ mol^-1"),
    "bond_dr":            ("bonded",             "product", "nm"),
    "bond_energy":        ("bonded",             "product", "kJ mol^-1"),
    "angle_dtheta":       ("bonded",             "product", "rad"),
    "dih_term":           ("bonded",             "product", "kJ mol^-1"),
    "bonded_force_acc":   ("bonded",             "accumulator", "kJ mol^-1 nm^-1"),
    # LINCS constraint
    "constraint_length":  ("lincs_constraint",   "operand", "nm"),
    "lincs_invmass":      ("lincs_constraint",   "operand", "amu^-1"),
    "lincs_coupling":     ("lincs_constraint",   "product", "1"),
    "lincs_rhs":          ("lincs_constraint",   "product", "nm"),
    "lincs_correction":   ("lincs_constraint",   "accumulator", "nm"),
    # integration (leap-frog update)
    "mass":               ("integration",        "operand", "amu"),
    "update_invmass":     ("integration",        "operand", "amu^-1"),
    "update_accel":       ("integration",        "product", "nm ps^-2"),
    "update_vel":         ("integration",        "product", "nm ps^-1"),
    "update_pos":         ("integration",        "product", "nm"),
}

# kernel -> (motif, description, operation_role)
KERNEL_INFO = {
    "nonbonded_coulomb": (
        "n_body_methods",
        "Short-range real-space Coulomb interaction over Verlet neighbour pairs (nbnxm).",
        "per pair: E = f * q_i * q_j * rinv * erfc(beta*r); force summed onto each atom",
    ),
    "nonbonded_lj": (
        "n_body_methods",
        "Short-range Lennard-Jones (van der Waals) interaction over Verlet neighbour pairs.",
        "per pair: E = c12 * rinv^12 - c6 * rinv^6; force summed onto each atom",
    ),
    "nonbonded_lj": (
        "n_body_methods",
        "Short-range Lennard-Jones (van der Waals) interaction over Verlet neighbour pairs.",
        "per pair: E = c12*r^-12 - c6*r^-6; force summed onto each atom",
    ),
    "pme_reciprocal": (
        "spectral_methods",
        "PME reciprocal space: charge spreading, 3D FFT, influence-function solve, gather.",
        "grid energy = sum_m |F(q)(m)|^2 * influence(m)",
    ),
    "bonded": (
        "n_body_methods",
        "Listed bonded forces: bonds, angles, dihedrals (+ Urey-Bradley / CMAP in CHARMM).",
        "e.g. bond: E = 0.5*k*(b-b0)^2; force = k*(b-b0)",
    ),
    "lincs_constraint": (
        "sparse_linear_algebra",
        "LINCS: reset bond lengths via a truncated series of the constraint-coupling inverse.",
        "solve S*x = rhs approximately by matrix-vector iterations",
    ),
    "cmap": (
        "structured_grids",
        "CHARMM CMAP: 2D phi/psi backbone dihedral correction map (grid interpolation).",
        "E_corr = bicubic-interpolate(CMAP grid, phi, psi)",
    ),
    "integration": (
        "dense_linear_algebra",
        "Leap-frog Verlet time integration (per-atom vector update).",
        "a = F/m; v += a*dt; x += v*dt  (AXPY-like per-atom update)",
    ),
}

DEFAULT_SOURCE_URL = "https://www.mpinat.mpg.de/grubmueller/bench"


def slug(s):
    return re.sub(r"[^a-z0-9_\-]", "-", s.lower()).strip("-")


def resolve_summaries(patterns):
    out = []
    for p in patterns:
        if os.path.isdir(p):
            out += sorted(glob.glob(os.path.join(p, "**", "summary.json"), recursive=True))
        elif any(c in p for c in "*?[]"):
            out += sorted(glob.glob(p))
        else:
            out.append(p)
    # de-dup, keep existing
    seen, uniq = set(), []
    for f in out:
        if f in seen:
            continue
        seen.add(f)
        if os.path.exists(f):
            uniq.append(f)
        else:
            print(f"  ! missing: {f}")
    return uniq


def records_from_summary(summary, app_version, source_url, license_, app_name="GROMACS"):
    """Split one loaded summary.json dict into a list of schema records.
    app_name distinguishes force fields in the cross-app comparison
    (e.g. 'GROMACS' vs 'CHARMM' -- same engine, different force field / dataset)."""
    label = summary["label"]
    n_atoms = summary.get("n_atoms")
    tpr = summary.get("tpr", "")
    date = datetime.date.today().isoformat()

    # group operand stats by kernel
    by_kernel = {}  # kernel -> list of quantity dicts
    for op_name, stats in summary["operands"].items():
        if op_name not in OPERAND_MAP:
            print(f"    (skipping unmapped operand '{op_name}' in {label})")
            continue
        kernel, role, units = OPERAND_MAP[op_name]
        s = dict(stats)
        s.pop("name", None)  # 'name' is not part of the schema's stats block
        by_kernel.setdefault(kernel, []).append(
            {"name": op_name, "role": role, "units": units, "stats": s})

    records = []
    for kernel, quantities in by_kernel.items():
        motif, desc, op_role = KERNEL_INFO[kernel]
        rec = {
            "schema_version": "1.0",
            "id": f"{slug(app_name)}.{kernel}.{slug(label)}",
            "app": {"name": app_name, "version": app_version,
                    "url": "https://www.gromacs.org", "license": "LGPL-2.1"},
            "kernel": {"name": kernel, "description": desc, "operation_role": op_role},
            "motif": motif,
            "input_case": {"name": label,
                           "size": {"unit": "atoms", "value": n_atoms} if n_atoms else {"unit": "atoms", "value": 0},
                           "source_url": source_url, "license": license_},
            "provenance": {"method": "static_extract", "tool": "summary_to_records.py",
                           "tool_version": "1.0", "date": date, "seed": 0,
                           "notes": f"Derived from {os.path.basename(tpr) or label} operand summary; "
                                    "accumulation and format_fit layers require instrumented trace + replay."},
            "quantities": quantities,
            "notes": "Static operand-distribution record; accumulation/format_fit pending instrumentation.",
        }
        records.append(rec)
    return records


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="summary.json files, dirs, or globs")
    ap.add_argument("--app-name", default="GROMACS",
                    help="app/force-field label for the records (e.g. CHARMM) -> id prefix + app.name")
    ap.add_argument("--app-version", default="2025.3")
    ap.add_argument("--source-url", default=DEFAULT_SOURCE_URL)
    ap.add_argument("--license", default="CC-BY-4.0")
    ap.add_argument("--out", default=os.path.join(HERE, "..", "dataset", "records"))
    ap.add_argument("--no-validate", action="store_true")
    args = ap.parse_args()

    summaries = resolve_summaries(args.inputs)
    if not summaries:
        sys.exit("no summary.json files found")
    os.makedirs(args.out, exist_ok=True)

    validator = None
    if not args.no_validate:
        try:
            import jsonschema
            schema = json.load(open(SCHEMA_PATH))
            jsonschema.Draft7Validator.check_schema(schema)
            validator = jsonschema.Draft7Validator(schema)
        except ImportError:
            print("  (jsonschema not installed -> skipping validation)")

    n_records, n_fail = 0, 0
    for spath in summaries:
        summary = json.load(open(spath))
        recs = records_from_summary(summary, args.app_version, args.source_url,
                                    args.license, app_name=args.app_name)
        for rec in recs:
            out_path = os.path.join(args.out, rec["id"] + ".json")
            if validator is not None:
                errs = sorted(validator.iter_errors(rec), key=lambda e: list(e.path))
                if errs:
                    n_fail += 1
                    print(f"  FAIL {rec['id']}")
                    for e in errs:
                        loc = "/".join(str(x) for x in e.path) or "(root)"
                        print(f"        {loc}: {e.message}")
                    continue
            with open(out_path, "w") as f:
                json.dump(rec, f, indent=2)
            n_records += 1
            print(f"  OK   {rec['id']}  ({len(rec['quantities'])} quantities)")

    print(f"\nwrote {n_records} record(s) to {os.path.abspath(args.out)}"
          + (f"; {n_fail} failed validation" if n_fail else "; all valid"))


if __name__ == "__main__":
    main()
