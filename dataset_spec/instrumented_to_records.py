#!/usr/bin/env python3
"""
instrumented_to_records.py -- convert the RUNTIME instrumented nonbonded stats
(Method A: kernel_ref_inner.h hooks, see instrumentation/INSTRUMENTATION.md) into
kernel_characterization records (provenance.method = "instrumented_trace").

Input
-----
runs_instr/<APP>/<fp32|fp64>/<system>/{nga_gromacs_stats.json, nga_gromacs_stats.csv}

  nga_gromacs_stats.json: {"minexp":..,"maxexp":..,"operands":{name: {count,
    count_zero, count_nan, abs_min, abs_max, decades, signed_min, signed_max,
    log2_hist: [...]}}}   -- abs_min/abs_max/count/count_zero are the authoritative,
    UNCLAMPED extremes (tracked as running min/max, never binned).

  nga_gromacs_stats.csv: "operand,binade_exp,count" -- the same histogram, but as
    (exponent, count) pairs instead of a positional array, so it is unambiguous
    regardless of how the JSON array was offset/clamped. This is used as the
    histogram GROUND TRUTH.

Why not just trust the JSON's flat log2_hist array: cross-checking csv_sum against
(count - count_zero) shows the CSV is short by exactly the count of one extreme
high-exponent bin for `lj_rinv6` in every run (the CSV writer's loop range does not
reach exponents that far out). That single missing bin is recovered here from the
JSON's abs_max/count discrepancy and appended as one synthetic histogram bin, so no
data is silently dropped -- see instrumentation/INSTRUMENTATION.md and
docs/format_fit_runtime_results.md ("the rinv / lj_rinv6 clamp artifact") for the
full diagnosis: this bin, and `rinv`'s own isolated high bin, are a divide-by-near-
zero guard on masked/padding Verlet-cluster pairs, not real chemistry. This script
does not trim anything -- it emits the RAW measured record; trimming is a disclosed,
separate analysis step (dataset_spec/format_fit.py consumers decide whether to trim).

Output
------
One record per (app, precision, kernel, system): dataset/records_instrumented/
  <app>.<kernel>.<system>-<precision>.json
validated against dataset_spec/kernel_characterization.schema.json.

Usage
-----
    python3 dataset_spec/instrumented_to_records.py [--runs runs_instr] [--out dataset/records_instrumented]
"""
import argparse, csv, datetime, glob, json, math, os, re, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
SCHEMA_PATH = os.path.join(HERE, "kernel_characterization.schema.json")

# runtime operand name -> (kernel, role, units)
# nonbonded_coulomb / nonbonded_lj entries mirror summary_to_records.py's OPERAND_MAP
# (same kernel taxonomy); nb_fscal/nb_force_acc get their own kernel because they are
# the LJ+Coulomb *combined* force scalar / per-atom accumulator (INSTRUMENTATION.md
# section 3.3 tags them "nonbonded / n_body" and "nonbonded (accumulation)" rather
# than either nonbonded_coulomb or nonbonded_lj alone).
RUNTIME_OPERAND_MAP = {
    "coulomb_qq":    ("nonbonded_coulomb", "product", "kJ mol^-1 nm (Coulomb-prefactor-scaled q_i*q_j; see notes)"),
    "coulomb_term":  ("nonbonded_coulomb", "product", "kJ mol^-1"),
    "rinv":          ("nonbonded_coulomb", "product", "nm^-1"),
    "lj_c6":         ("nonbonded_lj",      "operand",  "kJ mol^-1 nm^6"),
    "lj_c12":        ("nonbonded_lj",      "operand",  "kJ mol^-1 nm^12"),
    "lj_rinv6":      ("nonbonded_lj",      "product",  "nm^-6"),
    "lj_term6":      ("nonbonded_lj",      "product",  "kJ mol^-1"),
    "lj_term12":     ("nonbonded_lj",      "product",  "kJ mol^-1"),
    "lj_force":      ("nonbonded_lj",      "product",  "kJ mol^-1 nm^-1"),
    "nb_fscal":      ("nonbonded_force",   "product",     "kJ mol^-1 nm^-1"),
    "nb_force_acc":  ("nonbonded_force",   "accumulator", "kJ mol^-1 nm^-1"),
}

KERNEL_INFO = {
    "nonbonded_coulomb": (
        "n_body_methods",
        "Short-range real-space Coulomb interaction over Verlet neighbour pairs (nbnxm); "
        "runtime products from the reference kernel's inner loop.",
        "per pair: coulomb_qq = f*q_i*q_j; coulomb_term = coulomb_qq * rinv",
    ),
    "nonbonded_lj": (
        "n_body_methods",
        "Short-range Lennard-Jones interaction over Verlet neighbour pairs; runtime "
        "products from the reference kernel's inner loop.",
        "per pair: lj_rinv6 = rinv^6; lj_term6 = c6*lj_rinv6; lj_term12 = c12*lj_rinv6^2; "
        "lj_force = 12*lj_term12 - 6*lj_term6",
    ),
    "nonbonded_force": (
        "n_body_methods",
        "Combined LJ+Coulomb per-pair force scalar and its running per-atom accumulation "
        "(leap-frog force sum), captured at the point of force combine in the reference kernel.",
        "nb_fscal = per-pair net force prefactor (LJ + Coulomb); "
        "nb_force_acc = running sum of nb_fscal-derived force onto atom i",
    ),
}

DEFAULT_SOURCE_URL = "https://www.mpinat.mpg.de/grubmueller/bench"


def slug(s):
    return re.sub(r"[^a-z0-9_\-]", "-", s.lower()).strip("-")


def load_csv_hist(csv_path):
    """-> {operand: [(binade_exp, count), ...]} sorted by exponent."""
    out = defaultdict(list)
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            try:
                e = int(row["binade_exp"]); c = int(row["count"])
            except (KeyError, ValueError):
                continue
            out[row["operand"]].append((e, c))
    return {op: sorted(xs) for op, xs in out.items()}


def build_stats(op_name, json_stats, csv_bins):
    """-> dynamic_range_stats dict (schema-shaped), recovering any bin the CSV
    writer dropped (verified against count - count_zero)."""
    count = json_stats["count"]
    count_zero = json_stats["count_zero"]
    nonzero_expected = count - count_zero
    bins = list(csv_bins)
    csv_sum = sum(c for _, c in bins)
    missing = nonzero_expected - csv_sum
    recovered_note = None
    if missing > 0:
        abs_max = json_stats["abs_max"]
        top_exp = math.floor(math.log2(abs_max)) if abs_max > 0 else 0
        bins.append((top_exp, missing))
        bins.sort()
        recovered_note = (
            f"recovered {missing} count(s) missing from the CSV histogram as one "
            f"synthetic bin at exponent {top_exp} (from abs_max={abs_max:.6g}); "
            f"the CSV writer's bin range does not reach this exponent -- see "
            f"docs/format_fit_runtime_results.md"
        )
    elif missing < 0:
        # more in the CSV than expected -- surface loudly rather than silently truncating
        recovered_note = (
            f"WARNING: csv histogram sum ({csv_sum}) exceeds count-count_zero "
            f"({nonzero_expected}) by {-missing}; investigate before trusting this record"
        )

    edges = [e for e, _ in bins] + [bins[-1][0] + 1] if bins else []
    counts = [c for _, c in bins]

    abs_min = json_stats["abs_min"]
    abs_max = json_stats["abs_max"]
    binades_log2 = (math.log2(abs_max / abs_min)
                     if (abs_min and abs_min > 0 and abs_max and abs_max > 0) else None)

    stats = {
        "count": count,
        "count_zero": count_zero,
        "count_nonzero": nonzero_expected,
        "count_nan": json_stats.get("count_nan", 0),
        "signed_min": json_stats.get("signed_min"),
        "signed_max": json_stats.get("signed_max"),
        "abs_min": abs_min,
        "abs_max": abs_max,
        "decades": json_stats.get("decades"),
        "binades_log2": binades_log2,
    }
    if edges:
        stats["log2_hist"] = {"edges": edges, "counts": counts}
    return stats, recovered_note


def records_from_run(app_name, precision, system, json_stats_path, csv_path,
                      app_version, source_url, license_):
    json_stats = json.load(open(json_stats_path))
    csv_bins = load_csv_hist(csv_path)
    date = datetime.date.today().isoformat()

    by_kernel = defaultdict(list)
    for op_name, ops in json_stats["operands"].items():
        if op_name not in RUNTIME_OPERAND_MAP:
            print(f"    (skipping unmapped runtime operand '{op_name}')")
            continue
        kernel, role, units = RUNTIME_OPERAND_MAP[op_name]
        stats, note = build_stats(op_name, ops, csv_bins.get(op_name, []))
        q = {"name": op_name, "role": role, "units": units, "stats": stats}
        if note:
            q["_recovery_note"] = note  # stripped before writing; folded into provenance.notes
        by_kernel[kernel].append(q)

    records = []
    for kernel, quantities in by_kernel.items():
        motif, desc, op_role = KERNEL_INFO[kernel]
        notes = [f"{q['name']}: {q.pop('_recovery_note')}" for q in quantities if "_recovery_note" in q]
        rec = {
            "schema_version": "1.0",
            "id": f"{slug(app_name)}.{kernel}.{slug(system)}-{precision}",
            "app": {"name": app_name, "version": app_version,
                    "url": "https://www.gromacs.org", "license": "LGPL-2.1"},
            "kernel": {"name": kernel, "description": desc, "operation_role": op_role},
            "motif": motif,
            "input_case": {"name": system,
                           "size": {"unit": "atoms", "value": 0},
                           "source_url": source_url, "license": license_},
            "provenance": {
                "method": "instrumented_trace",
                "tool": "instrumented_to_records.py",
                "tool_version": "1.0",
                "date": date,
                "seed": 0,
                "notes": (
                    f"GROMACS reference kernel (kernel_ref_1x1/4x4, GMX_SIMD=NONE), "
                    f"real={precision}, one deterministic force evaluation "
                    f"(GMX_NBNXN_REF=1, -nsteps 0, -ntmpi 1 -ntomp 1). "
                    f"See instrumentation/INSTRUMENTATION.md for hook locations. "
                    + (" | ".join(notes) if notes else "")
                ).strip(),
            },
            "quantities": quantities,
            "notes": ("Runtime (instrumented_trace) record: operands + products the nonbonded "
                      "reference kernel actually forms, not just static .tpr parameters. "
                      "accumulation/format_fit layers still pending the paired-trace replay step."),
        }
        records.append(rec)
    return records


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", default=os.path.join(HERE, "..", "runs_instr"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "dataset", "records_instrumented"))
    ap.add_argument("--app-version", default="2026.3")
    ap.add_argument("--source-url", default=DEFAULT_SOURCE_URL)
    ap.add_argument("--license", default="CC-BY-4.0")
    ap.add_argument("--no-validate", action="store_true")
    args = ap.parse_args()

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

    n_ok, n_fail = 0, 0
    stats_files = sorted(glob.glob(os.path.join(args.runs, "*", "*", "*", "nga_gromacs_stats.json")))
    for jpath in stats_files:
        cpath = jpath.replace(".json", ".csv")
        if not os.path.exists(cpath):
            print(f"  ! no csv for {jpath}, skipping")
            continue
        parts = jpath.split(os.sep)
        # .../runs_instr/<APP>/<precision>/<system>/nga_gromacs_stats.json
        system, precision, app = parts[-2], parts[-3], parts[-4]
        recs = records_from_run(app, precision, system, jpath, cpath,
                                args.app_version, args.source_url, args.license)
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
            n_ok += 1
            print(f"  OK   {rec['id']}  ({len(rec['quantities'])} quantities)")

    print(f"\nwrote {n_ok} record(s) to {os.path.abspath(args.out)}"
          + (f"; {n_fail} failed validation" if n_fail else "; all valid"))


if __name__ == "__main__":
    main()
