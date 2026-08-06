#!/usr/bin/env bash
# build_gromacs.sh -- regenerate the GROMACS operand records + histograms from the
# raw benchmark inputs in benchmarks/ (the GROMACS counterpart of build_charmm.sh).
#
# Unlike CHARMM there is no system to build: the Max Planck benchmark .tpr files are
# the inputs. This script (1) ensures each benchmarks/<name>.tpr has a matching
# <name>.dump.txt (generating it with `gmx dump` if a gmx binary is available and the
# dump is missing), then (2) runs the operand pipeline over them and writes the
# gromacs.* records (app.name=GROMACS) and the value-distribution histograms.
#
#   * charges/masses come from the .tpr via MDAnalysis (no gmx needed for these)
#   * LJ c6/c12 + bonded + constraint params come from the .dump.txt
#
# Inputs live in benchmarks/ (git-ignored). If empty, fetch them first:
#   BENCH_URLS="<benchMEM url> ..." bash build/fetch_benchmarks.sh
#
# Usage (run from anywhere):
#   bash build/build_gromacs.sh
#   TPRS="benchMEM benchRIB" bash build/build_gromacs.sh   # subset by stem
set -e

ROOT="$(cd "$(dirname "$0")/.." && pwd)"     # repo root (parent of build/)
BM="$ROOT/benchmarks"                         # GROMACS benchmark .tpr live here (top level)

# Collect the .tpr set: either the explicit stems in $TPRS, or every top-level
# benchmarks/*.tpr (the charmm/ and zip-files/ subdirs are intentionally excluded).
TPRLIST=()
if [ -n "${TPRS:-}" ]; then
  for s in $TPRS; do TPRLIST+=("$BM/$s.tpr"); done
else
  for f in "$BM"/*.tpr; do [ -e "$f" ] && TPRLIST+=("$f"); done
fi
if [ ${#TPRLIST[@]} -eq 0 ]; then
  echo "No .tpr files in $BM. Fetch them first:"
  echo "  BENCH_URLS=\"<benchMEM url> ...\" bash build/fetch_benchmarks.sh"
  exit 1
fi

# Ensure a matching .dump.txt for each .tpr (needed for LJ/bonded/constraint operands).
HAVE_GMX="$(command -v gmx || true)"
for tpr in "${TPRLIST[@]}"; do
  [ -e "$tpr" ] || { echo "  ! missing: $tpr (skipping)"; continue; }
  dump="${tpr%.tpr}.dump.txt"
  if [ ! -s "$dump" ]; then
    if [ -n "$HAVE_GMX" ]; then
      echo ">> gmx dump -> $(basename "$dump")"
      gmx dump -s "$tpr" > "$dump" 2>/dev/null
    else
      echo "  ! no $(basename "$dump") and no gmx binary -- LJ/bonded operands will be"
      echo "    skipped for $(basename "$tpr") (charges/masses still extracted from .tpr)"
    fi
  fi
done

echo; echo "======================================================================"
echo " Reproducing GROMACS records from ${#TPRLIST[@]} benchmark .tpr in benchmarks/"
echo "======================================================================"

# (1) operands -> per-system summaries    (2) summaries -> schema records (app=GROMACS)
python3 "$ROOT/scripts/nga_gromacs.py" "${TPRLIST[@]}" --gmx-dump-dir "$BM" --out "$ROOT/runs"
python3 "$ROOT/dataset_spec/summary_to_records.py" "$ROOT/runs" --app-name GROMACS --out "$ROOT/dataset/records"

# (3) value-distribution histograms across every record (GROMACS + any CHARMM present)
python3 "$ROOT/dataset_spec/histograms.py" "$ROOT"/dataset/records/*.json --out "$ROOT/docs"

echo ">>> DONE. GROMACS records in dataset/records/gromacs.* ; histograms in docs/"
