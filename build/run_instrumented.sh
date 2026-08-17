#!/usr/bin/env bash
# run_instrumented.sh -- Method A capture: run the instrumented GROMACS in BOTH precisions
# over BOTH app input sets (GROMACS benchmarks + CHARMM36 systems), one force evaluation
# each, and collect the runtime operand traces + dynamic-range stats into a labelled tree.
#
#   apps        : GROMACS  = benchmarks/*.tpr
#                 CHARMM   = benchmarks/charmm/charmm_*.tpr   (same binary, CHARMM36 inputs)
#   precisions  : fp32 = $GMX_S (default build, `gmx`)
#                 fp64 = $GMX_D (double build,  `gmx_d`)
#
# Output layout:
#   runs_instr/<APP>/<PREC>/<system>/{nga_force_trace.csv, nga_gromacs_stats.json, .csv}
#
# Usage:
#   GMX_S=/path/build-single/bin/gmx GMX_D=/path/build-double/bin/gmx_d \
#     bash build/run_instrumented.sh
#   # subset:  SYSTEMS_GMX="benchMEM" SYSTEMS_CHARMM="charmm_1AKI" bash build/run_instrumented.sh
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BM="$ROOT/benchmarks"
OUT="$ROOT/runs_instr"
GMX_S="${GMX_S:-gmx}"       # fp32 (single/mixed) binary
GMX_D="${GMX_D:-gmx_d}"     # fp64 (double) binary

command -v "$GMX_S" >/dev/null 2>&1 || [ -x "$GMX_S" ] || { echo "!! fp32 binary not found: $GMX_S"; exit 1; }
command -v "$GMX_D" >/dev/null 2>&1 || [ -x "$GMX_D" ] || { echo "!! fp64 binary not found: $GMX_D"; exit 1; }

# The per-interaction trace (nga_force_trace.csv) is one row PER PAIR and explodes to
# 100s of GB on large systems. By default we DISCARD it by pointing the file at
# /dev/null (the tracer just fopen()s that name), so even a trace-compiled binary is
# safe and we still get the aggregate nga_gromacs_stats.json we need for distributions.
# Set KEEP_TRACE=1 to actually keep the CSV -- do that ONLY on a small system, ideally
# with the binary's runtime cap (NGA_TRACE_MAX / NGA_TRACE_STRIDE).
KEEP_TRACE="${KEEP_TRACE:-0}"

# NSTEPS=1 by default so ALL instrumented kernels fire: nonbonded/bonded/PME run in the
# force evaluation (would fire even at nsteps=0), but LINCS (constraints) and the leap-frog
# integration only run when an actual step is taken. One step is enough and stays nearly
# deterministic. (Integration hooks fire only for the md integrator, not steepest descent.)
NSTEPS="${NSTEPS:-1}"

# One deterministic step: reference kernel, single rank/thread.
run_one() {  # $1 app  $2 prec  $3 gmxbin  $4 tpr
  local app="$1" prec="$2" bin="$3" tpr="$4"
  local label; label="$(basename "${tpr%.tpr}")"
  local d="$OUT/$app/$prec/$label"
  [ -e "$tpr" ] || { echo "   ! missing $tpr (skip)"; return; }
  echo ">> [$app/$prec] $label"
  rm -rf "$d"; mkdir -p "$d"
  # discard the trace unless explicitly kept
  [ "$KEEP_TRACE" = "1" ] || ln -sf /dev/null "$d/nga_force_trace.csv"
  ( cd "$d" && \
    GMX_NBNXN_REF=1 GMX_DISABLE_GPU_DETECTION=1 \
    "$bin" mdrun -s "$tpr" -nsteps "$NSTEPS" -ntmpi 1 -ntomp 1 \
        -o traj.trr -e ener.edr -g md.log >/dev/null 2>&1 || true )
  # the aggregate stats are the deliverable for distributions
  [ -s "$d/nga_gromacs_stats.json" ] || echo "   ! no nga_gromacs_stats.json (flush not called?)"
  [ "$KEEP_TRACE" = "1" ] && { [ -s "$d/nga_force_trace.csv" ] || echo "   ! no trace CSV"; } || true
}

# --- build the system lists ---------------------------------------------------------
GMX_TPRS=()
if [ -n "${SYSTEMS_GMX:-}" ]; then for s in $SYSTEMS_GMX; do GMX_TPRS+=("$BM/$s.tpr"); done
else for f in "$BM"/*.tpr; do [ -e "$f" ] && GMX_TPRS+=("$f"); done; fi

CHARMM_TPRS=()
if [ -n "${SYSTEMS_CHARMM:-}" ]; then for s in $SYSTEMS_CHARMM; do CHARMM_TPRS+=("$BM/charmm/$s.tpr"); done
else for f in "$BM"/charmm/charmm_*.tpr; do [ -e "$f" ] && CHARMM_TPRS+=("$f"); done; fi

echo "======================================================================"
echo " Method A capture  ->  $OUT"
echo "   GROMACS systems: ${#GMX_TPRS[@]}   CHARMM systems: ${#CHARMM_TPRS[@]}"
echo "   fp32: $GMX_S    fp64: $GMX_D"
echo "======================================================================"

for tpr in "${GMX_TPRS[@]:-}";    do [ -n "$tpr" ] || continue
  run_one GROMACS fp32 "$GMX_S" "$tpr"; run_one GROMACS fp64 "$GMX_D" "$tpr"; done
for tpr in "${CHARMM_TPRS[@]:-}"; do [ -n "$tpr" ] || continue
  run_one CHARMM  fp32 "$GMX_S" "$tpr"; run_one CHARMM  fp64 "$GMX_D" "$tpr"; done

echo; echo ">> DONE. Traces + stats under $OUT/<APP>/<fp32|fp64>/<system>/"
echo ">> Next: pair fp32 vs fp64 rows and compute error (see"
echo "         instrumentation/DESIGN_precision_experiment.md)."
