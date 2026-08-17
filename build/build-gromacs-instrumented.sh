#!/usr/bin/env bash
# build-gromacs-instrumented.sh -- Method A: build the SAME instrumented GROMACS source
# tree in BOTH precisions, so runtime operands can be captured in fp32 and fp64.
#
#   single (default, mixed precision) -> binary `gmx`     : real = float  (fp32)
#   double (-DGMX_DOUBLE=on)          -> binary `gmx_d`   : real = double (fp64)
#
# The precision comes entirely from the build; the NGA logging headers are unchanged
# (they cast to double only for output, so a fp32 `real` is logged as its fp32 value).
#
# PREREQUISITE -- the source tree must already be PATCHED with the kernel hooks from
#   instrumentation/gromacs_nbnxm_instrument.md  (NGA_TRACE_ROW / NGA_RECORD in
#   kernels_reference/kernel_ref_inner.h, the one-shot nbfp dump, and NGA_FLUSH /
#   NGA_TRACE_CLOSE in runner.cpp). This script copies the headers and verifies the
#   patch is present, but it cannot paste the grep-anchored snippets for you.
#
# Usage:
#   GMX_SRC=/path/to/gromacs-2026.3 bash build/build-gromacs-instrumented.sh
#   GMX_SRC=... JOBS=8 bash build/build-gromacs-instrumented.sh   # parallel make
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GMX_SRC="${GMX_SRC:?set GMX_SRC=/path/to/gromacs/source/tree}"
JOBS="${JOBS:-4}"
# The instrumentation is enabled purely through compiler defines. Keep the space-
# separated flags as ONE cmake argument (a bash array preserves the quoting; a bare
# string would word-split and cmake would read "-DNGA_TRACE" as a stray cache define).
#
# NGA_INSTRUMENT  -> aggregate per-quantity stats (nga_gromacs_stats.json). CONSTANT
#                    memory/disk; this is what the value-distribution records need.
# NGA_TRACE       -> per-interaction CSV (nga_force_trace.csv). One row PER PAIR, so it
#                    explodes on large systems (100s of GB). OPT-IN only, and even then
#                    it is capped at runtime (see NGA_TRACE_MAX / NGA_TRACE_STRIDE).
# Enable the trace with:  TRACE=1 GMX_SRC=... bash build/build-gromacs-instrumented.sh
NGA_FLAGS="-DNGA_INSTRUMENT"
if [ "${TRACE:-0}" = "1" ]; then
  NGA_FLAGS="$NGA_FLAGS -DNGA_TRACE"
  echo ">> NGA_TRACE enabled -- per-interaction CSV (capped; use only on a SMALL system)"
else
  echo ">> aggregate stats only (no per-interaction trace). Use TRACE=1 to enable the trace."
fi
# Reference scalar kernel (exposes named scalar operands), no SIMD/GPU/MPI -> deterministic,
# single-eval friendly, and the operands are actually readable.
# GMX_OPENMP=OFF: Apple's clang ships no OpenMP, and we run single-threaded
# (-ntomp 1) for the reference kernel anyway, so threading is not needed.
COMMON_CMAKE=(
  -DGMX_SIMD=NONE -DGMX_GPU=OFF -DGMX_MPI=OFF -DGMX_OPENMP=OFF -DGMX_BUILD_OWN_FFTW=ON
  -DCMAKE_C_FLAGS="$NGA_FLAGS" -DCMAKE_CXX_FLAGS="$NGA_FLAGS"
)

echo ">> GROMACS source : $GMX_SRC"
[ -f "$GMX_SRC/CMakeLists.txt" ] || { echo "!! not a GROMACS source tree: $GMX_SRC"; exit 1; }

# --- 1. make the NGA headers visible to the build -----------------------------------
NGA_DIR="$GMX_SRC/src/gromacs/nga"
mkdir -p "$NGA_DIR"
cp "$ROOT/instrumentation/nga_range_stats.hpp" "$NGA_DIR/"
cp "$ROOT/mixed_precision/nga_trace.hpp"        "$NGA_DIR/"
echo ">> copied NGA headers into src/gromacs/nga/"

# --- 2. verify the kernel is actually patched (else the build is a silent no-op) ----
KERN="$GMX_SRC/src/gromacs/nbnxm/kernels_reference/kernel_ref_inner.h"
if ! grep -q "NGA_TRACE_ROW\|NGA_RECORD" "$KERN" 2>/dev/null; then
  echo "!! No NGA hooks found in $(basename "$KERN")."
  echo "!! Paste the snippets from instrumentation/gromacs_nbnxm_instrument.md (§B) first,"
  echo "!! plus NGA_FLUSH/NGA_TRACE_CLOSE in src/gromacs/mdrun/runner.cpp (§C)."
  exit 2
fi
echo ">> kernel hooks present in kernel_ref_inner.h"

# --- 3. two out-of-tree builds ------------------------------------------------------
build_one() {  # $1 = label ; remaining args = extra cmake defines (precision)
  local label="$1"; shift
  local bdir="$GMX_SRC/build-$label"
  echo; echo ">> configuring $label build ($bdir)"
  cmake -S "$GMX_SRC" -B "$bdir" "${COMMON_CMAKE[@]}" "$@"
  echo ">> building $label (make -j$JOBS)"
  cmake --build "$bdir" --target gmx -j "$JOBS"
}

build_one single                    # -> gmx    (fp32)
build_one double -DGMX_DOUBLE=on     # -> gmx_d  (fp64)

echo; echo "======================================================================"
echo " DONE. Instrumented binaries:"
echo "   fp32 : $GMX_SRC/build-single/bin/gmx"
echo "   fp64 : $GMX_SRC/build-double/bin/gmx_d"
echo " Next: capture runtime operands with"
echo "   GMX_S=$GMX_SRC/build-single/bin/gmx GMX_D=$GMX_SRC/build-double/bin/gmx_d \\"
echo "     bash build/run_instrumented.sh"
echo "======================================================================"
