#!/usr/bin/env bash
# build-gromacs.sh -- reproducible GROMACS build for the NGA operand/product study.
#
# Two paths:
#   --spack     : let Spack resolve the toolchain + FFTW and build GROMACS (recommended
#                 on the HPC cluster; matches how the SPICE side was pinned).
#   (default)   : manual CMake build from a release tarball, single-node.
#
# Add --instrument to compile the reference kernels with -DNGA_INSTRUMENT so the
# force-loop product logging in instrumentation/ is active. A build WITHOUT that
# flag is a stock GROMACS (use it for the operand-only, gmx-dump path).
#
# Usage:
#   ./build-gromacs.sh                      # manual, stock (double precision)
#   ./build-gromacs.sh --instrument         # manual, with NGA force-loop logging
#   ./build-gromacs.sh --spack --instrument # Spack build + instrumentation
#
# After building, source the resulting GMXRC to get `gmx` on PATH.
set -euo pipefail

GMX_VERSION="${GMX_VERSION:-2025.3}"          # pin a version you can reproduce (latest stable as of 2026-07)
PREFIX="${PREFIX:-$PWD/gromacs-nga-install}"
BUILD_DIR="${BUILD_DIR:-$PWD/gromacs-build}"
NGA_HDR_SRC="$(cd "$(dirname "$0")/../instrumentation" && pwd)/nga_range_stats.hpp"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"

USE_SPACK=0
INSTRUMENT=0
for a in "$@"; do
  case "$a" in
    --spack) USE_SPACK=1 ;;
    --instrument) INSTRUMENT=1 ;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown arg: $a" >&2; exit 2 ;;
  esac
done

# Double precision is the right default here: it gives the reference distribution
# the mixed/reduced-precision formats are to be compared AGAINST. (GROMACS single
# is the existing default-mixed build; double is what we instrument for ground truth.)
GMX_DOUBLE="${GMX_DOUBLE:-ON}"

echo ">> GROMACS $GMX_VERSION  double=$GMX_DOUBLE  instrument=$INSTRUMENT  spack=$USE_SPACK"
echo ">> prefix: $PREFIX"

# ---------------------------------------------------------------------------
if [[ "$USE_SPACK" == "1" ]]; then
  command -v spack >/dev/null || { echo "spack not on PATH; run: . <spack>/share/spack/setup-env.sh" >&2; exit 1; }
  VAR_DOUBLE=$([[ "$GMX_DOUBLE" == "ON" ]] && echo "+double" || echo "~double")
  # Spack can't inject -DNGA_INSTRUMENT directly; use cxxflags for the define and
  # pre-place the header where the instrumented sources #include it.
  SPEC="gromacs@${GMX_VERSION} ${VAR_DOUBLE} +mpi~cuda build_type=RelWithDebInfo"
  if [[ "$INSTRUMENT" == "1" ]]; then
    echo ">> For an instrumented Spack build: apply the instrumentation patch to a"
    echo ">>   'spack develop' checkout, then add cxxflags=-DNGA_INSTRUMENT:"
    echo ">>   spack develop gromacs@${GMX_VERSION}"
    echo ">>   (copy nga_range_stats.hpp into the dev source tree first)"
    SPEC+=" cxxflags=-DNGA_INSTRUMENT"
  fi
  echo ">> spack install $SPEC"
  spack install $SPEC
  echo ">> done. Load with:  spack load gromacs@${GMX_VERSION}"
  exit 0
fi

# ---------------------------------------------------------------------------
# Manual CMake build
mkdir -p "$BUILD_DIR"; cd "$BUILD_DIR"
TARBALL="gromacs-${GMX_VERSION}.tar.gz"
if [[ ! -d "gromacs-${GMX_VERSION}" ]]; then
  [[ -f "$TARBALL" ]] || curl -fsSLO "https://ftp.gromacs.org/gromacs/${TARBALL}"
  tar xf "$TARBALL"
fi
SRC="$BUILD_DIR/gromacs-${GMX_VERSION}"

# place instrumentation header inside the tree
mkdir -p "$SRC/src/gromacs/nga"
cp "$NGA_HDR_SRC" "$SRC/src/gromacs/nga/"
echo ">> copied nga_range_stats.hpp into src/gromacs/nga/"

if [[ "$INSTRUMENT" == "1" ]]; then
  echo ">> NOTE: apply the kernel hooks from instrumentation/gromacs_nbnxm_instrument.md"
  echo ">>       (guarded by #ifdef NGA_INSTRUMENT) before this build if not already done."
  EXTRA_CXX="-DNGA_INSTRUMENT"
else
  EXTRA_CXX=""
fi

cmake -S "$SRC" -B "$SRC/build" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DGMX_DOUBLE="$GMX_DOUBLE" \
  -DGMX_MPI=OFF \
  -DGMX_GPU=OFF \
  -DGMX_SIMD=NONE \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_CXX_FLAGS="$EXTRA_CXX"

cmake --build "$SRC/build" -j "$JOBS"
cmake --install "$SRC/build"

echo ""
echo ">> installed to $PREFIX"
echo ">> source it:   source $PREFIX/bin/GMXRC"
echo ">> sanity:      gmx --version | head -20"
[[ "$INSTRUMENT" == "1" ]] && \
  echo ">> instrumented: a run will emit nga_gromacs_stats.{json,csv} in CWD"
# Note: -DGMX_SIMD=NONE keeps the scalar reference kernels active, which is what
# the per-interaction product hooks read. For the operand-only (gmx dump) path,
# SIMD can be left at native for speed.
