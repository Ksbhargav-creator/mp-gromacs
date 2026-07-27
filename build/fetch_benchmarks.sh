#!/usr/bin/env bash
# fetch_benchmarks.sh -- download the Max Planck (Grubmueller/Kutzner) free GROMACS
# benchmark .tpr files. CC-BY 4.0. This is the chosen dataset source (size ladder
# from ~82k to 12M atoms) for the operand/product distribution study.
#
#   benchMEM   ~82k atoms   membrane protein     <- small prototype (start here)
#   benchRIB   ~2M atoms    ribosome in water    <- mid scale
#   benchPEP   ~12M atoms   peptides in water    <- stress scale
#
# Landing page (verify current filenames/links there):
#   https://www.mpinat.mpg.de/grubmueller/bench
set -euo pipefail
DEST="${DEST:-$PWD/benchmarks}"
mkdir -p "$DEST"; cd "$DEST"

# The download links on the landing page are the source of truth. Set BENCH_URLS
# to the direct .zip/.tpr URLs you copy from that page, e.g.:
#   export BENCH_URLS="https://www.mpinat.mpg.de/.../benchMEM.zip"
: "${BENCH_URLS:=}"

if [[ -z "$BENCH_URLS" ]]; then
  cat <<'EOF'
No BENCH_URLS set. Open the benchmark page and copy the direct download link(s):
    https://www.mpinat.mpg.de/grubmueller/bench
Then:
    export BENCH_URLS="<benchMEM url> <benchRIB url> ..."
    ./fetch_benchmarks.sh
Start with benchMEM (82k atoms) as the small prototype case.
EOF
  exit 0
fi

for url in $BENCH_URLS; do
  echo ">> fetching $url"
  curl -fSLO "$url"
done
# unzip any archives
for z in *.zip; do [[ -f "$z" ]] && unzip -o "$z"; done
echo ">> benchmark .tpr files in $DEST:"
find "$DEST" -name '*.tpr' -print
