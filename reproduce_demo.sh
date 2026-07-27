#!/usr/bin/env bash
# reproduce_demo.sh -- reproduce the adenylate-kinase operand histogram on your
# machine, end to end, with no GROMACS build required.
#
# It (1) installs the pure-Python deps, (2) locates the bundled adk_oplsaa.tpr
# (adenylate kinase, OPLS-AA, 47,681 atoms) that ships inside MDAnalysisTests,
# and (3) runs the operand pipeline to regenerate the histogram.
#
#   ./reproduce_demo.sh
#
# Output figure: figs/operand_histograms_adk_oplsaa.png
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

echo ">> [1/3] installing Python deps (MDAnalysis, numpy, matplotlib + test data)"
# MDAnalysisTests carries the demo .tpr. Use --only-binary for MDAnalysis so it
# doesn't compile from source (slow); the Tests pkg is pure-python.
python3 -m pip install --upgrade "numpy<2" matplotlib >/dev/null
python3 -m pip install --only-binary=:all: "MDAnalysis>=2.7,<3" >/dev/null
python3 -m pip install "MDAnalysisTests>=2.7,<3" >/dev/null || \
  python3 -m pip install --no-build-isolation "MDAnalysisTests>=2.7,<3" >/dev/null

echo ">> [2/3] locating bundled adk_oplsaa.tpr"
TPR="$(python3 - <<'PY'
import importlib.util, os, sys
s = importlib.util.find_spec("MDAnalysisTests")   # find_spec does NOT run __init__
if s is None:
    sys.exit("MDAnalysisTests not importable")
p = os.path.join(list(s.submodule_search_locations)[0], "data", "adk_oplsaa.tpr")
assert os.path.exists(p), p
print(p)
PY
)"
echo "   using: $TPR"

echo ">> [3/3] running operand pipeline"
python3 "$HERE/scripts/extract_operands.py" --tpr "$TPR" --label adk_oplsaa --out "$HERE/out_adk"
python3 "$HERE/scripts/plot_histograms.py" "$HERE/out_adk/operands.npz" \
        --summary "$HERE/out_adk/summary.json" --out "$HERE/figs"

echo ""
echo ">> DONE. Figure: $HERE/figs/operand_histograms_adk_oplsaa.png"
echo ">> Dynamic-range numbers are in $HERE/out_adk/summary.json"
