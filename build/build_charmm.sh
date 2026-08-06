#!/usr/bin/env bash
# Route 1: build CHARMM36 systems with the GROMACS CHARMM36 force-field port, then run
# them through the value-distribution pipeline -- no CHARMM-GUI, no instrumented build.
#
# All CHARMM benchmark artifacts (downloaded PDBs, generated .tpr and .dump.txt, run
# outputs) live under  benchmarks/charmm/  -- the CHARMM counterpart of the GROMACS
# benchmarks in benchmarks/. Records go to dataset/records, histograms to docs/.
#
# Prereqs:
#   * a normal `gmx` install (the stock tool; not a custom build)
#   * the CHARMM36 GROMACS force field, e.g. charmm36-jul2022.ff/, placed inside
#     benchmarks/charmm/ (download charmm36-*.ff.tgz from
#     http://mackerell.umaryland.edu/charmm_ff.shtml -> GROMACS section, then `tar xf`
#     it there). pdb2gmx auto-detects a force field in its working directory.
#
# Usage (run from anywhere -- paths are resolved relative to the repo):
#   TER="0\n0\n" bash build/build_charmm.sh
#   FF=charmm36-jul2022 PDBS="1AKI 1UBQ 4AKE" bash build/build_charmm.sh
set -e

ROOT="$(cd "$(dirname "$0")/.." && pwd)"     # repo root (parent of build/)
BM="$ROOT/benchmarks/charmm"                 # all CHARMM benchmark files live here
MDP="$(cd "$(dirname "$0")" && pwd)/min.mdp"
FF="${FF:-charmm36-jul2022}"
PDBS="${PDBS:-1AKI 1UBQ 4AKE}"               # lysozyme, ubiquitin, adenylate kinase(=adk)
FIRST="$(echo $PDBS | awk '{print $1}')"

mkdir -p "$BM"
cd "$BM"                                      # do all gmx work inside benchmarks/charmm

for pdb in $PDBS; do
  [ -f "$pdb.pdb" ] || wget -q "https://files.rcsb.org/download/$pdb.pdb"
  # keep only protein backbone/side-chain atoms; drop crystal waters, ligands, ions
  # (HETATM) -- these have no hydrogens and make pdb2gmx fail ("missing atoms in Water").
  grep -E '^(ATOM|TER|END)' "$pdb.pdb" > "$pdb.clean.pdb"
  # -ter: choose termini MANUALLY. The CHARMM36 auto-terminus misfires here
  #   ("atom C1 not found in 1MET while combining tdb and rtp"). When prompted, pick
  #   the physiological standard: NH3+ for the START terminus and COO- for the END
  #   terminus (usually index 0 each time). These are complete caps -- do NOT pick
  #   "None" (that leaves dangling bonds: "dangling bond at terminal ends").
  #   For fully non-interactive runs, set TER="0\n0\n" (piped below).
  if [ -n "${TER:-}" ]; then
    printf "$TER" | gmx pdb2gmx -f "$pdb.clean.pdb" -o "$pdb.gro" -p "$pdb.top" -ff "$FF" -water tip3p -ignh -ter
  else
    gmx pdb2gmx -f "$pdb.clean.pdb" -o "$pdb.gro" -p "$pdb.top" -ff "$FF" -water tip3p -ignh -ter
  fi
  gmx editconf -f "$pdb.gro" -o "$pdb.box.gro" -c -d 1.0 -bt cubic
  gmx grompp  -f "$MDP" -c "$pdb.box.gro" -p "$pdb.top" -o "charmm_$pdb.tpr" -maxwarn 5
  gmx dump    -s "charmm_$pdb.tpr" > "charmm_$pdb.dump.txt" 2>/dev/null
done

# Build the records + histograms (CHARMM records are labelled app.name=CHARMM).
# --gmx-dump-dir "$BM" points the pipeline at the dumps we just wrote here.
python3 "$ROOT/scripts/nga_gromacs.py" "$BM"/charmm_*.tpr --gmx-dump-dir "$BM" --out "$BM/runs"
python3 "$ROOT/dataset_spec/summary_to_records.py" "$BM/runs" --app-name CHARMM --out "$ROOT/dataset/records"
python3 "$ROOT/dataset_spec/histograms.py" "$ROOT"/dataset/records/*.json --out "$ROOT/docs"
echo ">>> DONE. CHARMM artifacts in benchmarks/charmm/ ; cross-app figure: docs/hist_apps.png"
