# =============================================================================
#  mp-gromacs -- NGA/posit kernel-characterization pipeline
# =============================================================================
#  One image that covers both paths in the repo:
#
#   (A) LEAN REPRODUCTION  (default) -- regenerate the format-fit result and the
#       value-distribution histograms from the committed dataset/records/*.json.
#       Needs only numpy + matplotlib + jsonschema. Runs in seconds, no GROMACS.
#
#   (B) REPRODUCE GROMACS -- re-extract the gromacs.* records + histograms straight
#       from the Max Planck benchmark .tpr in benchmarks/ (build/build_gromacs.sh).
#       Charges/masses come from the .tpr via MDAnalysis; LJ/bonded/constraint
#       params from the .dump.txt (generated here with `gmx dump` if missing).
#       benchmarks/ is git-ignored, so mount it at run time.
#
#   (C) BUILD CHARMM -- build fresh CHARMM36 systems from PDB (build/build_charmm.sh).
#       Needs the `gmx` binary and MDAnalysis, both installed below. Stock apt
#       GROMACS is used; the parser reads charges/params from `gmx dump` text, so
#       the exact gmx version does not matter (independent of the .tpr tpx version).
#
#  Build:   docker build -t mp-gromacs .
#  Run (A): docker run --rm -v "$PWD/docs:/repo/docs" mp-gromacs
#  Run (B): docker run --rm -v "$PWD/benchmarks:/repo/benchmarks" \
#                          -v "$PWD/dataset:/repo/dataset" -v "$PWD/docs:/repo/docs" \
#                          mp-gromacs bash build/build_gromacs.sh
#  Run (C): docker run --rm -it -v "$PWD:/repo" mp-gromacs \
#               bash -lc 'TER="0\n0\n" bash build/build_charmm.sh'
# =============================================================================
FROM python:3.11-slim-bookworm

# --- system deps -------------------------------------------------------------
#  gromacs : provides `gmx` (pdb2gmx/grompp/dump) for the full CHARMM path
#  wget    : build_charmm.sh downloads PDBs from the RCSB
#  procps  : small conveniences for interactive debugging
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        gromacs \
        wget \
        ca-certificates \
        procps \
    && rm -rf /var/lib/apt/lists/*

# --- python deps (cached separately from the source for fast rebuilds) -------
WORKDIR /repo
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# --- repository source -------------------------------------------------------
COPY . .

# scripts import sibling modules (e.g. parse_gmx_dump) by directory; make the
# key package dirs importable from anywhere, and keep matplotlib headless.
ENV PYTHONPATH=/repo:/repo/scripts:/repo/dataset_spec \
    MPLBACKEND=Agg \
    PYTHONUNBUFFERED=1

RUN chmod +x reproduce_results.sh build/build_gromacs.sh build/build_charmm.sh 2>/dev/null || true

# Default: the lean, GROMACS-free reproduction a reviewer/mentor can run as-is.
CMD ["bash", "reproduce_results.sh"]
