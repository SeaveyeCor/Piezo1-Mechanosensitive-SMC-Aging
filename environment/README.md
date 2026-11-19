# Environments

## R (Seurat pipeline, post-processing, 04_cell_level scripts)
`R_sessionInfo.txt`: R 4.4.2 (Windows) session used for the analysis. Install the listed package versions (for example with `remotes::install_version`).

## Python (RNA velocity, `03_velocity/`)
`requirements_velocity.txt` holds the exact package set of the Python 3.12 environment used for the velocity analysis
(scvelo 0.3.4, scanpy 1.12, anndata 0.12.10, numpy 1.26.4, scikit-misc 0.5.2). The set has one declared-dependency
conflict (fast-array-utils 1.3.1 declares numpy>=2), so install it without dependency resolution:

    python3.12 -m venv venv
    venv/bin/pip install --no-deps -r environment/requirements_velocity.txt

STAR 2.7.11b performs the alignment (`03_velocity/run_starsolo.sh`); velocyto is not needed.
