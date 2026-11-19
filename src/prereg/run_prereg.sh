#!/usr/bin/env bash
# Order of the preregistered plan A runs (plan A, Section 7). Each stage is run only after the previous deposit.
# Required: DEPMAP_DIR (DepMap 26Q1: gene_effect.csv, Model.csv, ScreenSequenceMap.csv), RESULTS_DIR.
set -euo pipefail
: "${DEPMAP_DIR:?set DEPMAP_DIR}"; : "${RESULTS_DIR:?set RESULTS_DIR}"
cd "$(dirname "$0")"; PY=${PYTHON:-python3}; STAGE=${1:?usage: run_prereg.sh v1b|simulation|real}
case "$STAGE" in
  v1b)          # before v1b: implementation checks and enumeration of the simulation cells (simulated / permuted data only)
    $PY implementation_checks.py --data-dir "$DEPMAP_DIR" --out "$RESULTS_DIR/prereg_v1b"
    $PY simulation.py --enumerate --out "$RESULTS_DIR/prereg_v1b" ;;
  simulation)   # after v1b is deposited; outputs go into v2
    O="$RESULTS_DIR/prereg_v2"
    $PY sim_inputs.py --data-dir "$DEPMAP_DIR" --out "$O"
    $PY culture_coef_fit.py --data-dir "$DEPMAP_DIR" --out "$O"
    $PY simulation.py --inputs "$O/sim_inputs.npz" --out "$O" --cells 0-719 --workers "${WORKERS:-4}"
    $PY simulation_mixture.py --inputs "$O/sim_inputs.npz" --coef "$O/culture_coefficients.csv" --out "$O" --reps 0-19 --workers 2
    $PY simulation_decide.py --results "$O/sim_cells_results_0_719.csv" --out "$O" ;;
  real)         # after v2 is deposited (needs V2_DOI); outputs go into v3
    : "${V2_DOI:?set V2_DOI}"
    $PY real_data.py --data-dir "$DEPMAP_DIR" --primary-results "$RESULTS_DIR" --out "$RESULTS_DIR/prereg_v3" \
        --v2-doi "$V2_DOI" --decision "$RESULTS_DIR/prereg_v2/sim_decision.json" --n-calibration 1000 ;;
esac
