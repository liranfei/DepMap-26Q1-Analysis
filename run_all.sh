#!/usr/bin/env bash
# Reproduce all results, tables and figures.
# Usage: DEPMAP_DIR=<dir with DepMap files> TCGA_DIR=<dir with Xena files> bash run_all.sh   (optional: RESULTS_DIR, FIG_DIR)
set -euo pipefail
: "${DEPMAP_DIR:?set DEPMAP_DIR}"; : "${TCGA_DIR:?set TCGA_DIR}"
mkdir -p "${RESULTS_DIR:-results}" "${FIG_DIR:-figures}"
export RESULTS_DIR="$(cd "${RESULTS_DIR:-results}" && pwd)" FIG_DIR="$(cd "${FIG_DIR:-figures}" && pwd)" DEPMAP_DIR="$(cd "$DEPMAP_DIR" && pwd)" TCGA_DIR="$(cd "$TCGA_DIR" && pwd)"
cd "$(dirname "$0")/src"
PY=${PYTHON:-python3}
$PY run_pipeline.py --data-dir "$DEPMAP_DIR" --out "$RESULTS_DIR" --top-n 0 --exclude-noncancerous     # primary analysis (about 5 min)
$PY genome_perm.py  --data-dir "$DEPMAP_DIR" --out "$RESULTS_DIR" --n 1000 --exclude-noncancerous       # permutation test of the whole pipeline
$PY prepare_mutations.py
$PY sens_primary.py; $PY ttest_alt.py; $PY part2_depmap.py; $PY part3_tcga_enrich.py; $PY cross_check.py
$PY cohort_stats.py; $PY batch_composition.py; $PY make_tables.py; $PY figs_batch1.py; $PY figs_batch2.py; $PY figs_supp.py
if [ -n "${DEPMAP22_DIR:-}" ]; then export DEPMAP22_DIR="$(cd "$DEPMAP22_DIR" && pwd)"; $PY cross_release.py; $PY figs_crossrelease.py; fi   # optional: DepMap 22Q1 comparison
if [ -n "${SANGER_DIR:-}" ]; then export SANGER_DIR="$(cd "$SANGER_DIR" && pwd)"; $PY independent_sanger.py; $PY figs_sanger.py; fi   # optional: Sanger Project Score comparison (figshare 10.6084/m9.figshare.14461980)
$PY annotate_lit_drug.py; $PY make_tables_ext.py   # optional; needs internet (PubMed, DGIdb); results are date-dependent, the stored versions are in results/
if [ -n "${TCGA_CDR:-}" ]; then $PY clinical_extra.py; fi   # optional: TCGA-CDR file (Survival_SupplementalTable_S1_20171025_xena_sp) for PFI and adjusted Cox models
$PY benchmark_metrics.py   # comparison with other lineage-selectivity measures (Table 3, S9 Table)
$PY growth_confound.py; $PY sanger_matched.py   # culture-format checks (S10 Table); effect-size adjusted Sanger comparison (needs independent_sanger.py output)
echo "ALL DONE"
