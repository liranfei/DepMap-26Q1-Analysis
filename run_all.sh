#!/usr/bin/env bash
# Reproduce all results, tables and figures, in dependency order.
# Required: DEPMAP_DIR (DepMap 26Q1 files), TCGA_DIR (Xena TCGA-KIRC/PAAD star_counts and survival files, unzipped).
# Optional: DEPMAP22_DIR (22Q1 CRISPR_gene_effect.csv), SANGER_DIR (Project Score Chronos gene_effect.csv), RNAI_DIR (DEMETER2 D2_combined_gene_dep_scores.csv), TCGA_CDR (TCGA-CDR table),
#           PAAD_S1 (Table S1 of Raphael et al. 2017), RUN_ANNOTATION=1 (re-query PubMed/DGIdb; otherwise the stored results/annotation_* files are used),
#           RESULTS_DIR, FIG_DIR.  Internet access is needed for the Enrichr gene-set libraries.
set -euo pipefail
: "${DEPMAP_DIR:?set DEPMAP_DIR}"; : "${TCGA_DIR:?set TCGA_DIR}"
mkdir -p "${RESULTS_DIR:-results}" "${FIG_DIR:-figures}"
export RESULTS_DIR="$(cd "${RESULTS_DIR:-results}" && pwd)" FIG_DIR="$(cd "${FIG_DIR:-figures}" && pwd)" DEPMAP_DIR="$(cd "$DEPMAP_DIR" && pwd)" TCGA_DIR="$(cd "$TCGA_DIR" && pwd)"
cd "$(dirname "$0")/src"
PY=${PYTHON:-python3}

# 1. primary analysis and permutation tests
$PY run_pipeline.py --data-dir "$DEPMAP_DIR" --out "$RESULTS_DIR" --top-n 0 --exclude-noncancerous      # primary analysis (about 5 min)
$PY genome_perm.py  --data-dir "$DEPMAP_DIR" --out "$RESULTS_DIR" --n 1000 --exclude-noncancerous        # unrestricted permutation of the whole procedure
for ST in growth batch growth_batch patient; do                                                          # stratified and patient-block permutations
    $PY genome_perm_strat.py --data-dir "$DEPMAP_DIR" --out "$RESULTS_DIR" --n 1000 --exclude-noncancerous --strata $ST
done
$PY patient_one_model.py                                                                                 # one model per patient

# 2. sensitivity, genotype, TCGA, enrichment and other downstream analyses
$PY prepare_mutations.py
$PY sens_primary.py; $PY ttest_alt.py; $PY part2_depmap.py
$PY part3_tcga_enrich.py                                                                                 # TCGA (log2 CPM) and enrichment (Enrichr, internet)
$PY kirc_paired.py
$PY cross_check.py; $PY cohort_stats.py; $PY batch_composition.py; $PY data_description.py
$PY benchmark_metrics.py                                                                                 # Table 3, S9 Table
$PY revision2_checks.py                                                                                  # co-dependency, restricted enrichment background, RCC subsets, subsampling
$PY growth_confound.py                                                                                   # culture-format checks (S10 Table)
$PY revision3_checks.py                                                                                  # MDM2 without TP53 variant, co-dependency null, split-half reference, BH for growth models
if [ -n "${TCGA_CDR:-}" ]; then
    $PY clinical_extra.py; $PY ph_check.py                                                               # PFI, adjusted Cox, proportional hazards
    if [ -n "${PAAD_S1:-}" ]; then $PY paad_pdac_subset.py; fi                                           # curated PDAC samples
fi

# 3. comparisons with other data sets
if [ -n "${DEPMAP22_DIR:-}" ]; then
    export DEPMAP22_DIR="$(cd "$DEPMAP22_DIR" && pwd)"
    $PY cross_release.py; $PY cross_release_decomp.py; $PY figs_crossrelease.py                          # Fig 11, S6 Table
fi
if [ -n "${SANGER_DIR:-}" ]; then
    export SANGER_DIR="$(cd "$SANGER_DIR" && pwd)"
    $PY independent_sanger.py                                                                            # earlier, non-independent comparison (last sheet of S7 Table)
    $PY sanger_disjoint.py; $PY figs_sanger_disjoint.py                                                  # separated design: Fig 12, S7 Table
fi
if [ -n "${TCGA_CDR:-}" ] && [ -n "${SANGER_DIR:-}" ]; then
    $PY revision4_checks.py                                                                              # Welch, low-n permutation, batch CIs, spline, median-of-ratios
fi

if [ -n "${SANGER_DIR:-}" ] && [ -n "${DEPMAP22_DIR:-}" ]; then
    $PY revision5_checks.py; $PY revision6_checks.py; $PY confidence_tiers.py                            # need the 22Q1 and Project Score outputs above
    if [ -n "${RNAI_DIR:-}" ]; then export RNAI_DIR="$(cd "$RNAI_DIR" && pwd)"; $PY rnai_validation.py; fi   # DEMETER2 RNAi check
fi

# 4. literature annotation (date-dependent; stored versions are used unless RUN_ANNOTATION=1)
if [ "${RUN_ANNOTATION:-0}" = "1" ]; then $PY annotate_lit_drug.py; fi
$PY make_tables_ext.py                                                                                   # S8 Table

# 5. tables and figures
$PY make_tables.py; $PY make_tables_s6_s9_s10.py                                                         # S1-S6, S9, S10 Tables
$PY figs_batch1.py; $PY figs_batch2.py; $PY figs_supp.py                                                 # Figs 1-10, S1-S5 Figs
echo "ALL DONE"
