# Lineage-specific cancer dependencies from DepMap 26Q1 (analysis code, version 2)

Code, intermediate results and figures for the manuscript *"A selectivity-based analytical strategy for identifying lineage-specific cancer dependencies from DepMap data"*.

## Important note on versions

Version 1 of this repository (scripts `01_*` to `23_*`, now in `legacy_v1/`) produced the results of the first manuscript version. While revising the manuscript we found that the gene-effect matrix and the sample annotation had been merged with an inner join that silently dropped 348 of the 1,208 cell lines (the analysis used 860 lines), and that the multiple-testing correction had been applied only to pairs that had already passed an effect-size filter. **The scripts in `legacy_v1/` are superseded and must not be used to reproduce the current results; they are kept only for transparency.** Version 2 (`src/`) analyses all 1,208 models, excludes 16 models annotated as non-cancerous (1,192 cancer cell lines), tests every gene in every eligible lineage and applies Benjamini-Hochberg adjustment over the whole test family before the effect-size rule is applied. The merge in `src/run_pipeline.py` stops with an error if any cell line is lost or lacks a lineage annotation.

## Analysis in brief

1. Merge DepMap 26Q1 `CRISPRGeneEffect` (saved as `gene_effect.csv`) with `Model.csv` on ModelID (all 1,208 IDs must match); remove models with `OncotreePrimaryDisease == "Non-Cancerous"`.
2. Lineage = `OncotreeLineage`; lineages with at least 5 cancer cell lines are eligible (26).
3. For every gene (18,531) and eligible lineage: one-sided Mann-Whitney U test (target lineage lower than all other cancer lines; asymptotic, tie and continuity corrected); missing values removed per gene; at least 5 non-missing target values (467,091 tests).
4. Benjamini-Hochberg over all tests; candidate pairs have q < 0.05, target median Chronos < -1 and selectivity (target median minus median of all other lineages) < -0.5 (94 pairs, 58 genes, 21 lineages).
5. Permutation test (1,000 label shuffles, complete procedure repeated), sensitivity analyses, genotype, TCGA, enrichment, batch/library, co-dependency and other analyses (see `src/`).
6. Replication in Sanger Project Score screens with the two data sets separated (`src/sanger_disjoint.py`): the 26Q1 `CRISPRGeneEffect` matrix is estimated from the combined Broad (Avana) and Sanger (KY) screens, so for the 314 cancer lines with a KY screen the 26Q1 values are partly derived from the Project Score data. Candidates are therefore re-derived from the 878 lines without a KY screen (identified from `ScreenSequenceMap.csv`) and tested in the Project Score data, which then share no screens and no cell lines with the discovery data. The earlier direct comparison (`src/independent_sanger.py`) is not independent and is kept for transparency only.

Adjusted p-values refer to the whole family of significant pairs; the effect-size thresholds are a prioritisation rule and do not guarantee an FDR below 5% for the prioritised subset.

## Reproducing the results

Requirements: Python 3.13 (tested), packages in `requirements.txt` (numpy 2.4.4, pandas 2.3.3, scipy 1.17.1, statsmodels 0.14.6, lifelines 0.30.3, matplotlib 3.10.8).

Input data (not included; public): download from the DepMap portal (release DepMap Public 26Q1) the files `CRISPRGeneEffect.csv` (rename to `gene_effect.csv`), `Model.csv`, `ScreenSequenceMap.csv`, `CRISPRInferredCommonEssentials.csv` and `OmicsSomaticMutations.csv`, and from the UCSC Xena GDC hub the files `TCGA-KIRC.star_counts.tsv.gz`, `TCGA-KIRC.survival.tsv.gz`, `TCGA-PAAD.star_counts.tsv.gz` and `TCGA-PAAD.survival.tsv.gz` (unzipped). SHA-256 checksums of the files that were used are in `checksums/input_files_sha256.txt`; they could not be compared with official checksums. Optional (comparison with the earlier release): `CRISPR_gene_effect.csv` of DepMap 22Q1 from figshare (https://doi.org/10.6084/m9.figshare.19139906.v1, file id 34008491), placed in a folder given as `DEPMAP22_DIR`. Optional (comparison with independent Sanger screens): `gene_effect.csv` of "Project Score Chronos" from figshare (https://doi.org/10.6084/m9.figshare.14461980.v2, file id 28340607), placed in a folder given as `SANGER_DIR`; `annotate_lit_drug.py` queries PubMed and DGIdb online (results depend on the query date; the versions used in the paper are in `results/`). Optional (adjusted survival models, progression-free interval): `Survival_SupplementalTable_S1_20171025_xena_sp` from the UCSC Xena TCGA Pan-Cancer Atlas hub (https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp), path given as `TCGA_CDR`. The Xena STAR counts are log2(count + 1) without between-sample normalisation; `part3_tcga_enrich.py` converts them back to counts and analyses log2(CPM + 1).

```bash
pip install -r requirements.txt
DEPMAP_DIR=/path/to/depmap_files TCGA_DIR=/path/to/xena_files bash run_all.sh
```

Outputs are written to `results/` (tables, JSON summaries) and `figures/` (TIFF files and supporting tables). The main analysis takes roughly 5 minutes and the permutation test roughly 10 minutes on a laptop. `results/` in this repository contains the results reported in the manuscript.

## Verification

Downstream analyses were cross-checked with a second implementation (`src/cross_check.py`, `results/cross_check_results.csv`). In addition, the primary analysis (467,091 tests; 9,805 pairs with q < 0.05; 94 candidate pairs; 1,000 permutations) and the TCGA secondary analyses were re-implemented by ChatGPT (OpenAI) from a written specification, without our code, and gave identical numbers. Analysis code in this repository was written with the help of Claude (Anthropic) and reviewed by the authors.

## Files

| Path | Content |
|---|---|
| `src/run_pipeline.py` | merge, cohort, testing, BH, candidate selection |
| `src/genome_perm.py` | fast genome-wide permutation test (validated against `run_pipeline.py`) |
| `src/sens_primary.py`, `src/ttest_alt.py` | multiplicity, threshold, top-N, minimum-n sensitivity; Student t-test screen |
| `src/part2_depmap.py`, `src/part3_tcga_enrich.py` | genotype, TP53-MDM2, batch/library, co-dependency, subsampling, subtypes; TCGA and enrichment |
| `src/cross_check.py` | second implementation of the downstream analyses |
| `src/sanger_disjoint.py`, `src/figs_sanger_disjoint.py` | separated Broad/Sanger replication (Fig 12, S7 Table) |
| `src/patient_one_model.py` | primary procedure with one model per patient (41 patients contribute 87 models) |
| `src/revision4_checks.py` | Welch t-tests, permutation p for low-n pairs, batch-adjusted CIs, DCAF7 coverage, TCGA spline and median-of-ratios sensitivity |
| `src/kirc_paired.py` | HNF1B in TCGA-KIRC tumour versus normal samples of the same patients (Wilcoxon signed-rank) |
| `src/paad_pdac_subset.py` | KRAS expression and survival in the 150 curated PDAC samples (Table S1 of Raphael et al. 2017; set `PAAD_S1`) |
| `src/cross_release_decomp.py` | 22Q1 versus 26Q1 on the same 1,036 lines: composition versus data/processing |
| `src/revision3_checks.py`, `src/ph_check.py` | MDM2 in TP53 wild-type lines, matched null for co-dependency, split-half reference correlation, BH over pairs for growth-adjusted models; proportional-hazards tests of the TCGA Cox models |
| `src/revision2_checks.py` | within-lineage co-dependency, restricted enrichment background, HNF1B/PAX8 in renal cell carcinoma lines, subsampling at the BH threshold |
| `src/figlib.py`, `src/figs_*.py`, `src/make_tables.py` | figures (PLOS ONE format) and supporting tables |
| `results/` | result tables (`all_tests.csv.gz` contains all 467,091 tests) |
| `checksums/` | SHA-256 checksums of input files |

## License

MIT (see `LICENSE`).
