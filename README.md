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

Adjusted p-values refer to the whole family of significant pairs; the effect-size thresholds are a prioritisation rule and do not guarantee an FDR below 5% for the prioritised subset.

## Reproducing the results

Requirements: Python 3.13 (tested), packages in `requirements.txt` (numpy 2.4.4, pandas 2.3.3, scipy 1.17.1, statsmodels 0.14.6, lifelines 0.30.3, matplotlib 3.10.8).

Input data (not included; public): download from the DepMap portal (release DepMap Public 26Q1) the files `CRISPRGeneEffect.csv` (rename to `gene_effect.csv`), `Model.csv`, `ScreenSequenceMap.csv`, `CRISPRInferredCommonEssentials.csv` and `OmicsSomaticMutations.csv`, and from the UCSC Xena GDC hub the files `TCGA-KIRC.star_counts.tsv.gz`, `TCGA-KIRC.survival.tsv.gz`, `TCGA-PAAD.star_counts.tsv.gz` and `TCGA-PAAD.survival.tsv.gz` (unzipped). SHA-256 checksums of the files that were used are in `checksums/input_files_sha256.txt`; they could not be compared with official checksums. Optional (comparison with the earlier release): `CRISPR_gene_effect.csv` of DepMap 22Q1 from figshare (https://doi.org/10.6084/m9.figshare.19139906.v1, file id 34008491), placed in a folder given as `DEPMAP22_DIR`. Note that the Xena expression values are already log2(count + 1) and must not be transformed again.

```bash
pip install -r requirements.txt
DEPMAP_DIR=/path/to/depmap_files TCGA_DIR=/path/to/xena_files bash run_all.sh
```

Outputs are written to `results/` (tables, JSON summaries) and `figures/` (TIFF files and supporting tables). The main analysis takes roughly 5 minutes and the permutation test roughly 10 minutes on a laptop. `results/` in this repository contains the results reported in the manuscript.

## Verification

The primary analysis (467,091 tests; 9,805 pairs with q < 0.05; 94 candidate pairs; 1,000 permutations) was re-implemented independently from a written specification and gave identical numbers; downstream analyses were cross-checked with a second implementation (`src/cross_check.py`, `results/cross_check_results.csv`) and re-implemented independently from a written specification.

## Files

| Path | Content |
|---|---|
| `src/run_pipeline.py` | merge, cohort, testing, BH, candidate selection |
| `src/genome_perm.py` | fast genome-wide permutation test (validated against `run_pipeline.py`) |
| `src/sens_primary.py`, `src/ttest_alt.py` | multiplicity, threshold, top-N, minimum-n sensitivity; Student t-test screen |
| `src/part2_depmap.py`, `src/part3_tcga_enrich.py` | genotype, TP53-MDM2, batch/library, co-dependency, subsampling, subtypes; TCGA and enrichment |
| `src/cross_check.py` | second implementation of the downstream analyses |
| `src/figlib.py`, `src/figs_*.py`, `src/make_tables.py` | figures (PLOS ONE format) and supporting tables |
| `results/` | result tables (`all_tests.csv.gz` contains all 467,091 tests) |
| `checksums/` | SHA-256 checksums of input files |

## License

MIT (see `LICENSE`).
