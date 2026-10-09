# Pre-specified checks: FERMT2 (kindlin-2) as a candidate example

Written 2026-10-09, before any of checks 1-6 below were computed. The file's SHA-256 is recorded in the analysis output so that later edits can be detected.

## Why FERMT2, and what was already known (NOT pre-specified)

SNAP23 in bowel, the earlier example, was not reproduced in the 3D organoid screens of Neiswender et al. (Nature 2026; figshare 10.6084/m9.figshare.29472362). All 17 prioritised pairs in lineages covered by those screens were then examined there (post hoc). FERMT2 was chosen after seeing these results, because it is the least-studied pair that was reproduced in 3D. The following were therefore already known when this plan was written and are reported as such, not as tests:

- 26Q1 2D: FERMT2 prioritised in six lineages (kidney, skin, soft tissue, CNS/brain, pleura, ovary/fallopian tube); CNS/brain selectivity -0.621, q < 1e-4, bootstrap 95% upper bound -0.414 (does not lie entirely below -0.5); robustness tier low (CNS), high (skin), medium (kidney).
- Project Score (lines disjoint from the 856-line discovery cohort): CNS/brain q = 3e-6, selectivity -0.48; FERMT2 was significant (q = 0.037) but not prioritised in the discovery cohort.
- RNAi (DEMETER2): CNS/brain q = 0.10 (not supported); skin q = 1e-4 and kidney q = 3e-4 (supported).
- 3D NextGen screens: CNS/brain vs other NextGen models, mean difference -0.263, one-sided p = 2e-5, BH q = 4e-4 over 17 pairs. Note that NextGen CNS models are spheroids or coated 2D cultures whereas the other NextGen models are dome organoids, so this comparison is not matched for culture format.
- PubMed records for FERMT2 with the lineage terms: 0-1 per pair; no DGIdb drug-gene interaction.

## Primary pair

FERMT2 in CNS/brain (the only FERMT2 pair that can be examined in 3D). Skin and kidney are reported with the same checks as secondary pairs; they do not enter the decision.

## Checks (computed after this file is saved)

Data: DepMap 26Q1 CRISPRGeneEffect, Model.csv, OmicsSomaticMutations, OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded; cancer models only (Non-Cancerous excluded), as in the main analysis.

1. Culture format. FERMT2 mediates integrin adhesion, so the 2D comparison with all other lines, which include suspension cultures, could inflate selectivity. Compare CNS/brain lines with adherent lines of other lineages only (GrowthPattern = Adherent), one-sided Mann-Whitney.
   PASS: selectivity <= -0.31 (at least half of the primary -0.621) and p < 0.05.
2. Genotype. OLS of FERMT2 gene effect on a CNS/brain indicator, without and with TP53, PTEN, NF1 (high/moderate impact), IDH1 and EGFR (hotspot) as covariates, HC3 standard errors.
   PASS: adjusted coefficient retains >= 70% of the unadjusted coefficient and p < 0.05.
3. Disease subtypes. Selectivity of each CNS/brain OncotreePrimaryDisease with >= 5 lines vs all other cancer lines.
   PASS: selectivity < 0 with one-sided p < 0.05 in at least two thirds of these subtypes.
4. Carcinoma comparison. Because the 3D comparison was CNS vs carcinoma organoids, repeat the 2D comparison with carcinoma lineages only (bowel, pancreas, esophagus/stomach, breast, prostate, lung, head and neck, bladder/urinary tract, kidney excluded, liver, biliary tract, ampulla of Vater, uterus, ovary/fallopian tube, cervix, thyroid). This checks whether 2D and 3D agree on the same contrast.
   PASS: selectivity < 0 with one-sided p < 0.05.
5. Paralogue buffering (direction fixed in advance). FERMT1 (epithelial kindlin) and FERMT3 (haematopoietic kindlin) may substitute for FERMT2. Prediction: higher FERMT1 expression goes with weaker FERMT2 dependency (positive Spearman rho between FERMT1 expression and FERMT2 gene effect) across all cancer lines, and also after removing lineage means.
   PASS: rho > 0 with p < 0.05 in both analyses. Reported as mechanism, not a validation criterion.
6. Co-dependency with the integrin-adhesion module. Pre-specified partner set: ITGB1, ITGAV, ITGB5, ITGA3, ITGA5, ILK, LIMS1, PARVA, TLN1, PTK2, PXN, VCL, FERMT1 (13 genes). Pearson correlation of gene effects after removing lineage means, against all genes with data in >= 300 lines.
   PASS: at least one partner among the top 20 co-dependencies of FERMT2.

## Decision rule

FERMT2 (CNS/brain) is used as the new example in the manuscript only if checks 1, 2 and 4 all pass. Checks 3, 5 and 6 are reported whatever their outcome. If the decision rule fails, this is reported to the user and FERMT2 is not used as the example.
