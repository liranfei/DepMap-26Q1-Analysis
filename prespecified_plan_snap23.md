# Pre-specified dry-lab validation plan: SNAP23 dependency in the bowel lineage

Written and saved **before any of the analyses below were run** (2026-10-05).
Target chosen from the 94 prioritised pairs of the DepMap 26Q1 analysis because it is
(i) in the high robustness tier, (ii) sparsely covered in the literature (5 PubMed records,
no DGIdb drug-gene interaction), (iii) without a genotype explanation in S2 Table,
(iv) in an adherent lineage, so it is not subject to the suspension-culture confounding of the
blood lineages, and (v) replicated in the Sanger Project Score screens of non-overlapping lines.

Data: local DepMap 26Q1 files (gene_effect.csv = CRISPRGeneEffect, Model.csv,
OmicsSomaticMutations.csv, OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv),
Sanger Project Score Chronos (sanger_gene_effect.csv.gz), DEMETER2 RNAi
(D2_combined_gene_dep_scores.csv). Conventions follow the manuscript: non-cancerous models
excluded; lineage = OncotreeLineage (lower-cased); one-sided Mann-Whitney (target lower);
genotype from default-entry mutations (Hotspot == True; VepImpact HIGH/MODERATE for
protein-altering); OLS with HC3 robust standard errors and a t reference distribution.

## Checks and pass criteria (fixed in advance)

| # | Check | Pass criterion |
|---|---|---|
| V1 | Primary recomputation: bowel vs all other cancer lines | one-sided p < 0.05 after BH over the full family (q < 0.05), median Chronos < -1, selectivity < -0.5, and the bootstrap 95% percentile interval of selectivity entirely below -0.5 (1,000 resamples, seed 0) |
| V2 | Genotype independence: OLS gene effect ~ bowel + KRAS hotspot + APC protein-altering + TP53 protein-altering + BRAF hotspot + PIK3CA hotspot | bowel coefficient keeps >= 70% of its unadjusted value and p < 0.05 |
| V3 | Culture format: bowel vs adherent non-bowel lines only | selectivity < -0.5 and one-sided p < 0.05 |
| V4 | Disease subtype: each bowel subtype with >= 5 lines vs all non-bowel lines | selectivity < -0.5 in every tested subtype (not driven by one subtype) |
| V5 | Paralogue buffering (new mechanistic prediction, direction fixed in advance): higher expression of the paralogues SNAP25 and SNAP29 should make cells LESS dependent, i.e. a POSITIVE Spearman correlation between paralogue expression and the SNAP23 gene effect | positive rho with p < 0.05 for at least one paralogue, across all cancer lines |
| V6 | Co-dependency: Pearson correlation of the SNAP23 gene effect with every other gene after subtracting lineage means; top 20 partners inspected | at least one gene of the pre-specified SNARE/vesicle-fusion set appears in the top 20: STX1A, STX2, STX3, STX4, STX7, STX8, VAMP2, VAMP3, VAMP4, VAMP7, VAMP8, SNAP25, SNAP29, NAPA, NAPB, NSF, STXBP1, STXBP2, STXBP3, SEC22B |
| V7 | Orthogonal data: Sanger Project Score (bowel vs rest) and DEMETER2 RNAi (bowel vs rest) | reported for both; pass when one-sided p < 0.05 with negative selectivity |
| V8 | Expression sanity: SNAP23 expressed in bowel lines | median log2(TPM + 1) > 1 |

## Overall verdict rule (fixed in advance)

- **Validated**: V1, V2, V3 and V4 all pass, and at least two of V5, V6, V7 pass.
- **Partly validated**: V1-V4 pass but fewer than two of V5-V7 pass.
- **Not validated**: any of V1-V4 fails.

Every result is reported whether or not it passes, including failures. No check is added,
removed or re-thresholded after seeing the results; anything computed afterwards is labelled
post hoc.
