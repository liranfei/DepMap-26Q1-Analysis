# Pre-specified plan A: culture-format-stratified lineage test (main analysis of the BMC Bioinformatics version)

Status: DRAFT 9 (final), revised 2026-11-19 after nine reviews, all by an AI model (Claude, Anthropic) in sessions separate from the one that drafted the plans: two rounds in one session and seven rounds of audit in another; not independent human review. The content of plan A has not changed since draft 6, except the release-tag mapping in Section 7. Frozen on 2026-11-19 as Zenodo v1 (repository release v2.5.0-prereg); the SHA-256 of this file is in `SHA256SUMS.txt`. This plan and plan B are frozen together as Zenodo version 1 (with a DOI) before any analysis in Sections 4–6 is run. Its SHA-256 is written into every output file of this plan. All parameters below are fixed by this text; nothing is left to be chosen in the scripts.

## 0. What was already known when this plan was written (not pre-specified)

The idea of stratifying by culture format was proposed **after** the following results had been seen, so this is a post hoc improvement of the original analysis. Only the new runs described below are pre-specified.

- Primary analysis (DepMap 26Q1, 1,192 cancer models, 467,091 gene–lineage tests, one-sided Mann–Whitney vs all other cancer lines, BH): 9,805 pairs with q < 0.05, 94 prioritised pairs (q < 0.05, median < −1, selectivity < −0.5). With the Hodges–Lehmann shift in place of the median difference and the same thresholds, 68 pairs (S9 Table).
- Culture-format composition of the cohort: Adherent 923, Suspension 163, Mixed 65, Unknown 41. Lymphoid: 83 suspension, 4 adherent, 7 other; Myeloid: 39 suspension, 3 adherent, 1 other. 41 non-blood suspension lines (Lung 16, Esophagus/Stomach 7, Bone, Breast and CNS/Brain 3 each, Eye and Head and Neck 2 each, others 1 each).
- Permutation of lineage labels: unrestricted, mean 0.0 (max 0) prioritised pairs per permutation; within the 4-level GrowthPattern strata (Adherent/Suspension/Mixed/Unknown), mean 11.3 (max 25); within pDNA batch, mean 0.004 (max 1); 1,000 permutations each.
- OLS adjustment for the 4-level GrowthPattern (HC3) of the 94 pairs: median ratio of adjusted to unadjusted lineage coefficient 0.99 for the 52 non-blood pairs and 0.56 for the 42 blood-lineage pairs; 0.84 for the blood pairs when the suspension effect is estimated from the 41 non-blood suspension lines. Adjustment for batch/library: median ratio 1.00 (min 0.75).
- Model-level counts of lines in informative strata (Section 2) and the distribution of per-gene medians and within-lineage SDs of gene effect (Section 4) were computed before writing this plan; no lineage test statistic was computed. Per-gene pooled within-lineage SD: median 0.133 (5%–95%: 0.102–0.334); Spearman correlation with the per-gene median −0.64.

## 1. Question and estimand

The original procedure compares a lineage with all other cancer lines. When a lineage is (almost) only grown in one culture format, this contrast mixes lineage with culture format. We estimate the **within-culture-format lineage difference**: the lineage-vs-rest difference among lines grown in the same format, pooled over formats.

This is a different estimand, not the "true" lineage effect. Blood cells grow in suspension because of their biology, so stratification also removes part of the lineage signal. The Van Elteren weights are dominated by the stratum holding most target lines; for blood lineages this is the suspension stratum (about 79% of the weight for Lymphoid), where the other blood lineage makes up about half (Lymphoid: 39 of 80) to two thirds (Myeloid: 83 of 124) of the comparators, so the M2 estimand for blood lineages is partly a lymphoid-vs-myeloid contrast. This is stated next to the main results, and M4 is reported alongside. Results are reported as "adjusted for culture format", never as "deconfounded" or "true".

## 2. Data and strata

- Data: DepMap 26Q1 CRISPRGeneEffect (`gene_effect.csv`, SHA-256 e610a4cefb13a82b5b256b47eb08b63ff14843f8dbd0fb164bc0a32688e5b89e), `Model.csv` (ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8), `ScreenSequenceMap.csv` (cd053d229750d2e35abb1dcd71410ab336f3f9eef122a8168108cb57f98359d3); cancer models only (Non-Cancerous excluded); all genes; lineage = OncotreeLineage, as in the primary analysis.
- Main strata (3): `GrowthPattern` = Adherent; Suspension; Other (Mixed, Unknown, missing and any other value).
- Sensitivity strata: culture format × screening library. Library level per model = the set of libraries among its 2DS screens that pass QC and are not drug-treated, engineered or excluded from the combined dataset (same screen filter as `genome_perm_strat.py`): "Avana" (822 models), "KY" (116), "Avana+KY" (195), "other" (any set containing TKOv3, Brunello or Humagne-CD; 6 models; Humagne-CD does not occur among the filtered 26Q1 screens), "unknown" (no such screen; 53). Using the set avoids a tie rule.
- Per gene, only lines with a non-missing gene effect are used.
- Informative stratum (per gene × lineage): at least 3 target-lineage lines and at least 3 other lines in that stratum (sensitivity: at least 5 each). Non-informative strata are dropped for that test.
- A gene × lineage pair is testable if the informative strata together contain at least 5 target-lineage lines (the MIN_N of the primary analysis). The test family is all testable pairs; its size is reported and will differ from 467,091.
- For each pair, the fraction of the lineage's lines with data that fall in informative strata is reported; pairs below 0.80 are flagged "estimate covers only part of the lineage". Expected at model level: Eye flagged (0.70); Ampulla of Vater untestable (4 lines in informative strata).

## 3. Methods compared

M1. **Unstratified (original)**: one-sided Mann–Whitney vs all other cancer lines; existing `all_tests.csv`, not rerun. Changes are reported as a chain: M1 (median difference, 94 pairs) → M1 (HL shift, 68 pairs) → M1 (HL shift, restricted to the M2 test family, BH recomputed in that family; HL computed, with the exact method of `hodges_lehmann.py`, for any pair that becomes significant only in the smaller family) → M2.

M2. **Stratified (new procedure)**: one-sided Van Elteren test over informative strata s, with n_t,s target lines, n_o,s other lines, n_s = n_t,s + n_o,s, W_s the target rank sum within s (average ranks for ties), t_k the sizes of tied groups within s:
  - T = Σ_s W_s/(n_s + 1); E[T] = Σ_s n_t,s/2;
  - Var[T] = Σ_s n_t,s·n_o,s / [12 (n_s + 1)²] · [(n_s + 1) − Σ_k (t_k³ − t_k) / (n_s (n_s − 1))];
  - continuity correction c = 0.5 · max_s 1/(n_s + 1); p = Φ((T − E[T] + c) / √Var[T]) (alternative: target lineage lower).
  - With a single stratum this equals `scipy.stats.mannwhitneyu(..., alternative="less", method="asymptotic")`; the script asserts this on simulated data. Before running (step 2 of Section 7), the implementation is also checked against a separately written per-stratum implementation (agreement to 1e-10) and against exact within-stratum permutation p values on small simulated examples (the exact p is compared with the normal-approximation p computed with c = 0; acceptable: |Δp| < 0.01 or relative difference < 10%); results are archived. (R is not installed on the analysis machine, so `coin` is not used.) With several strata, c is set by the smallest informative stratum and is slightly conservative there; in the culture × library sensitivity analysis, where strata can have only 6 lines, this is noted with the results.
  - BH over the M2 test family.
  - Primary effect size: stratified Hodges–Lehmann shift, the shift Δ for which T equals its null expectation. Because each target–other pair in stratum s contributes 1/(n_s + 1) to T, this equals the weighted median of all within-stratum pairwise differences (target − other), each weighted 1/(n_s + 1). Computed with exact rational weights (integer weights L/(n_s + 1), L the least common multiple of the (n_s + 1)): the smallest difference at which the cumulative weight first reaches half the total; if it equals half exactly, the mean of that difference and the next distinct difference. With one stratum this is the ordinary HL shift. Computed only for pairs with q < 0.05 (as in the primary analysis).
  - Secondary effect size (used for exploratory rankings over all pairs): weighted mean of within-stratum median differences, weights n_t,s·n_o,s / (n_s + 1).
  - Heterogeneity check (secondary): where a lineage has ≥ 5 target lines in each of ≥ 2 strata, within-stratum HL shifts are compared and pairs with opposite signs are flagged (the test assumes a common shift across strata).
  - Absolute effect: median gene effect of all target-lineage lines with data, as in the primary analysis (version restricted to informative strata reported separately).
  - Prioritised: q < 0.05, absolute effect < −1, stratified HL shift < −0.5.

M3. **OLS with covariate**: gene effect ~ target-lineage indicator + culture-format indicators (3 levels), HC3 standard errors, one-sided p for the lineage coefficient < 0 from the t distribution with n − 4 degrees of freedom (as `use_t=True` in the existing scripts), BH over the same family as M2. Prioritised: q < 0.05, absolute effect < −1, coefficient < −0.5.

M4. **Blood vs non-blood suspension (secondary)**: blood (Lymphoid + Myeloid merged) vs the 41 non-blood suspension lines, one-sided Mann–Whitney, all genes (one test per gene), BH over this family. Reported: q, median difference, blood median; "supported" = q < 0.05, median difference < −0.5 and blood median < −1. The comparator group is heterogeneous (16 of 41 are lung lines); stated with the results. This addresses dependencies shared by the two blood lineages, which M2 tends to cancel.

## 4. Simulation with known truth (run first; does not use real lineage effects)

Purpose: show what each method does under culture-format confounding with a known answer. This is the main evidence for the procedure.

### 4.1 Generating model

Real lineage labels and real 3-level culture formats of the 1,192 models. For each simulated gene g:
y = μ_g + β·[target] + δ·[suspension] + σ_g·k·ε

- (μ_g, σ_g, missing-value mask): sampled **jointly**, with replacement, from the real genes in `gene_effect.csv` (per-gene median over all cancer lines; per-gene pooled within-lineage SD; the gene's pattern of missing lines), so that their correlation and real missingness are kept. No lineage test is used.
- Effects are in units of each gene's SD: β = c_β·σ_g, δ = c_δ·σ_g.
- ε: (i) t with 5 df scaled to SD 1; (ii) skewed: log-normal with sdlog = 0.5, standardised to **mean 0** and SD 1, sign-flipped so that the long tail is negative.
- k = 1 for non-target lines; for target lines k ∈ {1, 1.5}.
- Target lineages: Lymphoid; Myeloid; Skin; Lung. **Shared blood scenario**: β applied to Lymphoid and Myeloid lines at the same time; tests still run per lineage (Lymphoid vs rest, Myeloid vs rest) for M1–M3; M4 run as well.
- Misspecification scenario (reported only, not used for decisions): culture effect differs by lineage (δ/2 for suspension lines of Lung, the largest non-blood suspension group) and Mixed lines receive δ/2.

### 4.2 Grid and replication

- c_β ∈ {0, −0.1, −0.2, −0.35, −0.5}; c_δ ∈ {0, −0.25, −0.5, −1}.
- Genes per cell: 20,000 for null cells (c_β = 0); 2,000 for the other cells.
- Seeds: cell number i (in the order of the enumeration written in the script header, archived at v1b) uses seed 20261119 + i.
- Power = rejection rate at unadjusted one-sided α = 0.05 (not after BH). Monte Carlo intervals: two-sided 95% Clopper–Pearson.

### 4.3 Null hypothesis under unequal spread

With k = 1 the location-shift model holds, and β = 0 is a null for every method. With k = 1.5 and symmetric t noise (mean = median), β = 0 is still a null for every method, but rank tests are expected (known in advance) to exceed nominal type I error under unequal variance: by normal approximation, M1 about 0.08 for all four target lineages; M2 about 0.05 for Lymphoid (balanced groups in the suspension stratum) and about 0.08 for the others; OLS with HC3 is expected to be more robust.

Known in advance about efficiency: rank tests are more efficient than t-based tests under the noise used here (asymptotic relative efficiency of Wilcoxon vs t about 1.24 for t with 5 df and about 1.58 for the log-normal), and with negatively skewed noise the one-sided t test is somewhat conservative. M3 is therefore compared with its own unadjusted version, not with M1 (E3 below). With k = 1.5 and skewed noise, no single null holds for all methods; these cells are reported descriptively and enter no criterion.

### 4.4 Genome-like mixture

- Test family of the real size: 18,531 simulated genes (the real number) × all lineages, keeping every pair that is testable under the sampled missing-value masks by the rules of each method; BH over the whole family as in the real analysis; noise (i), k = 1.
- 2% of genes non-null: β drawn with equal probability from {−0.5, −1.0, −1.5} (absolute units, so that the absolute and selectivity thresholds are meaningful), in one lineage drawn in proportion to lineage size; δ = 0.
- 5% of genes culture-only: β = 0, δ drawn from the negative suspension coefficients of an OLS fit, per real gene, of gene effect on lineage fixed effects plus the 3-level culture indicators (absolute units); this fit is part of the v1b scripts and is run after v1b and before the simulation (estimation noise makes this distribution slightly too wide, i.e. conservative).
- All other genes: β = 0, δ = 0. Effect of the Other stratum: 0 for all genes.
- 20 replicates, seeds 30261119 + r (seed ranges in plans A and B do not overlap: cells 20261119 + i, mixture 30261119 + r, permutations 42 + b, plan B 40261120–40261122).
- Empirical FDR of the full prioritisation rule under M1, M2 and M3, under two definitions of truth:
  - (T1) within-culture estimand: a prioritised pair is true if β ≤ −0.5 and μ_g + β < −1. This definition favours the stratified methods by construction and is labelled so.
  - (T2) lineage-vs-rest estimand (the M1 estimand), defined on means: with f_target and f_rest the suspension fractions of the target lineage and of all other lines, computed over the lines with data for that gene, true if β + δ·(f_target − f_rest) < −0.5 and μ_g + β + δ·f_target < −1.

### 4.5 Pre-specified expectations (stated so that they can fail)

- E1: with c_δ ≠ 0 and a blood target (single-lineage scenarios), M1 type I error > 0.10. If E1 fails, the simulation has not reproduced the confounding; this is reported and discussed, and the decision rule is still applied.
- E2: M2 and M3 type I error: upper Clopper–Pearson bound ≤ 0.06 in every null cell with k = 1 (both noise types; single-lineage scenarios; the shared-blood scenario is excluded because its null depends on whether the other blood lineage is a comparator).
- E3 (power cost of the adjustment itself): for Skin and Lung, in **every** cell with c_β ∈ {−0.1, −0.2, −0.35}, k = 1, where the unadjusted comparator's power is between 0.2 and 0.9: M2 power ≥ 0.9 × M1 power, and M3 power ≥ 0.9 × the power of M3u (the same OLS with HC3 without culture indicators, run in the simulation only). If no such cell exists for a lineage, E3 is assessed in the cell where the comparator's power is closest to 0.5.
- Reported, not used for the decision: type I error with k = 1.5 and symmetric noise (Section 4.3), the shared-blood scenario, the misspecification scenario, bias of the effect sizes, and the mixture FDR.

### 4.6 Decision rule (applied to the simulation results only, before any real-data result exists)

- M2 meets E2 and E3 → M2 is the recommended adjusted procedure.
- M2 fails either and M3 meets both → M3 is recommended.
- Both meet both → M2 is recommended (no distributional assumption); M3 is reported as agreeing or not.
- Neither meets both → no adjusted method is recommended; both are reported. Wherever this plan or plan B refers to "the recommended method", M2 is then used and labelled "not recommended by the simulation".
- If M2 loses most power in the shared-blood scenario, M4 is presented as the necessary complement for blood lineages.
- The recommendation is not changed after the real-data results are seen.

## 5. Real-data analyses (run after the simulation and after plan B's list is frozen)

Primary outputs (reported whatever they show):
1. Number of prioritised pairs under M2; of the 94 M1 pairs, how many are retained and how many lost; new pairs gained; all split into blood and non-blood lineages; decomposed along the M1 chain of Section 3.
2. Agreement between M2 and M3 (pairs prioritised by both, by one only).
3. M4: results for all genes, and for the 12 genes prioritised in both blood lineages under M1 (12 M4 tests).

Secondary:
4. Sensitivity analyses for M2: strata culture × library (small strata such as Suspension × KY, 7 lines, are expected to be dropped; the number of lines dropped is reported); "at least 5 per stratum"; Mixed, Unknown and missing GrowthPattern excluded (annotation quality varies and may partly follow how lineages were curated; stated as a limitation).
5. M2 repeated with one model per patient (same rule as `patient_one_model.py`: the smallest ModelID is kept).
6. Calibration check: whole M2 procedure under unrestricted permutation of lineage labels (1,000 permutations, seeds 42 + b as in `genome_perm.py`). Expectation: mean number of prioritised pairs ≤ 1 per permutation (M1: 0.0). If exceeded, this is reported and investigated as a possible implementation error before any further use. The within-stratum permutation is the null of M2 itself, so low counts there are not presented as an advantage.
7. Lineages that become untestable and pairs flagged under Section 2 are listed.

## 6. Interpretation rules (fixed in advance)

- A pair is "retained" if it is prioritised by the recommended method (Section 4.6). "Lost" pairs are not called false; they are pairs whose selectivity depends on the comparison with lines of another culture format.
- No claim of "first" or "novel test": Van Elteren (1960) is a standard method. The contribution is a pre-specified, genome-wide, culture-aware procedure with a simulation benchmark and a comparison with the unstratified approach and with regression adjustment; the relation to linear-model lineage tests is discussed in the manuscript. A statement that no such procedure was found is made only if a literature search with archived queries, dates and hit counts is deposited with Zenodo v1; otherwise no such statement is made.
- Limitations stated in advance: lines within a lineage are not exchangeable (subtypes, patients of origin); GrowthPattern is a curated annotation; the simulation uses real missingness but a simplified generating model.
- Wording fixed in advance for the main possible outcomes:
  - Most blood-lineage pairs lost: "the selectivity of these pairs depends on comparison with lines grown in a different culture format".
  - Most pairs retained: "the prioritised pairs are largely robust to culture-format adjustment".
  - Pairs gained: "adjustment reveals pairs whose selectivity was masked by culture-format differences in the comparison group".
  - M2 and M3 disagree: both lists are reported, with the simulation result that favours one or neither.
- The recommended method is chosen from the simulation (Section 4), not from positive-control performance (plan B); M1, M2 and M3 are all reported.

## 7. Freezing order and deviations

1. **v1**: plans A and B; plan B's screening materials (Section 5 of plan B); the computing environment (Python 3.13.12 with numpy 2.4.4, pandas 2.3.3, scipy 1.17.1, statsmodels 0.14.6 and the remaining packages in a lock file, with its SHA-256). Optionally also registered on OSF for an independent registry timestamp.
2. **v1b** (before any simulation is run): the simulation and M2/M3/M4 scripts, the enumeration of simulation cells, the implementation checks of Section 3 (single-stratum equality, independent implementation, exact permutation) and their results. Scripts are tested only on simulated data and on data with permuted lineage labels.
3. Simulation run (Section 4).
4. **v2**: simulation outputs; the written statement of the recommended method derived from them under Section 4.6, with SHA-256 (time evidence that it preceded any real-data result); plan B's screening products (Section 5 of plan B).
5. Real-data analyses of Section 5; outputs and their SHA-256 deposited as **v3**.
6. Plan B rankings (R0–R2 use the existing `all_tests.csv`, whose SHA-256 is in v2; stratified comparisons use the v3 outputs).

Nobody (authors, screeners, adjudicator) opens `all_tests.csv`, `final_targets.csv`, S2 Table or the outputs of step 5 to look up positive controls before v2 is deposited.

Release tags: "Zenodo v1" in these plans is the repository release `v2.5.0-prereg`; v1b, v2 and v3 are the releases `v2.5.1-prereg-v1b`, `v2.5.2-prereg-v2` and `v2.5.3-prereg-v3`. Each has its own Zenodo version DOI (shown on the Zenodo record and the GitHub release page); all are versions of the same Zenodo concept record as the earlier analysis releases. A bug fix or deviation is released under a further patch number with the suffix of the stage it belongs to.

Deviations: a **bug** is code that does not do what this text says; it is fixed in a new Zenodo version and results before and after the fix are reported. Any situation this text does not determine, which requires a new choice, is a **deviation**: it is recorded with its reason in a new Zenodo version and the manuscript reports results with and without it. Use of AI tools in writing the scripts and this plan is documented in Methods by the authors.
