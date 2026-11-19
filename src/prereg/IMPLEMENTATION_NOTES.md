# Implementation notes for plan A (v1b)

These notes belong to the v1b deposit (release `v2.5.1-prereg-v1b`). They list the points that plan A
(`preregistration/plan_A_stratified_procedure.md`, frozen in v1 = `v2.5.0-prereg`) does not fix and that the code
therefore has to fix. All of them are fixed here **before** any simulation cell, mixture replicate or real-data
analysis is run. The code cites them as [N1] to [N22].

The notes do not change any parameter of the plan. If one of them is later found to contradict the plan, it is
treated as a bug (plan A, Section 7).

## Data and strata

- **N1 Lineage names.** Lineage names are converted to lower case, as in the primary analysis (`run_pipeline.py`), so that the pairs match `all_tests.csv`.
- **N2 The Other stratum.**
  - It contains every `GrowthPattern` that is not Adherent or Suspension.
  - In the 1,192 cancer models with gene effects these are Mixed (65) and Unknown (41). The values Dome and Spheroid in `Model.csv` occur only in models without gene effects.
  - The library level raises an error for any library set not listed in plan A, Section 2. None occurs.
- **N3 Zero variance.** If Var[T] = 0 (all values tied within every informative stratum), the pair is not testable and does not enter the family. The number of such pairs is returned.

## M3 and M3u

- **N4 Empty culture levels and leverage-1 cells.**
  - A culture level with no line with data for a gene is dropped from that gene's model; df = n − (number of remaining columns), so df = n − 4 when all three levels occur.
  - A gene is not testable under M3/M3u if any populated design cell has leverage 1 (1 − h < 1e-8), because the HC3 standard error is then undefined. This happens only for lineages with a single line; such pairs are outside the M2 family.

## Simulation (Sections 4.1–4.3)

- **N5 Test families in the simulation.** M3 is evaluated on the genes that are testable under M2 (in the real data its family is the M2 family). M3u and M1 are evaluated on the genes with ≥ 5 target lines with data and ≥ 1 other line.
- **N6 Denominator of rejection rates.** Rejection rates use the testable genes of the cell as denominator. The number of testable genes is written next to each rate.
- **N7 Effect-size bias.**
  - Computed on the first 200 simulated genes of each cell, restricted to the genes testable under M1.
  - Expressed as (estimate − β)/σ_g; the mean and the median are reported.
  - Estimators: M1 median difference, M1 HL, M2 stratified HL, M2 secondary effect, M3 coefficient, M3u coefficient.
  - Restricting to 200 genes keeps the run time of the exact HL computations manageable. Bias is reported only, not used for any decision.
- **N8 Pooled within-lineage SD.**
  - Definition: σ_g = sqrt( Σ_l Σ_{i∈l} (y_i − ȳ_l)² / (n_g − L_g) ), over the lines with data, with L_g the number of lineages that have at least one line with data.
  - This definition reproduces the values in plan A, Section 0: median 0.133 (0.1330), 5%–95% 0.102–0.334 (0.1015–0.3338), Spearman correlation with the median −0.64 (−0.635).
  - μ_g is the median over all cancer lines with data.
- **N9 Cells used for E1–E3.**
  - E1 uses the cells single:lymphoid and single:myeloid with c_β = 0, c_δ ≠ 0 and k = 1, for both noise types. E1 is met if the M1 rate is > 0.10 in every one of these cells.
  - E2 uses the single:* scenarios (the misspecification and shared-blood scenarios are excluded), with c_β = 0, k = 1, both noise types and every c_δ.
  - E3 fallback: if no cell is in the power range [0.2, 0.9], the cell with comparator power closest to 0.5 is used; ties are broken by the lower cell number.
- **N10 "M2 loses most power" (Section 4.6, last rule).**
  - Uses the shared-blood cells with k = 1 and c_β ≠ 0, for the lymphoid and the myeloid tests.
  - The rule applies if the median over these cells of (M2 rate / M1 rate) is < 0.5.
- **N11 Mixture.**
  - Random draws are made in a fixed order; see the header of `simulation_mixture.py`.
  - round(0.02 × 18,531) = 371 non-null genes and round(0.05 × 18,531) = 927 culture-only genes.
  - The lineage of a non-null gene is drawn among all lineages, in proportion to their numbers of lines (plan A, Section 4.4). A pair whose lineage is not testable under a method simply does not enter that method's family.
  - M1 in the mixture is computed as a single-stratum Van Elteren test with ≥ 1 line per group. This is identical to the Mann–Whitney test (check C1).
  - The empirical FDR of a replicate is 0 when nothing is prioritised; the number of replicates with nothing prioritised is reported next to the mean FDR.
  - The truth function (T1, T2) is checked against hand-computed values (check C10).
- **N13 Misspecification scenario.**
  - "Mixed lines" means `GrowthPattern` = Mixed (Unknown lines receive no culture effect).
  - In the shared-blood scenario, β and k apply to both blood lineages.
  - Cells are generated in blocks of 2,000 genes. The block size is part of the random-number stream.

## Real data (Section 5)

- **N12 The 12 shared blood genes.** These are the genes prioritised under M1 in both lymphoid and myeloid. Their M4 results are reported with q from the all-gene M4 family and with q from BH within these genes.
- **N14 Effect sizes and the heterogeneity check.**
  - The stratified HL shift, the absolute effect restricted to informative strata (`abs_inf`: the median gene effect of the target lines in informative strata; not HL-based), and the heterogeneity check are computed for pairs with q < 0.05 (as plan A prescribes for the HL shift). The secondary effect size is computed for all testable pairs.
  - Heterogeneity check: the within-stratum HL shifts of informative strata with ≥ 5 target lines are compared, and a pair is flagged if two of these shifts have opposite signs (zero counts as neither sign).
- **N15 Differences for the stratified HL.** The differences target − other are computed in float64. The exact rational weights apply to the weighting, not to the differences.
- **N16 One model per patient.** The smallest ModelID by string order is kept, as in `patient_one_model.py`.
- **N17 M4 testability.** A gene is testable under M4 if it has ≥ 5 blood lines with data and ≥ 1 comparison line with data.
- **N21 Size of the M3 family.** M3 is run on the M2 family, but a pair of the M2 family is dropped from M3 when its design has a leverage-1 cell or is rank deficient [N4]. The M3 family can therefore be slightly smaller than the M2 family; both sizes are reported. M3 testability is checked against an independent determination (check C11).
- **N22 Checks of the primary-analysis files.** With real labels, `real_data.py` stops unless `final_targets.csv` has 94 pairs, all present in `all_tests.csv`; the HL file gives 68 pairs with median < −1 and HL < −0.5; 12 genes are prioritised in both lymphoid and myeloid; M4 has 41 comparison lines; and the lineage names of the three files are lower case and occur in the cohort [N1].
- **N18 Guard.** With real labels, `real_data.py` runs only when the v2 DOI and the v2 decision file are given. The comparison with the primary analysis (`all_tests.csv`, `final_targets.csv`, HL file) is read only in that mode.

## Testing before v1b

- **N19 Seeds used for testing.** These seeds are outside every seed range of plans A and B:
  - implementation checks: 50261119;
  - permuted-label checks: permutation seed 7, pair sample seed 11;
  - code test of the simulation: 990000000 + i (200 genes per cell);
  - code test of the mixture: 990500000 + r;
  - code test of `real_data.py`: permutation seed 7.

  Every code test run before v1b, what was looked at, and every code change made after looking, are listed in `results/prereg_v1b/test_runs_TEST/README.md`. Those outputs are deposited there with the suffix `_TEST`. They are not results of plan A: they come from test seeds, from 200 genes per cell, from a synthetic δ pool, or from permuted lineage labels.
- **N20 Hash of `sim_inputs.npz`.** The SHA-256 of `sim_inputs.npz` is recorded in every simulation output. A compressed `.npz` written again can differ byte-wise (zip metadata) while its arrays are identical, so the arrays, not the file bytes, are the reproducible object. Every simulation output therefore also records the SHA-256 of the arrays (`arrays_sha256` in `strat_core.py`).
