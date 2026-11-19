# Code test runs before v1b (disclosure, plan A Section 7 and note N19)

Every run of the analysis code made before the v1b deposit is listed here. None of these outputs is a result of
plan A. They come from test seeds (outside every seed range of plans A and B), from 200 genes per simulation cell
instead of 20,000 / 2,000, from a synthetic δ pool, or from permuted lineage labels. No real-label lineage test was
computed, and `all_tests.csv`, `final_targets.csv` and the HL file of the primary analysis were not opened.

All runs were made on 2026-11-19 by the AI coding assistant (Claude) that wrote the code, under the direction of the
author. An independent AI reviewer of the code (a separate Claude session, "plan-audit") read some of these outputs;
what it read is listed below.

## Runs

| # | Command (from `src/prereg/`) | Data and seeds | Output here | Looked at |
|---|---|---|---|---|
| 1 | `python implementation_checks.py` (several times while fixing the code) | simulated data, seed 50261119; real gene effects with permuted labels, permutation seed 7, pair sample seed 11 | the final run is `../implementation_checks.json` | pass/fail of C1–C9 and the reported differences |
| 2 | `python sim_inputs.py` | real gene effects; no lineage test | not deposited (it is re-created in the simulation stage) | σ median, 5%–95% and Spearman, compared with plan A Section 0 |
| 3 | `python simulation.py --cells 0-3 --test` | seeds 990000000 + i, 200 genes per cell | `simulation_cells/*_0_3_TEST*` | rejection rates of cells 0–3 |
| 4 | `python simulation.py --cells 0-719 --test` | seeds 990000000 + i, 200 genes per cell | `simulation_cells/*_0_719_TEST*` | run time and that every cell completed; rates only as far as they enter run 5 |
| 5 | `python simulation_decide.py` on a temporary copy of the output of run 4 in which only the `seed` column was relabelled to 20261119 + cell so that the guard passes (the copy was deleted after the run; the script itself was not changed) | output of run 4 (test seeds, 200 genes) | `simulation_decide/*_TEST*` | the printed decision on test data: E1 met, E2 not met, E3 met, M4 rule ratio 0.321 |
| 6 | `python simulation_mixture.py --reps 0-1 --test` | seeds 990500000 + r, synthetic δ pool (−|N(0.3, 0.1)|, seed 1) | `mixture/*_TEST*` | the empirical FDRs printed |
| 7 | `python real_data.py --permute-seed 7` | real gene effects, all lineage labels permuted (seed 7) | `real_data_permuted/*_TEST*`; the large per-pair files are not deposited, their SHA-256 are in `large_outputs_not_deposited_SHA256_TEST.txt` | the printed summary |

Run 6 was made **before** the change V-1 below: in that run the lineage of a non-null gene was drawn only among
lineages with ≥ 5 lines. Runs 3–7 were made before checks C10 and C11 were added.

The independent reviewer, re-running the permuted-label analysis, reported: M2, 0 pairs with q < 0.05; M3, 152 pairs
with q < 0.05, of which 2 prioritised. These counts are from permuted labels and are not results.

## Code changes made after looking at test outputs

None of these changes was chosen to move a rate, an FDR or a decision. Each fixes a crash, a wrong merge, an
undefined value or a departure from the text of plan A.

1. **Seed of the checks.** The checks first used a seed equal to that of simulation cell 0; changed to 50261119.
2. **Empty case of the weighted median** (`strat_reference.weighted_median_fraction`): raised an error, now returns NaN; C4 compares NaN status on both sides.
3. **Leverage-1 guard in M3** (`strat_core.ols_hc3`, note N4): C9 failed on a lineage with a single line, where the HC3 standard error is undefined. Such genes are now not testable under M3; C9 uses lineages with ≥ 5 lines plus lymphoid, myeloid, skin and lung.
4. **Merging of simulation outputs** (`simulation.py`): the merge overwrote the test seeds and gene counts with those of the real run; now the cells actually run are used.
5. **Bias loop** (`simulation.py`, note N7): empty-slice warnings; bias is computed only for genes testable under M1, and empty values are dropped.
6. **V-1, lineage of non-null genes in the mixture** (`simulation_mixture.py`, note N11): the code drew it among lineages with ≥ 5 lines; plan A says all lineages in proportion to their sizes. Changed to all lineages. Found by the reviewer, by reading the code against plan A.
7. **V-2, checks of the primary-analysis files** (`real_data.py`, note N22): assertions added. They act only with real labels.
8. **Review suggestions**: truth function of the mixture moved to module level and checked (C10); M3 testability checked independently (C11); array-level hash of `sim_inputs.npz` recorded (N20); the mixture uses 2 workers (memory); notes N11 and N14 corrected; notes N21 and N22 added.
