"""S12 Table: formal test of the prioritisation rule (threshold_test.py).
Environment: RESULTS_DIR (reads threshold_test_pairs.csv and threshold_test_summary.json)."""
import os, json, pandas as pd
O = os.environ.get("RESULTS_DIR", "results")
P = pd.read_csv(O + "/threshold_test_pairs.csv"); S = json.load(open(O + "/threshold_test_summary.json"))
P = P[P.prioritised | (P.q_shift < 0.05)].copy()
P.insert(2, "Lineage_group", P.Lineage.map(lambda l: "blood" if l in ("lymphoid", "myeloid") else "non-blood"))
P = P.sort_values(["q_shift", "p_shift"]).rename(columns={
    "n_target": "n_target", "Chronos_median": "Median_Chronos_target", "Selectivity": "Selectivity_difference_in_medians",
    "prioritised": "Prioritised_in_primary_analysis", "tier": "Robustness_tier", "sel_ci_high": "Bootstrap_95CI_upper_selectivity",
    "p_shift": "p_selectivity_threshold", "q_shift": "q_selectivity_threshold_BH",
    "p_sign": "p_median_threshold_sign_test", "p_iut": "p_both_thresholds_IUT", "q_iut": "q_both_thresholds_BH",
    "q_iut_BY": "q_both_thresholds_BY", "formal_rule": "Both_thresholds_q_lt_0.05_BH", "formal_rule_BY": "Both_thresholds_q_lt_0.05_BY"})
defs = pd.DataFrame([
 ("Family", f"All {S['n_tests']:,} gene-lineage tests of the primary analysis; BH (and BY) adjustment over the whole family."),
 ("Selectivity threshold test", "H0: target lineage lower by no more than 0.5 (shift >= -0.5). One-sided asymptotic Mann-Whitney test "
  "(tie and continuity corrections) after adding 0.5 to the Chronos scores of the target lineage. Under a location-shift model this "
  "tests whether the target lineage is lower by more than 0.5; it addresses the location shift estimated by Hodges-Lehmann, "
  "not the difference in medians itself. Thresholded null hypotheses as in TREAT (McCarthy DJ, Smyth GK. Testing significance relative to a fold-change threshold is a TREAT. Bioinformatics. 2009;25(6):765-71)."),
 ("Median threshold test", "H0: median Chronos of the target lineage >= -1. One-sided exact sign test: probability of at least k of "
  "n target lines below -1 under Binomial(n, 0.5)."),
 ("Both thresholds", "Intersection-union test: p = max(p_selectivity_threshold, p_median_threshold); a pair is declared to meet "
  "the rule when the adjusted p is < 0.05."),
 ("Validation", f"With no shift the vectorised p-values reproduce all {S['validation']['n_ref']:,} primary p-values "
  f"(max |log10 p| difference {S['validation']['max_abs_log10p_diff']:.1e}); shifted p-values agree with "
  f"scipy.stats.mannwhitneyu for {S['validation']['shifted_checked_pairs']} pairs (max difference "
  f"{S['validation']['shifted_max_abs_log10p_diff']:.1e})."),
 ("Rows", "The 94 prioritised pairs and any further pair with q < 0.05 in the selectivity threshold test."),
], columns=["Item", "Definition"])
summ = pd.DataFrame([
 ("Pairs with q < 0.05, selectivity threshold test (BH)", S["shift_test_q_lt_0_05"]),
 ("  of which prioritised in the primary analysis", int(((P.q_selectivity_threshold_BH < 0.05) & P.Prioritised_in_primary_analysis).sum())),
 ("Pairs with q < 0.05, selectivity threshold test (BY)", S["shift_test_q_lt_0_05_BY"]),
 ("Pairs with q < 0.05 for both thresholds (BH)", S["formal_rule_BH"]),
 ("Pairs with q < 0.05 for both thresholds (BY)", S["formal_rule_BY"]),
 ("Smallest q for both thresholds (BH)", round(S["min_q_iut"], 3)),
 ("Prioritised pairs with unadjusted p < 0.05, selectivity threshold test", S["prioritised_unadjusted_p_shift_lt_0_05"]),
 ("Prioritised pairs with unadjusted p < 0.05, median threshold test", S["prioritised_unadjusted_p_sign_lt_0_05"]),
 ("Prioritised pairs with unadjusted p < 0.05, both thresholds", S["prioritised_unadjusted_p_iut_lt_0_05"]),
 ("Prioritised pairs with bootstrap 95% CI entirely below -0.5", S["boot_ci_below_threshold_among_prioritised"]),
], columns=["Quantity", "Value"], dtype=object)
with pd.ExcelWriter(O + "/S12_Table.xlsx") as w:
    summ.to_excel(w, sheet_name="Summary", index=False)
    P.to_excel(w, sheet_name="Pairs", index=False)
    defs.to_excel(w, sheet_name="Definitions", index=False)
print("S12 written", len(P), "rows")
