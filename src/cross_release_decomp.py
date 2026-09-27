"""Why do fewer 26Q1 candidates meet the effect-size rule in 22Q1?  Decomposition into cell-line composition and data/processing.

1. Whether the 22Q1 CRISPR_gene_effect matrix contains Sanger (KY-only) models, and how many of the shared lines received a KY screen in 26Q1.
2. The complete procedure on the 1,036 cancer lines present in both releases, once with 26Q1 and once with 22Q1 gene effects (restricted to genes present in
   both releases).  Same lines and genes: the difference between the two runs reflects data and processing; the difference between the primary 26Q1
   analysis and the 26Q1 run on shared lines reflects cell-line composition (and the gene set).
3. Scale of selectivity in the two releases on the shared lines (SD over all shared tests, slope of 22Q1 on 26Q1 selectivity), and the reverse direction
   (22Q1 candidates in 26Q1 on the same lines), as a check of regression towards the mean.
Output: cross_release_decomp.json, cross_release_decomp_pairs.csv."""
import os, json, numpy as np, pandas as pd
from scipy.stats import spearmanr
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); D22 = os.environ.get("DEPMAP22_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; R = {}
d26 = rp.load(D); d26 = d26[d26.primary_disease != "Non-Cancerous"].reset_index(drop=True)
ge = pd.read_csv(os.path.join(D22, "CRISPR_gene_effect.csv")); ge = ge.rename(columns={ge.columns[0]: "DepMap_ID"}); ge["DepMap_ID"] = ge.DepMap_ID.astype(str).str.strip()
meta = d26[["DepMap_ID", "lineage", "primary_disease"]]
d22 = ge.merge(meta, on="DepMap_ID", how="inner")                     # cancer lines of the 26Q1 cohort that are also in 22Q1
shared_ids = set(d22.DepMap_ID); R["shared_cancer_lines"] = len(shared_ids)

s = pd.read_csv(os.path.join(D, "ScreenSequenceMap.csv")); u = s[(s.PassesQC == True) & (s.ExcludeFromCRISPRCombined == False) & s.Library.notna()]
lib = u.groupby("ModelID").Library.agg(lambda x: set(x)); ky_only = {m for m, l in lib.items() if l == {"KY"}}; any_ky = {m for m, l in lib.items() if "KY" in l}
ids22 = set(ge.DepMap_ID)
R["KY_only_models_26Q1"] = len(ky_only); R["KY_only_models_in_22Q1_matrix"] = len(ids22 & ky_only)
R["any_KY_models_26Q1"] = len(any_ky); R["any_KY_models_in_22Q1_matrix"] = len(ids22 & any_ky)
R["shared_lines_with_KY_screen_in_26Q1"] = len(shared_ids & any_ky)

genes = sorted(set(c for c in ge.columns if c != "DepMap_ID") & set(c for c in d26.columns if c not in ("DepMap_ID", "lineage", "primary_disease")))
a26 = d26[d26.DepMap_ID.isin(shared_ids)].sort_values("DepMap_ID").reset_index(drop=True); a22 = d22.sort_values("DepMap_ID").reset_index(drop=True)
assert (a26.DepMap_ID.values == a22.DepMap_ID.values).all()
c26p, c22p = O + "all_tests_26Q1_sharedlines.csv.gz", O + "all_tests_22Q1_sharedlines.csv.gz"
if os.path.exists(c26p) and os.path.exists(c22p): t26, t22 = pd.read_csv(c26p), pd.read_csv(c22p)
else:
    t26 = rp.test_all(a26, genes, a26.lineage.values); t22 = rp.test_all(a22, genes, a22.lineage.values); t26.to_csv(c26p, index=False); t22.to_csv(c22p, index=False)
f26 = rp.select(t26); f22 = rp.select(t22); k26s = set(zip(f26.Gene, f26.Lineage)); k22s = set(zip(f22.Gene, f22.Lineage))
F = pd.read_csv(O + "final_targets.csv"); k94 = set(zip(F.Gene, F.Lineage))
meet = lambda t: t[(t.q_value < rp.Q_MAX) & (t.Chronos_median < rp.CHRONOS_MAX) & (t.Selectivity < rp.SEL_MAX)]
P = F[["Gene", "Lineage", "Selectivity", "Chronos_median"]].merge(t26[["Gene", "Lineage", "Selectivity", "Chronos_median", "q_value"]].rename(columns={"Selectivity": "Selectivity_26shared", "Chronos_median": "Chronos_median_26shared", "q_value": "q_value_26shared"}), on=["Gene", "Lineage"], how="left") \
     .merge(t22[["Gene", "Lineage", "Selectivity", "Chronos_median", "q_value"]].rename(columns={"Selectivity": "Selectivity_22shared", "Chronos_median": "Chronos_median_22shared", "q_value": "q_value_22shared"}), on=["Gene", "Lineage"], how="left")
P["meets_26shared"] = [(g, l) in k26s for g, l in zip(P.Gene, P.Lineage)]; P["meets_22shared"] = [(g, l) in k22s for g, l in zip(P.Gene, P.Lineage)]
P.to_csv(O + "cross_release_decomp_pairs.csv", index=False)
tp = P.dropna(subset=["q_value_26shared", "q_value_22shared"])
R["shared_run"] = {"genes": len(genes), "tests": len(t26), "candidates_26Q1_shared_lines": len(f26), "candidates_22Q1_shared_lines": len(f22), "overlap_candidates": len(k26s & k22s),
                   "primary94_testable_in_both": len(tp), "primary94_meet_rule_26Q1_shared": int(tp.meets_26shared.sum()), "primary94_meet_rule_22Q1_shared": int(tp.meets_22shared.sum()),
                   "primary94_median_sel_primary": float(tp.Selectivity.median()), "primary94_median_sel_26shared": float(tp.Selectivity_26shared.median()), "primary94_median_sel_22shared": float(tp.Selectivity_22shared.median()),
                   "primary94_median_chronos_26shared": float(tp.Chronos_median_26shared.median()), "primary94_median_chronos_22shared": float(tp.Chronos_median_22shared.median()),
                   "reverse_22cands_meeting_rule_in_26shared": len(k22s & k26s), "reverse_22cands_total": len(k22s),
                   "reverse_22cands_median_sel_22": float(f22.Selectivity.median()), "reverse_22cands_median_sel_26": float(t26.set_index(["Gene", "Lineage"]).loc[list(k22s)].Selectivity.median())}
M = t26.merge(t22, on=["Gene", "Lineage"], suffixes=("_26", "_22"))
R["scale"] = {"sd_selectivity_26": float(M.Selectivity_26.std()), "sd_selectivity_22": float(M.Selectivity_22.std()),
              "slope_22_on_26": float(np.polyfit(M.Selectivity_26, M.Selectivity_22, 1)[0]), "slope_26_on_22": float(np.polyfit(M.Selectivity_22, M.Selectivity_26, 1)[0]),
              "spearman": float(spearmanr(M.Selectivity_26, M.Selectivity_22)[0]),
              "sig26_sel_lt_-0.3_median_ratio_22_to_26": float((M[(M.q_value_26 < 0.05) & (M.Selectivity_26 < -0.3)].Selectivity_22 / M[(M.q_value_26 < 0.05) & (M.Selectivity_26 < -0.3)].Selectivity_26).median()),
              "sig22_sel_lt_-0.3_median_ratio_26_to_22": float((M[(M.q_value_22 < 0.05) & (M.Selectivity_22 < -0.3)].Selectivity_26 / M[(M.q_value_22 < 0.05) & (M.Selectivity_22 < -0.3)].Selectivity_22).median())}
json.dump(R, open(O + "cross_release_decomp.json", "w"), indent=1, default=float); print(json.dumps(R, indent=1, default=float))
