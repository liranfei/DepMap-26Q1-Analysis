"""Confidence tiers for the 94 prioritised pairs, combining checks that were already computed (rules fixed before the tiers were computed).

Five lines of evidence per pair (pass / fail / uninformative):
 1. threshold      bootstrap 95% interval of selectivity entirely below -0.5                         (revision6_pairs.csv)
 2. one model      retained with one model per patient                                               (all_tests_onepatient.csv.gz)
 3. release        22Q1: q < 0.05 with negative selectivity; not testable = uninformative              (cross_release_candidates.csv)
 4. culture        blood pairs: difference in medians against non-blood suspension lines < -0.5 (fewer than ten comparator lines = uninformative);
                   non-blood pairs: selectivity < -0.5 against all other lines except suspension lines
 5. independent    prioritised again in the discovery cohort without Project Score screens AND q < 0.05 with negative selectivity in the Project Score data;
                   not testable in the Project Score data (or in the discovery cohort) = uninformative
Tiers: high = no failure and at least three passes; low = two or more failures, or fewer than ten target lines, or the gene has data for fewer than half of
the cancer lines; medium = all others.
Output: confidence_tiers.csv, confidence_tiers_summary.json."""
import os, json, numpy as np, pandas as pd
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
gp = df.DepMap_ID.map(pd.read_csv(os.path.join(D, "Model.csv")).set_index("ModelID").GrowthPattern).fillna("Unknown").values
F = pd.read_csv(O + "final_targets.csv"); key = ["Gene", "Lineage"]
P6 = pd.read_csv(O + "revision6_pairs.csv"); T1 = pd.read_csv(O + "all_tests_onepatient.csv.gz"); k1 = set(map(tuple, rp.select(T1)[key].values))
X22 = pd.read_csv(O + "cross_release_candidates.csv"); G = pd.read_csv(O + "growth_confound_pairs.csv"); S = pd.read_csv(O + "sanger_disjoint_primary94.csv")
T = F.merge(P6[key + ["sel_ci_high", "coverage"]], on=key).merge(X22[key + ["q_value_22", "Selectivity_22"]], on=key, how="left") \
     .merge(G[key + ["median_diff_vs_suspension"]], on=key, how="left").merge(S[key + ["q_value_disc", "candidate_disc", "q_value_S", "Selectivity_S"]], on=key, how="left")
assert len(T) == len(F) == 94
nosusp = []
for _, r in T.iterrows():
    v = pd.to_numeric(df[r.Gene], errors="coerce"); m = (df.lineage == r.Lineage).values
    nosusp.append(float(v[m].median() - v[~m & (gp != "Suspension")].median()))
T["selectivity_without_suspension"] = nosusp
blood = T.Lineage.isin(["lymphoid", "myeloid"])
ev = pd.DataFrame(index=T.index)
ev["threshold"] = np.where(T.sel_ci_high < -0.5, "pass", "fail")
ev["one_model"] = np.where([tuple(x) in k1 for x in T[key].values], "pass", "fail")
ev["release"] = np.where(T.q_value_22.isna(), "uninformative", np.where((T.q_value_22 < 0.05) & (T.Selectivity_22 < 0), "pass", "fail"))
cul_b = np.where(T.median_diff_vs_suspension.isna(), "uninformative", np.where(T.median_diff_vs_suspension < -0.5, "pass", "fail"))
cul_n = np.where(T.selectivity_without_suspension < -0.5, "pass", "fail")
ev["culture"] = np.where(blood, cul_b, cul_n)
testable_S = T.q_value_S.notna() & T.q_value_disc.notna()
ev["independent"] = np.where(~testable_S, "uninformative", np.where((T.candidate_disc == True) & (T.q_value_S < 0.05) & (T.Selectivity_S < 0), "pass", "fail"))
npass = (ev == "pass").sum(1); nfail = (ev == "fail").sum(1)
caveat = (T.n_target < 10) | (T.coverage < 0.5)
T["tier"] = np.where((nfail >= 2) | caveat, "low", np.where((nfail == 0) & (npass >= 3), "high", "medium"))
T = pd.concat([T, ev.add_prefix("ev_")], axis=1); T["n_pass"] = npass; T["n_fail"] = nfail
T.to_csv(O + "confidence_tiers.csv", index=False)
R = {"tier_counts": T.tier.value_counts().to_dict(), "evidence_counts": {c: ev[c].value_counts().to_dict() for c in ev},
     "high": [f"{g.split(' (')[0]} {l}" for g, l in T[T.tier == "high"][key].values], "tier_by_blood": pd.crosstab(T.tier, blood.map({True: "blood", False: "non-blood"})).to_dict()}
json.dump(R, open(O + "confidence_tiers_summary.json", "w"), indent=1, default=str); print(json.dumps(R, indent=1, default=str))
