"""Is the lineage signal confounded by per-cell-line screen dynamic range?

Chronos scales each cell line so that the median of common-essential reference genes is about -1.
Any residual variation in that scaling (or in the spread of essential effects) could make essential
genes look more negative in some lineages. This script quantifies the residual variation, tests
whether it differs between lineages, and repeats the lineage coefficient of each prioritised pair
with the per-line quality measures as covariates. Results are reported separately for candidates
that are, and are not, on the DepMap common-essential list.
Environment: DEPMAP_DIR, RESULTS_DIR.
"""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import kruskal

D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"
ce = {x.split(" (")[0] for x in pd.read_csv(os.path.join(D, "CRISPRInferredCommonEssentials.csv")).iloc[:, 0].astype(str)}
ge = pd.read_csv(os.path.join(D, "gene_effect.csv"), low_memory=False)
ge = ge.rename(columns={ge.columns[0]: "DepMap_ID"}).set_index("DepMap_ID")
syms = np.array([c.split(" (")[0] for c in ge.columns]); isce = np.isin(syms, list(ce))
X = ge.values.astype("float32")

Q = pd.DataFrame({"DepMap_ID": ge.index,
                  "med_CE": np.nanmedian(X[:, isce], axis=1),      # scaling anchor (~ -1 by construction)
                  "sd_CE": np.nanstd(X[:, isce], axis=1),          # spread of essential effects
                  "q10": np.nanpercentile(X, 10, axis=1)})         # dynamic range of the strong tail
M = pd.read_csv(os.path.join(D, "Model.csv"))[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "GrowthPattern"]]
M.columns = ["DepMap_ID", "lineage", "disease", "growth"]; M["lineage"] = M.lineage.str.strip().str.lower()
Q = Q.merge(M, on="DepMap_ID").query('disease != "Non-Cancerous"').reset_index(drop=True)
Q.to_csv(O + "screen_quality_per_line.csv", index=False)

R = {"residual_scaling": {k: {"mean": float(Q[k].mean()), "sd": float(Q[k].sd() if hasattr(Q[k], "sd") else Q[k].std()),
                              "min": float(Q[k].min()), "max": float(Q[k].max())} for k in ("med_CE", "sd_CE", "q10")}}
R["differs_between_lineages_kruskal_p"] = {k: float(kruskal(*[g[k].values for _, g in Q.groupby("lineage") if len(g) >= 5])[1])
                                           for k in ("med_CE", "sd_CE", "q10")}
blood = Q.lineage.isin(["lymphoid", "myeloid"])
R["blood_vs_other"] = {k: {"blood": float(Q[blood][k].median()), "other": float(Q[~blood][k].median()),
                           "difference": float(Q[blood][k].median() - Q[~blood][k].median())}
                       for k in ("med_CE", "sd_CE", "q10")}

F = pd.read_csv(O + "final_targets.csv")
rows = []
for _, r in F.iterrows():
    g = r.Gene
    y = pd.to_numeric(ge[g], errors="coerce").reindex(Q.DepMap_ID).values
    d = Q.assign(y=y).dropna(subset=["y"])
    d["tgt"] = (d.lineage == r.Lineage).astype(int)
    if d.tgt.sum() < 5: continue
    m0 = smf.ols("y ~ tgt", d).fit(cov_type="HC3", use_t=True)
    m1 = smf.ols("y ~ tgt + med_CE + sd_CE + q10", d).fit(cov_type="HC3", use_t=True)
    rows.append((g.split(" (")[0], r.Lineage, float(m0.params.tgt), float(m1.params.tgt),
                 float(m1.params.tgt / m0.params.tgt), float(m1.pvalues.tgt),
                 g.split(" (")[0] in ce, r.Lineage in ("lymphoid", "myeloid")))
A = pd.DataFrame(rows, columns=["Gene_symbol", "Lineage", "coef_unadjusted", "coef_adjusted", "ratio",
                                "p_adjusted", "common_essential", "blood_lineage"])
A.to_csv(O + "screen_quality_adjusted_pairs.csv", index=False)
grp = lambda m: {"n": int(m.sum()), "median_ratio": float(A[m].ratio.median()),
                 "min_ratio": float(A[m].ratio.min()), "n_still_p_lt_0.05": int((A[m].p_adjusted < 0.05).sum())}
R["lineage_coefficient_adjusted_for_screen_quality"] = {
    "all_pairs": grp(A.index == A.index), "common_essential_genes": grp(A.common_essential),
    "other_genes": grp(~A.common_essential), "blood_lineages": grp(A.blood_lineage),
    "non_blood_lineages": grp(~A.blood_lineage)}
json.dump(R, open(O + "screen_quality_summary.json", "w"), indent=1)
print(json.dumps(R, indent=1))
