"""Hodges-Lehmann sensitivity analysis: the estimator that corresponds to the Mann-Whitney test.

For every pair that is significant in the primary analysis (q < 0.05 over the full family), the
Hodges-Lehmann estimator (median of all target-minus-other differences) is computed and compared
with the difference in medians used by the prioritisation rule. Because the rule is applied only to
significant pairs, this is enough to determine the candidate list that an HL-based rule would give.
Environment: DEPMAP_DIR, RESULTS_DIR.
"""
import os, json, numpy as np, pandas as pd

D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results")
CHRONOS_MAX, SEL_MAX, Q_MAX = -1.0, -0.5, 0.05

T = pd.read_csv(os.path.join(O, "all_tests.csv"))
sig = T[T.q_value < Q_MAX].copy()
genes = sorted(sig.Gene.unique())

ge = pd.read_csv(os.path.join(D, "gene_effect.csv"), usecols=["Unnamed: 0"] + genes)
ge = ge.rename(columns={"Unnamed: 0": "DepMap_ID"})
M = pd.read_csv(os.path.join(D, "Model.csv"))[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"]]
M.columns = ["DepMap_ID", "lineage", "disease"]
M["lineage"] = M.lineage.str.strip().str.lower()
df = ge.merge(M, on="DepMap_ID", how="left", validate="one_to_one")
assert df.lineage.notna().all()
df = df[df.disease != "Non-Cancerous"].reset_index(drop=True)

def hl(a, b, cap=4_000_000):
    """Median of all target-minus-other pairwise differences."""
    if a.size * b.size <= cap:
        return float(np.median(a[:, None] - b[None, :]))
    raise ValueError(f"Exact HL calculation exceeds cap: {a.size} x {b.size}")

rows = []
for g, grp in sig.groupby("Gene", sort=False):
    v = pd.to_numeric(df[g], errors="coerce").values
    lin = df.lineage.values
    ok = ~np.isnan(v)
    vv, ll = v[ok], lin[ok]
    for lineage, sel, med, q, n in zip(grp.Lineage, grp.Selectivity, grp.Chronos_median, grp.q_value, grp.n_target):
        m = ll == lineage
        rows.append((g, lineage, int(m.sum()), float(med), float(sel), hl(vv[m], vv[~m]), float(q)))
H = pd.DataFrame(rows, columns=["Gene", "Lineage", "n_target", "Chronos_median", "Selectivity_median_diff",
                                "Selectivity_HL", "q_value"])
H.to_csv(os.path.join(O, "hodges_lehmann_significant_pairs.csv"), index=False)

prim = (H.Chronos_median < CHRONOS_MAX) & (H.Selectivity_median_diff < SEL_MAX)
hlr  = (H.Chronos_median < CHRONOS_MAX) & (H.Selectivity_HL < SEL_MAX)
q_eq = float(np.quantile(H.Selectivity_HL, (H.Selectivity_median_diff < SEL_MAX).mean()))
matched = (H.Chronos_median < CHRONOS_MAX) & (H.Selectivity_HL < q_eq)
R = {"n_significant_pairs": int(len(H)),
     "n_prioritised_median_rule": int(prim.sum()),
     "n_prioritised_HL_rule_same_threshold": int(hlr.sum()),
     "overlap": int((prim & hlr).sum()),
     "only_median_rule": int((prim & ~hlr).sum()),
     "only_HL_rule": int((~prim & hlr).sum()),
     "spearman_median_vs_HL": float(H[["Selectivity_median_diff", "Selectivity_HL"]].corr(method="spearman").iloc[0, 1]),
     "pearson_median_vs_HL": float(H[["Selectivity_median_diff", "Selectivity_HL"]].corr().iloc[0, 1]),
     "median_abs_difference": float((H.Selectivity_HL - H.Selectivity_median_diff).abs().median()),
     "HL_posthoc_quantile_matched_threshold": q_eq}
json.dump(R, open(os.path.join(O, "hodges_lehmann_summary.json"), "w"), indent=1)
matched_summary = {"matched_threshold": q_eq,
                   "n_HL_matched": int(matched.sum()),
                   "overlap_with_94": int((prim & matched).sum()),
                   "lost": int((prim & ~matched).sum()),
                   "new": int((~prim & matched).sum())}
json.dump(matched_summary, open(os.path.join(O, "hodges_lehmann_matched.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
