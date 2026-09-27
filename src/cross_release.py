"""Comparison with the earlier DepMap 22Q1 release.  The same procedure is applied to the 22Q1 CRISPR gene effect (CRISPR_gene_effect.csv);
lineage and the non-cancerous flag are taken from the 26Q1 Model.csv through the (persistent) model identifiers, so that only the release differs."""
import os, json, gzip, numpy as np, pandas as pd
from scipy.stats import spearmanr
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); D22 = os.environ.get("DEPMAP22_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"
ge = pd.read_csv(os.path.join(D22, "CRISPR_gene_effect.csv")); ge = ge.rename(columns={ge.columns[0]: "DepMap_ID"}); ge["DepMap_ID"] = ge.DepMap_ID.astype(str).str.strip()
meta = pd.read_csv(os.path.join(D, "Model.csv"))[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"]]; meta.columns = ["DepMap_ID", "lineage", "primary_disease"]
n_raw = len(ge); df = ge.merge(meta, on="DepMap_ID", how="left"); n_unmatched = int(df.lineage.isna().sum()); df = df[df.lineage.notna()]
n_nc = int((df.primary_disease == "Non-Cancerous").sum()); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True); df["lineage"] = df.lineage.str.strip().str.lower()
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
t = rp.test_all(df, genes, df.lineage.values); f = rp.select(t)
T26 = pd.read_csv(O + ("all_tests.csv" if os.path.exists(O + "all_tests.csv") else "all_tests.csv.gz")); F26 = pd.read_csv(O + "final_targets.csv"); ids26 = set(rp.load(D).query("primary_disease != 'Non-Cancerous'").DepMap_ID)
lc = df.lineage.value_counts(); M = T26.merge(t, on=["Gene", "Lineage"], suffixes=("_26", "_22"))
c26 = F26.merge(t, on=["Gene", "Lineage"], how="left", suffixes=("_26", "_22")); testable = c26.q_value_22.notna(); k22 = set(zip(f.Gene, f.Lineage)); k26 = set(zip(F26.Gene, F26.Lineage))
R = dict(models_22Q1_raw=n_raw, models_unmatched_in_Model_csv=n_unmatched, non_cancerous_removed=n_nc, cancer_models_22Q1=len(df), models_shared_with_26Q1_cohort=int(df.DepMap_ID.isin(ids26).sum()),
         genes_22Q1=len(genes), eligible_lineages_22Q1=int((lc >= rp.MIN_N).sum()), n_tests_22Q1=len(t), n_q_lt_0_05_22Q1=int((t.q_value < rp.Q_MAX).sum()), n_candidates_22Q1=len(f), genes_22Q1_candidates=int(f.Gene.nunique()), lineages_22Q1_candidates=int(f.Lineage.nunique()),
         shared_tests=len(M), spearman_selectivity_all_shared=float(spearmanr(M.Selectivity_26, M.Selectivity_22)[0]),
         c26_testable_in_22=int(testable.sum()), c26_q_lt_0_05_negative_in_22=int(((c26.q_value_22 < rp.Q_MAX) & (c26.Selectivity_22 < 0)).sum()), c26_also_candidates_in_22=len(k26 & k22),
         c26_median_selectivity_26=float(c26.Selectivity_26.median()), c26_median_selectivity_22=float(c26.Selectivity_22.median()), spearman_selectivity_candidates=float(spearmanr(c26[testable].Selectivity_26, c26[testable].Selectivity_22)[0]),
         c22_also_candidates_in_26=len(k22 & k26), c22_total=len(k22))
missing = c26[~testable][["Gene", "Lineage"]].values.tolist(); R["c26_not_testable_in_22"] = missing
json.dump(R, open(O + "cross_release_summary.json", "w"), indent=1); t.to_csv(O + "all_tests_22Q1.csv.gz", index=False)
c26.to_csv(O + "cross_release_candidates.csv", index=False)
_c=c26[c26.q_value_22.notna() & (c26.q_value_22 < rp.Q_MAX) & (c26.Selectivity_22 < 0)]; _a=_c.Chronos_median_22 >= -1; _b=_c.Selectivity_22 >= -0.5
R.update(c26_signeg_fail_only_median_rule=int((_a & ~_b).sum()), c26_signeg_fail_only_selectivity=int((~_a & _b).sum()), c26_signeg_fail_both=int((_a & _b).sum()), median_chronos_c26_26Q1=float(c26.Chronos_median_26.median()), median_chronos_c26_22Q1=float(c26.Chronos_median_22.median()))
k26s = set(zip(T26.Gene[(T26.q_value < 0.05) & (T26.Chronos_median < -1) & (T26.Selectivity < -0.5)], T26.Lineage[(T26.q_value < 0.05) & (T26.Chronos_median < -1) & (T26.Selectivity < -0.5)]))
_o=f[[ (g, l) not in k26 for g, l in zip(f.Gene, f.Lineage)]].merge(T26, on=["Gene", "Lineage"], how="left", suffixes=("_22", "_26"))
R.update(c22_only_count=len(_o), c22_only_in26_sig_negative=int(((_o.q_value_26 < 0.05) & (_o.Selectivity_26 < 0)).sum()), shared_models_fraction=R["models_shared_with_26Q1_cohort"] / R["cancer_models_22Q1"])
_st = ["not testable" if pd.isna(q) else ("meets rule" if (q < rp.Q_MAX and m < rp.CHRONOS_MAX and s_ < rp.SEL_MAX) else ("significant, below rule" if (q < rp.Q_MAX and s_ < 0) else "not significant")) for q, m, s_ in zip(c26.q_value_22, c26.Chronos_median_22, c26.Selectivity_22)]
R["c26_categories"] = pd.Series(_st).value_counts().to_dict()   # status in 22Q1 of the 26Q1 candidates (Fig 11B)
json.dump(R, open(O + "cross_release_summary.json", "w"), indent=1); print(json.dumps(R, indent=1))
