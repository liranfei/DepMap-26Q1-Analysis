"""Independent-screen check: the same procedure applied to the Sanger Project Score CRISPR screens (Chronos-processed, figshare 10.6084/m9.figshare.14461980).
Project Score screens were generated at a different institute, with a different library and protocol; cell lines overlap with the Broad screens, so this is a check of
reproducibility across independently generated screens, not a test in new cell lines.  Lineage / non-cancerous flag are taken from the 26Q1 Model.csv."""
import os, json, numpy as np, pandas as pd
from scipy.stats import spearmanr
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); DS = os.environ.get("SANGER_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"
ge = pd.read_csv(os.path.join(DS, "gene_effect.csv")); ge = ge.rename(columns={ge.columns[0]: "DepMap_ID"}); ge["DepMap_ID"] = ge.DepMap_ID.astype(str).str.strip()
meta = pd.read_csv(os.path.join(D, "Model.csv"))[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"]]; meta.columns = ["DepMap_ID", "lineage", "primary_disease"]
n_raw = len(ge); df = ge.merge(meta, on="DepMap_ID", how="left"); n_unmatched = int(df.lineage.isna().sum()); df = df[df.lineage.notna()]
n_nc = int((df.primary_disease == "Non-Cancerous").sum()); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True); df["lineage"] = df.lineage.str.strip().str.lower()
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
t = rp.test_all(df, genes, df.lineage.values); f = rp.select(t)
T26 = pd.read_csv(O + "all_tests.csv"); F26 = pd.read_csv(O + "final_targets.csv"); ids26 = set(rp.load(D).query("primary_disease != 'Non-Cancerous'").DepMap_ID)
lc = df.lineage.value_counts(); M = T26.merge(t, on=["Gene", "Lineage"], suffixes=("_26", "_S"))
c = F26.merge(t, on=["Gene", "Lineage"], how="left", suffixes=("_26", "_S")); testable = c.q_value_S.notna(); k26 = set(zip(F26.Gene, F26.Lineage)); kS = set(zip(f.Gene, f.Lineage))
signeg = (c.q_value_S < rp.Q_MAX) & (c.Selectivity_S < 0); direc = testable & (c.Selectivity_S < 0); nom = testable & (c.p_value_S < 0.05) & (c.Selectivity_S < 0)
# baseline: the same replication measure for gene-lineage pairs with q<0.05 in 26Q1 and selectivity in [-0.5,-0.3) that are NOT candidates, to show that candidates replicate more often than merely-significant pairs
Ms = M[(M.q_value_26 < rp.Q_MAX) & (M.Selectivity_26 < 0)]; cand_keys = pd.MultiIndex.from_frame(F26[["Gene", "Lineage"]])
Ms = Ms[~pd.MultiIndex.from_frame(Ms[["Gene", "Lineage"]]).isin(cand_keys)]
R = dict(models_sanger_raw=n_raw, models_unmatched_in_Model_csv=n_unmatched, non_cancerous_removed=n_nc, cancer_models_sanger=len(df), models_shared_with_26Q1_cohort=int(df.DepMap_ID.isin(ids26).sum()),
         genes_sanger=len(genes), eligible_lineages_sanger=int((lc >= rp.MIN_N).sum()), lineages_sanger=lc[lc >= rp.MIN_N].to_dict(), n_tests_sanger=len(t), n_q_lt_0_05_sanger=int((t.q_value < rp.Q_MAX).sum()), n_candidates_sanger=len(f),
         shared_tests=len(M), spearman_selectivity_all_shared=float(spearmanr(M.Selectivity_26, M.Selectivity_S)[0]),
         c26_testable=int(testable.sum()), c26_negative_direction=int(direc.sum()), c26_nominal_p_lt_0_05_negative=int(nom.sum()), c26_q_lt_0_05_negative=int(signeg.sum()), c26_also_candidates=len(k26 & kS),
         c26_median_selectivity_26=float(c[testable].Selectivity_26.median()), c26_median_selectivity_S=float(c[testable].Selectivity_S.median()), spearman_selectivity_candidates=float(spearmanr(c[testable].Selectivity_26, c[testable].Selectivity_S)[0]),
         baseline_noncandidate_sig_pairs=int(len(Ms)), baseline_n_testable_in_sanger=int(Ms.q_value_S.notna().sum()), baseline_frac_q_lt_0_05_negative=float(((Ms.q_value_S < rp.Q_MAX) & (Ms.Selectivity_S < 0)).sum() / max(1, Ms.q_value_S.notna().sum())),
         baseline_frac_nominal_negative=float(((Ms.p_value_S < 0.05) & (Ms.Selectivity_S < 0)).sum() / max(1, Ms.p_value_S.notna().sum())),
         cand_frac_q_lt_0_05_negative_of_testable=float(signeg.sum() / max(1, testable.sum())), cand_frac_nominal_negative_of_testable=float(nom.sum() / max(1, testable.sum())))
R["c26_not_testable"] = c[~testable][["Gene", "Lineage"]].values.tolist(); R["c26_replicating_rule"] = sorted(map(list, k26 & kS))
json.dump(R, open(O + "independent_sanger_summary.json", "w"), indent=1, default=str); t.to_csv(O + "all_tests_sanger.csv.gz", index=False); c.to_csv(O + "independent_sanger_candidates.csv", index=False); print(json.dumps({k: v for k, v in R.items() if k not in ('c26_not_testable', 'c26_replicating_rule')}, indent=1, default=str))
