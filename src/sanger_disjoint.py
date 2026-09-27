"""Separation of the DepMap 26Q1 data from the Sanger Project Score screens.

The 26Q1 CRISPRGeneEffect matrix is estimated from all screens that pass QC and are not flagged ExcludeFromCRISPRCombined in ScreenSequenceMap.csv;
these include the Sanger Project Score screens (library "KY").  For every cell line with such a KY screen, the 26Q1 gene effect is therefore partly or
entirely derived from the Sanger screen, and a comparison of 26Q1 with the Project Score data is not independent for these lines.

Here the complete selection procedure is repeated on the 26Q1 cancer cell lines WITHOUT any KY screen (discovery cohort), and the resulting candidates are
tested in the Sanger Project Score data.  The two data sets then share neither screens nor cell lines (every Sanger model has a KY screen and is therefore
not in the discovery cohort).  Remaining link: the 26Q1 combined data set is harmonised across libraries, so a small amount of shared processing cannot be excluded.
Outputs: sanger_disjoint_summary.json, sanger_disjoint_candidates.csv (discovery candidates and their Sanger results), all_tests_noKY.csv.gz,
sanger_disjoint_primary94.csv (the 94 primary candidates in the discovery cohort and in the Sanger data), sanger_disjoint_bands.csv."""
import os, json, numpy as np, pandas as pd
from scipy.stats import spearmanr
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); DS = os.environ.get("SANGER_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"

# ---- 26Q1 discovery cohort: cancer lines without any KY screen used in the combined data set
s = pd.read_csv(os.path.join(D, "ScreenSequenceMap.csv"))
used = s[(s.PassesQC == True) & (s.ExcludeFromCRISPRCombined == False) & s.Library.notna()]
ky_models = set(used.loc[used.Library == "KY", "ModelID"])
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
assert df.DepMap_ID.isin(set(used.ModelID)).all(), "cohort model without screen annotation"
n_cohort = len(df); dfd = df[~df.DepMap_ID.isin(ky_models)].reset_index(drop=True)
# also exclude discovery models derived from the same patient as any Project Score model (PatientID in Model.csv), so that the two sets share no patients
_pm = pd.read_csv(os.path.join(D, "Model.csv"))[["ModelID", "PatientID"]].set_index("ModelID").PatientID
_sg_ids = pd.read_csv(os.path.join(DS, "gene_effect.csv"), usecols=[0]).iloc[:, 0].astype(str).str.strip()
_sg_patients = set(_pm.reindex(_sg_ids).dropna()); _shared = dfd.DepMap_ID.map(_pm).isin(_sg_patients)
n_excluded_shared_patient = int(_shared.sum()); dfd = dfd[~_shared].reset_index(drop=True)
genes = [c for c in dfd.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
td = rp.test_all(dfd, genes, dfd.lineage.values); fd = rp.select(td).sort_values("Selectivity")

# ---- Sanger Project Score (Chronos), same procedure
ge = pd.read_csv(os.path.join(DS, "gene_effect.csv")); ge = ge.rename(columns={ge.columns[0]: "DepMap_ID"}); ge["DepMap_ID"] = ge.DepMap_ID.astype(str).str.strip()
meta = pd.read_csv(os.path.join(D, "Model.csv"))[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"]]; meta.columns = ["DepMap_ID", "lineage", "primary_disease"]
ds = ge.merge(meta, on="DepMap_ID", how="left"); ds = ds[ds.lineage.notna() & (ds.primary_disease != "Non-Cancerous")].reset_index(drop=True); ds["lineage"] = ds.lineage.str.strip().str.lower()
shared_lines = int(ds.DepMap_ID.isin(set(dfd.DepMap_ID)).sum()); assert shared_lines == 0, "discovery and replication cohorts share cell lines"
shared_patients = len(set(dfd.DepMap_ID.map(_pm).dropna()) & set(ds.DepMap_ID.map(_pm).dropna())); assert shared_patients == 0, "discovery and replication cohorts share patients"
gs = [c for c in ds.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
ts = rp.test_all(ds, gs, ds.lineage.values); fs = rp.select(ts)

# ---- replication of the discovery candidates
c = fd.merge(ts, on=["Gene", "Lineage"], how="left", suffixes=("_disc", "_S")); testable = c.q_value_S.notna()
c["replicated"] = np.where(testable, (c.q_value_S < rp.Q_MAX) & (c.Selectivity_S < 0), np.nan)
c["meets_rule_S"] = np.where(testable, (c.q_value_S < rp.Q_MAX) & (c.Chronos_median_S < rp.CHRONOS_MAX) & (c.Selectivity_S < rp.SEL_MAX), np.nan)
M = td.merge(ts, on=["Gene", "Lineage"], suffixes=("_disc", "_S"))
keys = pd.MultiIndex.from_frame(fd[["Gene", "Lineage"]])
Ms = M[(M.q_value_disc < rp.Q_MAX) & (M.Selectivity_disc < 0)].copy(); Ms["candidate"] = pd.MultiIndex.from_frame(Ms[["Gene", "Lineage"]]).isin(keys)
Ms["rep"] = (Ms.q_value_S < rp.Q_MAX) & (Ms.Selectivity_S < 0)
edges = [-np.inf, -0.7, -0.5, -0.4, -0.3, -0.2, 0]
Ms["band"] = pd.cut(Ms.Selectivity_disc, edges, right=False)
bands = Ms.groupby("band", observed=True).agg(n=("rep", "size"), replicated=("rep", "sum"), n_candidates=("candidate", "sum")).reset_index()
bands["frac"] = bands.replicated / bands.n; bands["band"] = bands.band.astype(str)
base = Ms[~Ms.candidate]

# ---- the 94 primary candidates in the discovery cohort and in the Sanger data
F26 = pd.read_csv(O + "final_targets.csv")
p = F26[["Gene", "Lineage", "Selectivity", "Chronos_median", "n_target", "q_value"]].merge(td[["Gene", "Lineage", "Selectivity", "Chronos_median", "n_target", "q_value"]], on=["Gene", "Lineage"], how="left", suffixes=("", "_disc"))
p = p.merge(ts[["Gene", "Lineage", "Selectivity", "Chronos_median", "n_target", "q_value"]].rename(columns=lambda x: x if x in ("Gene", "Lineage") else x + "_S"), on=["Gene", "Lineage"], how="left")
p["candidate_disc"] = pd.MultiIndex.from_frame(p[["Gene", "Lineage"]]).isin(keys)
tS = p.q_value_S.notna(); pc = p[p.candidate_disc & tS]

lcd = dfd.lineage.value_counts(); lcs = ds.lineage.value_counts()
R = dict(discovery_models_excluded_shared_patient=n_excluded_shared_patient, cohort_26Q1=n_cohort, cohort_models_with_KY_screen=int(df.DepMap_ID.isin(ky_models).sum()), discovery_models=len(dfd), discovery_eligible_lineages=int((lcd >= rp.MIN_N).sum()),
         discovery_tests=len(td), discovery_q_lt_0_05=int((td.q_value < rp.Q_MAX).sum()), discovery_candidates=len(fd), discovery_candidate_genes=int(fd.Gene.nunique()), discovery_candidate_lineages=int(fd.Lineage.nunique()),
         sanger_models=len(ds), sanger_models_shared_with_discovery=shared_lines, sanger_eligible_lineages=int((lcs >= rp.MIN_N).sum()), sanger_tests=len(ts), sanger_q_lt_0_05=int((ts.q_value < rp.Q_MAX).sum()),
         disc_testable_in_sanger=int(testable.sum()), disc_negative_direction=int((testable & (c.Selectivity_S < 0)).sum()), disc_nominal=int((testable & (c.p_value_S < 0.05) & (c.Selectivity_S < 0)).sum()),
         disc_replicated=int((c.replicated == 1).sum()), disc_meet_rule_in_sanger=int((c.meets_rule_S == 1).sum()),
         spearman_selectivity_all_shared=float(spearmanr(M.Selectivity_disc, M.Selectivity_S)[0]), shared_tests=len(M),
         spearman_selectivity_candidates=float(spearmanr(c[testable].Selectivity_disc, c[testable].Selectivity_S)[0]) if testable.sum() > 2 else None,
         median_sel_candidates_disc=float(c[testable].Selectivity_disc.median()), median_sel_candidates_S=float(c[testable].Selectivity_S.median()),
         baseline_noncandidate_sig_testable=int(len(base)), baseline_frac_replicated=float(base.rep.mean()),
         primary94_in_discovery_tested=int(p.q_value_disc.notna().sum()), primary94_candidates_in_discovery=int(p.candidate_disc.sum()),
         primary94_q05_neg_in_discovery=int(((p.q_value_disc < rp.Q_MAX) & (p.Selectivity_disc < 0)).sum()),
         primary94_disc_candidates_testable_in_sanger=len(pc), primary94_disc_candidates_replicated=int(((pc.q_value_S < rp.Q_MAX) & (pc.Selectivity_S < 0)).sum()),
         primary94_disc_candidates_meet_rule_in_sanger=int(((pc.q_value_S < rp.Q_MAX) & (pc.Chronos_median_S < rp.CHRONOS_MAX) & (pc.Selectivity_S < rp.SEL_MAX)).sum()),
         sanger_lineages=lcs[lcs >= rp.MIN_N].to_dict())
R["replicated_pairs"] = c[c.replicated == 1][["Gene", "Lineage"]].values.tolist(); R["not_replicated_testable"] = c[testable & (c.replicated == 0)][["Gene", "Lineage"]].values.tolist()
json.dump(R, open(O + "sanger_disjoint_summary.json", "w"), indent=1, default=str)
c.to_csv(O + "sanger_disjoint_candidates.csv", index=False); p.to_csv(O + "sanger_disjoint_primary94.csv", index=False); bands.to_csv(O + "sanger_disjoint_bands.csv", index=False)
td.to_csv(O + "all_tests_noKY.csv.gz", index=False)
print(json.dumps({k: v for k, v in R.items() if k not in ("replicated_pairs", "not_replicated_testable")}, indent=1, default=str)); print(bands.to_string())
