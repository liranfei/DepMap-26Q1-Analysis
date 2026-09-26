"""Descriptive statistics of the cancer cohort and of the test family (results/cohort_stats.json)."""
import os, json, numpy as np, pandas as pd
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"
T = pd.read_csv(O + "all_tests.csv"); F = pd.read_csv(O + "final_targets.csv"); s = T.Selectivity
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"]; genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
X = df[genes].apply(pd.to_numeric, errors="coerce"); miss = int(X.isna().sum().sum())
C = pd.read_csv(O + "codependency_corr.csv", index_col=0); sym = {g: g.split(" (")[0] for g in C.index}; C.index = [sym[g] for g in C.index]; C.columns = [sym[g] for g in C.columns]
TF = [g for g in ["SOX10", "HNF1B", "PAX8", "IRF4", "MYB", "CBFB", "FOXA1"] if g in C.index]; MET = [g for g in ["ADSL", "PAICS", "UMPS", "CTPS1", "DHFR", "NMNAT1", "NAMPT", "FPGS", "TYMS", "SDHB", "UROD"] if g in C.index]
cross = np.abs(C.loc[TF, MET].values).ravel(); within = np.abs(C.loc[MET, MET].values[np.triu_indices(len(MET), 1)])
ce = {g.split(" (")[0] for g in pd.read_csv(os.path.join(D, "CRISPRInferredCommonEssentials.csv")).iloc[:, 0]}
med = X.median().dropna().sort_values(); top = [g.split(" (")[0] for g in med.index[:len(set(F.Gene))]]; cand = {g.split(" (")[0] for g in F.Gene}
abs_rank = dict(n=len(top), frac_common_essential=float(np.mean([g in ce for g in top])), overlap_with_candidates=len(set(top) & cand), median_of_median_chronos=float(med.iloc[:len(top)].median()))
rec = F.assign(sym=F.Gene.str.split(" (", regex=False).str[0]).groupby("sym").Lineage.nunique().sort_values(ascending=False)
H = F[F.Lineage.isin(["lymphoid", "myeloid"])]
R = dict(absolute_rank_baseline=abs_rank, recurrent_genes={k: int(v) for k, v in rec[rec >= 3].items()}, genes_single_lineage=int((rec == 1).sum()), missing_values=miss, missing_fraction=miss / X.size, genes_with_missing=int(X.isna().any().sum()), genes_tested=int(T.Gene.nunique()), genes_total=len(genes),
         frac_selectivity_within_0_2=float(((s >= -0.2) & (s <= 0.2)).mean()), frac_selectivity_lt_minus0_5=float((s < -0.5).mean()), median_selectivity=float(s.median()),
         heme_pairs=len(H), heme_genes=int(H.Gene.nunique()), n_pairs=len(F), tf_genes=TF, metabolic_genes=MET, median_abs_r_tf_vs_metabolic=float(np.median(cross)), median_abs_r_within_metabolic=float(np.median(within)))
json.dump(R, open(O + "cohort_stats.json", "w"), indent=1); print(json.dumps(R, indent=1))
