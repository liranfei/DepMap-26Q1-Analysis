"""S13 Table: prioritised pairs in the next-generation (NextGen) model screens of Neiswender et al. (Nature 2026).
Uses next_gen_dependency_lineage_tests.csv and model_metadata.csv from figshare 10.6084/m9.figshare.29472362
(results/nextgen_figshare/), as published, and the 26Q1 2D data for the matched contrast.
Environment: DEPMAP_DIR (26Q1 files), RESULTS_DIR."""
import os, json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

DATA = os.environ.get("DEPMAP_DIR", "."); OUT = os.environ.get("RESULTS_DIR", "results")
NG = os.path.join(OUT, "nextgen_figshare")
ng = pd.read_csv(os.path.join(NG, "next_gen_dependency_lineage_tests.csv"))
mm = pd.read_csv(os.path.join(NG, "model_metadata.csv"))
fmt = mm[mm.CRISPRScreenType.isin(["2DO", "3DO", "2DN", "3DN"])].groupby(["OncotreeLineage", "CRISPRScreenType"]).size().unstack(fill_value=0)
LMAP = {"bowel": ("Colorectal", "Bowel"), "pancreas": ("Pancreas", "Pancreas"), "esophagus/stomach": ("Esophagus/Stomach", "Esophagus/Stomach"),
        "breast": ("Breast", "Breast"), "prostate": ("Prostate", "Prostate"), "cns/brain": ("CNS/Brain", "CNS/Brain")}
F = pd.read_csv(os.path.join(OUT, "final_targets.csv"))
F = F[F.Lineage.isin(LMAP)]
# 2D matched contrast: target lineage vs the lineages represented among NextGen models
ng_lins = [l for l in fmt.index]
ge_hdr = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), nrows=0).columns
M = pd.read_csv(os.path.join(DATA, "Model.csv"), usecols=["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"])
genes = sorted(set(F.Gene))
ge = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), usecols=[ge_hdr[0]] + genes).rename(columns={ge_hdr[0]: "ModelID"})
D = ge.merge(M, on="ModelID"); D = D[D.OncotreePrimaryDisease != "Non-Cancerous"]
rows = []
for _, r in F.iterrows():
    nl, ml = LMAP[r.Lineage]
    x = ng[(ng.Dependency == r.Gene) & (ng.OncotreeLineage == nl)]
    if x.empty:
        continue
    x = x.iloc[0]
    vi, ni, vo, no = x["Ingroup Var"], x["Ingroup n"], x["Outgroup Var"], x["Outgroup n"]
    se = np.sqrt(vi / ni + vo / no); t = x["Mean Difference"] / se
    df = (vi / ni + vo / no) ** 2 / ((vi / ni) ** 2 / (ni - 1) + (vo / no) ** 2 / (no - 1))
    y = D[r.Gene]; tgt = y[(D.OncotreeLineage == ml) & y.notna()]
    oth = y[(D.OncotreeLineage != ml) & D.OncotreeLineage.isin([l for l in ng_lins if l != ml]) & y.notna()]
    rows.append(dict(Gene=r.Gene, Lineage=r.Lineage, NextGen_lineage=nl,
                     NextGen_target_models=int(ni), NextGen_target_3D=int(fmt.loc[ml, ["3DO", "3DN"]].sum()), NextGen_target_2D_coated=int(fmt.loc[ml, ["2DO", "2DN"]].sum()),
                     NextGen_other_models=int(no), NextGen_target_mean=x["Ingroup Mean"], NextGen_other_mean=x["Outgroup Mean"],
                     NextGen_mean_difference=x["Mean Difference"], NextGen_p_one_sided=float(stats.t.sf(-t, df)),
                     TwoD_selectivity_primary=r.Selectivity, TwoD_mean_difference_matched=float(tgt.mean() - oth.mean()), TwoD_n_target=int(len(tgt)), TwoD_n_other=int(len(oth))))
T = pd.DataFrame(rows)
T["NextGen_q_BH"] = multipletests(T.NextGen_p_one_sided, method="fdr_bh")[1]
T["Reproduced_q_lt_0.05"] = T.NextGen_q_BH < 0.05
T = T.sort_values("NextGen_p_one_sided")
T.to_csv(os.path.join(OUT, "nextgen_check_pairs.csv"), index=False)
# CFLAR in pancreas: 2D lineage coefficient with and without KRAS hotspot status (OLS, HC3)
import statsmodels.formula.api as smf
mu = pd.read_csv(os.path.join(DATA, "OmicsSomaticMutations.csv"), usecols=["ModelID", "HugoSymbol", "Hotspot", "IsDefaultEntryForModel"], low_memory=False)
mu = mu[mu.IsDefaultEntryForModel.astype(str).isin(["Yes", "True"])]
kr = set(mu[(mu.HugoSymbol == "KRAS") & (mu.Hotspot == True)].ModelID)
cg = [g for g in genes if g.startswith("CFLAR ")][0]
g = D[["ModelID", "OncotreeLineage", cg]].dropna().rename(columns={cg: "y"})
g["L"] = (g.OncotreeLineage == "Pancreas").astype(int); g["KRAS"] = g.ModelID.isin(kr).astype(int)
m0 = smf.ols("y ~ L", g).fit(cov_type="HC3", use_t=True); m1 = smf.ols("y ~ L + KRAS", g).fit(cov_type="HC3", use_t=True)
cflar = dict(coef_unadj=float(m0.params.L), coef_adj_KRAS=float(m1.params.L), ratio=float(m1.params.L / m0.params.L), p_adj=float(m1.pvalues.L), KRAS_coef=float(m1.params.KRAS), KRAS_p=float(m1.pvalues.KRAS))
summ = dict(n_pairs=int(len(T)), n_reproduced=int(T["Reproduced_q_lt_0.05"].sum()), reproduced=T[T["Reproduced_q_lt_0.05"]][["Gene", "Lineage"]].values.tolist(),
            screen_formats=fmt.to_dict("index"), cflar_pancreas_kras_adjustment=cflar)
json.dump(summ, open(os.path.join(OUT, "nextgen_check_summary.json"), "w"), indent=1)
defs = pd.DataFrame([
    ("Source", "Neiswender JV et al., Nature 2026; figshare 10.6084/m9.figshare.29472362, file next_gen_dependency_lineage_tests.csv (means, variances and numbers of NextGen models per lineage vs the remaining NextGen models), used as published."),
    ("Test", "One-sided Welch t-test (target lineage lower) computed from the published summaries; Benjamini-Hochberg over the 17 pairs."),
    ("Culture format", "NextGen screens are 3D (3DO dome organoids, 3DN spheroids) or on coated plates (2DO, 2DN); counts per lineage are given. CNS/brain NextGen models are diffuse gliomas (25 2DN, 14 3DN)."),
    ("2D matched contrast", "DepMap 26Q1: difference in means between the target lineage and the cancer lines of the other lineages represented among NextGen models."),
    ("Timing", "Added after the revised analysis was complete (post hoc)."),
], columns=["Item", "Definition"])
with pd.ExcelWriter(os.path.join(OUT, "S13_Table.xlsx")) as w:
    T.to_excel(w, sheet_name="Pairs", index=False)
    defs.to_excel(w, sheet_name="Definitions", index=False)
print(json.dumps(summ, indent=1)); print(T[["Gene", "Lineage", "NextGen_mean_difference", "NextGen_q_BH", "TwoD_mean_difference_matched"]].round(3).to_string(index=False))
