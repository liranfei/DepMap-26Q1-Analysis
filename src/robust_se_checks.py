"""Heteroscedasticity-robust (HC3) standard errors for every inference that rests on ordinary least squares, because gene-effect variances differ
between lineages (for example, between blood and solid lineages) and the default OLS standard errors assume equal variances.
Re-computes, with HC3: (1) batch/library-adjusted lineage coefficients (95% CI below 0), (2) growth-pattern-adjusted lineage coefficients (BH over the
94 pairs; blood pairs also without the other blood lineage, BH over 42), (3) lineage-adjusted genotype differences (KRAS, CTNNB1, MDM2/TP53, APC, BRAF),
(4) genotype-adjusted lineage coefficients and lineage-by-genotype interactions. Output: robust_se_checks.json."""
import os, json, numpy as np, pandas as pd, statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; R = {}
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True); lin = df.lineage.values
F = pd.read_csv(O + "final_targets.csv"); M = pd.read_csv(os.path.join(D, "Model.csv")).set_index("ModelID")
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]; col = lambda s: [c for c in genes if c.startswith(s + " (")][0]

# ---- 1. batch / library (same covariates as revision4_checks.py)
s = pd.read_csv(os.path.join(D, "ScreenSequenceMap.csv")); s = s[(s.ScreenType == "2DS") & (s.PassesQC == True) & (s.DrugTreated == False) & (s.IsEngineered == False) & (s.ExcludeFromCRISPRCombined == False)]
mode = lambda x: x.value_counts().index[0]
dfb = pd.DataFrame({"batch": df.DepMap_ID.map(s.groupby("ModelID").pDNABatch.agg(mode)).fillna("unknown"), "library": df.DepMap_ID.map(s.groupby("ModelID").Library.agg(mode)).fillna("unknown")})
hi_ols, hi_hc3 = [], []
for r in F.itertuples():
    d = pd.DataFrame({"y": df[r.Gene].astype(float), "tgt": (lin == r.Lineage).astype(int), "batch": dfb.batch, "library": dfb.library}).dropna()
    f = smf.ols("y~tgt+C(batch)+C(library)", d).fit(); hi_ols.append(f.conf_int().loc["tgt", 1]); hi_hc3.append(f.get_robustcov_results("HC3").conf_int()[list(f.params.index).index("tgt"), 1])
R["batch_ci_below_0"] = {"ols": int((np.array(hi_ols) < 0).sum()), "hc3": int((np.array(hi_hc3) < 0).sum()), "n": len(F)}

# ---- 2. growth pattern
gp = pd.get_dummies(df.DepMap_ID.map(M.GrowthPattern).fillna("Unknown")).astype(float); gp = gp[[c for c in gp.columns if c != "Adherent"]]
rows = []
for _, r in F.iterrows():
    y = pd.to_numeric(df[r.Gene], errors="coerce"); ok = y.notna().values; t = (lin == r.Lineage).astype(float)
    f = sm.OLS(y[ok].values, sm.add_constant(np.column_stack([t[ok], gp.values[ok]]))).fit(); row = dict(Gene=r.Gene, Lineage=r.Lineage, p_ols=f.pvalues[1], p_hc3=f.get_robustcov_results("HC3").pvalues[1], coef=f.params[1])
    if r.Lineage in ("lymphoid", "myeloid"):
        k = ok & (lin != {"lymphoid": "myeloid", "myeloid": "lymphoid"}[r.Lineage]); Xg = gp.values[k]; Xg = Xg[:, Xg.std(0) > 0]
        f2 = sm.OLS(y[k].values, sm.add_constant(np.column_stack([t[k], Xg]))).fit(); row.update(p_ols_ex=f2.pvalues[1], p_hc3_ex=f2.get_robustcov_results("HC3").pvalues[1])
    rows.append(row)
G = pd.DataFrame(rows); heme = G.Lineage.isin(["lymphoid", "myeloid"])
for c in ("p_ols", "p_hc3"): G["q" + c[1:]] = multipletests(G[c], method="fdr_bh")[1]
H = G[heme].copy()
for c in ("p_ols_ex", "p_hc3_ex"): H["q" + c[1:]] = multipletests(H[c], method="fdr_bh")[1]
sig = lambda q, g: int(((q < 0.05) & (g.coef < 0)).sum())
R["growth"] = {"all94_ols": sig(G.q_ols, G), "all94_hc3": sig(G.q_hc3, G), "heme_ols": sig(G[heme].q_ols, G[heme]), "heme_hc3": sig(G[heme].q_hc3, G[heme]),
               "nonheme_ols": sig(G[~heme].q_ols, G[~heme]), "nonheme_hc3": sig(G[~heme].q_hc3, G[~heme]),
               "heme_excl_other_blood_ols": int((H.q_ols_ex < 0.05).sum()), "heme_excl_other_blood_hc3": int((H.q_hc3_ex < 0.05).sum())}

# ---- 3. genotype differences adjusted for lineage, and 4. lineage coefficients adjusted for genotype
mu = pd.read_csv(os.path.join(D, "OmicsSomaticMutations.csv"), usecols=["ModelID", "HugoSymbol", "VepImpact", "ProteinChange", "Hotspot", "IsDefaultEntryForModel"], low_memory=False)
mu = mu[(mu.IsDefaultEntryForModel == "Yes") & mu.HugoSymbol.isin(["APC", "BRAF", "CTNNB1", "KRAS", "TP53"])]
pa = lambda g: set(mu[(mu.HugoSymbol == g) & mu.VepImpact.isin(["HIGH", "MODERATE"])].ModelID); hot = lambda g: set(mu[(mu.HugoSymbol == g) & (mu.Hotspot == True)].ModelID)
apc = set(mu[(mu.HugoSymbol == "APC") & (mu.VepImpact == "HIGH")].ModelID); v600 = set(mu[(mu.HugoSymbol == "BRAF") & mu.ProteinChange.astype(str).str.startswith("p.V600")].ModelID)
def geno(gene, carriers, excluded):
    d = pd.DataFrame({"y": df[col(gene)].astype(float), "g": df.DepMap_ID.isin(carriers).astype(int), "lin": lin, "id": df.DepMap_ID}).dropna(subset=["y"])
    d = d[(d.g == 1) | ~d.id.isin(excluded)]; f = smf.ols("y ~ g + C(lin)", d).fit()
    return {"coef": float(f.params["g"]), "p_ols": float(f.pvalues["g"]), "p_hc3": float(f.get_robustcov_results("HC3").pvalues[list(f.params.index).index("g")])}
R["genotype_lineage_adjusted"] = {"KRAS_hotspot": geno("KRAS", hot("KRAS"), pa("KRAS") - hot("KRAS")), "CTNNB1_hotspot": geno("CTNNB1", hot("CTNNB1"), pa("CTNNB1") - hot("CTNNB1")),
    "MDM2_by_TP53": geno("MDM2", set(df.DepMap_ID) - pa("TP53"), set()), "CTNNB1_by_APC": geno("CTNNB1", apc, (pa("APC") - apc) | hot("CTNNB1")), "BRAF_V600": geno("BRAF", v600, pa("BRAF") - v600)}
def lg(gene, lineage, carriers, excluded):
    d = pd.DataFrame({"y": df[col(gene)].astype(float), "t": (lin == lineage).astype(int), "g": df.DepMap_ID.isin(carriers).astype(int), "id": df.DepMap_ID}).dropna(subset=["y"])
    d = d[(d.g == 1) | ~d.id.isin(excluded)]; f1 = smf.ols("y ~ t + g", d).fit(); f2 = smf.ols("y ~ t * g", d).fit()
    i1 = list(f1.params.index); i2 = list(f2.params.index); r1 = f1.get_robustcov_results("HC3"); r2 = f2.get_robustcov_results("HC3")
    return {"lineage_p_ols": float(f1.pvalues["t"]), "lineage_p_hc3": float(r1.pvalues[i1.index("t")]), "noncarrier_lineage_p_ols": float(f2.pvalues["t"]), "noncarrier_lineage_p_hc3": float(r2.pvalues[i2.index("t")]),
            "interaction_p_ols": float(f2.pvalues["t:g"]), "interaction_p_hc3": float(r2.pvalues[i2.index("t:g")])}
R["lineage_given_genotype"] = {"KRAS_pancreas": lg("KRAS", "pancreas", hot("KRAS"), pa("KRAS") - hot("KRAS")),
                               "CTNNB1_bowel": lg("CTNNB1", "bowel", apc | hot("CTNNB1"), (pa("APC") - apc) | (pa("CTNNB1") - hot("CTNNB1"))), "BRAF_skin": lg("BRAF", "skin", v600, pa("BRAF") - v600)}
json.dump(R, open(O + "robust_se_checks.json", "w"), indent=1, default=float); print(json.dumps(R, indent=1, default=float))
