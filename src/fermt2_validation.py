"""Pre-specified checks of FERMT2 (prespecified_plan_fermt2.md). Primary pair: CNS/brain; secondary: skin, kidney.
Environment: DEPMAP_DIR (26Q1 files), RESULTS_DIR."""
import os, json, hashlib
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
import statsmodels.formula.api as smf

DATA = os.environ.get("DEPMAP_DIR", "."); OUT = os.environ.get("RESULTS_DIR", "results")
PLAN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "prespecified_plan_fermt2.md")
GENE = "FERMT2"
hdr = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), nrows=0).columns
idc = hdr[0]
ge = pd.read_csv(os.path.join(DATA, "gene_effect.csv")).rename(columns={idc: "ModelID"}).set_index("ModelID")
M = pd.read_csv(os.path.join(DATA, "Model.csv"), usecols=["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "GrowthPattern"]).set_index("ModelID")
M = M.loc[M.index.intersection(ge.index)]
M = M[M.OncotreePrimaryDisease != "Non-Cancerous"]
ge = ge.loc[M.index]
col = [c for c in ge.columns if c.split(" ")[0] == GENE][0]
y = ge[col]; ok = y.notna()
lin = M.OncotreeLineage
mw = lambda a, b: float(mannwhitneyu(a, b, alternative="less", method="asymptotic")[1])
sel = lambda a, b: float(np.median(a) - np.median(b))
R = {"plan_sha256": hashlib.sha256(open(PLAN, "rb").read()).hexdigest()}

# mutations
mu = pd.read_csv(os.path.join(DATA, "OmicsSomaticMutations.csv"), usecols=["ModelID", "HugoSymbol", "Hotspot", "VepImpact", "IsDefaultEntryForModel"], low_memory=False)
mu = mu[mu.IsDefaultEntryForModel.astype(str).isin(["Yes", "True"])]
def mut(sym, hotspot):
    s = mu[mu.HugoSymbol == sym]
    return set(s[s.Hotspot == True].ModelID) if hotspot else set(s[s.VepImpact.isin(["HIGH", "MODERATE"])].ModelID)

def checks(target):
    t = (lin == target) & ok; o = (lin != target) & ok
    res = {"n_target": int(t.sum()), "selectivity_all": sel(y[t], y[o]), "p_all": mw(y[t], y[o])}
    # 1 culture format: other lineages, adherent only
    oa = o & (M.GrowthPattern == "Adherent")
    res["c1_adherent"] = {"n_other": int(oa.sum()), "selectivity": sel(y[t], y[oa]), "p": mw(y[t], y[oa])}
    res["c1_pass"] = bool(res["c1_adherent"]["selectivity"] <= res["selectivity_all"] / 2 and res["c1_adherent"]["p"] < 0.05)
    # 2 genotype
    g = pd.DataFrame({"y": y[ok], "L": (lin[ok] == target).astype(int)})
    for s_, hs in [("TP53", False), ("PTEN", False), ("NF1", False), ("IDH1", True), ("EGFR", True)]:
        g[s_] = g.index.isin(mut(s_, hs)).astype(int)
    m0 = smf.ols("y ~ L", g).fit(cov_type="HC3", use_t=True); m1 = smf.ols("y ~ L + TP53 + PTEN + NF1 + IDH1 + EGFR", g).fit(cov_type="HC3", use_t=True)
    res["c2_genotype"] = {"coef_unadj": float(m0.params.L), "coef_adj": float(m1.params.L), "ratio": float(m1.params.L / m0.params.L), "p_adj": float(m1.pvalues.L)}
    res["c2_pass"] = bool(res["c2_genotype"]["ratio"] >= 0.70 and res["c2_genotype"]["p_adj"] < 0.05)
    # 3 subtypes
    sub = {}
    for dz, idx in M[t].groupby("OncotreePrimaryDisease").groups.items():
        if len(idx) >= 5:
            a = y[idx]; sub[dz] = {"n": len(idx), "selectivity": sel(a, y[o]), "p": mw(a, y[o])}
    res["c3_subtypes"] = sub
    npass = sum(v["selectivity"] < 0 and v["p"] < 0.05 for v in sub.values())
    res["c3_pass"] = bool(len(sub) > 0 and npass >= 2 / 3 * len(sub))
    # 4 carcinoma comparison
    carc = ["Bowel", "Pancreas", "Esophagus/Stomach", "Breast", "Prostate", "Lung", "Head and Neck", "Bladder/Urinary Tract", "Liver",
            "Biliary Tract", "Ampulla of Vater", "Uterus", "Ovary/Fallopian Tube", "Cervix", "Thyroid"]
    oc = o & lin.isin([c for c in carc if c != target])
    res["c4_carcinoma"] = {"n_other": int(oc.sum()), "selectivity": sel(y[t], y[oc]), "p": mw(y[t], y[oc])}
    res["c4_pass"] = bool(res["c4_carcinoma"]["selectivity"] < 0 and res["c4_carcinoma"]["p"] < 0.05)
    return res

for tgt in ["CNS/Brain", "Skin", "Kidney"]:
    R[tgt] = checks(tgt)

# 5 paralogue buffering
ex = pd.read_csv(os.path.join(DATA, "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv"), low_memory=False)
if "IsDefaultEntryForModel" in ex.columns:
    ex = ex[ex.IsDefaultEntryForModel.astype(str).isin(["Yes", "True"])]
ex = ex.set_index("ModelID")
pcol = {p: [c for c in ex.columns if str(c).split(" ")[0] == p][0] for p in ["FERMT1", "FERMT3"]}
e = pd.DataFrame({"y": y[ok]}).join(ex[[pcol["FERMT1"], pcol["FERMT3"]]].rename(columns={pcol["FERMT1"]: "FERMT1", pcol["FERMT3"]: "FERMT3"}), how="inner").join(lin)
P = {}
for p in ["FERMT1", "FERMT3"]:
    r1, p1 = spearmanr(e[p], e.y)
    c = e.copy(); c[p] = c[p] - c.groupby("OncotreeLineage")[p].transform("mean"); c["y"] = c.y - c.groupby("OncotreeLineage").y.transform("mean")
    r2, p2 = spearmanr(c[p], c.y)
    P[p] = {"n": int(len(e)), "rho": float(r1), "p": float(p1), "rho_within_lineage": float(r2), "p_within_lineage": float(p2)}
R["c5_paralogue"] = P
R["c5_pass_FERMT1"] = bool(P["FERMT1"]["rho"] > 0 and P["FERMT1"]["p"] < 0.05 and P["FERMT1"]["rho_within_lineage"] > 0 and P["FERMT1"]["p_within_lineage"] < 0.05)
# carcinoma lines' FERMT1 vs CNS
R["c5_FERMT1_median"] = {k: float(e[e.OncotreeLineage == k].FERMT1.median()) for k in ["CNS/Brain", "Skin", "Kidney", "Bowel", "Lung", "Breast"]}

# 6 co-dependency (lineage-centred), genes with >= 300 lines
X = ge.loc[:, ge.notna().sum() >= 300]
Xc = X - X.groupby(lin).transform("mean")
yc = Xc[col]
cors = Xc.drop(columns=[col]).corrwith(yc).dropna().sort_values(ascending=False)
partners = ["ITGB1", "ITGAV", "ITGB5", "ITGA3", "ITGA5", "ILK", "LIMS1", "PARVA", "TLN1", "PTK2", "PXN", "VCL", "FERMT1"]
top20 = cors.head(20)
R["c6_top20"] = {k.split(" ")[0]: round(float(v), 3) for k, v in top20.items()}
R["c6_partner_ranks"] = {p: int(np.where(cors.index.str.split(" ").str[0] == p)[0][0] + 1) for p in partners if (cors.index.str.split(" ").str[0] == p).any()}
R["c6_n_genes"] = int(len(cors) + 1)
R["c6_pass"] = bool(any(k in partners for k in R["c6_top20"]))

prim = R["CNS/Brain"]
R["decision_use_as_example"] = bool(prim["c1_pass"] and prim["c2_pass"] and prim["c4_pass"])
json.dump(R, open(os.path.join(OUT, "fermt2_validation_summary.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
