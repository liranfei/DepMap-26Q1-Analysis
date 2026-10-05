"""Validation of a single candidate, SNAP23 in the bowel lineage (DepMap 26Q1).
Environment: DEPMAP_DIR (gene_effect.csv, Model.csv, OmicsSomaticMutations.csv,
OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv), SANGER_DIR, RNAI_DIR, RESULTS_DIR.
Original header: Dry-lab validation of SNAP23 dependency in the bowel lineage (DepMap 26Q1).
Checks and pass criteria are fixed in 00_PREREGISTERED_PLAN.md, written before this script was run.
Run: python3 snap23_validation.py     (expects the DepMap files in DATA)
"""
import json, os, re, gzip
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

D = os.environ.get("DEPMAP_DIR", "."); DS = os.environ.get("SANGER_DIR", "")
RD = os.environ.get("RNAI_DIR", ""); OUT = os.environ.get("RESULTS_DIR", "results")
DATA = D
SANGER = os.path.join(DS, "gene_effect.csv") if DS else ""
GENE, LIN = "SNAP23", "bowel"
R = {}

def mw(a, b):
    return float(mannwhitneyu(a, b, alternative="less", method="asymptotic")[1])

def sel(a, b):
    return float(np.median(a) - np.median(b))

# ---------- load 26Q1 ----------
hdr = pd.read_csv(DATA + "/gene_effect.csv", nrows=0).columns
need = [hdr[0]] + [c for c in hdr if c.split(" ")[0] in
        (GENE, "SNAP25", "SNAP29", "STX4", "VAMP3", "VAMP8", "RPA3")]
ge = pd.read_csv(DATA + "/gene_effect.csv", usecols=need)
ge = ge.rename(columns={ge.columns[0]: "DepMap_ID"})
M = pd.read_csv(DATA + "/Model.csv")
meta = M[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "GrowthPattern"]]
meta.columns = ["DepMap_ID", "lineage", "primary_disease", "subtype", "growth"]
meta = meta.copy(); meta["lineage"] = meta.lineage.str.strip().str.lower()
df = ge.merge(meta, on="DepMap_ID", how="left", validate="one_to_one")
assert df.lineage.notna().all(), "merge failure"
df["lineage"] = df.lineage.str.strip().str.lower()
df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
col = [c for c in df.columns if c.split(" ")[0] == GENE][0]
df["y"] = pd.to_numeric(df[col], errors="coerce")
d = df[df.y.notna()].copy()
tgt, oth = d[d.lineage == LIN], d[d.lineage != LIN]
R["cohort"] = {"n_cancer_lines": int(len(df)), "n_with_SNAP23_data": int(len(d)),
               "n_bowel": int(len(tgt)), "n_other": int(len(oth))}

# ---------- V1 primary ----------
a, b = tgt.y.values, oth.y.values
rng = np.random.default_rng(0)
boot = [sel(rng.choice(a, len(a), replace=True), rng.choice(b, len(b), replace=True)) for _ in range(1000)]
lo, hi = np.percentile(boot, [2.5, 97.5])
v1 = {"median_bowel": float(np.median(a)), "median_other": float(np.median(b)),
      "selectivity": sel(a, b), "p_one_sided": mw(a, b),
      "bootstrap_CI": [float(lo), float(hi)],
      "published_q_from_S2_Table": 4.971552796994737e-08,
      "published_bootstrap_CI_S2_Table": [-0.9467674591839504, -0.6161600962252058]}
v1["pass"] = bool(v1["median_bowel"] < -1 and v1["selectivity"] < -0.5 and hi < -0.5)
R["V1_primary"] = v1

# ---------- genotype ----------
mu = pd.read_csv(DATA + "/OmicsSomaticMutations.csv",
                 usecols=["ModelID", "HugoSymbol", "VepImpact", "Hotspot", "IsDefaultEntryForModel"],
                 low_memory=False)
mu = mu[(mu.IsDefaultEntryForModel == "Yes") &
        mu.HugoSymbol.isin(["KRAS", "APC", "TP53", "BRAF", "PIK3CA"])]
def hot(sym):  return set(mu[(mu.HugoSymbol == sym) & (mu.Hotspot == True)].ModelID)
def alt(sym):  return set(mu[(mu.HugoSymbol == sym) & mu.VepImpact.isin(["HIGH", "MODERATE"])].ModelID)
g = d.copy()
g["bowel"] = (g.lineage == LIN).astype(int)
for sym, f in [("KRAS", hot), ("BRAF", hot), ("PIK3CA", hot), ("APC", alt), ("TP53", alt)]:
    g[sym] = g.DepMap_ID.isin(f(sym)).astype(int)
m0 = smf.ols("y ~ bowel", g).fit(cov_type="HC3", use_t=True)
m1 = smf.ols("y ~ bowel + KRAS + APC + TP53 + BRAF + PIK3CA", g).fit(cov_type="HC3", use_t=True)
v2 = {"coef_unadjusted": float(m0.params.bowel), "p_unadjusted": float(m0.pvalues.bowel),
      "coef_adjusted": float(m1.params.bowel), "p_adjusted": float(m1.pvalues.bowel),
      "ratio_adj_over_unadj": float(m1.params.bowel / m0.params.bowel),
      "genotype_coefficients": {k: round(float(v), 3) for k, v in m1.params.items() if k not in ("Intercept", "bowel")},
      "n_bowel_mutated": {s: int(g[g.bowel == 1][s].sum()) for s in ["KRAS", "APC", "TP53", "BRAF", "PIK3CA"]}}
v2["pass"] = bool(v2["ratio_adj_over_unadj"] >= 0.70 and v2["p_adjusted"] < 0.05)
R["V2_genotype"] = v2

# ---------- V3 culture format ----------
adh = oth[oth.growth == "Adherent"]
v3 = {"n_adherent_other": int(len(adh)), "median_adherent_other": float(adh.y.median()),
      "selectivity_vs_adherent": sel(a, adh.y.values), "p_one_sided": mw(a, adh.y.values)}
v3["pass"] = bool(v3["selectivity_vs_adherent"] < -0.5 and v3["p_one_sided"] < 0.05)
R["V3_culture_format"] = v3

# ---------- V4 subtypes ----------
sub = {}
for s, grp in tgt.groupby("subtype"):
    if len(grp) >= 5:
        sub[s] = {"n": int(len(grp)), "median": float(grp.y.median()),
                  "selectivity": sel(grp.y.values, b), "p_one_sided": mw(grp.y.values, b)}
v4 = {"subtypes": sub, "all_subtypes_below_-0.5": bool(sub) and all(v["selectivity"] < -0.5 for v in sub.values())}
v4["pass"] = v4["all_subtypes_below_-0.5"]
R["V4_subtype"] = v4

# ---------- expression (V5, V8) ----------
ehdr = pd.read_csv(DATA + "/OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv", nrows=0).columns
ecols = ["ModelID", "IsDefaultEntryForModel"] + [c for c in ehdr if c.split(" ")[0] in (GENE, "SNAP25", "SNAP29", "STX4", "VAMP3", "VAMP8")]
ex = pd.read_csv(DATA + "/OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv", usecols=ecols)
ex = ex[ex.IsDefaultEntryForModel == "Yes"].drop(columns=["IsDefaultEntryForModel"])
ex = ex.rename(columns={"ModelID": "DepMap_ID"})
ex.columns = [c.split(" ")[0] for c in ex.columns]
e = d[["DepMap_ID", "y", "lineage"]].merge(ex, on="DepMap_ID", how="inner")
R["V8_expression"] = {"n_lines_with_expression": int(len(e)),
                      "median_log2TPM1_bowel": float(e[e.lineage == LIN][GENE].median()),
                      "median_log2TPM1_other": float(e[e.lineage != LIN][GENE].median())}
R["V8_expression"]["pass"] = bool(R["V8_expression"]["median_log2TPM1_bowel"] > 1)
v5 = {}
for p in ["SNAP25", "SNAP29"]:
    ok = e[[p, "y"]].dropna()
    rho, pv = spearmanr(ok[p], ok.y)
    v5[p] = {"n": int(len(ok)), "spearman_rho": float(rho), "p_two_sided": float(pv),
             "direction_as_predicted(positive)": bool(rho > 0)}
# post hoc: same correlations after removing lineage means from both variables (controls for lineage composition)
el = e.copy()
for cc in [GENE, "SNAP25", "SNAP29", "y"]:
    el[cc] = el[cc] - el.groupby("lineage")[cc].transform("mean")
for p in ["SNAP25", "SNAP29"]:
    okl = el[[p, "y"]].dropna(); rl, pl = spearmanr(okl[p], okl.y)
    v5[p + "_lineage_centred(post_hoc)"] = {"n": int(len(okl)), "spearman_rho": float(rl), "p_two_sided": float(pl)}
ok = e[[GENE, "y"]].dropna(); rho, pv = spearmanr(ok[GENE], ok.y)
v5["SNAP23_own_expression(post_hoc)"] = {"n": int(len(ok)), "spearman_rho": float(rho), "p_two_sided": float(pv)}
v5["pass"] = any(v["spearman_rho"] > 0 and v["p_two_sided"] < 0.05 for k, v in v5.items() if k in ("SNAP25", "SNAP29"))
R["V5_paralogue_buffering"] = v5

# ---------- V6 co-dependency (all genes) ----------
full = pd.read_csv(DATA + "/gene_effect.csv", low_memory=False)
full = full.rename(columns={full.columns[0]: "DepMap_ID"})
full = full[full.DepMap_ID.isin(d.DepMap_ID)].set_index("DepMap_ID")
lin = d.set_index("DepMap_ID").lineage.reindex(full.index)
X = full.astype("float32")
Xc = X - X.groupby(lin.values).transform("mean")          # remove lineage means
y = Xc[col]
keep = Xc.notna().sum() >= 300
Xk = Xc.loc[:, keep].drop(columns=[col])
r = Xk.corrwith(y)
top = r.sort_values(ascending=False).head(20)
PART = {"STX1A","STX2","STX3","STX4","STX7","STX8","VAMP2","VAMP3","VAMP4","VAMP7","VAMP8",
        "SNAP25","SNAP29","NAPA","NAPB","NSF","STXBP1","STXBP2","STXBP3","SEC22B"}
top_syms = [t.split(" ")[0] for t in top.index]
v6 = {"top20_codependencies_positive_r": {k.split(" ")[0]: round(float(v), 3) for k, v in top.items()},
      "prespecified_SNARE_partners_in_top20": sorted(set(top_syms) & PART)}
v6["top20_anticorrelated(post_hoc)"] = {k.split(" ")[0]: round(float(v), 3) for k, v in r.sort_values().head(20).items()}
v6["rank_of_prespecified_partners(post_hoc)"] = {g: (int((r > r[c]).sum()) + 1, round(float(r[c]), 3))
    for g in sorted(PART) for c in [next((x for x in r.index if x.split(" ")[0] == g), None)] if c is not None}
v6["n_genes_tested"] = int(len(r))
v6["pass"] = bool(v6["prespecified_SNARE_partners_in_top20"])
R["V6_codependency"] = v6
r.rename("pearson_r_after_lineage_centering").sort_values().to_csv(OUT + "/snap23_codependency_all_genes.csv")

# ---------- V7 orthogonal ----------
v7 = {}
sg = pd.read_csv(SANGER, low_memory=False)
sg = sg.rename(columns={sg.columns[0]: "DepMap_ID"})
scol = [c for c in sg.columns if c.split(" ")[0] == GENE]
if scol:
    s = sg[["DepMap_ID", scol[0]]].copy(); s["y"] = pd.to_numeric(s[scol[0]], errors="coerce")
    s = s.merge(meta, on="DepMap_ID", how="left")
    s = s[s.y.notna() & s.lineage.notna() & (s.primary_disease != "Non-Cancerous")]
    sa, sb = s[s.lineage == LIN].y.values, s[s.lineage != LIN].y.values
    v7["sanger_project_score"] = {"n_bowel": int(len(sa)), "n_other": int(len(sb)),
        "median_bowel": float(np.median(sa)), "selectivity": sel(sa, sb), "p_one_sided": mw(sa, sb),
        "pass": bool(sel(sa, sb) < 0 and mw(sa, sb) < 0.05)}
R2 = pd.read_csv(os.path.join(RD, "D2_combined_gene_dep_scores.csv"), index_col=0).T
by_ccle = M.dropna(subset=["CCLEName"]).drop_duplicates("CCLEName").set_index("CCLEName").ModelID
strip = M.dropna(subset=["StrippedCellLineName"])
strip = strip[~strip.StrippedCellLineName.duplicated(keep=False)].set_index("StrippedCellLineName").ModelID
ids = [by_ccle[c] if c in by_ccle.index else strip.get(c.split("_")[0]) for c in R2.index]
R2["DepMap_ID"] = ids
rcol = [c for c in R2.columns if isinstance(c, str) and c.split(" ")[0] == GENE]
if rcol:
    rr = R2[["DepMap_ID", rcol[0]]].copy(); rr["y"] = pd.to_numeric(rr[rcol[0]], errors="coerce")
    rr = rr.dropna(subset=["DepMap_ID"]).drop_duplicates("DepMap_ID").merge(meta, on="DepMap_ID", how="left")
    rr = rr[rr.y.notna() & rr.lineage.notna() & (rr.primary_disease != "Non-Cancerous")]
    ra, rb = rr[rr.lineage == LIN].y.values, rr[rr.lineage != LIN].y.values
    v7["rnai_demeter2"] = {"n_bowel": int(len(ra)), "n_other": int(len(rb)),
        "median_bowel": float(np.median(ra)), "selectivity": sel(ra, rb), "p_one_sided": mw(ra, rb),
        "pass": bool(sel(ra, rb) < 0 and mw(ra, rb) < 0.05)}
v7["pass"] = any(v.get("pass") for v in v7.values() if isinstance(v, dict))
R["V7_orthogonal"] = v7

# ---------- verdict ----------
core = all(R[k]["pass"] for k in ["V1_primary", "V2_genotype", "V3_culture_format", "V4_subtype"])
extra = sum(bool(R[k]["pass"]) for k in ["V5_paralogue_buffering", "V6_codependency", "V7_orthogonal"])
R["verdict"] = {"core_V1_V4_all_pass": bool(core), "n_passed_of_V5_V7": int(extra),
                "verdict": "validated" if core and extra >= 2 else ("partly validated" if core else "not validated")}
json.dump(R, open(OUT + "/snap23_validation_summary.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(R, indent=1, ensure_ascii=False))
