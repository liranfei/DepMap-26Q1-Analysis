"""Third-round checks (revision 3).

1. MDM2 lineage candidates restricted to TP53 wild-type cell lines (TP53 status: high/moderate-impact variant = mutant).
2. Null distribution for the within-lineage co-dependency of the metabolic candidates: sets of common-essential, non-candidate genes matched to each
   metabolic candidate on pan-cancer median Chronos, residualised on lineage in the same way (1,000 draws).
3. Split-half reference for the Sanger comparison: the discovery cohort (lines without a KY screen) is split at random, within lineage, into two disjoint
   halves; selectivity is computed in each half and correlated across all tests (20 splits).  This gives the correlation expected between two disjoint sets
   of cell lines screened on the same platform.
4. Growth-pattern adjusted models: BH adjustment over the 94 pairs (the batch-adjusted models already use BH over the 94 pairs).
5. Consistency of the two mean-based rules (target mean < -1 vs target median < -1, each with difference in means < -0.5).
6. Distinct genes behind the KEGG pathways with q < 0.05 (restricted background); HNF1B in clear cell RCC lines only.
Output: revision3_summary.json, mdm2_tp53wt.csv, codep_null.csv, splithalf_selectivity.csv."""
import os, json, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; R = {}
mw = lambda x, y: float(mannwhitneyu(x, y, alternative="less", method="asymptotic")[1])
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
col = lambda s: [c for c in genes if c.startswith(s + " (")][0]
F = pd.read_csv(O + "final_targets.csv"); T = pd.read_csv(O + "all_tests.csv"); tg = sorted(F.Gene.unique())
X = df[genes].apply(pd.to_numeric, errors="coerce")

# ---- 1. MDM2 in TP53 wild-type lines
mut = pd.read_csv(O + "mut_3genes.csv"); mut = mut[mut.IsDefaultEntryForModel == "Yes"]
t53 = set(mut[(mut.HugoSymbol == "TP53") & mut.VepImpact.isin(["HIGH", "MODERATE"])].ModelID); wt = ~df.DepMap_ID.isin(t53)
y = X[col("MDM2")]; rows = []
for L in F[F.Gene == col("MDM2")].Lineage:
    tl = df.lineage == L
    for name, mask in (("all lines", pd.Series(True, index=df.index)), ("TP53 wild-type only", wt)):
        a = y[tl & mask].dropna(); b = y[~tl & mask].dropna()
        rows.append(dict(Lineage=L, cohort=name, n_target=len(a), n_target_TP53_mut=int((tl & ~wt).sum()) if name == "all lines" else 0, n_other=len(b),
                         median_target=float(a.median()), median_other=float(b.median()), selectivity=float(a.median() - b.median()), p=mw(a, b) if len(a) >= 5 else np.nan))
M = pd.DataFrame(rows); M.to_csv(O + "mdm2_tp53wt.csv", index=False); R["MDM2_TP53wt"] = M.round(4).to_dict("records")

# ---- 2. co-dependency null (matched common-essential genes)
res = lambda Z: Z - Z.groupby(df.lineage).transform("mean")
MET = [col(g) for g in ["NAMPT", "NMNAT1", "ADSL", "PAICS", "UMPS", "CTPS1", "DHFR", "TYMS"]]
iu = np.triu_indices(len(MET), 1); obs = float(np.median(np.abs(res(X[MET]).corr().values[iu])))
ce = set(pd.read_csv(D + "/CRISPRInferredCommonEssentials.csv").iloc[:, 0]); med = X.median()
pool = [g for g in genes if g in ce and g not in tg and X[g].notna().mean() > 0.95]; pm = med[pool]
near = {g: pm.iloc[np.argsort(np.abs(pm.values - med[g]))[:30]].index.tolist() for g in MET}
rng = np.random.default_rng(42); XR = res(X[sorted(set(sum(near.values(), [])))]); null = []
for _ in range(1000):
    pick = []
    for g in MET:
        c = [x for x in near[g] if x not in pick]; pick.append(c[rng.integers(len(c))])
    null.append(float(np.median(np.abs(XR[pick].corr().values[iu]))))
null = np.array(null); pd.DataFrame({"median_abs_r_null": null}).to_csv(O + "codep_null.csv", index=False)
R["codep_null"] = {"observed_median_abs_r_residual": obs, "null_mean": float(null.mean()), "null_95pct": float(np.percentile(null, 95)), "null_max": float(null.max()),
                   "empirical_p": float((1 + (null >= obs).sum()) / (1 + len(null))), "pool_size": len(pool), "matched_neighbours": 30}

# ---- 3. split-half reference for disjoint cell-line sets
s = pd.read_csv(D + "/ScreenSequenceMap.csv"); u = s[(s.PassesQC == True) & (s.ExcludeFromCRISPRCombined == False) & s.Library.notna()]
ky = set(u.loc[u.Library == "KY", "ModelID"]); dd = df[~df.DepMap_ID.isin(ky)].reset_index(drop=True); XD = dd[genes].apply(pd.to_numeric, errors="coerce").values
def sel(idx):
    lin = dd.lineage.values[idx]; V = XD[idx]; out = {}
    for L in np.unique(lin):
        m = lin == L
        if m.sum() < rp.MIN_N: continue
        a = V[m]; ok = (~np.isnan(a)).sum(0) >= rp.MIN_N
        s_ = np.nanmedian(a, 0) - np.nanmedian(V[~m], 0); s_[~ok] = np.nan; out[L] = s_
    return out
rows = []
for k in range(20):
    rng = np.random.default_rng(100 + k); h1 = []
    for L, idx in dd.groupby("lineage").groups.items():
        idx = np.array(list(idx)); rng.shuffle(idx); h1 += list(idx[: len(idx) // 2])
    h1 = np.array(sorted(h1)); h2 = np.setdiff1d(np.arange(len(dd)), h1); A, B = sel(h1), sel(h2)
    va, vb = [], []
    for L in set(A) & set(B): ok = ~np.isnan(A[L]) & ~np.isnan(B[L]); va.append(A[L][ok]); vb.append(B[L][ok])
    va, vb = np.concatenate(va), np.concatenate(vb); rows.append(dict(split=k, n_tests=len(va), n_lineages=len(set(A) & set(B)), spearman=float(spearmanr(va, vb)[0])))
SH = pd.DataFrame(rows); SH.to_csv(O + "splithalf_selectivity.csv", index=False)
R["splithalf"] = {"n_splits": len(SH), "half_size": int(len(dd) // 2), "spearman_median": float(SH.spearman.median()), "spearman_min": float(SH.spearman.min()), "spearman_max": float(SH.spearman.max()), "median_n_tests": int(SH.n_tests.median())}

# ---- 4. growth-pattern adjusted models: BH over the 94 pairs
G = pd.read_csv(O + "growth_confound_pairs.csv"); G["q_growth_adjusted"] = multipletests(G.p_growth_adjusted_two_sided, method="fdr_bh")[1]
G["ratio"] = G.coef_growth_adjusted / G.coef_unadjusted; G["heme"] = G.Lineage.isin(["lymphoid", "myeloid"]); G.to_csv(O + "growth_confound_pairs.csv", index=False)
R["growth_bh"] = {k: {"n": int(len(g)), "q_lt_0.05_negative": int(((g.q_growth_adjusted < 0.05) & (g.coef_growth_adjusted < 0)).sum()), "median_ratio": float(g.ratio.median()),
                      "ratio_lt_0.5": int((g.ratio < 0.5).sum()), "min_ratio": float(g.ratio.min())} for k, g in (("heme", G[G.heme]), ("nonheme", G[~G.heme]), ("all", G))}
hs = G[G.heme].dropna(subset=["p_vs_suspension_one_sided"]); qs = multipletests(hs.p_vs_suspension_one_sided, method="fdr_bh")[1]
R["suspension_bh"] = {"testable": len(hs), "not_testable_lt10_comparators": int(G.heme.sum() - len(hs)), "negative": int((hs.median_diff_vs_suspension < 0).sum()), "lt_-0.5": int((hs.median_diff_vs_suspension < -0.5).sum()), "q_lt_0.05": int((qs < 0.05).sum()), "median_diff": float(hs.median_diff_vs_suspension.median())}
qe = multipletests(G.p_excl_engineered_one_sided, method="fdr_bh")[1]; R["excl_engineered_bh"] = {"q_lt_0.05": int((qe < 0.05).sum()), "sel_lt_-0.5": int((G.median_diff_excl_engineered < -0.5).sum())}
B = pd.read_csv(O + "batch_library_adjusted_pairs.csv"); R["batch_ratio"] = {"min": float(B.ratio.min()), "median": float(B.ratio.median()), "q_lt_0.05": int((B.q_adj < 0.05).sum())}

# ---- 5. the two mean-based rules
TT = pd.read_csv(O + "ttest_alt_all_tests.csv.gz"); BM = pd.read_csv(O + "benchmark_metrics_all_tests.csv.gz")
R["mean_rule_columns"] = {"ttest_alt": TT.columns.tolist(), "benchmark": BM.columns.tolist()}

# ---- 6. KEGG distinct genes; clear cell RCC
E = pd.read_csv(O + "enrichment_restricted_background.csv"); k = E[(E.library == "KEGG_2021_Human") & (E.q < 0.05)]
gs = pd.Series(";".join(k.genes).split(";")).value_counts(); R["kegg_restricted"] = {"n_pathways": len(k), "n_distinct_genes": int(len(gs)), "genes_in_ge_10_pathways": gs[gs >= 10].to_dict()}
md = pd.read_csv(D + "/Model.csv")[["ModelID", "OncotreeSubtype"]].set_index("ModelID").OncotreeSubtype; df["subtype"] = df.DepMap_ID.map(md)
R["kidney_subtypes_in_cohort"] = df[df.lineage == "kidney"].subtype.value_counts().to_dict()
ccr = df.subtype == "Renal Clear Cell Carcinoma"; o = df.lineage != "kidney"
for g in ["HNF1B", "PAX8"]:
    a = X.loc[ccr, col(g)].dropna(); b = X.loc[o, col(g)].dropna(); R[f"{g}_clear_cell"] = {"n": len(a), "median": float(a.median()), "selectivity": float(a.median() - b.median()), "p": mw(a, b)}
json.dump(R, open(O + "revision3_summary.json", "w"), indent=1, default=float); print(json.dumps(R, indent=1, default=float))
