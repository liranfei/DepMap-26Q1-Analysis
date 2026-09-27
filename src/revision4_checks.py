"""Fourth-round checks (revision 4).

1. Welch (unequal-variance) t-tests for every gene x eligible lineage, BH over the whole family (sensitivity to the equal-variance Student t-test).
2. Permutation p-values for the three prioritised pairs with fewer than ten target lines (10,000 shuffles of target membership among non-missing lines).
3. Batch/library-adjusted lineage coefficients of the 94 pairs with 95% confidence intervals (same models in both fits).
4. DCAF7 coverage by screening library; target-median criterion after removing drug-adapted or engineered models.
5. TCGA: restricted cubic spline (3 knots at the 10th/50th/90th percentiles) versus linear expression in Cox models (likelihood-ratio test of the
   non-linear term), and a DESeq2-type median-of-ratios normalisation as a sensitivity analysis to the CPM normalisation.
6. Clopper-Pearson intervals for the Sanger replication proportions.
Output: revision4_summary.json, batch_coef_ci.csv."""
import os, re, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import t as tdist, mannwhitneyu, wilcoxon, chi2
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint
from lifelines import CoxPHFitter
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; TD = os.environ.get("TCGA_DIR", "."); CDR = os.environ.get("TCGA_CDR", ""); R = {}
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
X = df[genes].apply(pd.to_numeric, errors="coerce").values; lin = df.lineage.values
F = pd.read_csv(O + "final_targets.csv"); k94 = set(zip(F.Gene, F.Lineage))

# ---- 1. Welch t-tests (vectorised)
rows = []
for L in pd.Series(lin).value_counts()[lambda s: s >= rp.MIN_N].index:
    m = lin == L; A, B = X[m], X[~m]
    na, nb = (~np.isnan(A)).sum(0), (~np.isnan(B)).sum(0)
    ma, mb = np.nanmean(A, 0), np.nanmean(B, 0); va, vb = np.nanvar(A, 0, ddof=1), np.nanvar(B, 0, ddof=1)
    se = np.sqrt(va / na + vb / nb); tt = (ma - mb) / se; dfw = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    p = 2 * tdist.sf(np.abs(tt), dfw); ok = (na >= rp.MIN_N) & np.isfinite(p)
    rows.append(pd.DataFrame({"Gene": np.array(genes)[ok], "Lineage": L, "mean_diff": (ma - mb)[ok], "p_welch": p[ok]}))
W = pd.concat(rows, ignore_index=True); W["q_welch"] = multipletests(W.p_welch, method="fdr_bh")[1]
w94 = W[[k in k94 for k in zip(W.Gene, W.Lineage)]]
R["welch"] = {"n_tests": len(W), "q_lt_0_05_negative": int(((W.q_welch < 0.05) & (W.mean_diff < 0)).sum()), "primary94_q_lt_0_05_negative": int(((w94.q_welch < 0.05) & (w94.mean_diff < 0)).sum()), "primary94_max_q": float(w94.q_welch.max())}

# ---- 2. permutation p-values for low-n pairs
rng = np.random.default_rng(42); low = []
for r in F[F.n_target < 10].itertuples():
    y = df[r.Gene].astype(float); ok = y.notna().values; v = y.values[ok]; tgt = (lin[ok] == r.Lineage)
    u_obs = mannwhitneyu(v[tgt], v[~tgt], alternative="less", method="asymptotic").statistic
    null = np.array([mannwhitneyu(v[s], v[~s], alternative="less", method="asymptotic").statistic for s in (rng.permutation(tgt) for _ in range(10000))])
    low.append(dict(Gene=r.Gene, Lineage=r.Lineage, n_target=int(tgt.sum()), p_asymptotic=float(r.p_value), p_permutation=float((1 + (null <= u_obs).sum()) / (1 + len(null)))))
R["low_n_permutation"] = low

# ---- 3. batch/library-adjusted coefficients with CI
s = pd.read_csv(D + "/ScreenSequenceMap.csv"); s = s[(s.ScreenType == "2DS") & (s.PassesQC == True) & (s.DrugTreated == False) & (s.IsEngineered == False) & (s.ExcludeFromCRISPRCombined == False)]
mode = lambda x: x.value_counts().index[0]
bt = s.groupby("ModelID").pDNABatch.agg(mode); lb = s.groupby("ModelID").Library.agg(mode)
dfb = pd.DataFrame({"batch": df.DepMap_ID.map(bt).fillna("unknown"), "library": df.DepMap_ID.map(lb).fillna("unknown")})
rows = []
for r in F.itertuples():
    d = pd.DataFrame({"y": df[r.Gene].astype(float), "tgt": (lin == r.Lineage).astype(int), "batch": dfb.batch, "library": dfb.library}).dropna()
    m0 = smf.ols("y~tgt", d).fit(); m1 = smf.ols("y~tgt+C(batch)+C(library)", d).fit()
    rows.append(dict(Gene=r.Gene, Lineage=r.Lineage, n=len(d), coef_unadj=m0.params.tgt, ci_unadj_low=m0.conf_int().loc["tgt", 0], ci_unadj_high=m0.conf_int().loc["tgt", 1],
                     coef_adj=m1.params.tgt, ci_adj_low=m1.conf_int().loc["tgt", 0], ci_adj_high=m1.conf_int().loc["tgt", 1]))
BC = pd.DataFrame(rows); BC.to_csv(O + "batch_coef_ci.csv", index=False)
R["batch_ci"] = {"n_pairs": len(BC), "adj_ci_excludes_0": int((BC.ci_adj_high < 0).sum()), "min_ratio": float((BC.coef_adj / BC.coef_unadj).min()), "models_per_fit": int(BC.n.max()), "batch_unknown_models": int((dfb.batch == "unknown").sum())}

# ---- 4. DCAF7 coverage; engineered removal and the median criterion
lib_all = s.groupby("ModelID").Library.agg(lambda x: "+".join(sorted(set(x)))); g = [c for c in genes if c.startswith("DCAF7 (")][0]
cov = pd.DataFrame({"lib": df.DepMap_ID.map(lib_all).fillna("none"), "has": df[g].notna(), "lineage": lin})
R["DCAF7_coverage"] = {"non_missing": int(cov.has.sum()), "by_library_non_missing": cov[cov.has].lib.value_counts().to_dict(), "myeloid_non_missing": int((cov.has & (cov.lineage == "myeloid")).sum()), "myeloid_total": int((cov.lineage == "myeloid").sum())}
eng = df.DepMap_ID.map(pd.read_csv(D + "/Model.csv").set_index("ModelID").EngineeredModel).notna().values   # same definition as growth_confound.py
med = []
for r in F.itertuples():
    y = df[r.Gene].astype(float).values; ok = ~np.isnan(y) & ~eng; med.append(np.median(y[ok & (lin == r.Lineage)]))
R["engineered_median_check"] = {"n_removed": int(eng.sum()), "target_median_lt_-1_after_removal": int((np.array(med) < -1).sum()), "n_pairs": len(F)}

# ---- 5. TCGA spline and median-of-ratios sensitivity
def rcs(x, k):
    c = lambda u: np.clip(u, 0, None) ** 3
    return (c(x - k[0]) - c(x - k[1]) * (k[2] - k[0]) / (k[2] - k[1]) + c(x - k[2]) * (k[1] - k[0]) / (k[2] - k[1])) / (k[2] - k[0]) ** 2
def spline_lr(d, T, E):
    d = d[[T, E, "z"]].dropna(); d = d[d[T] > 0].copy(); k = np.percentile(d.z, [10, 50, 90]); d["s1"] = rcs(d.z.values, k)
    l0 = CoxPHFitter().fit(d[[T, E, "z"]], T, E).log_likelihood_; l1 = CoxPHFitter().fit(d[[T, E, "z", "s1"]], T, E).log_likelihood_
    return {"n": len(d), "events": int(d[E].sum()), "LR_nonlinear_p": float(chi2.sf(2 * (l1 - l0), 1))}
cdr = pd.read_csv(CDR, sep="\t"); cdr = cdr[cdr["sample"].str.split("-").str[3].str[:2] == "01"].drop_duplicates("_PATIENT").rename(columns={"_PATIENT": "pid"})
def med_ratio(path, ensg):
    Xc = pd.read_csv(path, sep="\t", index_col=0); C = np.power(2.0, Xc) - 1
    pos = (C > 0).all(axis=1); lg = np.log(C[pos]); sf = np.exp((lg.sub(lg.mean(axis=1), axis=0)).median(axis=0))
    row = [i for i in Xc.index if ensg in str(i)][0]; return np.log2(C.loc[row] / sf + 1)
for label, ct, ensg in [("KIRC_HNF1B", "KIRC", "ENSG00000275410"), ("PAAD_KRAS", "PAAD", "ENSG00000133703")]:
    m = pd.read_csv(O + f"tcga_{label}_survival_table.csv").merge(cdr[cdr["cancer type abbreviation"] == ct][["pid", "PFI", "PFI.time"]], on="pid", how="left")
    out = {"OS_spline": spline_lr(m, "OS.time", "OS"), "PFI_spline": spline_lr(m, "PFI.time", "PFI")}
    e = med_ratio(os.path.join(TD, f"TCGA-{ct}.star_counts.tsv"), ensg); code = pd.Series(e.index.str.split("-").str[3].str[:2], index=e.index)
    tum = e[code == "01"].groupby(e[code == "01"].index.str[:12]).mean().rename("x_mr").reset_index().rename(columns={"index": "pid"})
    mm = m.merge(tum, on="pid"); mm["z_mr"] = (mm.x_mr - mm.x_mr.mean()) / mm.x_mr.std(); c = CoxPHFitter().fit(mm[["OS.time", "OS", "z_mr"]], "OS.time", "OS").summary.loc["z_mr"]
    out["median_of_ratios_OS"] = {"HR": float(np.exp(c["coef"])), "CI": [float(np.exp(c["coef lower 95%"])), float(np.exp(c["coef upper 95%"]))], "p": float(c["p"]), "spearman_with_cpm": float(mm[["x", "x_mr"]].corr("spearman").iloc[0, 1])}
    if ct == "KIRC":
        nor = e[code == "11"].groupby(e[code == "11"].index.str[:12]).mean(); tt = e[code == "01"].groupby(e[code == "01"].index.str[:12]).mean(); b = tt.index.intersection(nor.index)
        out["median_of_ratios_paired"] = {"n": len(b), "median_difference": float((tt[b] - nor[b]).median()), "p": float(wilcoxon(tt[b], nor[b]).pvalue)}
    R[f"TCGA_{label}"] = out

# ---- 6. Clopper-Pearson intervals for replication
S = json.load(open(O + "sanger_disjoint_summary.json")); k, n = S["disc_replicated"], S["disc_testable_in_sanger"]
kb = round(S["baseline_frac_replicated"] * S["baseline_noncandidate_sig_testable"]); nb = S["baseline_noncandidate_sig_testable"]
R["replication_ci"] = {"candidates": [k, n, *proportion_confint(k, n, method="beta")], "reference": [kb, nb, *proportion_confint(kb, nb, method="beta")]}
json.dump(R, open(O + "revision4_summary.json", "w"), indent=1, default=float); print(json.dumps(R, indent=1, default=float))
