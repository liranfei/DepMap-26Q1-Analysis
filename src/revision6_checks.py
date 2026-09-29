"""Sixth-round checks (revision 6).

1. Co-dependency null drawn from all non-candidate genes (not only common-essential genes), matched on median and standard deviation of the gene effect.
2. Bootstrap 95% intervals of selectivity for the 94 prioritised pairs (1,000 resamples of target and comparison lines separately, seed 42).
3. Median of the comparison group and number of other eligible lineages with median Chronos < -1, per prioritised pair.
4. Coverage of the candidate genes (fraction of cancer lines with data).
5. Genotype: CTNNB1 dependency by APC truncating (HIGH impact) variant; BRAF dependency by BRAF V600 variant (default entries).
6. 22Q1 cross-release comparison with a quantile-matched selectivity threshold (same tail fraction as -0.5 in 26Q1, shared lines).
7. Non-blood pairs: selectivity against adherent lines of other lineages only.
8. Number of non-blood suspension lines in the Project Score data.
Output: revision6_summary.json, revision6_pairs.csv."""
import os, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; R = {}
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]; col = lambda s: [c for c in genes if c.startswith(s + " (")][0]
X = df[genes].apply(pd.to_numeric, errors="coerce"); F = pd.read_csv(O + "final_targets.csv"); tg = set(F.Gene)
M = pd.read_csv(os.path.join(D, "Model.csv")).set_index("ModelID"); gp = df.DepMap_ID.map(M.GrowthPattern).fillna("Unknown").values
lc = df.lineage.value_counts(); elig = list(lc[lc >= rp.MIN_N].index)

# ---- 1. co-dependency null from all non-candidate genes, matched on median and SD
res = lambda Z: Z - Z.groupby(df.lineage).transform("mean")
MET = [col(g) for g in ["NAMPT", "NMNAT1", "ADSL", "PAICS", "UMPS", "CTPS1", "DHFR", "TYMS"]]; iu = np.triu_indices(len(MET), 1)
obs = float(np.median(np.abs(res(X[MET]).corr().values[iu])))
pool = [g for g in genes if g not in tg and X[g].notna().mean() > 0.95]; md, sd = X[pool].median(), X[pool].std()
z = lambda v, s: (v - s.mean()) / s.std()
near = {g: md.index[np.argsort(np.hypot(z(md, md) - (X[g].median() - md.mean()) / md.std(), z(sd, sd) - (X[g].std() - sd.mean()) / sd.std()).values)[:30]].tolist() for g in MET}
rng = np.random.default_rng(42); XR = res(X[sorted(set(sum(near.values(), [])))]); null = []
for _ in range(1000):
    pick = []
    for g in MET:
        c = [x for x in near[g] if x not in pick]; pick.append(c[rng.integers(len(c))])
    null.append(float(np.median(np.abs(XR[pick].corr().values[iu]))))
null = np.array(null)
R["codep_null_all_genes_median_sd"] = {"observed": obs, "null_mean": float(null.mean()), "null_95pct": float(np.percentile(null, 95)), "null_max": float(null.max()),
                                       "empirical_p": float((1 + (null >= obs).sum()) / (1 + len(null))), "pool_size": len(pool)}

# ---- 2-4. per-pair bootstrap, comparison-group median, other lineages below -1, coverage
rng = np.random.default_rng(42); rows = []
for _, r in F.iterrows():
    v = X[r.Gene]; m = (df.lineage == r.Lineage).values; a = v[m].dropna().values; b = v[~m].dropna().values
    bs = np.array([np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(1000)])
    others = [L for L in elig if L != r.Lineage and v[(df.lineage == L).values].notna().sum() >= rp.MIN_N and v[(df.lineage == L).values].median() < -1]
    rows.append(dict(Gene=r.Gene, Lineage=r.Lineage, selectivity=r.Selectivity, sel_ci_low=float(np.percentile(bs, 2.5)), sel_ci_high=float(np.percentile(bs, 97.5)),
                     p_boot_sel_lt_minus0_5=float((bs < -0.5).mean()), median_comparison=float(np.median(b)), n_other_lineages_median_lt_minus1=len(others),
                     coverage=float(v.notna().mean())))
P = pd.DataFrame(rows)
R["bootstrap"] = {"ci_high_lt_minus0_5": int((P.sel_ci_high < -0.5).sum()), "ci_high_lt_0": int((P.sel_ci_high < 0).sum()), "median_p_boot": float(P.p_boot_sel_lt_minus0_5.median())}
R["comparison_median"] = {"gt_minus0_5": int((P.median_comparison > -0.5).sum()), "between_minus1_and_minus0_5": int(((P.median_comparison <= -0.5) & (P.median_comparison > -1)).sum()),
                          "le_minus1": int((P.median_comparison <= -1).sum()), "pairs_with_other_lineage_below_minus1": int((P.n_other_lineages_median_lt_minus1 > 0).sum())}
R["coverage_lt_0_5"] = P[P.coverage < 0.5][["Gene", "Lineage", "coverage"]].to_dict("records")

# ---- 5. genotype: APC truncating vs CTNNB1; BRAF V600 vs BRAF
mu = pd.read_csv(os.path.join(D, "OmicsSomaticMutations.csv"), usecols=["ModelID", "HugoSymbol", "VepImpact", "ProteinChange", "IsDefaultEntryForModel"], low_memory=False)
mu = mu[(mu.IsDefaultEntryForModel == "Yes") & mu.HugoSymbol.isin(["APC", "BRAF"])]
apc = set(mu[(mu.HugoSymbol == "APC") & (mu.VepImpact == "HIGH")].ModelID); v600 = set(mu[(mu.HugoSymbol == "BRAF") & mu.ProteinChange.astype(str).str.startswith("p.V600")].ModelID)
def geno(gene, carriers, lineage):
    d = pd.DataFrame({"y": X[col(gene)], "g": df.DepMap_ID.isin(carriers).astype(int), "lin": df.lineage}).dropna(); fit = smf.ols("y ~ g + C(lin)", d).fit()
    w = d[d.lin == lineage]
    return {"n_carrier": int(d.g.sum()), "n_other": int((1 - d.g).sum()), "median_carrier": float(d[d.g == 1].y.median()), "median_other": float(d[d.g == 0].y.median()),
            "p_one_sided": float(mannwhitneyu(d[d.g == 1].y, d[d.g == 0].y, alternative="less")[1]), "lineage_adjusted_diff": float(fit.params["g"]), "lineage_adjusted_p": float(fit.pvalues["g"]),
            f"in_{lineage}": {"n_carrier": int(w.g.sum()), "n_other": int((1 - w.g).sum()), "median_carrier": float(w[w.g == 1].y.median()) if w.g.sum() else None,
                              "median_other": float(w[w.g == 0].y.median()) if (1 - w.g).sum() else None,
                              "p_one_sided": float(mannwhitneyu(w[w.g == 1].y, w[w.g == 0].y, alternative="less")[1]) if w.g.sum() and (1 - w.g).sum() else None}}
R["APC_trunc_CTNNB1"] = geno("CTNNB1", apc, "bowel"); R["BRAF_V600_BRAF"] = geno("BRAF", v600, "skin")

# ---- 6. 22Q1 with a quantile-matched selectivity threshold (shared lines and genes)
s26 = pd.read_csv(O + "all_tests_26Q1_sharedlines.csv.gz"); s22 = pd.read_csv(O + "all_tests_22Q1_sharedlines.csv.gz")
frac = float((s26.Selectivity < -0.5).mean()); thr = float(s22.Selectivity.quantile(frac))
k = F[["Gene", "Lineage"]].merge(s22, on=["Gene", "Lineage"])
R["cross_release_quantile"] = {"tail_fraction_26Q1": frac, "matched_threshold_22Q1": thr, "testable": len(k),
                               "meet_rule_matched_threshold": int(((k.q_value < 0.05) & (k.Chronos_median < -1) & (k.Selectivity < thr)).sum()),
                               "meet_rule_minus0_5": int(((k.q_value < 0.05) & (k.Chronos_median < -1) & (k.Selectivity < -0.5)).sum())}

# ---- 7. non-blood pairs against adherent lines only
rows = []
for _, r in F[~F.Lineage.isin(["lymphoid", "myeloid"])].iterrows():
    v = X[r.Gene]; m = (df.lineage == r.Lineage).values; o = ~m & (gp == "Adherent")
    rows.append(dict(Gene=r.Gene, Lineage=r.Lineage, selectivity=r.Selectivity, selectivity_vs_adherent=float(v[m].median() - v[o].median())))
A = pd.DataFrame(rows); P = P.merge(A[["Gene", "Lineage", "selectivity_vs_adherent"]], on=["Gene", "Lineage"], how="left")
R["nonblood_vs_adherent"] = {"n": len(A), "lt_minus0_5": int((A.selectivity_vs_adherent < -0.5).sum()), "median_change": float((A.selectivity_vs_adherent - A.selectivity).median()),
                             "max_change": float((A.selectivity_vs_adherent - A.selectivity).max()),
                             "FERMT2": A[A.Gene.str.startswith("FERMT2")][["Lineage", "selectivity", "selectivity_vs_adherent"]].round(3).to_dict("records"),
                             "not_below_minus0_5": A[A.selectivity_vs_adherent >= -0.5][["Gene", "Lineage", "selectivity", "selectivity_vs_adherent"]].round(3).to_dict("records")}

# ---- 8. non-blood suspension lines in the Project Score data
sd_ = os.environ.get("SANGER_DIR")
if sd_:
    ids = pd.read_csv(os.path.join(sd_, "gene_effect.csv"), usecols=[0]).iloc[:, 0]; mm = M.reindex(ids.values)
    ok = (mm.OncotreePrimaryDisease != "Non-Cancerous").values
    R["sanger_nonblood_suspension"] = int((ok & (mm.GrowthPattern == "Suspension").values & ~mm.OncotreeLineage.isin(["Lymphoid", "Myeloid"]).values).sum())

P.to_csv(O + "revision6_pairs.csv", index=False); json.dump(R, open(O + "revision6_summary.json", "w"), indent=1, default=float); print(json.dumps(R, indent=1, default=float))
