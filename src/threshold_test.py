"""Formal test of the prioritisation rule (TREAT-style thresholded hypotheses), S12 Table.

The primary analysis tests H0: no lineage difference and then applies the thresholds (selectivity < -0.5, target median
< -1) as a filter, so the prioritised list has no false discovery rate guarantee. Here the thresholds themselves are
the null hypotheses (cf. McCarthy & Smyth 2009, TREAT):

  H0a: shift >= -0.5   one-sided Mann-Whitney test after adding 0.5 to the target-lineage values (target still lower?);
                       under a location-shift model this tests whether the target lineage is lower by more than 0.5.
  H0b: target median >= -1   one-sided sign test (exact binomial on the number of target lines with Chronos < -1).

A pair meets the rule formally when both nulls are rejected; the intersection-union p-value max(pa, pb) is a valid
p-value for the union null, and Benjamini-Hochberg (and Benjamini-Yekutieli) is applied over the same 467,091 tests as
the primary analysis. p-values replicate scipy.stats.mannwhitneyu(alternative='less', method='asymptotic'); with a
shift of 0 they are checked against all_tests.csv, and with the shift against scipy on random pairs.
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy.stats import norm, binom, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import run_pipeline as rp

ap = argparse.ArgumentParser(); ap.add_argument("--data-dir", default="."); ap.add_argument("--out", default="out")
ap.add_argument("--exclude-noncancerous", action="store_true"); a = ap.parse_args()
SHIFT = -rp.SEL_MAX                                                    # 0.5
df = rp.load(a.data_dir)
if a.exclude_noncancerous:
    df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
X = df[genes].apply(pd.to_numeric, errors="coerce").values.astype(np.float64)          # lines x genes
ok = ~np.isnan(X); N = ok.sum(0).astype(np.float64)
lin = df["lineage"].values; lc = pd.Series(lin).value_counts(); elig = sorted(lc[lc >= rp.MIN_N].index)


def tie_term(Y):
    """sum(t^3 - t) over tied groups, per column, ignoring NaN."""
    S = np.sort(Y, axis=0); eq = S[1:] == S[:-1]; out = np.zeros(Y.shape[1])
    for j in np.where(eq.any(0))[0]:
        v = Y[ok[:, j], j]; _, c = np.unique(v, return_counts=True); out[j] = ((c ** 3) - c).sum()
    return out


def mw_less(m, shift):
    """One-sided (target < rest) asymptotic Mann-Whitney p for all genes; target values shifted up by `shift`."""
    Y = X.copy(); Y[m] += shift
    R = pd.DataFrame(Y).rank(axis=0, method="average", na_option="keep").values
    okm = ok & m[:, None]
    n1 = okm.sum(0).astype(np.float64); n2 = N - n1
    R1 = np.where(okm, R, 0.0).sum(0)
    U1 = R1 - n1 * (n1 + 1) / 2; U2 = n1 * n2 - U1; mu = n1 * n2 / 2
    den = N * (N - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.sqrt(n1 * n2 / 12 * ((N + 1) - tie_term(Y) / den))
        p = norm.sf((U2 - mu - 0.5) / s)
    return p, n1, n2


rows = []
for e in elig:
    m = lin == e
    p0, n1, n2 = mw_less(m, 0.0)
    ps, _, _ = mw_less(m, SHIFT)
    Xt = np.where(ok & m[:, None], X, np.nan); Xo = np.where(ok & ~m[:, None], X, np.nan)
    with np.errstate(all="ignore"):
        med_t = np.nanmedian(Xt, axis=0); med_o = np.nanmedian(Xo, axis=0)
    k = (Xt < rp.CHRONOS_MAX).sum(0)                                   # NaN compares False
    psign = binom.sf(k - 1, n1.astype(int), 0.5)                       # P(K >= k) under median = -1
    valid = (n1 >= rp.MIN_N) & (n2 > 0)
    for j in np.where(valid)[0]:
        rows.append((genes[j], e, med_t[j], med_t[j] - med_o[j], int(n1[j]), p0[j], ps[j], psign[j]))
    print(e, int(valid.sum()), flush=True)
T = pd.DataFrame(rows, columns=["Gene", "Lineage", "Chronos_median", "Selectivity", "n_target", "p_value", "p_shift", "p_sign"])

# ---- validation 1: shift 0 reproduces the primary analysis
ref = pd.read_csv(os.path.join(a.out, "all_tests.csv"))
v = ref.merge(T, on=["Gene", "Lineage"], suffixes=("_ref", ""), validate="one_to_one")
lp = lambda x: np.log10(np.clip(x, 1e-300, 1))
val = dict(n_tests=len(T), n_ref=len(ref), n_matched=len(v),
           max_abs_log10p_diff=float(np.max(np.abs(lp(v.p_value) - lp(v.p_value_ref)))),
           max_abs_selectivity_diff=float(np.max(np.abs(v.Selectivity - v.Selectivity_ref))),
           max_abs_median_diff=float(np.max(np.abs(v.Chronos_median - v.Chronos_median_ref))))
assert val["n_tests"] == val["n_ref"] == val["n_matched"], val
assert val["max_abs_log10p_diff"] < 1e-6 and val["max_abs_selectivity_diff"] < 1e-9, val

# ---- validation 2: shifted p-values against scipy on random pairs (plus the 94 prioritised pairs)
F = pd.read_csv(os.path.join(a.out, "final_targets.csv"))
chk = pd.concat([T.sample(300, random_state=1), T.merge(F[["Gene", "Lineage"]], on=["Gene", "Lineage"])])
gi = {g: j for j, g in enumerate(genes)}; diffs = []
for g, e, ps in zip(chk.Gene, chk.Lineage, chk.p_shift):
    x = X[:, gi[g]]; m = (lin == e) & ~np.isnan(x); o = (lin != e) & ~np.isnan(x)
    pr = mannwhitneyu(x[m] + SHIFT, x[o], alternative="less", method="asymptotic")[1]
    diffs.append(abs(lp(np.array([ps]))[0] - lp(np.array([pr]))[0]))
val["shifted_checked_pairs"] = len(chk); val["shifted_max_abs_log10p_diff"] = float(max(diffs))
assert val["shifted_max_abs_log10p_diff"] < 1e-6, val
print("VALIDATION", json.dumps(val), flush=True)

# ---- multiplicity over the same family as the primary analysis
T["p_iut"] = np.maximum(T.p_shift, T.p_sign)
for col, src in [("q_shift", "p_shift"), ("q_iut", "p_iut")]:
    T[col] = multipletests(T[src], method="fdr_bh")[1]
T["q_iut_BY"] = multipletests(T.p_iut, method="fdr_by")[1]
T["prioritised"] = T.set_index(["Gene", "Lineage"]).index.isin(F.set_index(["Gene", "Lineage"]).index)
T["formal_rule"] = T.q_iut < rp.Q_MAX
T["formal_rule_BY"] = T.q_iut_BY < rp.Q_MAX
T.to_csv(os.path.join(a.out, "threshold_test_all.csv.gz"), index=False)

tiers = pd.read_csv(os.path.join(a.out, "confidence_tiers.csv"))[["Gene", "Lineage", "tier", "sel_ci_high"]]
P = T[T.prioritised | T.formal_rule | (T.q_shift < rp.Q_MAX)].merge(tiers, on=["Gene", "Lineage"], how="left")
P["boot_ci_below_threshold"] = P.sel_ci_high < rp.SEL_MAX
P.loc[~P.prioritised, "boot_ci_below_threshold"] = np.nan
blood = {"lymphoid", "myeloid"}
fr, pr_ = P[P.formal_rule], P[P.prioritised]
summ = dict(
    validation=val, shift=SHIFT, n_tests=len(T),
    shift_test_q_lt_0_05=int((T.q_shift < 0.05).sum()),
    shift_test_q_lt_0_05_and_median_lt_minus1=int(((T.q_shift < 0.05) & (T.Chronos_median < rp.CHRONOS_MAX)).sum()),
    formal_rule_BH=int(T.formal_rule.sum()), formal_rule_BY=int(T.formal_rule_BY.sum()),
    formal_rule_genes=int(fr.Gene.nunique()), formal_rule_lineages=int(fr.Lineage.nunique()),
    formal_rule_in_prioritised=int((fr.prioritised).sum()), formal_rule_not_prioritised=int((~fr.prioritised).sum()),
    formal_rule_not_prioritised_pairs=fr[~fr.prioritised][["Gene", "Lineage", "Chronos_median", "Selectivity", "n_target", "q_iut"]].round(4).to_dict("records"),
    prioritised_passing_shift_test=int((pr_.q_shift < 0.05).sum()),
    prioritised_not_formal_binding_test=dict(
        selectivity_threshold=int(((~pr_.formal_rule) & (pr_.p_shift >= pr_.p_sign)).sum()),
        median_threshold=int(((~pr_.formal_rule) & (pr_.p_sign > pr_.p_shift)).sum())),
    formal_rule_by_lineage=fr.Lineage.value_counts().to_dict(),
    formal_rule_blood=int(fr.Lineage.isin(blood).sum()),
    prioritised_by_tier_formal=pr_.groupby("tier").formal_rule.agg(["sum", "count"]).to_dict("index"),
    boot_ci_below_threshold_among_prioritised=int(pr_.boot_ci_below_threshold.sum()),
    formal_vs_bootstrap_crosstab=pd.crosstab(pr_.formal_rule, pr_.boot_ci_below_threshold).to_dict(),
    min_n_target_formal=int(fr.n_target.min()) if len(fr) else None,
    bonferroni_level_first_rejection=0.05 / len(T),
    shift_test_q_lt_0_05_BY=int((multipletests(T.p_shift, method="fdr_by")[1] < 0.05).sum()),
    shift_test_pairs=T[T.q_shift < 0.05].merge(tiers, on=["Gene", "Lineage"], how="left")
        [["Gene", "Lineage", "n_target", "Chronos_median", "Selectivity", "q_shift", "prioritised", "tier", "sel_ci_high"]].round(4).to_dict("records"),
    prioritised_unadjusted_p_shift_lt_0_05=int((pr_.p_shift < 0.05).sum()),
    prioritised_unadjusted_p_sign_lt_0_05=int((pr_.p_sign < 0.05).sum()),
    prioritised_unadjusted_p_iut_lt_0_05=int((pr_.p_iut < 0.05).sum()),
    prioritised_unadjusted_shift_vs_bootstrap=pd.crosstab(pr_.p_shift < 0.05, pr_.boot_ci_below_threshold).to_dict(),
    min_q_iut=float(T.q_iut.min()), min_q_iut_pairs=T.nsmallest(2, "p_iut")[["Gene", "Lineage", "p_iut", "q_iut"]].to_dict("records"))
json.dump(summ, open(os.path.join(a.out, "threshold_test_summary.json"), "w"), indent=1, default=int)
cols = ["Gene", "Lineage", "n_target", "Chronos_median", "Selectivity", "prioritised", "tier", "sel_ci_high",
        "p_shift", "q_shift", "p_sign", "p_iut", "q_iut", "q_iut_BY", "formal_rule", "formal_rule_BY"]
P.sort_values(["formal_rule", "q_iut"], ascending=[False, True])[cols].to_csv(os.path.join(a.out, "threshold_test_pairs.csv"), index=False)
print(json.dumps({k: v for k, v in summ.items() if k != "formal_rule_not_prioritised_pairs"}, indent=1, default=int))
