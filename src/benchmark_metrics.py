"""Benchmark of the selectivity/effect-size rule against published lineage-versus-rest metrics on the same test family.
Metrics per (gene, lineage): (1) Meyers-style difference in means with two-tailed pooled-variance t-test, BH over the family (Meyers 2017);
(2) Cohen's d (pooled SD); (3) AUC (probability that a target-lineage value is lower than a value of another lineage) from the Mann-Whitney U;
(4) our selectivity (difference in medians).  For each metric the top-K pairs (K = number of prioritised pairs) among pairs with q<0.05 and target median<-1 are compared with the prioritised list."""
import os, json, numpy as np, pandas as pd
from scipy.stats import t as tdist, rankdata
from statsmodels.stats.multitest import multipletests
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "primary") + "/"
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]; X = df[genes].apply(pd.to_numeric, errors="coerce").values.astype(float); ok = ~np.isnan(X); X0 = np.where(ok, X, 0.0)
lc = df.lineage.value_counts(); elig = list(lc[lc >= rp.MIN_N].index); N = ok.sum(0).astype(float); S = X0.sum(0); SS = (X0 ** 2).sum(0)
R = pd.DataFrame(X).rank(axis=0, method="average", na_option="keep").values; R0 = np.where(ok, R, 0.0)
T = pd.read_csv(O + "all_tests.csv"); rows = []
for lin in elig:
    m = (df.lineage == lin).values; n1 = ok[m].sum(0).astype(float); s1 = X0[m].sum(0); q1 = (X0[m] ** 2).sum(0); r1 = R0[m].sum(0); n2 = N - n1
    valid = (n1 >= rp.MIN_N) & (n2 > 0); m1 = s1 / n1; m2 = (S - s1) / n2
    v1 = (q1 - n1 * m1 ** 2) / (n1 - 1); v2 = ((SS - q1) - n2 * m2 ** 2) / (n2 - 1); sp2 = ((n1 - 1) * v1 + (n2 - 1) * v2) / (n1 + n2 - 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        tt = (m1 - m2) / np.sqrt(sp2 * (1 / n1 + 1 / n2)); p2 = 2 * tdist.sf(np.abs(tt), n1 + n2 - 2); d = (m1 - m2) / np.sqrt(sp2)
        U1 = r1 - n1 * (n1 + 1) / 2; auc_lower = 1 - U1 / (n1 * n2)          # P(target < other) + 0.5 P(tie)
    for j in np.where(valid)[0]: rows.append((genes[j], lin, m1[j] - m2[j], p2[j], d[j], auc_lower[j]))
B = pd.DataFrame(rows, columns=["Gene", "Lineage", "mean_diff", "p_t2", "cohens_d", "auc_lower"]); B["q_t"] = multipletests(B.p_t2, method="fdr_bh")[1]
M = T.merge(B, on=["Gene", "Lineage"]); assert len(M) == len(T)
F = pd.read_csv(O + "final_targets.csv"); K = len(F); key = lambda d: set(zip(d.Gene, d.Lineage)); ours = key(F)
base = M[(M.Chronos_median < rp.CHRONOS_MAX)]                          # same median-Chronos filter for all metrics
out = dict(n_prioritised=K, n_tests=len(M))
def topk(df_, col, asc, extra=None, k=K):
    d = df_ if extra is None else df_[extra]; return key(d.sort_values(col, ascending=asc).head(k))
# published-style rules
meyers = M[(M.q_t < rp.Q_MAX) & (M.mean_diff < 0)]; out["meyers_style_q<0.05_negative_meandiff"] = int(len(meyers))
out["meyers_style_with_effect_rule(median<-1, meandiff<-0.5)"] = int(((M.q_t < rp.Q_MAX) & (M.Chronos_median < -1) & (M.mean_diff < -0.5)).sum())
rk = {"selectivity (difference in medians)": ("Selectivity", True), "difference in means": ("mean_diff", True), "Cohen's d": ("cohens_d", True), "AUC (target lower)": ("auc_lower", False)}
sig = base[(base.q_value < rp.Q_MAX)]; rr = {}
for name, (col, asc) in rk.items():
    tk = topk(sig, col, asc); rr[name] = dict(overlap_with_prioritised=len(tk & ours), jaccard=len(tk & ours) / len(tk | ours))
    # rank of prioritised pairs under this metric among all sig pairs
    s2 = sig.assign(rank=sig[col].rank(ascending=asc, method="min")); s2 = s2[[k in ours for k in zip(s2.Gene, s2.Lineage)]]; rr[name]["median_rank_of_prioritised_among_sig"] = float(s2["rank"].median()); rr[name]["n_sig_pool"] = int(len(sig))
out["topK_comparison"] = rr
from scipy.stats import spearmanr
out["spearman_all_tests"] = {c: float(spearmanr(M.Selectivity, M[c])[0]) for c in ["mean_diff", "cohens_d", "auc_lower"]}
# thresholds analogous to Cohen's d and AUC: how many prioritised pairs have |d|>=0.8 / AUC>=0.8 (conventional large effects)
Fm = F.merge(B, on=["Gene", "Lineage"]); out["prioritised_cohens_d_median"] = float(Fm.cohens_d.median()); out["prioritised_frac_d_le_-0.8"] = float((Fm.cohens_d <= -0.8).mean()); out["prioritised_auc_median"] = float(Fm.auc_lower.median()); out["prioritised_frac_auc_ge_0.8"] = float((Fm.auc_lower >= 0.8).mean())
out["prioritised_in_meyers_style"] = int(sum(k in key(meyers) for k in ours))
B.to_csv(O + "benchmark_metrics_all_tests.csv.gz", index=False); json.dump(out, open(O + "benchmark_metrics_summary.json", "w"), indent=1); print(json.dumps(out, indent=1))
