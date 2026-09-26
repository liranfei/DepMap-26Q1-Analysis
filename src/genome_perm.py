"""Fast genome-wide label-permutation test of the FULL selection pipeline.
Ranks of each gene across cell lines do not depend on lineage labels, so they are computed once; each permutation only
re-sums ranks per (permuted) lineage. p-values replicate scipy.stats.mannwhitneyu(alternative='less', method='asymptotic')
(tie correction + continuity correction). Observed result is validated against the slow reference implementation first.
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy.stats import norm
import run_pipeline as rp

ap = argparse.ArgumentParser(); ap.add_argument("--data-dir", default="."); ap.add_argument("--out", default="out")
ap.add_argument("--n", type=int, default=1000); ap.add_argument("--exclude-noncancerous", action="store_true"); a = ap.parse_args()
df = rp.load(a.data_dir)
if a.exclude_noncancerous:
    df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
X = df[genes].apply(pd.to_numeric, errors="coerce").values.astype(np.float64)          # lines x genes
ok = ~np.isnan(X); N = ok.sum(0).astype(np.float64)
R = pd.DataFrame(X).rank(axis=0, method="average", na_option="keep").values; Rz = np.where(ok, R, 0.0)
tie = np.zeros(X.shape[1])
for j in range(X.shape[1]):
    v = X[ok[:, j], j]
    if len(v): _, c = np.unique(v, return_counts=True); tie[j] = ((c ** 3) - c).sum()
lin = df["lineage"].values; lc = pd.Series(lin).value_counts(); elig = list(lc[lc >= rp.MIN_N].index)
den = N * (N - 1); den[den == 0] = np.nan

def pvals(labels):
    L = np.stack([(labels == e).astype(np.float64) for e in elig])                  # lineages x lines
    n1 = L @ ok.astype(np.float64); R1 = L @ Rz; n2 = N[None, :] - n1
    U1 = R1 - n1 * (n1 + 1) / 2; U2 = n1 * n2 - U1; mu = n1 * n2 / 2
    s = np.sqrt(n1 * n2 / 12 * ((N[None, :] + 1) - tie[None, :] / den[None, :]))
    with np.errstate(divide="ignore", invalid="ignore"):
        p = norm.sf((U2 - mu - 0.5) / s)
    valid = (n1 >= rp.MIN_N) & (n2 > 0) & np.isfinite(p)
    return p, valid, n1

def bh(p):
    m = len(p); o = np.argsort(p); q = np.empty(m); r = p[o] * m / (np.arange(m) + 1)
    q[o] = np.minimum.accumulate(r[::-1])[::-1]; return np.minimum(q, 1.0)

def final_count(labels, p, valid):
    q = np.full(p.shape, np.nan); q[valid] = bh(p[valid])
    li, gj = np.where(valid & (q < rp.Q_MAX)); cnt = 0
    for i, j in zip(li, gj):
        m = labels == elig[i]; x = X[:, j]
        t = x[m & ok[:, j]]; o = x[~m & ok[:, j]]
        if len(t) and np.median(t) < rp.CHRONOS_MAX and np.median(t) - np.median(o) < rp.SEL_MAX: cnt += 1
    return cnt, int((valid & (q < rp.Q_MAX)).sum()), q

# ---- validation against the reference implementation
p, valid, n1 = pvals(lin); cnt, nq, q = final_count(lin, p, valid)
ref = pd.read_csv(os.path.join(a.out, "all_tests.csv")); li = {e: i for i, e in enumerate(elig)}; gi = {g: j for j, g in enumerate(genes)}
pv = np.array([p[li[l], gi[g]] for g, l in zip(ref.Gene.values, ref.Lineage.values)])
val = dict(n_valid=int(valid.sum()), n_ref=len(ref), max_abs_log10p_diff=float(np.nanmax(np.abs(np.log10(np.clip(pv, 1e-300, 1)) - np.log10(np.clip(ref.p_value.values, 1e-300, 1))))),
           n_q_lt_0_05=nq, n_q_ref=int((ref.q_value < rp.Q_MAX).sum()), n_final=cnt, n_final_ref=int(len(pd.read_csv(os.path.join(a.out, "final_targets.csv")))))
print("VALIDATION", json.dumps(val), flush=True)
assert val["n_valid"] == val["n_ref"] and val["n_q_lt_0_05"] == val["n_q_ref"] and val["n_final"] == val["n_final_ref"], "vectorized != reference"

# ---- permutations (seeds 42..)
res = []
for b in range(a.n):
    lab = np.random.default_rng(42 + b).permutation(lin); pp, vv, _ = pvals(lab); c, nq_, _ = final_count(lab, pp, vv); res.append((c, nq_))
    if b % 100 == 0: print(b, res[-1], flush=True)
r = pd.DataFrame(res, columns=["n_final_pairs", "n_q_lt_0_05"]); r.to_csv(os.path.join(a.out, "permutation_counts_genome.csv"), index=False)
x = r.n_final_pairs.values; out = dict(validation=val, observed_final=cnt, observed_q_lt_0_05=nq, n_perm=a.n, null_final_mean=float(x.mean()), null_final_max=int(x.max()),
        null_q_mean=float(r.n_q_lt_0_05.mean()), null_q_max=int(r.n_q_lt_0_05.max()), empirical_p_final=float((1 + (x >= cnt).sum()) / (1 + len(x))),
        empirical_p_q=float((1 + (r.n_q_lt_0_05.values >= nq).sum()) / (1 + len(r))))
json.dump(out, open(os.path.join(a.out, "permutation_summary_genome.json"), "w"), indent=1); print(json.dumps(out, indent=1))
