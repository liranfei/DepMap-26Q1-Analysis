"""Genome-like mixture of plan A, Section 4.4. Run only after v1b is deposited (and after culture_coef_fit.py).

Per replicate r = 0..19 (seed 30261119 + r), in this order of random draws:
    1. 18,531 genes sampled with replacement from the real genes: (mu, sigma, missing mask) jointly;
    2. a random permutation of the 18,531 simulated genes: the first round(0.02 n) = 371 are non-null, the next
       round(0.05 n) = 927 are culture-only, the rest null [N11];
    3. non-null genes: beta from {-0.5, -1.0, -1.5} with equal probability, then one lineage drawn in proportion to its
       number of lines among all lineages [N11];
    4. culture-only genes: delta drawn with replacement from the negative suspension coefficients of
       culture_coef_fit.py;
    5. noise t5 / sqrt(5/3), k = 1: y = mu + beta [lineage] + delta [suspension] + sigma eps; Other stratum effect 0.
Tests: every lineage, every gene; each method keeps the pairs testable under its rule (M1: >= 5 target lines with
data and >= 1 other; M2: Section 2; M3: the M2 family); BH over the whole family of each method.
Prioritised: M1 q < 0.05, median < -1, median difference < -0.5; M2 q < 0.05, median < -1, stratified HL < -0.5;
M3 q < 0.05, median < -1, coefficient < -0.5 (median = median of all target-lineage lines with data).
Truth: T1 beta <= -0.5 and mu + beta < -1; T2 beta + delta (f_t - f_r) < -0.5 and mu + beta + delta f_t < -1, with
f_t, f_r the suspension fractions of the target-lineage lines and of all other lines with data for the gene.
Empirical FDR = false prioritised / prioritised (0 if none prioritised; the count of replicates with none is reported).

Usage: python simulation_mixture.py --inputs DIR/sim_inputs.npz --coef DIR/culture_coefficients.csv --out DIR
       [--reps 0-19] [--workers 4] [--test]
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat_core as C  # noqa: E402

SEED0 = 30261119
BETAS = (-0.5, -1.0, -1.5)


def truth(mu_j, beta_j, beta_lin_j, delta_j, lin, lineage, susp, ok_j):
    """T1 and T2 of plan A, Section 4.4 for the pair (gene j, lineage lin)."""
    b = beta_j if beta_lin_j == lin else 0.0
    tg = lineage == lin
    ft = susp[ok_j & tg].mean()
    fr = susp[ok_j & ~tg].mean()
    t1 = bool((b <= -0.5) and (mu_j + b < -1))
    t2 = bool((b + delta_j * (ft - fr) < -0.5) and (mu_j + b + delta_j * ft < -1))
    return t1, t2


def run_rep(r, inp, delta_pool, test=False):
    t0 = time.time()
    rng = np.random.default_rng((990500000 if test else SEED0) + r)
    lineage, culture = inp["lineage"], inp["culture"]
    n = inp["mu"].shape[0]
    idx = rng.integers(0, n, size=n)
    mu, sd, mask = inp["mu"][idx], inp["sigma"][idx], inp["mask"][:, idx]
    perm = rng.permutation(n)
    n_nn, n_cu = int(round(0.02 * n)), int(round(0.05 * n))
    nn, cu = perm[:n_nn], perm[n_nn:n_nn + n_cu]
    lins_all = np.unique(lineage)
    sizes = pd.Series(lineage).value_counts()          # all lineages, in proportion to their numbers of lines
    elig = sizes
    beta = np.zeros(n)
    beta_lin = np.full(n, "", dtype=object)
    beta[nn] = rng.choice(BETAS, size=n_nn)
    beta_lin[nn] = rng.choice(elig.index.to_numpy(object), size=n_nn, p=(elig / elig.sum()).to_numpy())
    delta = np.zeros(n)
    delta[cu] = rng.choice(delta_pool, size=n_cu, replace=True)
    eps = rng.standard_t(5, size=(len(lineage), n)) / np.sqrt(5.0 / 3.0)
    susp = (culture == 1).astype(float)
    B = (lineage[:, None] == beta_lin[None, :]) * beta[None, :]
    Y = mu[None, :] + B + susp[:, None] * delta[None, :] + sd[None, :] * eps
    Y[~mask] = np.nan
    ok = ~np.isnan(Y)
    sr = C.StrataRanks(Y, culture)
    sr1 = C.StrataRanks(Y, np.zeros(len(lineage), int))
    recs = {m: [] for m in ("M1", "M2", "M3")}
    for lin in lins_all:
        tg = lineage == lin
        m1 = C.van_elteren(sr1, tg, min_per=1, min_target=C.MIN_TARGET)          # == Mann-Whitney (check C1)
        m2 = C.van_elteren(sr, tg)
        m3 = C.ols_hc3(Y, tg, culture)
        recs["M1"].append((lin, m1["testable"], m1["p"], None))
        recs["M2"].append((lin, m2["testable"], m2["p"], None))
        recs["M3"].append((lin, m2["testable"] & m3["ok"], m3["p"], m3["coef"]))
    out = {}
    for meth, lst in recs.items():
        L = [(lin, np.where(t)[0], p[t], None if c is None else c[t]) for lin, t, p, c in lst]
        allp = np.concatenate([x[2] for x in L])
        q = C.bh(allp)
        pos = 0
        prior, n_fam, n_sig = [], len(allp), int((q < C.Q_MAX).sum())
        for lin, js, p, coef in L:
            qq = q[pos:pos + len(js)]
            pos += len(js)
            tg = lineage == lin
            for kk in np.where(qq < C.Q_MAX)[0]:
                j = js[kk]
                y = Y[:, j]
                a_, o_ = y[ok[:, j] & tg], y[ok[:, j] & ~tg]
                med = float(np.median(a_))
                if not med < C.ABS_MAX:
                    continue
                if meth == "M1":
                    shift = med - float(np.median(o_))
                elif meth == "M2":
                    shift = C.m2_effects(y, tg, culture)["hl"]
                else:
                    shift = float(coef[kk])
                if shift < C.SHIFT_MAX:
                    prior.append((j, lin))
        tv = [truth(mu[j], beta[j], beta_lin[j], delta[j], lin, lineage, susp, ok[:, j]) for j, lin in prior]
        nf1 = sum(not a for a, _ in tv)
        nf2 = sum(not b for _, b in tv)
        # number of true pairs in the family (T1: non-null pairs with mu + beta < -1 and testable)
        true_t1 = 0
        for lin, js, p, coef in L:
            sel = (beta_lin[js] == lin) & (beta[js] <= -0.5) & (mu[js] + beta[js] < -1)
            true_t1 += int(sel.sum())
        out[meth] = dict(n_family=n_fam, n_q05=n_sig, n_prioritised=len(prior), n_false_T1=nf1, n_false_T2=nf2,
                         FDR_T1=nf1 / len(prior) if prior else 0.0, FDR_T2=nf2 / len(prior) if prior else 0.0,
                         n_true_T1_in_family=true_t1, n_true_T1_found=sum(a for a, _ in tv))
    return dict(rep=r, seed=(990500000 if test else SEED0) + r, n_nonnull=n_nn, n_culture_only=n_cu,
                runtime_s=round(time.time() - t0, 1), methods=out)


def _worker(args):
    r, path, delta_pool, test = args
    inp = dict(np.load(path, allow_pickle=True))
    return run_rep(r, inp, delta_pool, test)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True)
    ap.add_argument("--coef")
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", default="0-19")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--test", action="store_true", help="code test: seeds 990500000 + r; synthetic delta pool if no --coef")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.coef:
        cc = pd.read_csv(a.coef, comment="#").suspension_coefficient.to_numpy()
        pool = cc[cc < 0]
    else:
        assert a.test, "--coef is required outside test mode"
        pool = -np.abs(np.random.default_rng(1).normal(0.3, 0.1, 5000))
    lo, hi = map(int, a.reps.split("-"))
    from multiprocessing import Pool
    with Pool(a.workers) as p:
        res = p.map(_worker, [(r, a.inputs, pool, a.test) for r in range(lo, hi + 1)], chunksize=1)
    rows = [dict(rep=x["rep"], seed=x["seed"], method=m, **v) for x in res for m, v in x["methods"].items()]
    D = pd.DataFrame(rows)
    tag = a.reps.replace("-", "_") + ("_TEST" if a.test else "")
    meta = dict(plan_A_sha256=C.plan_a_sha256(), inputs_sha256=C.sha256_file(a.inputs), inputs_arrays_sha256=C.arrays_sha256(a.inputs),
                coef_sha256=C.sha256_file(a.coef) if a.coef else None, test_mode=a.test)
    with open(os.path.join(a.out, f"sim_mixture_{tag}.csv"), "w") as f:
        f.write(f"# plan_A_sha256={meta['plan_A_sha256']} inputs_sha256={meta['inputs_sha256']} coef_sha256={meta['coef_sha256']}\n")
        D.to_csv(f, index=False)
    S = D.groupby("method").agg(FDR_T1_mean=("FDR_T1", "mean"), FDR_T2_mean=("FDR_T2", "mean"),
                                n_prioritised_mean=("n_prioritised", "mean"),
                                reps_with_none=("n_prioritised", lambda x: int((x == 0).sum()))).reset_index()
    meta["summary"] = S.to_dict("records")
    meta["runtimes_s"] = [x["runtime_s"] for x in res]
    json.dump(meta, open(os.path.join(a.out, f"sim_mixture_summary_{tag}.json"), "w"), indent=1)
    print(S.to_string())


if __name__ == "__main__":
    main()
