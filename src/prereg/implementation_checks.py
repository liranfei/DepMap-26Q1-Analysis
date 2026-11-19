"""Implementation checks of plan A, Section 3 (archived at v1b, before any simulation is run).

C1  single stratum: M2 p == scipy.stats.mannwhitneyu(alternative="less", method="asymptotic"); also M1 and M4 code
C2  M2 (vectorised) == separately written per-stratum implementation, agreement < 1e-10
C3  exact within-stratum permutation p vs normal approximation with c = 0: |dp| < 0.01 or relative difference < 10%
C4  stratified HL: one stratum == ordinary HL; several strata == brute force with Fraction weights (exact)
C5  M3 / M3u == statsmodels OLS HC3 (use_t), coefficient, SE and p
C6  secondary effect size, vectorised == per gene
C7  cohort, culture-format and library counts at model level (annotations only)
C8  coverage expectations at model level (Section 2: Eye 0.70, Ampulla of Vater untestable)
C9  real gene effects with unrestricted permuted lineage labels (one permutation, seed 7, not one of the seeds
    42 + b of the calibration check): whole M2/M3 code path runs; vectorised == reference on 300 random pairs
C10 truth() of the mixture (Section 4.4) == hand-computed T1/T2 on a 10-line example: non-null gene in its own and in
    another lineage, culture-only gene in two lineages, culture-only gene with a missing line
C11 M3 testability == independent determination per gene (design with the populated culture levels has full column
    rank, n > p and maximal leverage < 1 - 1e-8), on simulated data with small target groups and missing values

Only simulated data and data with permuted lineage labels are used. Usage:
    python implementation_checks.py --data-dir $DEPMAP_DIR --out $RESULTS_DIR/prereg_v1b
"""
import argparse
import json
import os
import sys
import time

import numpy as np
from scipy.stats import mannwhitneyu, norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat_core as C          # noqa: E402
import strat_reference as R     # noqa: E402
from env_check import environment  # noqa: E402


def sim_matrix(rng, n, g, n_strata, tie_round=None, miss=0.0):
    Y = rng.normal(size=(n, g))
    if tie_round:
        Y = np.round(Y * tie_round) / tie_round
    if miss:
        Y[rng.random((n, g)) < miss] = np.nan
    strata = rng.integers(0, n_strata, size=n)
    return Y, strata


def c1(rng):
    worst_m2 = worst_m1 = worst_m4 = 0.0
    n_cmp = 0
    for rep in range(40):
        n = int(rng.integers(12, 300))
        Y, _ = sim_matrix(rng, n, 25, 1, tie_round=[None, 10, 3][rep % 3], miss=[0, 0.1][rep % 2])
        tgt = rng.random(n) < rng.uniform(0.05, 0.5)
        if tgt.sum() < 5:
            tgt[:5] = True
        Y[tgt] -= rng.uniform(0, 0.8)
        sr = C.StrataRanks(Y, np.zeros(n, int))
        v = C.van_elteren(sr, tgt, min_per=1, min_target=5)
        m1 = C.mann_whitney_less(Y, tgt)
        cmp_ = rng.random(n) < 0.6
        m4 = C.mann_whitney_less(Y, tgt, comparator=cmp_)
        for j in range(Y.shape[1]):
            ok = ~np.isnan(Y[:, j])
            a, b = Y[ok & tgt, j], Y[ok & ~tgt, j]
            if len(a) < 5 or len(b) < 1 or np.unique(np.r_[a, b]).size < 2:
                continue
            ps = mannwhitneyu(a, b, alternative="less", method="asymptotic").pvalue
            worst_m2 = max(worst_m2, abs(v["p"][j] - ps))
            worst_m1 = max(worst_m1, abs(m1["p"][j] - ps))
            b4 = Y[ok & ~tgt & cmp_, j]
            if len(b4) >= 1 and np.unique(np.r_[a, b4]).size >= 2:
                worst_m4 = max(worst_m4, abs(m4["p"][j] - mannwhitneyu(a, b4, alternative="less", method="asymptotic").pvalue))
            n_cmp += 1
    ok = max(worst_m2, worst_m1, worst_m4) < 1e-12
    return dict(n_genes_compared=n_cmp, max_abs_diff_M2=worst_m2, max_abs_diff_M1=worst_m1, max_abs_diff_M4=worst_m4,
                criterion="< 1e-12", passed=bool(ok))


def c2(rng):
    worst, n_cmp, n_testable_mismatch = 0.0, 0, 0
    for rep in range(30):
        n = int(rng.integers(60, 400))
        S = [3, 3, 15][rep % 3]
        Y, strata = sim_matrix(rng, n, 30, S, tie_round=[None, 20, 4][rep % 3], miss=[0, 0.05, 0.2][rep % 3])
        tgt = rng.random(n) < rng.uniform(0.05, 0.4)
        Y[tgt] -= rng.uniform(0, 0.5)
        for mp in (3, 5):
            sr = C.StrataRanks(Y, strata)
            v = C.van_elteren(sr, tgt, min_per=mp)
            for j in range(Y.shape[1]):
                p_ref, t_ref, *_ = R.van_elteren_one(Y[:, j], tgt, strata, min_per=mp)
                if bool(v["testable"][j]) != bool(t_ref):
                    n_testable_mismatch += 1
                    continue
                if t_ref:
                    worst = max(worst, abs(v["p"][j] - p_ref))
                    n_cmp += 1
    return dict(n_pairs_compared=n_cmp, testable_mismatches=n_testable_mismatch, max_abs_diff_p=worst,
                criterion="< 1e-10 and no testability mismatch", passed=bool(worst < 1e-10 and n_testable_mismatch == 0))


def c3(rng):
    rows = []
    designs = [[(10, 4), (12, 5), (9, 3)], [(14, 5), (11, 4)], [(8, 3), (8, 3), (10, 4), (9, 3)], [(16, 6), (12, 4)]]
    for ex in range(24):
        des = designs[ex % len(designs)]
        y, tg, st = [], [], []
        shift = [0.0, -0.4, -0.8, -1.2][ex % 4]
        rnd = [None, 4][ex % 2]
        for k, (n, n1) in enumerate(des):
            v = rng.normal(size=n)
            t_ = np.zeros(n, bool)
            t_[rng.choice(n, n1, replace=False)] = True
            v[t_] += shift
            if rnd:
                v = np.round(v * rnd) / rnd
            y.append(v); tg.append(t_); st.append(np.full(n, k))
        y, tg, st = np.concatenate(y), np.concatenate(tg), np.concatenate(st)
        _, _, T, E, V = R.van_elteren_one(y, tg, st, min_per=1, min_target=1)
        p_norm0 = float(norm.cdf((T - E) / np.sqrt(V)))
        p_exact = float(R.exact_van_elteren_p(y, tg, st))
        sr = C.StrataRanks(y[:, None], st)
        pv = C.van_elteren(sr, tg, min_per=1, min_target=1, return_parts=True)
        p_core0 = float(norm.cdf((pv["T"][0] - pv["E"][0]) / np.sqrt(pv["V"][0])))
        d = abs(p_exact - p_norm0)
        rel = d / p_exact if p_exact > 0 else np.inf
        rows.append(dict(design=des, shift=shift, ties=bool(rnd), p_exact=p_exact, p_normal_c0=p_norm0,
                         p_normal_c0_core=p_core0, abs_diff=d, rel_diff=rel, ok=bool(d < 0.01 or rel < 0.10)))
    core_match = max(abs(r["p_normal_c0"] - r["p_normal_c0_core"]) for r in rows)
    return dict(examples=rows, n_ok=sum(r["ok"] for r in rows), n=len(rows), core_vs_reference_max_diff=core_match,
                criterion="every example |dp| < 0.01 or relative < 10%", passed=bool(all(r["ok"] for r in rows) and core_match < 1e-12))


def c4(rng):
    single_worst, multi_mismatch, n_single, n_multi = 0.0, 0, 0, 0
    for rep in range(60):
        a = rng.normal(size=int(rng.integers(3, 40)))
        b = rng.normal(size=int(rng.integers(3, 60)))
        if rep % 2:
            a, b = np.round(a * 2) / 2, np.round(b * 2) / 2
        hl1 = C.stratified_hl([(a[:, None] - b[None, :]).ravel()], [len(a) + len(b)])
        single_worst = max(single_worst, abs(hl1 - C.hl_shift(a, b)))
        n_single += 1
    for rep in range(60):
        n = int(rng.integers(12, 45))
        y = rng.normal(size=n)
        if rep % 3 == 0:
            y = np.round(y * 2) / 2                     # many ties, exact halves occur
        st = rng.integers(0, 3, size=n)
        tg = rng.random(n) < 0.4
        ref = R.weighted_median_fraction(y, tg, st, min_per=3)
        e = C.m2_effects(y, tg, st, min_per=3)
        if not np.isfinite(e["hl"]) or not np.isfinite(ref):
            multi_mismatch += int(np.isfinite(e["hl"]) != np.isfinite(ref))
            continue
        n_multi += 1
        if e["hl"] != ref:
            multi_mismatch += 1
    # explicit exact-half example: one stratum, 2 x 2 differences {0, 1, 2, 3}: median (1 + 2) / 2 = 1.5
    half = C.stratified_hl([np.array([0.0, 1.0, 2.0, 3.0])], [3])
    return dict(single_stratum_max_abs_diff=single_worst, n_single=n_single, multi_stratum_examples=n_multi,
                multi_stratum_mismatches=multi_mismatch, exact_half_example=half,
                criterion="single == ordinary HL (0 difference); multi == Fraction brute force (exact); half example 1.5",
                passed=bool(single_worst == 0.0 and multi_mismatch == 0 and half == 1.5))


def c5(rng):
    worst = dict(coef=0.0, se=0.0, p=0.0)
    n_cmp = 0
    for rep in range(25):
        n = int(rng.integers(40, 300))
        Y = rng.normal(size=(n, 12)) * rng.uniform(0.5, 2)
        Y[rng.random(Y.shape) < 0.1] = np.nan
        tgt = rng.random(n) < 0.2
        cult = rng.choice(3, size=n, p=[0.7, 0.2, 0.1])
        if rep % 5 == 0:
            Y[cult == 2, 0] = np.nan                    # a culture level without data for gene 0 [N4]
        Y[tgt] -= 0.3
        Y[cult == 1] -= 0.4
        for culture in (cult, None):
            f = C.ols_hc3(Y, tgt, culture)
            for j in range(Y.shape[1]):
                if not f["ok"][j]:
                    continue
                co, se, p, df = R.ols_hc3_statsmodels(Y[:, j], tgt, culture)
                assert df == f["df"][j], (df, f["df"][j])
                worst["coef"] = max(worst["coef"], abs(co - f["coef"][j]))
                worst["se"] = max(worst["se"], abs(se - f["se"][j]) / se)
                worst["p"] = max(worst["p"], abs(p - f["p"][j]))
                n_cmp += 1
    return dict(n_compared=n_cmp, max_abs_diff_coef=worst["coef"], max_rel_diff_se=worst["se"], max_abs_diff_p=worst["p"],
                criterion="all < 1e-8, df identical", passed=bool(max(worst.values()) < 1e-8))


def c6(rng):
    worst, n_cmp = 0.0, 0
    for rep in range(10):
        n = int(rng.integers(80, 300))
        Y, st = sim_matrix(rng, n, 20, 3, miss=0.1)
        tg = rng.random(n) < 0.25
        sec = C.secondary_effect_all(Y, tg, st)
        for j in range(Y.shape[1]):
            e = C.m2_effects(Y[:, j], tg, st)
            if np.isfinite(e["secondary"]) or np.isfinite(sec[j]):
                worst = max(worst, abs(e["secondary"] - sec[j]))
                n_cmp += 1
    return dict(n_compared=n_cmp, max_abs_diff=worst, criterion="< 1e-12", passed=bool(worst < 1e-12))


def c7(coh):
    cult = {lev: int((coh["culture"] == k).sum()) for k, lev in enumerate(C.CULTURE_LEVELS)}
    lib = {k: int((coh["library"] == k).sum()) for k in ("Avana", "KY", "Avana+KY", "other", "unknown")}
    exp_c = {"Adherent": 923, "Suspension": 163, "Other": 106}
    exp_l = {"Avana": 822, "KY": 116, "Avana+KY": 195, "other": 6, "unknown": 53}
    return dict(n_lines=int(len(coh["lineage"])), n_genes=int(len(coh["genes"])), culture=cult, library=lib,
                expected_culture=exp_c, expected_library=exp_l,
                passed=bool(cult == exp_c and lib == exp_l and len(coh["lineage"]) == 1192 and len(coh["genes"]) == 18531))


def c8(coh):
    """Model-level coverage (all lines treated as having data): fraction of each lineage's lines in informative
    strata, and lineages untestable at model level."""
    out = {}
    for lin in sorted(set(coh["lineage"])):
        tg = coh["lineage"] == lin
        n_inf = 0
        for k in range(3):
            m = coh["culture"] == k
            if (m & tg).sum() >= 3 and (m & ~tg).sum() >= 3:
                n_inf += int((m & tg).sum())
        out[lin] = dict(n=int(tg.sum()), n_in_informative=n_inf, frac=n_inf / tg.sum(), testable=n_inf >= 5)
    eye = out.get("eye", {})
    amp = out.get("ampulla of vater", {})
    return dict(per_lineage=out, eye_frac=eye.get("frac"), ampulla_testable=amp.get("testable"),
                passed=bool(abs(eye.get("frac", 0) - 0.70) < 1e-9 and amp.get("testable") is False))


def c9(coh, rng_pairs):
    Y = coh["Y"]
    perm = np.random.default_rng(7).permutation(len(coh["lineage"]))
    lab = coh["lineage"][perm]
    lc = {l: int((lab == l).sum()) for l in set(lab)}
    sr = C.StrataRanks(Y, coh["culture"])
    worst_p = worst_m3 = 0.0
    n_testable_total, mism = 0, 0
    pvals = []
    sample = []
    for lin in sorted(lc):
        tg = lab == lin
        v = C.van_elteren(sr, tg)
        n_testable_total += int(v["testable"].sum())
        pvals.append(v["p"][v["testable"]])
        js = np.where(v["testable"])[0]
        if len(js):
            for j in rng_pairs.choice(js, size=min(12, len(js)), replace=False):
                sample.append((lin, int(j), float(v["p"][j])))
    m3_lin = [l for l in sorted(lc) if lc[l] >= C.MIN_TARGET][:3] + ['lymphoid', 'myeloid', 'skin', 'lung']   # M3 is run on the M2 family only
    for lin, j, p in sample[:300]:
        pr, tr, *_ = R.van_elteren_one(Y[:, j], lab == lin, coh["culture"])
        if not tr:
            mism += 1
            continue
        worst_p = max(worst_p, abs(p - pr))
    n_m3 = 0
    for lin in m3_lin:
        tg = lab == lin
        f = C.ols_hc3(Y, tg, coh["culture"])
        for j in rng_pairs.choice(np.where(f["ok"])[0], size=20, replace=False):
            co, se, p, df = R.ols_hc3_statsmodels(Y[:, j], tg, coh["culture"])
            worst_m3 = max(worst_m3, abs(p - f["p"][j]))
            n_m3 += 1
    return dict(permutation_seed=7, n_testable_pairs=n_testable_total, n_pairs_compared=min(300, len(sample)),
                testability_mismatches=mism, max_abs_diff_p_M2=worst_p, n_M3_compared=n_m3, max_abs_diff_p_M3=worst_m3,
                criterion="M2 < 1e-10, M3 < 1e-8, no mismatch",
                passed=bool(worst_p < 1e-10 and worst_m3 < 1e-8 and mism == 0))


def c10():
    from simulation_mixture import truth
    lineage = np.array(["a"] * 4 + ["b"] * 6, dtype=object)
    susp = np.array([1, 1, 0, 0, 1, 0, 0, 0, 0, 0], float)           # a: 2/4 suspension, b: 1/6
    allok = np.ones(10, bool)
    miss = allok.copy(); miss[0] = False                                # a: 1/3 suspension with data
    # (mu, beta, beta_lineage, delta, tested lineage, ok) -> hand-computed (T1, T2)
    cases = [
        ((-0.2, -1.0, "a", 0.0, "a", allok), (True, True)),     # beta -1 <= -0.5, mu + beta = -1.2 < -1
        ((-0.2, -1.0, "a", 0.0, "b", allok), (False, False)),   # gene not selective for b: b = 0
        ((-0.8, 0.0, "", -2.0, "a", allok), (False, True)),     # -2 (1/2 - 1/6) = -0.667 < -0.5; -0.8 - 1 = -1.8 < -1
        ((-0.8, 0.0, "", -2.0, "b", allok), (False, False)),    # -2 (1/6 - 1/2) = +0.667
        ((-0.8, 0.0, "", -2.0, "a", miss), (False, False)),     # -2 (1/3 - 1/6) = -0.333 >= -0.5
        ((-0.4, -0.5, "a", 0.0, "a", allok), (False, False)),   # mu + beta = -0.9 >= -1
    ]
    got = [truth(mu, b, bl, d, lin, lineage, susp, ok) for (mu, b, bl, d, lin, ok), _ in cases]
    exp = [e for _, e in cases]
    return dict(n_cases=len(cases), got=got, expected=exp, passed=bool(got == exp))


def c11(rng):
    n_cmp = mism = n_untestable = 0
    for rep in range(20):
        n = int(rng.integers(30, 120))
        Y, _ = sim_matrix(rng, n, 30, 1, miss=0.15)
        cult = rng.choice(3, size=n, p=[0.75, 0.2, 0.05])
        tg = np.zeros(n, bool)
        tg[rng.choice(n, size=int(rng.choice([1, 2, 3, 6, 15])), replace=False)] = True
        for cu in (cult, None):
            f = C.ols_hc3(Y, tg, cu)
            for j in range(Y.shape[1]):
                ok = ~np.isnan(Y[:, j])
                c = np.zeros(n, int) if cu is None else cu
                cols = [np.ones(n), tg.astype(float)] + [(c == k).astype(float) for k in (1, 2) if (ok & (c == k)).any()]
                X = np.column_stack(cols)[ok]
                full = np.linalg.matrix_rank(X) == X.shape[1] and X.shape[0] > X.shape[1]
                if full:
                    Q, _ = np.linalg.qr(X)
                    full = (Q * Q).sum(1).max() < 1 - 1e-8
                mism += int(bool(full) != bool(f["ok"][j]))
                n_untestable += int(not full)
                n_cmp += 1
    return dict(n_compared=n_cmp, n_untestable_independent=n_untestable, mismatches=mism, passed=bool(mism == 0 and n_untestable > 0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=os.environ.get("DEPMAP_DIR", "."))
    ap.add_argument("--out", default=os.path.join(os.environ.get("RESULTS_DIR", "results"), "prereg_v1b"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    rng = np.random.default_rng(50261119)              # seed of the checks (simulated data only; outside all plan seed ranges)
    res = dict(plan_A_sha256=C.plan_a_sha256(), environment=environment(), seed_checks=50261119)
    for name, fn in (("C1_single_stratum", c1), ("C2_independent_implementation", c2), ("C3_exact_permutation", c3),
                     ("C4_stratified_HL", c4), ("C5_OLS_HC3", c5), ("C6_secondary_effect", c6)):
        res[name] = fn(rng)
        print(name, res[name]["passed"], flush=True)
    coh = C.load_cohort(a.data_dir)
    res["C7_cohort_counts"] = c7(coh); print("C7", res["C7_cohort_counts"]["passed"], flush=True)
    res["C8_coverage_model_level"] = c8(coh); print("C8", res["C8_coverage_model_level"]["passed"], flush=True)
    res["C9_permuted_real_data"] = c9(coh, np.random.default_rng(11)); print("C9", res["C9_permuted_real_data"]["passed"], flush=True)
    res["C10_mixture_truth"] = c10(); print("C10", res["C10_mixture_truth"]["passed"], flush=True)
    res["C11_M3_testability"] = c11(rng); print("C11", res["C11_M3_testability"]["passed"], flush=True)
    res["all_passed"] = all(v["passed"] for k, v in res.items() if k.startswith("C"))
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(os.path.join(a.out, "implementation_checks.json"), "w") as f:
        json.dump(res, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    print("ALL PASSED" if res["all_passed"] else "SOME CHECKS FAILED", res["runtime_s"], "s")


if __name__ == "__main__":
    main()
