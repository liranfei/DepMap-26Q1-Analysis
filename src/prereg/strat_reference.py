"""Separately written reference implementations used only by implementation_checks.py (plan A, Section 3).

Written without reusing any code of strat_core.py: plain per-gene, per-stratum loops, scipy.stats.rankdata,
explicit formulas, Fractions for the weighted median, statsmodels for OLS. Slow but direct.
"""
from fractions import Fraction
import itertools

import numpy as np
import statsmodels.api as sm
from scipy.stats import norm, rankdata


def van_elteren_one(y, target, strata, min_per=3, min_target=5):
    """Van Elteren one-sided p (target lower) for one gene; returns (p, testable, T, E, Var)."""
    y = np.asarray(y, float)
    target = np.asarray(target, bool)
    strata = np.asarray(strata)
    keep = ~np.isnan(y)
    T = E = V = 0.0
    inv_sizes = []
    n_target_informative = 0
    for s in sorted(set(strata.tolist())):
        sel = keep & (strata == s)
        ys, ts = y[sel], target[sel]
        n1, n0 = int(ts.sum()), int((~ts).sum())
        if n1 < min_per or n0 < min_per:
            continue
        n = n1 + n0
        r = rankdata(ys)                              # average ranks
        W = r[ts].sum()
        _, counts = np.unique(ys, return_counts=True)
        tie_sum = sum(int(c) ** 3 - int(c) for c in counts)
        T += W / (n + 1)
        E += n1 / 2
        V += n1 * n0 / (12 * (n + 1) ** 2) * ((n + 1) - tie_sum / (n * (n - 1)))
        inv_sizes.append(1 / (n + 1))
        n_target_informative += n1
    if n_target_informative < min_target or V <= 0:
        return np.nan, False, T, E, V
    c = 0.5 * max(inv_sizes)
    return float(norm.cdf((T - E + c) / np.sqrt(V))), True, T, E, V


def weighted_median_fraction(y, target, strata, min_per=3):
    """Stratified HL by brute force with Fraction weights 1/(n_s + 1)."""
    y = np.asarray(y, float)
    target = np.asarray(target, bool)
    strata = np.asarray(strata)
    keep = ~np.isnan(y)
    items = []
    for s in sorted(set(strata.tolist())):
        sel = keep & (strata == s)
        a, b = y[sel & target], y[sel & ~target]
        if len(a) < min_per or len(b) < min_per:
            continue
        w = Fraction(1, len(a) + len(b) + 1)
        items += [(float(x - z), w) for x in a for z in b]
    if not items:
        return float('nan')
    items.sort(key=lambda t: t[0])
    total = sum(w for _, w in items)
    values = sorted(set(v for v, _ in items))
    acc = Fraction(0)
    byval = {}
    for v, w in items:
        byval[v] = byval.get(v, Fraction(0)) + w
    for i, v in enumerate(values):
        acc += byval[v]
        if acc >= total / 2:
            if acc == total / 2 and i + 1 < len(values):
                return (v + values[i + 1]) / 2
            return v
    raise RuntimeError


def exact_van_elteren_p(y, target, strata):
    """Exact within-stratum permutation p value P(T <= T_obs) by complete enumeration of all target assignments
    within each stratum (all strata are used; for small examples only)."""
    y = np.asarray(y, float)
    target = np.asarray(target, bool)
    strata = np.asarray(strata)
    dists, t_obs = [], Fraction(0)
    for s in sorted(set(strata.tolist())):
        sel = strata == s
        ys, ts = y[sel], target[sel]
        n, n1 = len(ys), int(ts.sum())
        r = rankdata(ys)
        r2 = [int(round(2 * x)) for x in r]           # average ranks are multiples of 1/2
        t_obs += Fraction(sum(x for x, t in zip(r2, ts) if t), 2 * (n + 1))
        d = {}
        for comb in itertools.combinations(range(n), n1):
            key = Fraction(sum(r2[i] for i in comb), 2 * (n + 1))
            d[key] = d.get(key, 0) + 1
        dists.append(d)
    tot = {Fraction(0): 1}
    for d in dists:
        new = {}
        for a, ca in tot.items():
            for b, cb in d.items():
                new[a + b] = new.get(a + b, 0) + ca * cb
        tot = new
    allc = sum(tot.values())
    return sum(c for v, c in tot.items() if v <= t_obs) / allc


def ols_hc3_statsmodels(y, target, culture=None):
    """M3 (culture given) or M3u (culture None) with statsmodels; one-sided p from t with model df."""
    y = np.asarray(y, float)
    keep = ~np.isnan(y)
    cols = [np.ones(keep.sum()), np.asarray(target, float)[keep]]
    if culture is not None:
        c = np.asarray(culture)[keep]
        for lev in (1, 2):
            if (c == lev).any():
                cols.append((c == lev).astype(float))
    X = np.column_stack(cols)
    fit = sm.OLS(y[keep], X).fit(cov_type="HC3", use_t=True)
    coef, se, df = fit.params[1], fit.bse[1], fit.df_resid
    from scipy.stats import t
    return float(coef), float(se), float(t.cdf(coef / se, df)), float(df)
