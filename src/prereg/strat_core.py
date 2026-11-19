"""Core methods of preregistered plan A (preregistration/plan_A_stratified_procedure.md, Zenodo v1 = v2.5.0-prereg).

M1  one-sided Mann-Whitney, target lineage vs all other cancer lines (vectorised; equals scipy asymptotic)
M2  one-sided Van Elteren test over informative culture-format strata, stratified HL shift (exact rational weights),
    secondary effect size, heterogeneity flag
M3  OLS gene effect ~ target + culture-format indicators, HC3, one-sided t with n - p df (p = 4 when all three formats occur)
M3u the same OLS without culture indicators (simulation only)
M4  blood (Lymphoid + Myeloid) vs non-blood suspension lines, one-sided Mann-Whitney

All functions work on a matrix Y (lines x genes, NaN = missing gene effect). Nothing here looks at any real lineage
test result; the loaders only read gene effects and annotations. Implementation choices that the plan does not fix
are listed in IMPLEMENTATION_NOTES.md and referred to as [N1], [N2], ...
"""
from __future__ import annotations

import hashlib
import math
import os
from functools import reduce

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata, t as tdist

Q_MAX, ABS_MAX, SHIFT_MAX = 0.05, -1.0, -0.5      # prioritisation thresholds (plan A, Section 3)
MIN_TARGET = 5                                     # MIN_N of the primary analysis
MIN_PER_STRATUM = 3                                # informative stratum: >= 3 target and >= 3 other lines (sensitivity: 5)
FRAC_FLAG = 0.80                                   # coverage flag
CULTURE_LEVELS = ("Adherent", "Suspension", "Other")
BLOOD = ("lymphoid", "myeloid")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PLAN_A = os.path.join(REPO, "preregistration", "plan_A_stratified_procedure.md")


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def arrays_sha256(npz_path, keys=("mu", "sigma", "mask", "lineage", "culture", "growth_pattern")):
    """SHA-256 of the arrays of sim_inputs.npz (reproducible, unlike the bytes of a compressed .npz) [N20]."""
    z = np.load(npz_path, allow_pickle=True)
    h = hashlib.sha256()
    for k in keys:
        a = z[k]
        h.update(k.encode())
        h.update(a.astype(str).tobytes() if a.dtype == object else np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def plan_a_sha256():
    """SHA-256 of plan A; asserts that it equals the value frozen in SHA256SUMS.txt (v1)."""
    sums = os.path.join(REPO, "preregistration", "SHA256SUMS.txt")
    frozen = [l.split()[0] for l in open(sums) if l.strip().endswith("plan_A_stratified_procedure.md")][0]
    actual = sha256_file(PLAN_A)
    assert actual == frozen, "plan A differs from the frozen version"
    return actual


# ----------------------------------------------------------------------------------------------------------------
# data
# ----------------------------------------------------------------------------------------------------------------
EXPECTED_SHA256 = {
    "gene_effect.csv": "e610a4cefb13a82b5b256b47eb08b63ff14843f8dbd0fb164bc0a32688e5b89e",
    "Model.csv": "ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8",
    "ScreenSequenceMap.csv": "cd053d229750d2e35abb1dcd71410ab336f3f9eef122a8168108cb57f98359d3",
}


def culture3(growth_pattern):
    """3-level culture format: Adherent, Suspension, Other (Mixed, Unknown, missing and any other value)."""
    gp = pd.Series(growth_pattern, dtype=object)
    return np.where(gp == "Adherent", 0, np.where(gp == "Suspension", 1, 2)).astype(np.int64)


def library_level(model_ids, ssm):
    """Set of libraries among the model's 2DS screens that pass QC and are not drug-treated, engineered or excluded
    from the combined dataset (screen filter of genome_perm_strat.py) -> Avana / KY / Avana+KY / other / unknown."""
    s = ssm[(ssm.ScreenType == "2DS") & (ssm.PassesQC == True) & (ssm.DrugTreated == False)
            & (ssm.IsEngineered == False) & (ssm.ExcludeFromCRISPRCombined == False)]
    sets = s.groupby("ModelID").Library.agg(lambda x: frozenset(x.dropna()))
    out = []
    for m in model_ids:
        st = sets.get(m, frozenset())
        if not st:
            out.append("unknown")
        elif st & {"TKOv3", "Brunello", "Humagne-CD"}:
            out.append("other")
        elif st == {"Avana"}:
            out.append("Avana")
        elif st == {"KY"}:
            out.append("KY")
        elif st == {"Avana", "KY"}:
            out.append("Avana+KY")
        else:                                           # any other library name would be a new case
            raise ValueError(f"unexpected library set for {m}: {sorted(st)}")
    return np.array(out, dtype=object)


def load_cohort(data_dir, genes=None, check_sha=True):
    """Cancer models of DepMap 26Q1 (Non-Cancerous excluded), all genes, lineage = OncotreeLineage (lower case, as in
    the primary analysis), 3-level culture format, library level, patient. Reads annotations and gene effects only."""
    if check_sha:
        for f, h in EXPECTED_SHA256.items():
            assert sha256_file(os.path.join(data_dir, f)) == h, f"{f}: SHA-256 differs from plan A"
    usecols = None if genes is None else (lambda c: c == "" or c.startswith("Unnamed") or c in set(genes))
    ge = pd.read_csv(os.path.join(data_dir, "gene_effect.csv"), usecols=usecols)
    ge = ge.rename(columns={ge.columns[0]: "ModelID"})
    ge["ModelID"] = ge.ModelID.astype(str).str.strip()
    M = pd.read_csv(os.path.join(data_dir, "Model.csv"))
    M["ModelID"] = M.ModelID.astype(str).str.strip()
    meta = M[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "GrowthPattern", "PatientID"]]
    df = ge.merge(meta, on="ModelID", how="left", validate="one_to_one")
    assert len(df) == len(ge) and df.OncotreeLineage.notna().all()
    df = df[df.OncotreePrimaryDisease != "Non-Cancerous"].reset_index(drop=True)
    gcols = [c for c in ge.columns if c != "ModelID"]
    ssm = pd.read_csv(os.path.join(data_dir, "ScreenSequenceMap.csv"))
    return dict(
        Y=df[gcols].apply(pd.to_numeric, errors="coerce").to_numpy(np.float64),
        genes=np.array(gcols, dtype=object),
        model_ids=df.ModelID.to_numpy(object),
        lineage=df.OncotreeLineage.str.strip().str.lower().to_numpy(object),
        culture=culture3(df.GrowthPattern.to_numpy(object)),
        growth_pattern=df.GrowthPattern.fillna("missing").to_numpy(object),
        library=library_level(df.ModelID, ssm),
        patient=df.PatientID.to_numpy(object),
    )


# ----------------------------------------------------------------------------------------------------------------
# ranks within strata
# ----------------------------------------------------------------------------------------------------------------
def _tie_term(sub):
    """sum over tie groups of (t^3 - t), per column of sub (rows = lines), NaN ignored."""
    S = np.sort(sub, axis=0)                           # NaN sorted last; NaN == NaN is False
    eq = S[1:] == S[:-1]
    out = np.zeros(sub.shape[1])
    for j in np.where(eq.any(axis=0))[0]:
        v = sub[:, j]
        _, c = np.unique(v[~np.isnan(v)], return_counts=True)
        out[j] = float(((c.astype(np.float64) ** 3) - c).sum())
    return out


class StrataRanks:
    """Average ranks of every gene within each stratum (lines with data only), stratum sizes and tie terms.
    Ranks do not depend on lineage labels, so they are computed once per (data, strata)."""

    def __init__(self, Y, strata):
        strata = np.asarray(strata)
        self.levels = np.unique(strata)
        self.n_lines, self.n_genes = Y.shape
        self.ok = ~np.isnan(Y)
        self.Rz = np.zeros_like(Y)
        S = len(self.levels)
        self.ind = np.zeros((S, self.n_lines))         # stratum indicator, strata x lines
        self.N = np.zeros((S, self.n_genes))
        self.tie = np.zeros((S, self.n_genes))
        for k, lev in enumerate(self.levels):
            idx = np.where(strata == lev)[0]
            self.ind[k, idx] = 1.0
            sub = Y[idx]
            r = rankdata(sub, axis=0, nan_policy="omit")
            self.Rz[idx] = np.where(np.isnan(sub), 0.0, r)
            self.N[k] = (~np.isnan(sub)).sum(axis=0)
            self.tie[k] = _tie_term(sub)
        self.okf = self.ok.astype(np.float64)


def van_elteren(sr: StrataRanks, target, min_per=MIN_PER_STRATUM, min_target=MIN_TARGET, return_parts=False):
    """M2 for one target-lineage indicator over all genes.
    p = Phi((T - E[T] + c) / sqrt(Var[T])), c = 0.5 * max over informative strata of 1/(n_s + 1)."""
    tgt = np.asarray(target, dtype=np.float64)
    A = sr.ind * tgt[None, :]                          # strata x lines, target lines of each stratum
    nt = A @ sr.okf                                    # strata x genes
    W = A @ sr.Rz
    N = sr.N
    no = N - nt
    inf = (nt >= min_per) & (no >= min_per)
    n_t_inf = np.where(inf, nt, 0).sum(0)
    n_t_all = nt.sum(0)
    testable = n_t_inf >= min_target
    with np.errstate(divide="ignore", invalid="ignore"):
        Np1 = N + 1.0
        T = np.where(inf, W / Np1, 0).sum(0)
        E = np.where(inf, nt / 2.0, 0).sum(0)
        V = np.where(inf, nt * no / (12.0 * Np1 ** 2) * (Np1 - sr.tie / (N * (N - 1.0))), 0).sum(0)
        c = 0.5 * np.where(inf, 1.0 / Np1, 0).max(0)
        z = (T - E + c) / np.sqrt(V)
        p = norm.cdf(z)
    bad = testable & ~(V > 0)                          # all values tied in every informative stratum [N3]
    p = np.where(testable & (V > 0), p, np.nan)
    out = dict(p=p, testable=testable & (V > 0), n_target_inf=n_t_inf, n_target_data=n_t_all,
               frac_inf=np.divide(n_t_inf, n_t_all, out=np.zeros_like(n_t_inf), where=n_t_all > 0),
               n_zero_variance=int(bad.sum()))
    if return_parts:
        out.update(T=T, E=E, V=V, c=c, z=z, inf=inf, nt=nt, no=no, N=N)
    return out


def mann_whitney_less(Y, target, comparator=None, min_target=MIN_TARGET):
    """One-sided Mann-Whitney (target lower), asymptotic with tie and continuity correction, all genes at once.
    comparator: boolean mask of comparison lines (default: all non-target lines). Equals scipy mannwhitneyu."""
    tgt = np.asarray(target, bool)
    cmp_ = ~tgt if comparator is None else np.asarray(comparator, bool) & ~tgt
    use = tgt | cmp_
    st = np.where(use, 0, -1)
    sr = StrataRanks(np.where(use[:, None], Y, np.nan), st)
    k = int(np.where(sr.levels == 0)[0][0])
    sub = dict(ind=sr.ind[[k]], okf=sr.okf, Rz=sr.Rz, N=sr.N[[k]], tie=sr.tie[[k]])
    A = sub["ind"] * tgt[None, :]
    nt = (A @ sub["okf"])[0]
    W = (A @ sub["Rz"])[0]
    N = sub["N"][0]
    no = N - nt
    with np.errstate(divide="ignore", invalid="ignore"):
        U1 = W - nt * (nt + 1) / 2.0
        mu = nt * no / 2.0
        s = np.sqrt(nt * no / 12.0 * ((N + 1) - sub["tie"][0] / (N * (N - 1))))
        p = norm.cdf((U1 - mu + 0.5) / s)
    testable = (nt >= min_target) & (no >= 1) & (s > 0)
    return dict(p=np.where(testable, p, np.nan), testable=testable, n_target=nt, n_other=no)


# ----------------------------------------------------------------------------------------------------------------
# effect sizes
# ----------------------------------------------------------------------------------------------------------------
def hl_shift(a, b):
    """Ordinary Hodges-Lehmann shift: median of all a_i - b_j (as hodges_lehmann.py)."""
    return float(np.median(np.asarray(a)[:, None] - np.asarray(b)[None, :]))


def stratified_hl(diff_by_stratum, n_by_stratum):
    """Weighted median of within-stratum pairwise differences, weight 1/(n_s + 1), with exact integer weights
    L/(n_s + 1), L = lcm of the (n_s + 1). Smallest difference at which the cumulative weight first reaches half the
    total; if it equals half exactly, the mean of that difference and the next distinct difference.
    diff_by_stratum: list of 1-d arrays (target - other differences of each informative stratum)
    n_by_stratum: list of n_s (lines with data in the stratum)."""
    L = reduce(math.lcm, [int(n) + 1 for n in n_by_stratum])
    w = [L // (int(n) + 1) for n in n_by_stratum]       # Python ints, exact
    d = np.concatenate([np.asarray(x, np.float64) for x in diff_by_stratum])
    sid = np.concatenate([np.full(len(x), k) for k, x in enumerate(diff_by_stratum)])
    o = np.argsort(d, kind="mergesort")
    d, sid = d[o], sid[o]
    ends = np.r_[np.where(d[1:] != d[:-1])[0], len(d) - 1]     # last index of each distinct value
    cnt = np.stack([np.cumsum(sid == k) for k in range(len(w))])  # strata x n, int64 counts
    total2 = sum(wk * len(x) for wk, x in zip(w, diff_by_stratum))   # total weight (we compare 2*cum with total)

    def cum2(j):                                       # 2 * cumulative weight up to distinct value j (exact)
        e = ends[j]
        return 2 * sum(wk * int(cnt[k, e]) for k, wk in enumerate(w))

    lo, hi = 0, len(ends) - 1
    while lo < hi:                                     # smallest j with 2*cum >= total
        mid = (lo + hi) // 2
        if cum2(mid) >= total2:
            hi = mid
        else:
            lo = mid + 1
    v = d[ends[lo]]
    if cum2(lo) == total2 and lo + 1 < len(ends):
        return float((v + d[ends[lo + 1]]) / 2.0)
    return float(v)


def m2_effects(y, target, strata, min_per=MIN_PER_STRATUM):
    """For one gene: stratified HL shift, secondary effect (weighted mean of within-stratum median differences,
    weights n_t n_o/(n_s + 1)), within-stratum HL shifts of strata with >= 5 target lines (heterogeneity check),
    absolute effects. y: gene effects of all lines (NaN allowed)."""
    y = np.asarray(y, np.float64)
    tgt = np.asarray(target, bool)
    ok = ~np.isnan(y)
    diffs, ns, sec_num, sec_den, within = [], [], 0.0, 0.0, {}
    for lev in np.unique(strata):
        m = (np.asarray(strata) == lev) & ok
        a, b = y[m & tgt], y[m & ~tgt]
        if len(a) >= min_per and len(b) >= min_per:
            n = len(a) + len(b)
            diffs.append((a[:, None] - b[None, :]).ravel())
            ns.append(n)
            wgt = len(a) * len(b) / (n + 1.0)
            sec_num += wgt * (np.median(a) - np.median(b))
            sec_den += wgt
            if len(a) >= 5:
                within[lev] = hl_shift(a, b)
    st = np.asarray(strata)
    inf_levels = [lev for lev in np.unique(st)
                  if ((st == lev) & ok & tgt).sum() >= min_per and ((st == lev) & ok & ~tgt).sum() >= min_per]
    inf_t = y[np.isin(st, inf_levels) & ok & tgt]
    hl = stratified_hl(diffs, ns) if diffs else np.nan
    signs = {np.sign(v) for v in within.values() if v != 0}
    return dict(hl=hl, secondary=sec_num / sec_den if sec_den > 0 else np.nan,
                abs_all=float(np.median(y[ok & tgt])) if (ok & tgt).any() else np.nan,
                abs_inf=float(np.median(inf_t)) if len(inf_t) else np.nan,
                n_strata_ge5=len(within), heterogeneous=bool(len(within) >= 2 and len(signs) > 1),
                within_hl=within)


def secondary_effect_all(Y, target, strata, min_per=MIN_PER_STRATUM):
    """Secondary effect size of M2 for all genes (vectorised): weighted mean of within-stratum median differences,
    weights n_t n_o/(n_s + 1), informative strata only."""
    tgt = np.asarray(target, bool)
    num = np.zeros(Y.shape[1])
    den = np.zeros(Y.shape[1])
    for lev in np.unique(strata):
        m = np.asarray(strata) == lev
        A, B = Y[m & tgt], Y[m & ~tgt]
        na, nb = (~np.isnan(A)).sum(0), (~np.isnan(B)).sum(0)
        inf = (na >= min_per) & (nb >= min_per)
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                md = np.nanmedian(A, axis=0) - np.nanmedian(B, axis=0) if len(A) and len(B) else np.zeros(Y.shape[1])
        w = np.where(inf, na * nb / (na + nb + 1.0), 0.0)
        num += np.where(inf, w * np.nan_to_num(md), 0.0)
        den += w
    return np.divide(num, den, out=np.full(Y.shape[1], np.nan), where=den > 0)


# ----------------------------------------------------------------------------------------------------------------
# M3 / M3u: OLS with HC3, vectorised over genes through cell sufficient statistics
# ----------------------------------------------------------------------------------------------------------------
def ols_hc3(Y, target, culture=None, cells=None):
    """OLS of each gene's effect on [1, target, (suspension, other)] with HC3 SEs, one-sided p for the target
    coefficient < 0 from t with n - p df. With culture=None this is M3u ([1, target]).
    Exact and vectorised: the design has only 6 distinct rows (target x culture), so X'X, X'y, residual sums and
    leverages are computed from per-cell counts, sums and sums of squares. A culture level without any line with
    data for a gene is dropped for that gene and df = n - (number of remaining columns) [N4]."""
    tgt = np.asarray(target, bool)
    G = Y.shape[1]
    cult = np.zeros(len(tgt), np.int64) if culture is None else np.asarray(culture)
    ok = ~np.isnan(Y)
    Y0 = np.where(ok, Y, 0.0)
    cell_ids = [(t_, c_) for t_ in (0, 1) for c_ in (0, 1, 2)]
    n_c, s_c, ss_c = [], [], []
    for t_, c_ in cell_ids:
        m = ((tgt == bool(t_)) & (cult == c_)).astype(np.float64)
        n_c.append(m @ ok)
        s_c.append(m @ Y0)
        ss_c.append(m @ (Y0 * Y0))
    n_c, s_c, ss_c = np.stack(n_c, 1), np.stack(s_c, 1), np.stack(ss_c, 1)          # genes x 6
    # full design rows: [1, target, susp, other]
    Xc = np.array([[1, t_, int(c_ == 1), int(c_ == 2)] for t_, c_ in cell_ids], np.float64)  # 6 x 4
    has_s = (n_c[:, [1, 4]].sum(1) > 0)
    has_o = (n_c[:, [2, 5]].sum(1) > 0)
    if culture is None:
        has_s[:] = False
        has_o[:] = False
    coef = np.full(G, np.nan)
    se = np.full(G, np.nan)
    df = np.full(G, np.nan)
    n_tot = n_c.sum(1)
    for hs in (False, True):
        for ho in (False, True):
            g = np.where((has_s == hs) & (has_o == ho))[0]
            if not len(g):
                continue
            cols = [0, 1] + ([2] if hs else []) + ([3] if ho else [])
            X = Xc[:, cols]                                                     # 6 x p
            p_ = len(cols)
            XtX = np.einsum("gc,ci,cj->gij", n_c[g], X, X)
            Xty = np.einsum("gc,ci->gi", s_c[g], X)
            det_ok = np.abs(np.linalg.det(XtX)) > 1e-9
            gi = g[det_ok]
            if not len(gi):
                continue
            XtX, Xty = XtX[det_ok], Xty[det_ok]
            inv = np.linalg.inv(XtX)
            b = np.einsum("gij,gj->gi", inv, Xty)
            f = np.einsum("ci,gi->gc", X, b)                                    # fitted value of each cell
            rss_c = ss_c[gi] - 2 * f * s_c[gi] + n_c[gi] * f * f                # residual sum of squares per cell
            rss_c = np.maximum(rss_c, 0.0)
            h = np.einsum("ci,gij,cj->gc", X, inv, X)                           # leverage of each cell's rows
            lev_one = ((n_c[gi] > 0) & (1.0 - h < 1e-8)).any(1)                  # a cell with leverage 1 [N4]
            with np.errstate(divide="ignore", invalid="ignore"):
                wgt = np.where(n_c[gi] > 0, rss_c / (1.0 - h) ** 2, 0.0)
            meat = np.einsum("gc,ci,cj->gij", wgt, X, X)
            cov = inv @ meat @ inv
            keep = ~lev_one
            gi, b, cov = gi[keep], b[keep], cov[keep]
            coef[gi] = b[:, 1]
            se[gi] = np.sqrt(cov[:, 1, 1])
            df[gi] = n_tot[gi] - p_
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat = coef / se
        p = tdist.cdf(tstat, df)
    good = np.isfinite(p) & (df > 0)
    return dict(coef=coef, se=se, df=df, p=np.where(good, p, np.nan), ok=good, n=n_tot)


# ----------------------------------------------------------------------------------------------------------------
# multiple testing and intervals
# ----------------------------------------------------------------------------------------------------------------
def bh(p):
    """Benjamini-Hochberg q values (as statsmodels fdr_bh), NaN kept as NaN."""
    p = np.asarray(p, np.float64)
    q = np.full(p.shape, np.nan)
    v = np.isfinite(p)
    pv = p[v]
    m = len(pv)
    if m:
        o = np.argsort(pv, kind="mergesort")
        r = pv[o] * m / (np.arange(m) + 1.0)
        qq = np.minimum.accumulate(r[::-1])[::-1]
        out = np.empty(m)
        out[o] = np.minimum(qq, 1.0)
        q[v] = out
    return q


def clopper_pearson(k, n, alpha=0.05):
    from scipy.stats import beta
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return lo, hi
