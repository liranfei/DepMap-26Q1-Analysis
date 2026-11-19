"""Real-data analyses of plan A, Section 5. Run only after v2 is deposited (simulation outputs, recommendation and
plan B screening products), as step 5 of Section 7; outputs and their SHA-256 are deposited as v3.

Guard: with real lineage labels the script refuses to run unless --v2-doi and --decision (sim_decision.json of v2)
are given. With --permute-seed S the lineage labels are permuted (unrestricted) before anything is computed; this is
the only mode used before v2 (code test), and in this mode no result file of the primary analysis is read.

Outputs (Section 5):
  1. M2 for all testable pairs; prioritised pairs; vs the 94 M1 pairs: retained / lost / gained under the recommended
     method, split blood / non-blood; the M1 chain (median diff -> HL -> HL in the M2 family -> M2)
  2. M3 on the M2 family; agreement M2 vs M3
  3. M4 for all genes; and for the genes prioritised in both blood lineages under M1 (q from the all-gene family
     and BH within these genes [N12])
  4. sensitivity analyses: culture x library strata; >= 5 per stratum; Other (Mixed / Unknown / missing) excluded
  5. one model per patient (smallest ModelID kept)
  6. calibration: unrestricted permutation of lineage labels, 1,000 permutations, seeds 42 + b
  7. untestable lineages and flagged pairs
Usage:
  python real_data.py --data-dir $DEPMAP_DIR --primary-results $RESULTS_DIR --out DIR --v2-doi 10.5281/... \
                      --decision v2/sim_decision.json [--n-calibration 1000]
  python real_data.py --data-dir $DEPMAP_DIR --out DIR --permute-seed 7 [--n-calibration 2]      (code test)
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
from env_check import environment  # noqa: E402


def m2_family(Y, genes, lineage, strata, min_per=C.MIN_PER_STRATUM, effects=True, sr=None):
    """Whole M2 procedure: all lineages x genes, testable pairs, BH over the family, effects for q < 0.05,
    prioritisation. Returns a DataFrame with one row per testable pair."""
    sr = sr or C.StrataRanks(Y, strata)
    parts = []
    for lin in np.unique(lineage):
        tg = lineage == lin
        v = C.van_elteren(sr, tg, min_per=min_per)
        js = np.where(v["testable"])[0]
        if not len(js):
            continue
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                med = np.nanmedian(Y[tg][:, js], axis=0)
        d = dict(Gene=genes[js], Lineage=lin, p_M2=v["p"][js], n_target_data=v["n_target_data"][js],
                 n_target_inf=v["n_target_inf"][js], frac_inf=v["frac_inf"][js], abs_all=med, gene_index=js)
        if effects:
            d["secondary_M2"] = C.secondary_effect_all(Y[:, js], tg, strata, min_per=min_per)
        parts.append(pd.DataFrame(d))
    T = pd.concat(parts, ignore_index=True)
    T["q_M2"] = C.bh(T.p_M2.to_numpy())
    T["flag_partial_lineage"] = T.frac_inf < C.FRAC_FLAG
    T["HL_M2"] = np.nan
    T["abs_inf"] = np.nan
    T["heterogeneous"] = pd.Series([pd.NA] * len(T), dtype="boolean")
    sig = np.where(T.q_M2 < C.Q_MAX)[0]
    for i in sig:
        j, lin = int(T.gene_index.iat[i]), T.Lineage.iat[i]
        e = C.m2_effects(Y[:, j], lineage == lin, strata, min_per=min_per)
        T.iat[i, T.columns.get_loc("HL_M2")] = e["hl"]
        T.iat[i, T.columns.get_loc("abs_inf")] = e["abs_inf"]
        T.iat[i, T.columns.get_loc("heterogeneous")] = e["heterogeneous"] if e["n_strata_ge5"] >= 2 else pd.NA
    T["prioritised_M2"] = (T.q_M2 < C.Q_MAX) & (T.abs_all < C.ABS_MAX) & (T.HL_M2 < C.SHIFT_MAX)
    return T


def m3_on_family(Y, lineage, culture, T):
    T = T.copy()
    T["coef_M3"] = np.nan
    T["p_M3"] = np.nan
    for lin, idx in T.groupby("Lineage").groups.items():
        js = T.loc[idx, "gene_index"].to_numpy()
        f = C.ols_hc3(Y[:, js], lineage == lin, culture)
        T.loc[idx, "coef_M3"] = np.where(f["ok"], f["coef"], np.nan)
        T.loc[idx, "p_M3"] = np.where(f["ok"], f["p"], np.nan)
    T["q_M3"] = C.bh(T.p_M3.to_numpy())
    T["prioritised_M3"] = (T.q_M3 < C.Q_MAX) & (T.abs_all < C.ABS_MAX) & (T.coef_M3 < C.SHIFT_MAX)
    return T


def m4_all(Y, genes, lineage, culture):
    blood = np.isin(lineage, C.BLOOD)
    comp = (culture == 1) & ~blood
    r = C.mann_whitney_less(Y, blood, comparator=comp)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mb = np.nanmedian(Y[blood], axis=0)
        mc = np.nanmedian(Y[comp], axis=0)
    D = pd.DataFrame(dict(Gene=genes, p_M4=r["p"], n_blood=r["n_target"], n_comparator=r["n_other"],
                          blood_median=mb, median_difference=mb - mc))
    D = D[r["testable"]].reset_index(drop=True)
    D["q_M4"] = C.bh(D.p_M4.to_numpy())
    D["supported_M4"] = (D.q_M4 < C.Q_MAX) & (D.median_difference < C.SHIFT_MAX) & (D.blood_median < C.ABS_MAX)
    return D, int(comp.sum())


def count_prioritised(Y, genes, lineage, strata, sr):
    T = m2_family(Y, genes, lineage, strata, effects=False, sr=sr)
    return int(T.prioritised_M2.sum()), int((T.q_M2 < C.Q_MAX).sum())


def write_csv(df, path, sha):
    with open(path, "w") as f:
        f.write(f"# plan_A_sha256={sha}\n")
        df.to_csv(f, index=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=os.environ.get("DEPMAP_DIR", "."))
    ap.add_argument("--primary-results", default=os.environ.get("RESULTS_DIR", "results"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--v2-doi")
    ap.add_argument("--decision")
    ap.add_argument("--permute-seed", type=int)
    ap.add_argument("--n-calibration", type=int, default=1000)
    a = ap.parse_args()
    test = a.permute_seed is not None
    if not test:
        assert a.v2_doi and a.decision, "real labels: --v2-doi and --decision (v2) are required (plan A, Section 7)"
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    sha = C.plan_a_sha256()
    coh = C.load_cohort(a.data_dir)
    Y, genes, cult = coh["Y"], coh["genes"], coh["culture"]
    lineage = coh["lineage"].copy()
    if test:
        lineage = lineage[np.random.default_rng(a.permute_seed).permutation(len(lineage))]
    S = dict(plan_A_sha256=sha, environment=environment(), test_mode_permute_seed=a.permute_seed, v2_doi=a.v2_doi)

    # 1-2. M2 and M3
    sr = C.StrataRanks(Y, cult)
    T = m2_family(Y, genes, lineage, cult, sr=sr)
    T = m3_on_family(Y, lineage, cult, T)
    T["blood"] = T.Lineage.isin(C.BLOOD)
    write_csv(T.drop(columns="gene_index"), os.path.join(a.out, "m2_m3_all_pairs.csv"), sha)
    S["M2"] = dict(n_family=len(T), n_q05=int((T.q_M2 < C.Q_MAX).sum()), n_prioritised=int(T.prioritised_M2.sum()),
                   n_prioritised_blood=int((T.prioritised_M2 & T.blood).sum()),
                   n_prioritised_heterogeneous=int((T.prioritised_M2 & (T.heterogeneous == True)).sum()),
                   n_flagged_partial_lineage_pairs=int(T.flag_partial_lineage.sum()))
    S["M3"] = dict(n_family=int(T.p_M3.notna().sum()), n_q05=int((T.q_M3 < C.Q_MAX).sum()),
                   n_prioritised=int(T.prioritised_M3.sum()))
    S["M2_vs_M3"] = dict(both=int((T.prioritised_M2 & T.prioritised_M3).sum()),
                         M2_only=int((T.prioritised_M2 & ~T.prioritised_M3).sum()),
                         M3_only=int((~T.prioritised_M2 & T.prioritised_M3).sum()))
    lins_all = sorted(set(lineage))
    S["untestable_lineages"] = [l for l in lins_all if l not in set(T.Lineage)]

    # 3. M4
    M4, n_comp = m4_all(Y, genes, lineage, cult)
    write_csv(M4, os.path.join(a.out, "m4_all_genes.csv"), sha)
    S["M4"] = dict(n_comparator_lines=n_comp, n_genes=len(M4), n_q05=int((M4.q_M4 < C.Q_MAX).sum()),
                   n_supported=int(M4.supported_M4.sum()))

    # 1 (continued) and 3 (12 genes): comparison with the primary analysis - real labels only
    if not test:
        rec = json.load(open(a.decision))["decision"]["recommended"]
        S["recommended_method"] = rec
        S["decision_file_sha256"] = C.sha256_file(a.decision)
        P = a.primary_results
        A = pd.read_csv(os.path.join(P, "all_tests.csv"))
        F = pd.read_csv(os.path.join(P, "final_targets.csv"))
        H = pd.read_csv(os.path.join(P, "hodges_lehmann_significant_pairs.csv"))
        # the files must be the primary-analysis files described in plan A, Section 5 [N22]
        lin_set = set(lineage)
        for name, D_ in (("final_targets", F), ("all_tests", A), ("hodges_lehmann", H)):
            assert (D_.Lineage == D_.Lineage.str.lower()).all(), f"{name}: lineage names not lower case [N1]"
            assert set(D_.Lineage) <= lin_set, f"{name}: lineage not in the cohort [N1]"
        k94 = set(zip(F.Gene, F.Lineage))
        assert len(k94) == 94, len(k94)
        assert k94 <= set(zip(A.Gene, A.Lineage)), "a final target is missing from all_tests.csv"
        assert n_comp == 41, n_comp
        col = "prioritised_M2" if rec == "M2" else "prioritised_M3"
        kre = set(zip(T.Gene[T[col]], T.Lineage[T[col]]))
        def split(s):
            return dict(total=len(s), blood=sum(l in C.BLOOD for _, l in s), non_blood=sum(l not in C.BLOOD for _, l in s))
        S["vs_M1_94"] = dict(method=rec, retained=split(k94 & kre), lost=split(k94 - kre), gained=split(kre - k94),
                             lost_pairs=sorted(map(list, k94 - kre)), gained_pairs=sorted(map(list, kre - k94)))
        # chain: M1 median diff -> M1 HL -> M1 HL in the M2 family (BH recomputed) -> M2
        hl = H.set_index(["Gene", "Lineage"]).Selectivity_HL
        c_hl = set(H[(H.Chronos_median < C.ABS_MAX) & (H.Selectivity_HL < C.SHIFT_MAX)][["Gene", "Lineage"]].itertuples(index=False, name=None))
        assert len(c_hl) == 68, len(c_hl)
        fam = set(zip(T.Gene, T.Lineage))
        Ar = A[[k in fam for k in zip(A.Gene, A.Lineage)]].copy()
        Ar["q_restricted"] = C.bh(Ar.p_value.to_numpy())
        sig = Ar[Ar.q_restricted < C.Q_MAX]
        gi = {g: j for j, g in enumerate(genes)}
        vals = []
        for g, l, m in zip(sig.Gene, sig.Lineage, sig.Chronos_median):
            if (g, l) in hl.index:
                h = float(hl.loc[(g, l)])
            else:                                            # newly significant in the smaller family: exact HL
                y = Y[:, gi[g]]
                ok = ~np.isnan(y)
                h = C.hl_shift(y[ok & (lineage == l)], y[ok & (lineage != l)])
            vals.append(((g, l), m, h))
        c_fam = {k for k, m, h in vals if m < C.ABS_MAX and h < C.SHIFT_MAX}
        c_m2 = set(zip(T.Gene[T.prioritised_M2], T.Lineage[T.prioritised_M2]))
        S["M1_chain"] = dict(M1_median_diff=split(k94), M1_HL=split(c_hl), M1_HL_in_M2_family=split(c_fam),
                             M2=split(c_m2), n_pairs_in_M2_family_from_all_tests=len(Ar),
                             n_newly_significant_HL_computed=sum((k not in hl.index) for k, _, _ in vals))
        both = sorted(set(F.Gene[F.Lineage == "lymphoid"]) & set(F.Gene[F.Lineage == "myeloid"]))
        assert len(both) == 12, len(both)
        m12 = M4[M4.Gene.isin(both)].copy()
        m12["q_within_these_genes"] = C.bh(m12.p_M4.to_numpy())
        write_csv(m12, os.path.join(a.out, "m4_shared_blood_genes.csv"), sha)
        S["M4_shared_blood_genes"] = dict(n_genes=len(both), n_tested=len(m12),
                                          n_supported_allgene_q=int(m12.supported_M4.sum()))

    # 4. sensitivity analyses
    strata_lib = np.array([f"{c}|{l}" for c, l in zip(cult, coh["library"])], dtype=object)
    T_lib = m2_family(Y, genes, lineage, strata_lib, effects=False)
    T_5 = m2_family(Y, genes, lineage, cult, min_per=5, effects=False, sr=sr)
    keep = cult != 2
    T_ex = m2_family(Y[keep], genes, lineage[keep], cult[keep], effects=False)
    main_set = set(zip(T.Gene[T.prioritised_M2], T.Lineage[T.prioritised_M2]))
    sens = {}
    for name, X in (("culture_x_library", T_lib), ("min5_per_stratum", T_5), ("other_excluded", T_ex)):
        s_ = set(zip(X.Gene[X.prioritised_M2], X.Lineage[X.prioritised_M2]))
        sens[name] = dict(n_family=len(X), n_prioritised=len(s_), overlap_main=len(s_ & main_set),
                          only_here=len(s_ - main_set), only_main=len(main_set - s_))
        write_csv(X.drop(columns="gene_index"), os.path.join(a.out, f"m2_sensitivity_{name}.csv"), sha)
    # lines in non-informative strata at model level (culture x library), per lineage
    dropped = {}
    for lin in lins_all:
        tg = lineage == lin
        n_drop = 0
        for lev in np.unique(strata_lib):
            m = strata_lib == lev
            if not ((m & tg).sum() >= 3 and (m & ~tg).sum() >= 3):
                n_drop += int((m & tg).sum())
        dropped[lin] = n_drop
    sens["culture_x_library_target_lines_dropped_model_level"] = dropped
    S["sensitivity"] = sens

    # 5. one model per patient (smallest ModelID kept)
    order = np.argsort(coh["model_ids"].astype(str), kind="mergesort")
    seen, keep1 = set(), np.zeros(len(lineage), bool)
    for i in order:
        if coh["patient"][i] not in seen:
            seen.add(coh["patient"][i])
            keep1[i] = True
    T_p = m2_family(Y[keep1], genes, lineage[keep1], cult[keep1], effects=False)
    s_p = set(zip(T_p.Gene[T_p.prioritised_M2], T_p.Lineage[T_p.prioritised_M2]))
    S["one_model_per_patient"] = dict(n_models=int(keep1.sum()), n_prioritised=len(s_p), overlap_main=len(s_p & main_set))

    # 6. calibration: unrestricted permutations, seeds 42 + b
    cal = []
    for b in range(a.n_calibration):
        lab = lineage[np.random.default_rng(42 + b).permutation(len(lineage))]
        cal.append(count_prioritised(Y, genes, lab, cult, sr))
    if cal:
        x = np.array([c[0] for c in cal])
        S["calibration"] = dict(n_permutations=len(cal), mean_prioritised=float(x.mean()), max_prioritised=int(x.max()),
                                mean_q05=float(np.mean([c[1] for c in cal])), expectation="mean <= 1",
                                exceeded=bool(x.mean() > 1))
        pd.DataFrame(cal, columns=["n_prioritised", "n_q05"]).to_csv(os.path.join(a.out, "m2_calibration_permutations.csv"), index=False)

    # 7. flagged pairs
    write_csv(T[T.flag_partial_lineage].drop(columns="gene_index"), os.path.join(a.out, "m2_flagged_partial_lineage.csv"), sha)
    S["runtime_s"] = round(time.time() - t0, 1)
    json.dump(S, open(os.path.join(a.out, "real_data_summary.json"), "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in S.items() if k not in ("environment",)}, indent=1, default=str)[:3000])


if __name__ == "__main__":
    main()
