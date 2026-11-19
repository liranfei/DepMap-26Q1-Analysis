"""Simulation of plan A, Sections 4.1-4.3 (cells of the grid). Run only after v1b is deposited.

ENUMERATION OF CELLS (fixed at v1b; cell number i is 0-based and uses seed 20261119 + i)
    for scenario in SCENARIOS:                      # 9 scenarios, order below
      for c_beta in (0, -0.1, -0.2, -0.35, -0.5):
        for c_delta in (0, -0.25, -0.5, -1):
          for noise in ("t5", "lognormal"):
            for k in (1, 1.5):
                i += 1
    -> 9 x 5 x 4 x 2 x 2 = 720 cells; the full table is written to sim_cells_enumeration.csv by --enumerate.

SCENARIOS
    single:lymphoid, single:myeloid, single:skin, single:lung        target lineage gets beta; tests for that lineage
    shared_blood                                                     beta (and k) for lymphoid AND myeloid lines;
                                                                     M1-M3 run for lymphoid and for myeloid; M4 run
    misspec:lymphoid, misspec:myeloid, misspec:skin, misspec:lung    as single:, but suspension lines of Lung get
                                                                     delta/2 and Mixed lines get delta/2 (reported only)

GENERATING MODEL (Section 4.1): y = mu_g + beta * target + delta * d_line + sigma_g * k_line * eps
    (mu_g, sigma_g, missing mask) sampled jointly with replacement from the real genes (sim_inputs.npz);
    beta = c_beta * sigma_g, delta = c_delta * sigma_g; d_line = 1 for suspension lines, 0 otherwise (misspec: see
    above); eps = t5 / sqrt(5/3), or log-normal(sdlog 0.5) standardised to mean 0, SD 1 and sign-flipped;
    k_line = k for target lines, 1 otherwise. Genes per cell: 20,000 if c_beta = 0, else 2,000, generated in blocks
    of 2,000 genes in a fixed order (block size is part of the random-number stream and is fixed here).

METHODS (strat_core.py): M1, M2, M3 (on the M2-testable genes, as the M2 family [N5]), M3u (on the M1-testable genes
    [N5]), M4 (shared_blood only). Rejection = p < 0.05 (one-sided, unadjusted), among testable genes [N6].
    Effect-size bias on the first 200 genes of each cell [N7].

Usage:
    python simulation.py --enumerate --out DIR                     # writes the enumeration only
    python simulation.py --inputs DIR/sim_inputs.npz --out DIR [--cells 0-719] [--workers 4]
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

SCENARIOS = ["single:lymphoid", "single:myeloid", "single:skin", "single:lung", "shared_blood",
             "misspec:lymphoid", "misspec:myeloid", "misspec:skin", "misspec:lung"]
C_BETA = (0.0, -0.1, -0.2, -0.35, -0.5)
C_DELTA = (0.0, -0.25, -0.5, -1.0)
NOISE = ("t5", "lognormal")
K = (1.0, 1.5)
SEED0 = 20261119
BLOCK = 2000
N_BIAS = 200
ALPHA = 0.05


def enumerate_cells():
    rows, i = [], 0
    for sc in SCENARIOS:
        for cb in C_BETA:
            for cd in C_DELTA:
                for nz in NOISE:
                    for k in K:
                        rows.append(dict(cell=i, seed=SEED0 + i, scenario=sc, c_beta=cb, c_delta=cd, noise=nz, k=k,
                                         n_genes=20000 if cb == 0 else 2000))
                        i += 1
    return pd.DataFrame(rows)


def noise(rng, kind, shape):
    if kind == "t5":
        return rng.standard_t(5, size=shape) / np.sqrt(5.0 / 3.0)
    x = np.exp(0.5 * rng.standard_normal(size=shape))
    m, v = np.exp(0.125), (np.exp(0.25) - 1.0) * np.exp(0.25)
    return -(x - m) / np.sqrt(v)


def scenario_setup(sc, lineage, culture, growth_pattern):
    """target lines (receive beta and k), tested lineages, per-line multiplier of delta."""
    kind, _, lin = sc.partition(":")
    susp = (culture == 1).astype(np.float64)
    d = susp.copy()
    if kind == "shared_blood":
        target = np.isin(lineage, C.BLOOD)
        tested = list(C.BLOOD)
    else:
        target = lineage == lin
        tested = [lin]
    if kind == "misspec":
        d[(lineage == "lung") & (culture == 1)] = 0.5
        d[growth_pattern == "Mixed"] = 0.5
    return target, tested, d


def run_cell(cell, inp):
    t0 = time.time()
    rng = np.random.default_rng(int(cell["seed"]))
    lineage, culture, gp = inp["lineage"], inp["culture"], inp["growth_pattern"]
    target, tested, dmult = scenario_setup(cell["scenario"], lineage, culture, gp)
    kline = np.where(target, cell["k"], 1.0)
    G = int(cell["n_genes"])
    acc = {}                                                     # (method, lineage) -> [n_testable, n_reject]
    bias = {}                                                    # (estimator, lineage) -> list of (est - beta)/sigma
    m4_comp = (culture == 1) & ~np.isin(lineage, C.BLOOD)
    for b0 in range(0, G, BLOCK):
        g = min(BLOCK, G - b0)
        idx = rng.integers(0, inp["mu"].shape[0], size=g)
        mu, sd, mask = inp["mu"][idx], inp["sigma"][idx], inp["mask"][:, idx]
        eps = noise(rng, cell["noise"], (len(lineage), g))
        beta = cell["c_beta"] * sd
        delta = cell["c_delta"] * sd
        Y = mu[None, :] + target[:, None] * beta[None, :] + dmult[:, None] * delta[None, :] + kline[:, None] * sd[None, :] * eps
        Y[~mask] = np.nan
        sr = C.StrataRanks(Y, culture)
        for lin in tested:
            tg = lineage == lin
            m1 = C.mann_whitney_less(Y, tg)
            m2 = C.van_elteren(sr, tg)
            m3 = C.ols_hc3(Y, tg, culture)
            m3u = C.ols_hc3(Y, tg, None)
            res = {"M1": (m1["testable"], m1["p"]), "M2": (m2["testable"], m2["p"]),
                   "M3": (m2["testable"] & m3["ok"], m3["p"]), "M3u": (m1["testable"] & m3u["ok"], m3u["p"])}
            for meth, (tst, p) in res.items():
                a = acc.setdefault((meth, lin), [0, 0])
                a[0] += int(tst.sum())
                a[1] += int((tst & (p < ALPHA)).sum())
            if b0 == 0:                                          # effect-size bias, first N_BIAS genes [N7]
                for j in range(min(N_BIAS, g)):
                    y = Y[:, j]
                    ok = ~np.isnan(y)
                    a_, o_ = y[ok & tg], y[ok & ~tg]
                    if not m1["testable"][j]:                   # bias only for genes testable under M1 [N7]
                        continue
                    e = C.m2_effects(y, tg, culture)
                    ests = {"M1_median_diff": float(np.median(a_) - np.median(o_)), "M1_HL": C.hl_shift(a_, o_),
                            "M2_HL": e["hl"], "M2_secondary": e["secondary"],
                            "M3_coef": m3["coef"][j] if m3["ok"][j] else np.nan,
                            "M3u_coef": m3u["coef"][j] if m3u["ok"][j] else np.nan}
                    for est, v in ests.items():
                        bias.setdefault((est, lin), []).append((v - beta[j]) / sd[j])
        if cell["scenario"] == "shared_blood":
            m4 = C.mann_whitney_less(Y, np.isin(lineage, C.BLOOD), comparator=m4_comp)
            a = acc.setdefault(("M4", "blood"), [0, 0])
            a[0] += int(m4["testable"].sum())
            a[1] += int((m4["testable"] & (m4["p"] < ALPHA)).sum())
    rows = []
    for (meth, lin), (n, r) in acc.items():
        lo, hi = C.clopper_pearson(r, n) if n else (np.nan, np.nan)
        rows.append(dict(cell=int(cell["cell"]), method=meth, tested_lineage=lin, n_testable=n, n_reject=r,
                         rate=r / n if n else np.nan, cp_lo=lo, cp_hi=hi))
    brows = []
    for (est, lin), v in bias.items():
        v = np.array(v, float)
        v = v[np.isfinite(v)]
        brows.append(dict(cell=int(cell["cell"]), estimator=est, tested_lineage=lin, n=int(len(v)),
                          mean_bias_sd_units=float(v.mean()) if len(v) else np.nan,
                          median_bias_sd_units=float(np.median(v)) if len(v) else np.nan))
    return rows, brows, round(time.time() - t0, 1)


def _worker(args):
    cell, path = args
    inp = dict(np.load(path, allow_pickle=True))
    return run_cell(cell, inp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cells", default="0-719")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--enumerate", action="store_true")
    ap.add_argument("--test", action="store_true",
                    help="code test before v1b: seeds 990000000 + i and 200 genes per cell; results are not simulation results")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    E = enumerate_cells()
    E.to_csv(os.path.join(a.out, "sim_cells_enumeration.csv"), index=False)
    if a.enumerate:
        print(len(E), "cells written")
        return
    lo, hi = map(int, a.cells.split("-"))
    cells = [r for r in E.to_dict("records") if lo <= r["cell"] <= hi]
    if a.test:
        for c in cells:
            c["seed"], c["n_genes"] = 990000000 + c["cell"], 200
    meta = dict(plan_A_sha256=C.plan_a_sha256(), inputs_sha256=C.sha256_file(a.inputs), inputs_arrays_sha256=C.arrays_sha256(a.inputs), cells=a.cells, test_mode=a.test)
    from multiprocessing import Pool
    t0 = time.time()
    with Pool(a.workers) as pool:
        out = pool.map(_worker, [(c, a.inputs) for c in cells], chunksize=1)
    Eu = pd.DataFrame(cells)                                     # the cells actually run (test mode: test seeds)
    R = pd.DataFrame([r for rows, _, _ in out for r in rows]).merge(Eu, on="cell")
    B = pd.DataFrame([r for _, br, _ in out for r in br]).merge(Eu, on="cell")
    tag = a.cells.replace("-", "_") + ("_TEST" if a.test else "")
    for name, df in (("sim_cells_results", R), ("sim_cells_bias", B)):
        path = os.path.join(a.out, f"{name}_{tag}.csv")
        with open(path, "w") as f:
            f.write(f"# plan_A_sha256={meta['plan_A_sha256']} inputs_sha256={meta['inputs_sha256']}\n")
            df.to_csv(f, index=False)
    meta.update(runtime_s=round(time.time() - t0, 1), cell_runtimes_s=[t for _, _, t in out])
    json.dump(meta, open(os.path.join(a.out, f"sim_cells_meta_{tag}.json"), "w"), indent=1)
    print("done", meta["runtime_s"], "s")


if __name__ == "__main__":
    main()
