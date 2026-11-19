"""Expectations E1-E3 and decision rule of plan A, Sections 4.5-4.6, applied to the simulation cells only.

Input: sim_cells_results_0_719.csv (simulation.py, all 720 cells, not test mode).
Output: sim_decision.json and sim_decision_statement.md (with SHA-256; deposited in v2 before any real-data result).

E1  single:lymphoid and single:myeloid, c_beta = 0, c_delta != 0, k = 1, both noise types:
    M1 rejection rate > 0.10 in every such cell [N9].
E2  for M2 and for M3 separately: single-lineage scenarios (single:*, not misspec:, not shared_blood), c_beta = 0,
    k = 1, both noise types, every c_delta: upper Clopper-Pearson bound <= 0.06 in every such cell.
E3  for skin and lung (single:), c_beta in {-0.1, -0.2, -0.35}, k = 1, both noise types, every c_delta:
    in every cell whose comparator power is in [0.2, 0.9]: M2 >= 0.9 x M1 and M3 >= 0.9 x M3u.
    If no such cell exists for a lineage and method pair, the cell with comparator power closest to 0.5
    (first in enumeration order if tied) is used [N9].
Decision (4.6): M2 meets E2 and E3 -> M2; else M3 meets both -> M3; both -> M2; neither -> M2 labelled
    "not recommended by the simulation".
M4 complement (4.6, last rule) [N10]: in shared_blood, cells with k = 1 and c_beta != 0, for lymphoid and myeloid,
    the median over cells of (M2 rate / M1 rate) < 0.5 -> "M2 loses most power", M4 presented as the necessary
    complement for blood lineages.
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat_core as C  # noqa: E402


def rate(R, cell, method, lin):
    r = R[(R.cell == cell) & (R.method == method) & (R.tested_lineage == lin)]
    assert len(r) == 1, (cell, method, lin, len(r))
    return r.iloc[0]


def e1(R, E):
    cells = E[E.scenario.isin(["single:lymphoid", "single:myeloid"]) & (E.c_beta == 0) & (E.c_delta != 0) & (E.k == 1)]
    rows = []
    for c in cells.itertuples():
        r = rate(R, c.cell, "M1", c.scenario.split(":")[1])
        rows.append(dict(cell=c.cell, scenario=c.scenario, c_delta=c.c_delta, noise=c.noise, M1_rate=r.rate, ok=bool(r.rate > 0.10)))
    return dict(cells=rows, met=bool(all(r["ok"] for r in rows)))


def e2(R, E, method):
    cells = E[E.scenario.str.startswith("single:") & (E.c_beta == 0) & (E.k == 1)]
    rows = []
    for c in cells.itertuples():
        r = rate(R, c.cell, method, c.scenario.split(":")[1])
        rows.append(dict(cell=c.cell, scenario=c.scenario, c_delta=c.c_delta, noise=c.noise, rate=r.rate,
                         cp_hi=r.cp_hi, n=int(r.n_testable), ok=bool(r.cp_hi <= 0.06)))
    return dict(cells=rows, met=bool(all(r["ok"] for r in rows)))


def e3(R, E, method, comparator):
    out = {}
    for lin in ("skin", "lung"):
        cells = E[(E.scenario == f"single:{lin}") & E.c_beta.isin([-0.1, -0.2, -0.35]) & (E.k == 1)].sort_values("cell")
        rows = []
        for c in cells.itertuples():
            a, b = rate(R, c.cell, method, lin), rate(R, c.cell, comparator, lin)
            rows.append(dict(cell=c.cell, c_beta=c.c_beta, c_delta=c.c_delta, noise=c.noise, power=a.rate,
                             comparator_power=b.rate, in_range=bool(0.2 <= b.rate <= 0.9),
                             ok=bool(a.rate >= 0.9 * b.rate)))
        used = [r for r in rows if r["in_range"]]
        fallback = False
        if not used:
            fallback = True
            used = [min(rows, key=lambda r: (abs(r["comparator_power"] - 0.5), r["cell"]))]
        out[lin] = dict(cells=rows, cells_used=[r["cell"] for r in used], fallback=fallback,
                        met=bool(all(r["ok"] for r in used)))
    return dict(per_lineage=out, met=bool(all(v["met"] for v in out.values())))


def m4_rule(R, E):
    cells = E[(E.scenario == "shared_blood") & (E.k == 1) & (E.c_beta != 0)]
    ratios = []
    for c in cells.itertuples():
        for lin in C.BLOOD:
            m1, m2 = rate(R, c.cell, "M1", lin).rate, rate(R, c.cell, "M2", lin).rate
            ratios.append(dict(cell=c.cell, lineage=lin, M1=m1, M2=m2, ratio=m2 / m1 if m1 > 0 else np.nan))
    med = float(np.nanmedian([r["ratio"] for r in ratios]))
    return dict(cells=ratios, median_ratio=med, m2_loses_most_power=bool(med < 0.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    R = pd.read_csv(a.results, comment="#")
    E = R[["cell", "scenario", "c_beta", "c_delta", "noise", "k", "seed", "n_genes"]].drop_duplicates()
    assert len(E) == 720 and (E.seed == 20261119 + E.cell).all(), "needs all 720 cells of the real (not test) run"
    res = dict(plan_A_sha256=C.plan_a_sha256(), results_sha256=C.sha256_file(a.results))
    res["E1"] = e1(R, E)
    res["E2"] = {m: e2(R, E, m) for m in ("M2", "M3")}
    res["E3"] = {"M2_vs_M1": e3(R, E, "M2", "M1"), "M3_vs_M3u": e3(R, E, "M3", "M3u")}
    m2ok = res["E2"]["M2"]["met"] and res["E3"]["M2_vs_M1"]["met"]
    m3ok = res["E2"]["M3"]["met"] and res["E3"]["M3_vs_M3u"]["met"]
    if m2ok:
        rec, label = "M2", "recommended by the simulation"
    elif m3ok:
        rec, label = "M3", "recommended by the simulation"
    else:
        rec, label = "M2", "not recommended by the simulation"
    res["decision"] = dict(M2_meets_E2_E3=m2ok, M3_meets_E2_E3=m3ok, recommended=rec, label=label,
                           M3_agrees=(m3ok if rec == "M2" and m2ok else None))
    res["M4_complement"] = m4_rule(R, E)
    os.makedirs(a.out, exist_ok=True)
    jp = os.path.join(a.out, "sim_decision.json")
    json.dump(res, open(jp, "w"), indent=1)
    txt = (f"# Recommended adjusted procedure (plan A, Section 4.6)\n\n"
           f"Derived from the simulation results only (`{os.path.basename(a.results)}`, SHA-256 {res['results_sha256']}), "
           f"before any real-data result of plan A exists.\n\n"
           f"- E1 (M1 type I error > 0.10 with a blood target and culture effect): {'met' if res['E1']['met'] else 'NOT met'}\n"
           f"- E2: M2 {'met' if res['E2']['M2']['met'] else 'NOT met'}; M3 {'met' if res['E2']['M3']['met'] else 'NOT met'}\n"
           f"- E3: M2 vs M1 {'met' if res['E3']['M2_vs_M1']['met'] else 'NOT met'}; M3 vs M3u {'met' if res['E3']['M3_vs_M3u']['met'] else 'NOT met'}\n"
           f"- **Recommended method: {rec} ({label})**\n"
           f"- M4 presented as the necessary complement for blood lineages: {'yes' if res['M4_complement']['m2_loses_most_power'] else 'no'} "
           f"(median M2/M1 power ratio in the shared-blood scenario {res['M4_complement']['median_ratio']:.3f})\n\n"
           f"plan A SHA-256: {res['plan_A_sha256']}\n")
    sp = os.path.join(a.out, "sim_decision_statement.md")
    open(sp, "w").write(txt)
    print(txt)
    print("statement SHA-256:", C.sha256_file(sp))


if __name__ == "__main__":
    main()
