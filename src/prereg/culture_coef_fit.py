"""Distribution of culture-only effects for the genome-like mixture (plan A, Section 4.4).

Per real gene: OLS of gene effect on lineage fixed effects (every lineage with at least one line with data for the
gene) plus the 3-level culture indicators (suspension, other; adherent = reference), lines with data only. The
suspension coefficients < 0 (absolute units) form the distribution from which delta of the culture-only genes is
drawn. Part of the v1b scripts; run after v1b and before the simulation. No lineage coefficient is saved or used.

Usage: python culture_coef_fit.py --data-dir $DEPMAP_DIR --out DIR     -> DIR/culture_coefficients.csv, summary json
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat_core as C  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=os.environ.get("DEPMAP_DIR", "."))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    coh = C.load_cohort(a.data_dir)
    Y, lin, cult = coh["Y"], coh["lineage"], coh["culture"]
    levels = np.unique(lin)
    L = (lin[:, None] == levels[None, :]).astype(float)                 # lines x lineages
    susp, oth = (cult == 1).astype(float), (cult == 2).astype(float)
    coef = np.full(Y.shape[1], np.nan)
    for j in range(Y.shape[1]):
        ok = ~np.isnan(Y[:, j])
        Lj = L[ok][:, L[ok].sum(0) > 0]
        X = np.column_stack([Lj, susp[ok], oth[ok]])
        b, *_ = np.linalg.lstsq(X, Y[ok, j], rcond=None)
        coef[j] = b[-2]
    D = pd.DataFrame(dict(Gene=coh["genes"], suspension_coefficient=coef))
    path = os.path.join(a.out, "culture_coefficients.csv")
    with open(path, "w") as f:
        f.write(f"# plan_A_sha256={C.plan_a_sha256()}\n")
        D.to_csv(f, index=False)
    neg = coef[coef < 0]
    s = dict(n_genes=int(len(coef)), n_negative=int(len(neg)), median_negative=float(np.median(neg)),
             q05_negative=float(np.quantile(neg, 0.05)), sha256=C.sha256_file(path))
    json.dump(s, open(os.path.join(a.out, "culture_coefficients_summary.json"), "w"), indent=1)
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
