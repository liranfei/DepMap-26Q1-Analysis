"""Inputs of the simulation (plan A, Section 4.1): per real gene the median gene effect over all cancer lines, the
pooled within-lineage SD and the missing-value mask, plus the real lineage labels and 3-level culture formats.
No lineage test is computed.

Pooled within-lineage SD of gene g [N8]: sqrt( sum over lineages l, lines i in l with data of (y_i - mean_l)^2
/ (n_g - L_g) ), n_g lines with data, L_g lineages with at least one line with data.

Usage: python sim_inputs.py --data-dir $DEPMAP_DIR --out DIR      -> DIR/sim_inputs.npz, DIR/sim_inputs_summary.json
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import strat_core as C  # noqa: E402


def pooled_within_sd(Y, lineage):
    ok = ~np.isnan(Y)
    ss = np.zeros(Y.shape[1])
    n = ok.sum(0).astype(float)
    L = np.zeros(Y.shape[1])
    for lin in np.unique(lineage):
        m = lineage == lin
        sub = Y[m]
        cnt = (~np.isnan(sub)).sum(0)
        mean = np.divide(np.nansum(sub, 0), cnt, out=np.zeros(Y.shape[1]), where=cnt > 0)
        ss += np.nansum((sub - mean[None, :]) ** 2, 0)
        L += cnt > 0
    return np.sqrt(ss / (n - L))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=os.environ.get("DEPMAP_DIR", "."))
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    coh = C.load_cohort(a.data_dir)
    Y = coh["Y"]
    mu = np.nanmedian(Y, axis=0)
    sd = pooled_within_sd(Y, coh["lineage"])
    assert np.isfinite(mu).all() and np.isfinite(sd).all() and (sd > 0).all()
    path = os.path.join(a.out, "sim_inputs.npz")
    np.savez_compressed(path, mu=mu, sigma=sd, mask=~np.isnan(Y), lineage=coh["lineage"], culture=coh["culture"],
                        growth_pattern=coh["growth_pattern"], genes=coh["genes"])
    q = np.quantile(sd, [0.05, 0.5, 0.95])
    s = dict(plan_A_sha256=C.plan_a_sha256(), n_genes=int(len(mu)), n_lines=int(Y.shape[0]),
             sigma_median=float(q[1]), sigma_q05=float(q[0]), sigma_q95=float(q[2]),
             spearman_sigma_vs_median=float(spearmanr(sd, mu).statistic),
             fraction_missing=float(np.isnan(Y).mean()), sim_inputs_sha256=C.sha256_file(path), sim_inputs_arrays_sha256=C.arrays_sha256(path),
             plan_section0_values=dict(sigma_median=0.133, sigma_q05=0.102, sigma_q95=0.334, spearman=-0.64))
    json.dump(s, open(os.path.join(a.out, "sim_inputs_summary.json"), "w"), indent=1)
    print(json.dumps(s, indent=1))


if __name__ == "__main__":
    main()
