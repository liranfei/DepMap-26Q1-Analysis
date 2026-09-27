"""Effect-size-matched reference for the Sanger comparison.  Candidates (prioritised 26Q1 pairs testable in the Sanger data) are compared with non-prioritised 26Q1 pairs
(q<0.05 and negative selectivity) that are testable in the same Sanger lineages, matched 1:5 without replacement on 26Q1 selectivity and 26Q1 median Chronos (nearest neighbour on standardised values),
and with a logistic regression of replication on prioritisation adjusted for these effect sizes and log target-lineage size.  Replication = q<0.05 (Sanger family) and negative selectivity in the Sanger data."""
import os, json, numpy as np, pandas as pd, statsmodels.api as sm
O = os.environ.get("RESULTS_DIR", "primary") + "/"
T26 = pd.read_csv(O + "all_tests.csv"); TS = pd.read_csv(O + "all_tests_sanger.csv.gz"); F = pd.read_csv(O + "final_targets.csv")
M = T26.merge(TS, on=["Gene", "Lineage"], suffixes=("_26", "_S")); M["cand"] = pd.MultiIndex.from_frame(M[["Gene", "Lineage"]]).isin(pd.MultiIndex.from_frame(F[["Gene", "Lineage"]]))
M["rep"] = ((M.q_value_S < 0.05) & (M.Selectivity_S < 0)).astype(int); pool = M[(M.q_value_26 < 0.05) & (M.Selectivity_26 < 0)].copy()
C = pool[pool.cand]; N = pool[~pool.cand].copy(); out = dict(n_candidates_testable=len(C), n_noncand_pool=len(N), cand_rep=int(C.rep.sum()), cand_rate=float(C.rep.mean()))
cols = ["Selectivity_26", "Chronos_median_26"]; sd = pool[cols].std(); k = 5; used = set(); rows = []
for _, r in C.sort_values("Selectivity_26").iterrows():
    d = (((N[cols].astype(float) - r[cols].astype(float)) / sd) ** 2).sum(axis=1).astype(float); d = d[~d.index.isin(used)].nsmallest(k); used |= set(d.index); rows.append(N.loc[d.index])
Mt = pd.concat(rows); out.update(n_matched=len(Mt), matched_rep=int(Mt.rep.sum()), matched_rate=float(Mt.rep.mean()),
    smd_selectivity=float((C.Selectivity_26.mean() - Mt.Selectivity_26.mean()) / sd["Selectivity_26"]), smd_median=float((C.Chronos_median_26.mean() - Mt.Chronos_median_26.mean()) / sd["Chronos_median_26"]),
    cand_mean_sel=float(C.Selectivity_26.mean()), matched_mean_sel=float(Mt.Selectivity_26.mean()), cand_mean_med=float(C.Chronos_median_26.mean()), matched_mean_med=float(Mt.Chronos_median_26.mean()))
# bootstrap CI for the difference in replication rate (resample candidates and matched sets)
rng = np.random.default_rng(1); ci = []
groups = [Mt.iloc[i * k:(i + 1) * k].rep.values for i in range(len(C))]; cr = C.sort_values("Selectivity_26").rep.values
for _ in range(5000):
    ix = rng.integers(0, len(C), len(C)); ci.append(cr[ix].mean() - np.concatenate([groups[i] for i in ix]).mean())
out["diff_rate"] = out["cand_rate"] - out["matched_rate"]; out["diff_CI95"] = [float(np.percentile(ci, 2.5)), float(np.percentile(ci, 97.5))]
# logistic regression on the whole pool
pool["logn"] = np.log(pool.n_target_26); X = sm.add_constant(pool[["cand", "Selectivity_26", "Chronos_median_26", "logn"]].astype(float)); lr = sm.Logit(pool.rep, X).fit(disp=0); s = lr.summary2().tables[1]
out["logit_OR_prioritised"] = float(np.exp(s.loc["cand", "Coef."])); out["logit_OR_CI"] = [float(np.exp(s.loc["cand", "[0.025"])), float(np.exp(s.loc["cand", "0.975]"]))]; out["logit_p"] = float(s.loc["cand", "P>|z|"]); out["logit_n"] = int(len(pool))
X2 = sm.add_constant(pool[["cand"]].astype(float)); l2 = sm.Logit(pool.rep, X2).fit(disp=0); out["unadjusted_OR"] = float(np.exp(l2.params["cand"]))
json.dump(out, open(O + "sanger_matched_summary.json", "w"), indent=1); print(json.dumps(out, indent=1))
# replication rate by 26Q1 selectivity band among all 26Q1 pairs with q<0.05 and negative selectivity that are testable in Sanger
bands = [(-9, -0.7), (-0.7, -0.5), (-0.5, -0.4), (-0.4, -0.3), (-0.3, -0.2), (-0.2, 0)]; bt = []
for lo, hi in bands:
    b = pool[(pool.Selectivity_26 >= lo) & (pool.Selectivity_26 < hi)]; bt.append(dict(band=f"[{lo},{hi})", n=int(len(b)), n_prioritised=int(b.cand.sum()), replication_rate=float(b.rep.mean()) if len(b) else None))
out["bands"] = bt; json.dump(out, open(O + "sanger_matched_summary.json", "w"), indent=1); print(json.dumps(bt, indent=1))
