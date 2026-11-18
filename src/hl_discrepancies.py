"""Largest differences between the Hodges-Lehmann estimate and the difference in medians
among the 94 prioritised pairs (descriptive; reads the output of hodges_lehmann.py).

Input:  results/hodges_lehmann_significant_pairs.csv
Output: results/hodges_lehmann_top_discrepancies.json
"""
import os, json
import pandas as pd

OUT = os.environ.get("RESULTS_DIR", "results")
h = pd.read_csv(os.path.join(OUT, "hodges_lehmann_significant_pairs.csv"))
p = h[(h.q_value < 0.05) & (h.Chronos_median < -1) & (h.Selectivity_median_diff < -0.5)].copy()
assert len(p) == 94, len(p)
p["abs_difference"] = (p.Selectivity_HL - p.Selectivity_median_diff).abs()
p = p.sort_values("abs_difference", ascending=False)
cols = ["Gene", "Lineage", "Selectivity_median_diff", "Selectivity_HL", "abs_difference"]
R = {"n_prioritised": int(len(p)),
     "median_abs_difference": float(p.abs_difference.median()),
     "top10": p[cols].head(10).to_dict(orient="records"),
     "FAM50A_all_lineages": p[p.Gene.str.startswith("FAM50A ")][cols].to_dict(orient="records")}
json.dump(R, open(os.path.join(OUT, "hodges_lehmann_top_discrepancies.json"), "w"), indent=1)
print(p[cols].head(10).round(3).to_string(index=False))
