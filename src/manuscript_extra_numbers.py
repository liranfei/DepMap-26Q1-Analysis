"""Numbers reported in the manuscript that were previously computed only interactively.

Writes them to one file so that each can be traced to code:
- cohort counts (models and genes in the CRISPR matrix, non-cancerous models, lineages, lineage sizes);
- growth pattern of the cancer cell lines with CRISPR data (suspension cultures and blood lineages);
- mean SNAP23 gene effect of the bowel lines in the 2D data;
- the reverse direction of the release comparison: the 22Q1 candidates evaluated in the complete 26Q1 analysis.

Environment: DEPMAP_DIR (gene_effect.csv, Model.csv), RESULTS_DIR (all_tests.csv or all_tests.csv.gz, all_tests_22Q1.csv.gz).
Output: results/manuscript_extra_numbers.json
"""
import os, json
import pandas as pd

DATA = os.environ.get("DEPMAP_DIR", "."); OUT = os.environ.get("RESULTS_DIR", "results")
R = {}

# cohort
hdr = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), nrows=0).columns
ids = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), usecols=[hdr[0]]).iloc[:, 0]
M = pd.read_csv(os.path.join(DATA, "Model.csv"))
Mi = M[M.ModelID.isin(set(ids))]
C = Mi[Mi.OncotreePrimaryDisease != "Non-Cancerous"]
vc = C.OncotreeLineage.value_counts()
R["cohort"] = {"models_in_crispr_matrix": int(len(ids)), "genes_in_crispr_matrix": int(len(hdr) - 1),
               "non_cancerous_removed": int(len(Mi) - len(C)), "cancer_lines": int(len(C)),
               "lineages": int(len(vc)), "lineages_with_ge5": int((vc >= 5).sum()),
               "models_in_lineages_with_ge5": int(vc[vc >= 5].sum()),
               "lineage_size_range_ge5": [int(vc[vc >= 5].min()), int(vc[vc >= 5].max())]}

# growth pattern
sus = C[C.GrowthPattern == "Suspension"]
R["growth_pattern"] = {"suspension": int(len(sus)),
                       "suspension_lymphoid_or_myeloid": int(sus.OncotreeLineage.isin(["Lymphoid", "Myeloid"]).sum()),
                       "suspension_other": int((~sus.OncotreeLineage.isin(["Lymphoid", "Myeloid"])).sum())}

# SNAP23 bowel mean in the 2D data
col = [c for c in hdr if c.split(" ")[0] == "SNAP23"][0]
g = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), usecols=[hdr[0], col]); g.columns = ["ModelID", "y"]
g = g.merge(C[["ModelID", "OncotreeLineage"]], on="ModelID").dropna(subset=["y"])
R["SNAP23_bowel_2D"] = {"n": int((g.OncotreeLineage == "Bowel").sum()),
                        "mean": float(g[g.OncotreeLineage == "Bowel"].y.mean())}

# reverse direction of the release comparison (complete analyses, not only the shared lines)
a = pd.read_csv(os.path.join(OUT, "all_tests_22Q1.csv.gz")); b = pd.read_csv(os.path.join(OUT, "all_tests.csv" if os.path.exists(os.path.join(OUT, "all_tests.csv")) else "all_tests.csv.gz"))
rule = lambda d: (d.q_value < 0.05) & (d.Chronos_median < -1) & (d.Selectivity < -0.5)
key = lambda d: d.Gene.str.split(" ").str[0] + "|" + d.Lineage
c22 = a[rule(a)].assign(k=lambda d: key(d))
b = b.assign(k=key(b)).set_index("k")
m = c22.join(b[["q_value", "Selectivity", "Chronos_median"]], on="k", rsuffix="_26")
R["reverse_22Q1_candidates_in_26Q1"] = {
    "n_22Q1_candidates": int(len(c22)), "testable_in_26Q1": int(m.q_value_26.notna().sum()),
    "q_lt_0_05_negative_in_26Q1": int(((m.q_value_26 < 0.05) & (m.Selectivity_26 < 0)).sum()),
    "meet_rule_in_26Q1": int(((m.q_value_26 < 0.05) & (m.Selectivity_26 < -0.5) & (m.Chronos_median_26 < -1)).sum()),
    "median_selectivity_22Q1": float(m.Selectivity.median()),
    "median_selectivity_26Q1": float(m.Selectivity_26.median())}

json.dump(R, open(os.path.join(OUT, "manuscript_extra_numbers.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
