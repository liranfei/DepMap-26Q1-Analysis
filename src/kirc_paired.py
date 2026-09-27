"""HNF1B expression in TCGA-KIRC tumours versus solid tissue normal samples of the same patients (paired, two-sided Wilcoxon signed-rank test).
Most normal samples come from patients who also have a tumour sample, so the unpaired Mann-Whitney comparison treats dependent samples as independent.
Input: results/tcga_KIRC_HNF1B_expression_samples.csv (library-size normalised, log2(CPM + 1)).  Output: kirc_paired.json."""
import os, json, numpy as np, pandas as pd
from scipy.stats import wilcoxon, mannwhitneyu
O = os.environ.get("RESULTS_DIR", "results") + "/"; TD = os.environ.get("TCGA_DIR", "."); ENSG = "ENSG00000275410"   # HNF1B
E = pd.read_csv(O + "tcga_KIRC_HNF1B_expression_samples.csv")   # written by part3_tcga_enrich.py (log2 CPM + 1)
s = pd.DataFrame({"sample": E["sample"], "x": E.log2_cpm}); s["pid"] = s["sample"].str[:12]; s["type"] = s["sample"].str.split("-").str[3].str[:2]
tum = s[s.type == "01"].groupby("pid").x.mean(); nor = s[s.type == "11"].groupby("pid").x.mean(); both = tum.index.intersection(nor.index)
d = tum[both] - nor[both]; w = wilcoxon(tum[both], nor[both])
R = dict(n_tumour_patients=len(tum), n_normal_patients=len(nor), n_paired=len(both), normals_without_tumour=int(len(nor.index.difference(tum.index))),
         median_tumour_paired=float(tum[both].median()), median_normal_paired=float(nor[both].median()), median_difference=float(d.median()), frac_tumour_lower=float((d < 0).mean()),
         wilcoxon_p_two_sided=float(w.pvalue), unpaired_mw_p_all_samples=float(mannwhitneyu(s[s.type == "01"].x, s[s.type == "11"].x)[1]))
json.dump(R, open(O + "kirc_paired.json", "w"), indent=1); print(json.dumps(R, indent=1))
