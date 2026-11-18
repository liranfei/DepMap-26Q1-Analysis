"""FERMT1 expression in CNS/brain lines and in the comparison groups of the FERMT2 checks (descriptive, post hoc).

fermt2_validation.py reports FERMT1 medians for six single lineages. This script adds the medians of the
comparison groups actually used in the FERMT2 checks: the carcinoma lineages of check 4 and all other
cancer cell lines. Cohort and lineage definitions are copied from fermt2_validation.py.

Environment: DEPMAP_DIR (gene_effect.csv, Model.csv, OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv), RESULTS_DIR.
Output: results/fermt2_fermt1_expression.json
"""
import os, json
import pandas as pd
from scipy.stats import mannwhitneyu

DATA = os.environ.get("DEPMAP_DIR", "."); OUT = os.environ.get("RESULTS_DIR", "results")
GENE, TARGET = "FERMT2", "CNS/Brain"
CARC = ["Bowel", "Pancreas", "Esophagus/Stomach", "Breast", "Prostate", "Lung", "Head and Neck", "Bladder/Urinary Tract", "Liver",
        "Biliary Tract", "Ampulla of Vater", "Uterus", "Ovary/Fallopian Tube", "Cervix", "Thyroid"]   # as in fermt2_validation.py, check 4

hdr = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), nrows=0).columns
col = [c for c in hdr if c.split(" ")[0] == GENE][0]
ge = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), usecols=[hdr[0], col]).rename(columns={hdr[0]: "ModelID", col: "y"}).set_index("ModelID")
M = pd.read_csv(os.path.join(DATA, "Model.csv"), usecols=["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"]).set_index("ModelID")
M = M.loc[M.index.intersection(ge.index)]
M = M[M.OncotreePrimaryDisease != "Non-Cancerous"]
y = ge.loc[M.index, "y"]; ok = y.notna()

EXPR = "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv"
eh = pd.read_csv(os.path.join(DATA, EXPR), nrows=0).columns
fc = [c for c in eh if str(c).split(" ")[0] == "FERMT1"][0]
ex = pd.read_csv(os.path.join(DATA, EXPR), usecols=["ModelID", "IsDefaultEntryForModel", fc])
ex = ex[ex.IsDefaultEntryForModel.astype(str).isin(["Yes", "True"])].set_index("ModelID")[fc].rename("FERMT1")
e = pd.DataFrame({"y": y[ok]}).join(ex, how="inner").join(M.OncotreeLineage)

t = e.OncotreeLineage == TARGET
groups = {"CNS/Brain": t,
          "carcinoma lineages (check 4)": (~t) & e.OncotreeLineage.isin(CARC),
          "all other cancer lines": ~t}
R = {"n_lines_with_CRISPR_and_expression": int(len(e)),
     "FERMT1_median_log2TPM1": {k: float(e.loc[m, "FERMT1"].median()) for k, m in groups.items()},
     "n": {k: int(m.sum()) for k, m in groups.items()},
     "p_one_sided_CNS_lower_than_carcinoma": float(mannwhitneyu(e.loc[t, "FERMT1"], e.loc[groups["carcinoma lineages (check 4)"], "FERMT1"], alternative="less")[1]),
     "FERMT1_median_by_carcinoma_lineage": {k: float(v) for k, v in e[e.OncotreeLineage.isin(CARC)].groupby("OncotreeLineage").FERMT1.median().sort_values().items()}}
json.dump(R, open(os.path.join(OUT, "fermt2_fermt1_expression.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
