"""SNAP23 in bowel: selectivity among cell lines with low paralogue expression (post hoc).

Added after the pre-specified checks (prespecified_plan_snap23.md) had been run; it is a
post hoc analysis and is labelled as such in the output.

Question: check 5 found that lines expressing more SNAP25 or SNAP29 are less dependent on
SNAP23. If bowel lines express less of these paralogues, part of the bowel selectivity could
reflect paralogue expression. The bowel-versus-rest contrast is therefore repeated among lines
with low expression of both paralogues.

Primary definition (chosen as the simplest definition after several definitions had been
computed, which are all reported below as a sensitivity analysis): a line has low paralogue
expression when its log2(TPM + 1) expression is below the median of all cancer cell lines with
SNAP23 CRISPR and expression data, for both SNAP25 and SNAP29; bowel lines and comparison lines
are restricted in the same way. Selectivity = median(bowel) - median(other), as in the primary
analysis; one-sided Mann-Whitney p and a bootstrap 95% interval (1,000 resamples, seed 0) are
reported. Alternative definitions are reported as a sensitivity analysis.

Environment: DEPMAP_DIR (gene_effect.csv, Model.csv,
OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv), RESULTS_DIR.
Output: results/snap23_paralogue_restriction.json
"""
import os, json
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

DATA = os.environ.get("DEPMAP_DIR", "."); OUT = os.environ.get("RESULTS_DIR", "results")
GENE, LIN, PAR = "SNAP23", "bowel", ("SNAP25", "SNAP29")
EXPR = "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv"

# same cohort construction as snap23_validation.py
hdr = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), nrows=0).columns
col = [c for c in hdr if c.split(" ")[0] == GENE][0]
ge = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), usecols=[hdr[0], col])
ge.columns = ["DepMap_ID", "y"]
M = pd.read_csv(os.path.join(DATA, "Model.csv"))[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"]]
M.columns = ["DepMap_ID", "lineage", "primary_disease"]
M["lineage"] = M.lineage.str.strip().str.lower()
d = ge.merge(M, on="DepMap_ID", how="left", validate="one_to_one")
d = d[d.primary_disease != "Non-Cancerous"]
d["y"] = pd.to_numeric(d.y, errors="coerce"); d = d[d.y.notna()]

eh = pd.read_csv(os.path.join(DATA, EXPR), nrows=0).columns
ec = ["ModelID", "IsDefaultEntryForModel"] + [c for c in eh if c.split(" ")[0] in PAR]
ex = pd.read_csv(os.path.join(DATA, EXPR), usecols=ec)
ex = ex[ex.IsDefaultEntryForModel == "Yes"].drop(columns="IsDefaultEntryForModel")
ex = ex.rename(columns={"ModelID": "DepMap_ID"}); ex.columns = [c.split(" ")[0] for c in ex.columns]
e = d.merge(ex, on="DepMap_ID", how="inner").dropna(subset=list(PAR))
is_t = e.lineage == LIN


def contrast(mask_t, mask_o, boot=True):
    a = e.loc[is_t & mask_t, "y"].values; b = e.loc[~is_t & mask_o, "y"].values
    r = {"n_bowel": int(len(a)), "n_other": int(len(b))}
    if len(a) < 5 or len(b) < 5:
        r["note"] = "fewer than five lines"; return r
    r.update(median_bowel=float(np.median(a)), median_other=float(np.median(b)),
             selectivity=float(np.median(a) - np.median(b)),
             p_one_sided=float(mannwhitneyu(a, b, alternative="less", method="asymptotic")[1]))
    if boot:
        rng = np.random.default_rng(0)
        s = [np.median(rng.choice(a, len(a))) - np.median(rng.choice(b, len(b))) for _ in range(1000)]
        r["bootstrap_95CI"] = [float(x) for x in np.percentile(s, [2.5, 97.5])]
    return r


ALL = pd.Series(True, index=e.index)
R = {"label": "post hoc (added after the pre-specified SNAP23 checks)",
     "n_lines_with_CRISPR_and_expression": int(len(e)),
     "check_V5_spearman_rho_all_lines": {p: float(spearmanr(e[p], e.y)[0]) for p in PAR},
     "paralogue_expression_median_bowel_vs_other": {
         p: [float(e.loc[is_t, p].median()), float(e.loc[~is_t, p].median())] for p in PAR},
     "unrestricted_same_lines": contrast(ALL, ALL)}

med = {p: e[p].median() for p in PAR}
low_both = (e.SNAP25 < med["SNAP25"]) & (e.SNAP29 < med["SNAP29"])
R["primary_definition"] = "both SNAP25 and SNAP29 below the cohort median; bowel and comparison lines restricted"
R["primary"] = contrast(low_both, low_both)

sens = {}
for q, name in [(0.5, "median"), (1 / 3, "lower third"), (0.25, "lower quarter")]:
    thr = {p: e[p].quantile(q) for p in PAR}
    lb = (e.SNAP25 < thr["SNAP25"]) & (e.SNAP29 < thr["SNAP29"])
    sens[f"both below {name}"] = contrast(lb, lb, boot=False)
    for p in PAR:
        m = e[p] < thr[p]; sens[f"{p} only below {name}"] = contrast(m, m, boot=False)
sens["both below median, comparison lines only restricted"] = contrast(ALL, low_both, boot=False)
R["sensitivity"] = sens

os.makedirs(OUT, exist_ok=True)
json.dump(R, open(os.path.join(OUT, "snap23_paralogue_restriction.json"), "w"), indent=1)
print(json.dumps(R, indent=1))
