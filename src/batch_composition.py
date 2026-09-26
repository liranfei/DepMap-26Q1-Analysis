"""Screening-batch composition of each lineage that has candidate pairs (chi-square test target lineage vs. all other cancer cell lines) and IntOGen driver-gene annotation."""
import os, json, numpy as np, pandas as pd
from scipy.stats import chi2_contingency
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True); F = pd.read_csv(O + "final_targets.csv")
s = pd.read_csv(os.path.join(D, "ScreenSequenceMap.csv")); s = s[(s.ScreenType == "2DS") & (s.PassesQC == True) & (s.DrugTreated == False) & (s.IsEngineered == False) & (s.ExcludeFromCRISPRCombined == False)]
bm = s.groupby("ModelID").pDNABatch.agg(lambda x: x.value_counts().index[0]); df["batch"] = df.DepMap_ID.map(bm).fillna("unknown")
rows = []
for L in sorted(F.Lineage.unique()):
    tgt = df.lineage == L; ct = pd.crosstab(tgt, df.batch); p = chi2_contingency(ct)[1]; top = df[tgt].batch.value_counts(normalize=True)
    rows.append(dict(Lineage=L, n=int(tgt.sum()), chi2_p=p, largest_batch=top.index[0], share_in_lineage=float(top.iloc[0]), share_in_all_cancer_lines=float((df.batch == top.index[0]).mean())))
B = pd.DataFrame(rows); B.to_csv(O + "batch_composition_by_lineage.csv", index=False)
R = dict(n_lineages=len(B), n_chi2_p_lt_0_05=int((B.chi2_p < 0.05).sum()), n_chi2_p_lt_0_05_bonferroni=int((B.chi2_p < 0.05 / len(B)).sum()),
         lymphoid=B[B.Lineage == "lymphoid"].iloc[0].to_dict(), myeloid=B[B.Lineage == "myeloid"].iloc[0].to_dict(), models_without_batch=int((df.batch == "unknown").sum()))
f = os.path.join(D, "IntOGen-DriverGenes.tsv")
if os.path.exists(f):
    g = pd.read_csv(f, sep="\t"); drivers = set(g.Symbol); cand = set(F.Gene.str.split(" (", regex=False).str[0]); hit = sorted(cand & drivers)
    R["intogen"] = dict(n_driver_genes=len(drivers), n_candidate_genes=len(cand), n_candidates_in_intogen=len(hit), genes=hit, pairs_with_driver_gene=int(F.Gene.str.split(" (", regex=False).str[0].isin(drivers).sum()))
json.dump(R, open(O + "batch_composition_summary.json", "w"), indent=1, default=float); print(json.dumps(R, indent=1, default=float))
