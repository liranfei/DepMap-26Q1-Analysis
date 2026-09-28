"""Fifth-round checks (revision 5).

1. Candidate genes measured in only part of the cohort: number of non-missing lines, target lines and the screening libraries of the lines with data
   (for such genes the lineage comparison is confined to the lines screened with the libraries that measured the gene).
2. Within-lineage co-dependency of the metabolic candidates after removing the means of lineage-by-growth-pattern groups instead of lineage means,
   and separately in blood and non-blood lines (checks whether growth pattern within lineages explains the residual correlation).
Output: revision5_summary.json."""
import os, json, numpy as np, pandas as pd
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; R = {}
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]; col = lambda s: [c for c in genes if c.startswith(s + " (")][0]
F = pd.read_csv(O + "final_targets.csv")
s = pd.read_csv(os.path.join(D, "ScreenSequenceMap.csv")); u = s[(s.PassesQC == True) & (s.ExcludeFromCRISPRCombined == False) & s.Library.notna()]
lib = df.DepMap_ID.map(u.groupby("ModelID").Library.agg(lambda x: "+".join(sorted(set(x)))))

# ---- 1. partially measured candidate genes
rows = []
for g in F.Gene.unique():
    ok = df[g].notna()
    if ok.sum() >= 1100: continue
    for L in F[F.Gene == g].Lineage:
        rows.append(dict(Gene=g, Lineage=L, non_missing=int(ok.sum()), target_non_missing=int((ok & (df.lineage == L)).sum()), target_total=int((df.lineage == L).sum()),
                         libraries_with_data=lib[ok].value_counts().to_dict(), libraries_missing=lib[~ok].value_counts().to_dict()))
R["partially_measured"] = rows

# ---- 2. co-dependency with lineage-by-growth-pattern means removed
gp = df.DepMap_ID.map(pd.read_csv(os.path.join(D, "Model.csv")).set_index("ModelID").GrowthPattern).fillna("Unknown")
MET = [col(g) for g in ["NAMPT", "NMNAT1", "ADSL", "PAICS", "UMPS", "CTPS1", "DHFR", "TYMS"]]; X = df[MET].astype(float); iu = np.triu_indices(len(MET), 1)
r1 = X - X.groupby(df.lineage).transform("mean"); r2 = X - X.groupby(df.lineage + "|" + gp).transform("mean"); blood = df.lineage.isin(["lymphoid", "myeloid"])
med = lambda Z: float(np.median(np.abs(Z.corr().values[iu])))
R["codep_growth"] = {"residual_lineage": med(r1), "residual_lineage_x_growth": med(r2), "residual_lineage_x_growth_nonblood": med(r2[~blood]), "residual_lineage_blood": med(r1[blood]),
                     "ADSL_PAICS_lineage_x_growth": float(r2.corr().loc[col("ADSL"), col("PAICS")])}
json.dump(R, open(O + "revision5_summary.json", "w"), indent=1, default=float); print(json.dumps(R, indent=1, default=float))
