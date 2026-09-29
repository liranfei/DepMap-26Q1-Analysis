"""Description of the 26Q1 cancer cohort by screening library and of the structure of missing gene effects.
Libraries per model: screens that passed QC and were used in the combined data set (ScreenSequenceMap.csv).
Structural missingness: share of missing values that lie in (library combination, gene) cells in which more than 90% of lines lack the gene.
Output: data_description.json."""
import os, json, pandas as pd
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
s = pd.read_csv(os.path.join(D, "ScreenSequenceMap.csv")); u = s[(s.PassesQC == True) & (s.ExcludeFromCRISPRCombined == False) & s.Library.notna()]
L = df.DepMap_ID.map(u.groupby("ModelID").Library.agg(lambda x: set(x))); combo = L.map(lambda x: "+".join(sorted(x)))
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]; M = df[genes].isna()
rate = M.groupby(combo.values).mean(); n = combo.value_counts()
structural = float(sum((row[row > 0.9] * n[g]).sum() for g, row in rate.iterrows()))
R = {"lines_by_library": {lib: int(L.map(lambda x: lib in x).sum()) for lib in sorted(set().union(*L))}, "library_combinations": combo.value_counts().to_dict(),
     "missing_values": int(M.values.sum()), "missing_share_structural": structural / M.values.sum(), "genes_with_missing": int((M.sum() > 0).sum())}
json.dump(R, open(O + "data_description.json", "w"), indent=1); print(json.dumps(R, indent=1))
