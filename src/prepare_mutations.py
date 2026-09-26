"""Extract the KRAS, CTNNB1 and TP53 rows of OmicsSomaticMutations.csv (all needed columns) to results/mut_3genes.csv."""
import os, pandas as pd
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results"); os.makedirs(O, exist_ok=True)
m = pd.read_csv(os.path.join(D, "OmicsSomaticMutations.csv"), usecols=["ModelID", "HugoSymbol", "VepImpact", "Hotspot", "LikelyLoF", "IsDefaultEntryForModel"], low_memory=False)
m[m.HugoSymbol.isin(["KRAS", "CTNNB1", "TP53"])].to_csv(os.path.join(O, "mut_3genes.csv"), index=False)
