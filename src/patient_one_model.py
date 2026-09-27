"""Sensitivity to non-independence of models from the same patient: 41 patients contribute 87 of the 1,192 cancer models (PatientID in Model.csv).
The complete primary procedure is repeated with one model per patient, chosen without reference to the results (the model with the smallest ModelID).
Output: patient_one_model.json, all_tests_onepatient.csv.gz."""
import os, json, pandas as pd
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
pm = pd.read_csv(os.path.join(D, "Model.csv"))[["ModelID", "PatientID"]].set_index("ModelID").PatientID; df["pid"] = df.DepMap_ID.map(pm)
assert df.pid.notna().all()
keep = df.sort_values("DepMap_ID").drop_duplicates("pid").index; d1 = df.loc[sorted(keep)].drop(columns="pid").reset_index(drop=True)
genes = [c for c in d1.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
t = rp.test_all(d1, genes, d1.lineage.values); f = rp.select(t)
F = pd.read_csv(O + "final_targets.csv"); k94 = set(zip(F.Gene, F.Lineage)); k1 = set(zip(f.Gene, f.Lineage))
c = F.merge(t, on=["Gene", "Lineage"], how="left", suffixes=("", "_1pp"))
R = dict(models=len(df), patients=int(df.pid.nunique()), multi_model_patients=int((df.pid.value_counts() > 1).sum()), models_in_multi=int(df.pid.value_counts()[df.pid.value_counts() > 1].sum()),
         models_after=len(d1), tests=len(t), q_lt_0_05=int((t.q_value < rp.Q_MAX).sum()), candidates=len(f), candidate_genes=int(f.Gene.nunique()), candidate_lineages=int(f.Lineage.nunique()),
         primary94_retained=len(k94 & k1), primary94_q05_negative=int(((c.q_value_1pp < rp.Q_MAX) & (c.Selectivity_1pp < 0)).sum()), new_candidates=len(k1 - k94),
         primary94_lost=sorted(map(list, k94 - k1)))
t.to_csv(O + "all_tests_onepatient.csv.gz", index=False); json.dump(R, open(O + "patient_one_model.json", "w"), indent=1); print(json.dumps(R, indent=1))
