"""KRAS expression and overall survival in TCGA-PAAD restricted to the curated PDAC samples of the TCGA marker paper (Raphael et al., Cancer Cell 2017,
Table S1, sheet FreezeSamples: 150 samples; 76 classed as high purity).  TCGA-PAAD also contains non-ductal and low-purity tumours.
Inputs: results/tcga_PAAD_KRAS_survival_table.csv (patient-level expression and OS), PAAD_S1 (path to Table S1 xlsx = NIHMS898701-supplement-2.xlsx),
TCGA_CDR (for age, sex, stage).  Expression is re-standardised within each subset.  Output: paad_pdac_subset.json."""
import os, re, json, numpy as np, pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
O = os.environ.get("RESULTS_DIR", "results") + "/"; S1 = os.environ.get("PAAD_S1", "TCGA_PAAD_2017_TableS1.xlsx"); CDR = os.environ.get("TCGA_CDR", "Survival_SupplementalTable_S1_20171025_xena_sp")
fs = pd.read_excel(S1, sheet_name="FreezeSamples", header=1); fs["pid"] = fs["Tumor Sample ID"].str[:12]
m = pd.read_csv(O + "tcga_PAAD_KRAS_survival_table.csv")
cdr = pd.read_csv(CDR, sep="\t"); cdr = cdr[cdr["sample"].str.split("-").str[3].str[:2] == "01"].drop_duplicates("_PATIENT").rename(columns={"_PATIENT": "pid"})
def stage(x):
    mm = re.match(r"^Stage (I|II|III|IV)[ABC]?$", str(x).strip()); return np.nan if mm is None else (1.0 if mm.group(1) in ("III", "IV") else 0.0)
m = m.merge(cdr[["pid", "age_at_initial_pathologic_diagnosis", "gender", "ajcc_pathologic_tumor_stage"]], on="pid", how="left", validate="one_to_one")
m["age10"] = pd.to_numeric(m.age_at_initial_pathologic_diagnosis, errors="coerce") / 10; m["male"] = (m.gender == "MALE").astype(float); m["adv"] = m.ajcc_pathologic_tumor_stage.map(stage)
R = {"freeze_samples": len(fs), "high_purity": int((fs["Purity Class (high or low)"] == "high").sum()), "all_patients": len(m)}
for name, ids in (("all_PAAD", set(m.pid)), ("curated_150", set(fs.pid)), ("high_purity_76", set(fs[fs["Purity Class (high or low)"] == "high"].pid))):
    s = m[m.pid.isin(ids)].copy(); s["z"] = (s.x - s.x.mean()) / s.x.std(); s["hi"] = (s.x >= s.x.median()).astype(int)
    lr = logrank_test(s[s.hi == 1]["OS.time"], s[s.hi == 0]["OS.time"], s[s.hi == 1].OS, s[s.hi == 0].OS)
    c = CoxPHFitter().fit(s[["OS.time", "OS", "z"]], "OS.time", "OS"); q = c.summary.loc["z"]; ph = proportional_hazard_test(c, s[["OS.time", "OS", "z"]], time_transform="rank").summary.loc["z", "p"]
    out = dict(n=len(s), events=int(s.OS.sum()), logrank_p=float(lr.p_value), HR=float(np.exp(q["coef"])), CI=[float(np.exp(q["coef lower 95%"])), float(np.exp(q["coef upper 95%"]))], p=float(q["p"]), ph_p=float(ph))
    a = s.dropna(subset=["age10", "adv"])[["OS.time", "OS", "z", "age10", "male", "adv"]]
    if a.OS.sum() >= 20:
        c2 = CoxPHFitter().fit(a, "OS.time", "OS"); q2 = c2.summary.loc["z"]
        out.update(adj_n=len(a), adj_HR=float(np.exp(q2["coef"])), adj_CI=[float(np.exp(q2["coef lower 95%"])), float(np.exp(q2["coef upper 95%"]))], adj_p=float(q2["p"]))
    R[name] = out
json.dump(R, open(O + "paad_pdac_subset.json", "w"), indent=1); print(json.dumps(R, indent=1))
