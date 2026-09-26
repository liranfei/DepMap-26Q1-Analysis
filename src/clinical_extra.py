"""Additional TCGA analyses (exploratory): median overall survival by expression group, progression-free interval (PFI), and Cox models adjusted for age, sex and stage.
Clinical covariates and PFI come from the TCGA Pan-Cancer Clinical Data Resource table distributed by UCSC Xena (Survival_SupplementalTable_S1_20171025_xena_sp)."""
import os, json, numpy as np, pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
O = os.environ.get("RESULTS_DIR", "primary") + "/"; CDR = os.environ.get("TCGA_CDR", "../tcga_cdr/Survival_SupplementalTable_S1_20171025_xena_sp")
cdr = pd.read_csv(CDR, sep="\t"); cdr = cdr[cdr["sample"].str.split("-").str[3].str[:2] == "01"].drop_duplicates("_PATIENT").copy(); cdr["pid"] = cdr["_PATIENT"]; assert cdr.pid.is_unique; DAYS = 30.4375; R = {}
def stage(x):
    x = str(x).upper()
    if "IV" in x: return 1
    if "III" in x: return 1
    if "II" in x or "I" in x: return 0
    return np.nan
for label, ct in [("KIRC_HNF1B", "KIRC"), ("PAAD_KRAS", "PAAD")]:
    m = pd.read_csv(O + f"tcga_{label}_survival_table.csv").merge(cdr[cdr["cancer type abbreviation"] == ct][["pid", "age_at_initial_pathologic_diagnosis", "gender", "ajcc_pathologic_tumor_stage", "PFI", "PFI.time", "OS", "OS.time"]], on="pid", how="left", suffixes=("", "_cdr"), validate="one_to_one"); assert len(m) == len(pd.read_csv(O + f"tcga_{label}_survival_table.csv"))
    out = dict(n=len(m), os_agrees_with_cdr=float(((m.OS == m.OS_cdr) & (m["OS.time"] == m["OS.time_cdr"])).mean()))
    hi, lo = m[m.hi == 1], m[m.hi == 0]
    for name, g in (("high", hi), ("low", lo)):
        k = KaplanMeierFitter().fit(g["OS.time"], g.OS); out[f"median_OS_months_{name}"] = float(k.median_survival_time_ / DAYS) if np.isfinite(k.median_survival_time_) else None; out[f"n_{name}"] = len(g); out[f"events_{name}"] = int(g.OS.sum())
    p = m.dropna(subset=["PFI", "PFI.time"]); p = p[p["PFI.time"] > 0]; p = p.assign(hiP=(p.x >= p.x.median()).astype(int), zP=lambda d: (d.x - d.x.mean()) / d.x.std())
    lr = logrank_test(p[p.hiP == 1]["PFI.time"], p[p.hiP == 0]["PFI.time"], p[p.hiP == 1].PFI, p[p.hiP == 0].PFI); c = CoxPHFitter().fit(p[["PFI.time", "PFI", "zP"]], "PFI.time", "PFI").summary.loc["zP"]
    out.update(n_pfi=len(p), events_pfi=int(p.PFI.sum()), logrank_p_pfi=float(lr.p_value), cox_pfi_HR=float(np.exp(c["coef"])), cox_pfi_CI=[float(np.exp(c["coef lower 95%"])), float(np.exp(c["coef upper 95%"]))], cox_pfi_p=float(c["p"]))
    a = m.assign(age=pd.to_numeric(m.age_at_initial_pathologic_diagnosis, errors="coerce"), male=(m.gender == "MALE").astype(float), adv=m.ajcc_pathologic_tumor_stage.map(stage)).dropna(subset=["age", "adv"])
    a = a[a["OS.time"] > 0]; a["age10"] = a.age / 10
    c2 = CoxPHFitter().fit(a[["OS.time", "OS", "z", "age10", "male", "adv"]], "OS.time", "OS"); s = c2.summary
    out.update(n_adjusted=len(a), events_adjusted=int(a.OS.sum()), adj_z_HR=float(np.exp(s.loc["z", "coef"])), adj_z_CI=[float(np.exp(s.loc["z", "coef lower 95%"])), float(np.exp(s.loc["z", "coef upper 95%"]))], adj_z_p=float(s.loc["z", "p"]),
               adj_stage_HR=float(np.exp(s.loc["adv", "coef"])), adj_age10_HR=float(np.exp(s.loc["age10", "coef"])), stage_definition="stage III-IV vs I-II (AJCC pathologic); patients without stage excluded")
    R[label] = out
json.dump(R, open(O + "tcga_clinical_extra.json", "w"), indent=1); print(json.dumps(R, indent=1))
