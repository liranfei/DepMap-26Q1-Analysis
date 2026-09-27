"""Proportional-hazards assumption of the TCGA Cox models (unadjusted continuous expression; age/sex/stage-adjusted), tested with scaled Schoenfeld residuals
(lifelines proportional_hazard_test, rank-transformed time).  Uses the patient tables written by part3_tcga_enrich.py and the TCGA-CDR file (TCGA_CDR).
A stage-stratified model is fitted as well (stage as strata).
Output: ph_check.json."""
import os, re, json, numpy as np, pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
O = os.environ.get("RESULTS_DIR", "results") + "/"; CDR = os.environ.get("TCGA_CDR", "Survival_SupplementalTable_S1_20171025_xena_sp")
cdr = pd.read_csv(CDR, sep="\t"); cdr = cdr[cdr["sample"].str.split("-").str[3].str[:2] == "01"].drop_duplicates("_PATIENT").copy(); cdr["pid"] = cdr["_PATIENT"]
def stage(x):
    mm = re.match(r"^Stage (I|II|III|IV)[ABC]?$", str(x).strip()); return np.nan if mm is None else (1.0 if mm.group(1) in ("III", "IV") else 0.0)
R = {}
for label, ct in [("KIRC_HNF1B", "KIRC"), ("PAAD_KRAS", "PAAD")]:
    m = pd.read_csv(O + f"tcga_{label}_survival_table.csv"); out = {}
    c = CoxPHFitter().fit(m[["OS.time", "OS", "z"]], "OS.time", "OS"); t = proportional_hazard_test(c, m[["OS.time", "OS", "z"]], time_transform="rank")
    out["unadjusted_z_p"] = float(t.summary.loc["z", "p"])
    a = m.merge(cdr[cdr["cancer type abbreviation"] == ct][["pid", "age_at_initial_pathologic_diagnosis", "gender", "ajcc_pathologic_tumor_stage"]], on="pid", how="left", validate="one_to_one")
    a = a.assign(age10=pd.to_numeric(a.age_at_initial_pathologic_diagnosis, errors="coerce") / 10, male=(a.gender == "MALE").astype(float), adv=a.ajcc_pathologic_tumor_stage.map(stage)).dropna(subset=["age10", "adv"])
    a = a[a["OS.time"] > 0][["OS.time", "OS", "z", "age10", "male", "adv"]]; c2 = CoxPHFitter().fit(a, "OS.time", "OS"); t2 = proportional_hazard_test(c2, a, time_transform="rank")
    out["adjusted"] = {k: float(v) for k, v in t2.summary.p.items()}; out["n_adjusted"] = len(a)
    # stage-stratified model (stage as strata instead of covariate), used when stage violates proportional hazards
    c3 = CoxPHFitter().fit(a, "OS.time", "OS", strata=["adv"]); s3 = c3.summary.loc["z"]; t3 = proportional_hazard_test(c3, a, time_transform="rank")
    out["stage_stratified_z"] = {"HR": float(np.exp(s3["coef"])), "CI": [float(np.exp(s3["coef lower 95%"])), float(np.exp(s3["coef upper 95%"]))], "p": float(s3["p"]), "ph_p": {k: float(v) for k, v in t3.summary.p.items()}}
    R[label] = out
json.dump(R, open(O + "ph_check.json", "w"), indent=1); print(json.dumps(R, indent=1))
