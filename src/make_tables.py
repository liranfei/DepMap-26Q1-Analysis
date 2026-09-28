"""Supporting tables S1-S5 (xlsx) from the result files."""
import os, json, pandas as pd
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; OUT = os.environ.get("FIG_DIR", "figures") + "/supporting_information/"; os.makedirs(OUT, exist_ok=True)
F0 = pd.read_csv(O + "final_targets.csv"); ce = set(pd.read_csv(os.path.join(D, "CRISPRInferredCommonEssentials.csv")).iloc[:, 0])
_df = rp.load(D); _df = _df[_df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
F0["n_other"] = [int((pd.to_numeric(_df[g], errors="coerce").notna() & (_df.lineage != l)).sum()) for g, l in zip(F0.Gene, F0.Lineage)]
F0["Gene_symbol"] = F0.Gene.str.extract(r"^(.+?)\s*\(")[0]; F0["in_official_common_essential_list"] = F0.Gene.isin(ce); F0["confidence_flag"] = ["low sample size (n<10)" if n < 10 else "n>=10" for n in F0.n_target]
F0 = F0.rename(columns={"Gene": "Gene (Entrez)", "Chronos_median": "Median Chronos, target lineage", "n_target": "n target (non-missing)", "n_other": "n other (non-missing)", "p_value": "p (one-sided MWU)", "q_value": "q (BH, full family)"})
F0[["Gene_symbol", "Gene (Entrez)", "Lineage", "n target (non-missing)", "n other (non-missing)", "Median Chronos, target lineage", "Selectivity", "p (one-sided MWU)", "q (BH, full family)", "in_official_common_essential_list", "confidence_flag"]].sort_values("Selectivity").to_csv(O + "Table_S1_candidates.csv", index=False)
_one = O + "all_tests_onepatient.csv.gz"
if os.path.exists(_one):
    _t1 = pd.read_csv(_one); _k1 = set(map(tuple, rp.select(_t1)[["Gene", "Lineage"]].values))
    _tc = pd.read_csv(O + "Table_S1_candidates.csv"); _tc["Retained with one model per patient"] = [(g, l) in _k1 for g, l in zip(_tc["Gene (Entrez)"], _tc.Lineage)]; _tc.to_csv(O + "Table_S1_candidates.csv", index=False)
F = pd.read_csv(O + "Table_S1_candidates.csv"); E = pd.read_csv(O + "enrichment_custom_background.csv").sort_values("q")
lc = pd.read_csv(O + "lineage_counts.csv"); lc.columns = ["Lineage", "Cancer cell lines"]; LT = lc.set_index("Lineage").join(F.groupby("Lineage").size().rename("Candidate pairs")).fillna(0).astype(int).reset_index(); LT["Eligible (n >= 5)"] = LT["Cancer cell lines"] >= 5
df = rp.load(D); nc = df[df.primary_disease == "Non-Cancerous"][["DepMap_ID", "lineage"]].rename(columns={"lineage": "Lineage"})
mo = pd.read_csv(os.path.join(D, "Model.csv"))[["ModelID", "CellLineName", "OncotreeSubtype"]].rename(columns={"ModelID": "DepMap_ID", "CellLineName": "Cell line", "OncotreeSubtype": "Subtype"}); nc = nc.merge(mo, on="DepMap_ID").drop(columns="DepMap_ID")
S = json.load(open(O + "sensitivity_summary.json")); rows = []
for k in [x for x in S if x.startswith(("multiplicity_", "selectivity<", "chronos_median<", "min_n_target"))]: rows.append({"Analysis": k.replace("multiplicity_", "Multiplicity procedure: "), **S[k]})
for r in S["topN"]: rows.append({"Analysis": f"Top-{r['top_n']} genes by variance", "n_significant": r["n_q05"], "n_final": r["n_final"], "n_genes": r["n_genes"], "n_lineages": r["n_lineages"]})
with pd.ExcelWriter(OUT + "S1_Table.xlsx") as w: LT.to_excel(w, index=False, sheet_name="Lineages"); nc.to_excel(w, index=False, sheet_name="Excluded non-cancerous")
with pd.ExcelWriter(OUT + "S2_Table.xlsx") as w: F.to_excel(w, index=False, sheet_name="Candidate pairs")
def _j(f):
    return json.load(open(O + f)) if os.path.exists(O + f) else {}
_add = []
_pp = _j("patient_one_model.json")
if _pp: _add.append({"Analysis": "One model per patient (smallest ModelID)", "Result": f"{_pp['models_after']} models; {_pp['candidates']} candidates; {_pp['primary94_retained']} of 94 retained; not retained: " + "; ".join(f"{g.split(' (')[0]} ({l})" for g, l in _pp["primary94_lost"])})
_pb = _j("permutation_summary_genome_strat_patient.json")
if _pb: _add.append({"Analysis": "Patient-block label permutation (1,000)", "Result": f"mean {_pb['null_final_mean']:.2f}, max {_pb['null_final_max']} prioritised pairs; max {_pb['null_q_max']} pairs with q < 0.05; empirical p = {_pb['empirical_p_final']:.4f}"})
_r4 = _j("revision4_summary.json")
if _r4:
    _add.append({"Analysis": "Welch t-test, BH over all tests", "Result": f"{_r4['welch']['primary94_q_lt_0_05_negative']} of 94 with q < 0.05 and negative mean difference"})
    for _x in _r4["low_n_permutation"]: _add.append({"Analysis": f"Permutation p (10,000), {_x['Gene'].split(' (')[0]} {_x['Lineage']}", "Result": f"n = {_x['n_target']}; asymptotic p = {_x['p_asymptotic']:.2g}; permutation p = {_x['p_permutation']:.2g}"})
    _add.append({"Analysis": "Batch/library-adjusted coefficients (94 pairs, post-selection)", "Result": f"95% CI below 0 for {_r4['batch_ci']['adj_ci_excludes_0']} pairs; minimum ratio adjusted/unadjusted {_r4['batch_ci']['min_ratio']:.2f} (batch_coef_ci.csv)"})
_r3 = _j("revision3_summary.json")
if _r3:
    _mp = O + "mdm2_tp53wt.csv"
    for _x in (pd.read_csv(_mp).to_dict("records") if os.path.exists(_mp) else _r3.get("MDM2_TP53wt", [])):  # unrounded p-values (the JSON values are rounded)
        if _x["cohort"] != "all lines": _add.append({"Analysis": f"MDM2 {_x['Lineage']}, no detected TP53 variant only", "Result": f"n = {_x['n_target']}; selectivity {_x['selectivity']:.2f}; p = {_x['p']:.2g}"})
    _add.append({"Analysis": "Within-lineage co-dependency, matched null (1,000 sets)", "Result": f"observed median |r| {_r3['codep_null']['observed_median_abs_r_residual']:.2f}; null mean {_r3['codep_null']['null_mean']:.2f}, max {_r3['codep_null']['null_max']:.2f}"})
_cd = _j("cross_release_decomp.json")
if _cd: _add.append({"Analysis": "22Q1 vs 26Q1 on the same 1,036 models", "Result": f"primary pairs meeting the rule: {_cd['shared_run']['primary94_meet_rule_26Q1_shared']} (26Q1 data) vs {_cd['shared_run']['primary94_meet_rule_22Q1_shared']} (22Q1 data)"})
with pd.ExcelWriter(OUT + "S3_Table.xlsx") as w:
    pd.DataFrame(rows).to_excel(w, index=False, sheet_name="Sensitivity analyses")
    if _add: pd.DataFrame(_add).to_excel(w, index=False, sheet_name="Additional checks")
with pd.ExcelWriter(OUT + "S4_Table.xlsx") as w:
    E.to_excel(w, index=False, sheet_name="Background all tested genes")
    if os.path.exists(O + "enrichment_restricted_background.csv"): pd.read_csv(O + "enrichment_restricted_background.csv").to_excel(w, index=False, sheet_name="Background median below -1")
S5COLS = {"pid": "Patient (TCGA barcode, 12 characters)", "x": "Tumour expression, log2(CPM+1) (mean per patient)", "OS.time": "Overall survival time (days)", "OS": "Death (1 = yes)", "hi": "High expression (>= median)=1", "z": "Standardised expression (z-score)"}
with pd.ExcelWriter(OUT + "S5_Table.xlsx") as w:
    for tag, sheet in [("KIRC_HNF1B", "TCGA-KIRC (HNF1B)"), ("PAAD_KRAS", "TCGA-PAAD (KRAS)")]: pd.read_csv(O + f"tcga_{tag}_survival_table.csv").rename(columns=S5COLS).to_excel(w, index=False, sheet_name=sheet)
