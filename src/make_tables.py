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
_r6 = O + "revision6_pairs.csv"
if os.path.exists(_r6):   # descriptive columns added in revision 6 (bootstrap interval, comparison-group median, coverage, adherent-only comparison)
    _p6 = pd.read_csv(_r6).rename(columns={"Gene": "Gene (Entrez)", "sel_ci_low": "Selectivity, bootstrap 95% CI low", "sel_ci_high": "Selectivity, bootstrap 95% CI high",
        "p_boot_sel_lt_minus0_5": "Bootstrap fraction with selectivity < -0.5", "median_comparison": "Median Chronos, all other lines",
        "n_other_lineages_median_lt_minus1": "Other eligible lineages with median < -1", "coverage": "Fraction of cancer lines with data for the gene",
        "selectivity_vs_adherent": "Selectivity against adherent lines of other lineages (non-blood pairs)",
        "q_vs_adherent": "q against adherent lines (BH over 52 non-blood pairs)"}).drop(columns="selectivity")
    F = F.merge(_p6, on=["Gene (Entrez)", "Lineage"], how="left", validate="one_to_one"); F["Partial gene coverage (< 50% of lines)"] = F["Fraction of cancer lines with data for the gene"] < 0.5
_ct = O + "confidence_tiers.csv"
if os.path.exists(_ct):   # confidence tiers (confidence_tiers.py)
    _t = pd.read_csv(_ct)[["Gene", "Lineage", "tier", "ev_threshold", "ev_one_model", "ev_release", "ev_culture", "ev_independent", "selectivity_without_suspension"]]
    _t.columns = ["Gene (Entrez)", "Lineage", "Robustness tier", "Evidence: bootstrap interval below -0.5", "Evidence: one model per patient", "Evidence: 22Q1", "Evidence: culture format",
                  "Evidence: re-prioritised in the discovery cohort and replicated in Project Score", "Selectivity without suspension lines in the comparison group"]
    F = F.merge(_t, on=["Gene (Entrez)", "Lineage"], how="left", validate="one_to_one")
_rv = O + "rnai_validation.csv"
if os.path.exists(_rv):   # RNAi check (rnai_validation.py)
    _r = pd.read_csv(_rv)[["Gene", "Lineage", "n_target_RNAi", "Chronos_median_RNAi", "Selectivity_RNAi", "q_value_RNAi", "confirmed"]]
    _r["confirmed"] = _r.confirmed.where(_r.q_value_RNAi.notna(), "not testable")
    _r.columns = ["Gene (Entrez)", "Lineage", "RNAi: n target", "RNAi: median DEMETER2 score, target lineage", "RNAi: selectivity", "RNAi: q (BH over all RNAi tests)", "RNAi: supported (q < 0.05, negative selectivity)"]
    F = F.merge(_r, on=["Gene (Entrez)", "Lineage"], how="left", validate="one_to_one")
F = F.rename(columns={"confidence_flag": "Low-n flag (fewer than 10 target lines)"})
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
_r6j = _j("revision6_summary.json"); _gcj = _j("growth_confound_summary.json"); _r5j = _j("revision5_summary.json")
if _gcj and "heme_median_ratio_excl_other_blood" in _gcj:
    _add.append({"Analysis": "Growth-pattern adjustment of the 42 blood pairs, with / without the other blood lineage", "Result": f"median ratio adjusted/unadjusted {_gcj['heme_median_ratio']:.2f} / {_gcj['heme_median_ratio_excl_other_blood']:.2f}; pairs with q < 0.05 {_r3['growth_bh']['heme']['q_lt_0.05_negative'] if _r3 else 'NA'} / {_gcj['heme_q_lt_0_05_excl_other_blood']} (growth_confound_pairs.csv)"})
if _r5j and "sanger_by_lineage" in _r5j:
    _L = _r5j["sanger_by_lineage"]; _add.append({"Analysis": "Sanger replication by lineage group", "Result": f"lymphoid {_L['lymphoid']['candidates_replicated']}/{_L['lymphoid']['candidates_testable']}, other {_L['other']['candidates_replicated']}/{_L['other']['candidates_testable']}; non-prioritised reference {_L['lymphoid']['reference_replicated']}/{_L['lymphoid']['reference_n']} and {_L['other']['reference_replicated']}/{_L['other']['reference_n']}; failing only the median criterion {_r5j['sanger_fail_median_only']['replicated']}/{_r5j['sanger_fail_median_only']['testable']}"})
if _r6j:
    _b = _r6j["bootstrap"]; _add.append({"Analysis": "Bootstrap 95% intervals of selectivity (1,000 resamples, seed 42; per pair in S2 Table)", "Result": f"entirely below -0.5: {_b['ci_high_lt_minus0_5']} pairs (seeds 43-45: {', '.join(str(v) for v in _b['ci_high_lt_minus0_5_other_seeds'].values())}); below 0: {_b['ci_high_lt_0']}"})
    _c = _r6j["comparison_median"]; _add.append({"Analysis": "Median of all other lines, prioritised pairs", "Result": f"above -0.5: {_c['gt_minus0_5']}; -1 to -0.5: {_c['between_minus1_and_minus0_5']}; -1 or below: {_c['le_minus1']}; another eligible lineage below -1: {_c['pairs_with_other_lineage_below_minus1']}"})
    _v = _r6j["nonblood_comparator_variants"]; _add.append({"Analysis": "Non-blood pairs with other comparison groups (selectivity < -0.5, of 52)", "Result": f"adherent lines only {_v['adherent_only']}; without suspension lines {_v['without_suspension']}; without blood lines {_v['without_blood']}"})
    for _k, _lab in (("APC_trunc_CTNNB1", "CTNNB1 by APC truncating variant"), ("BRAF_V600_BRAF", "BRAF by BRAF V600 variant")):
        _g = _r6j[_k]; _w = [v for k2, v in _g.items() if k2.startswith("in_")][0]
        _add.append({"Analysis": _lab, "Result": f"median {_g['median_carrier']:.2f} (n = {_g['n_carrier']}) vs {_g['median_other']:.2f} (n = {_g['n_other']}); lineage-adjusted difference {_g['lineage_adjusted_diff']:.2f}, p = {_g['lineage_adjusted_p']:.2g}; within lineage {_w['median_carrier']:.2f} (n = {_w['n_carrier']}) vs {_w['median_other']:.2f} (n = {_w['n_other']}), p = {_w['p_one_sided']:.2g}"})
    _bg = _r6j["bowel_CTNNB1_groups"]; _add.append({"Analysis": "CTNNB1 in bowel by APC/CTNNB1 group", "Result": f"APC truncating {_bg['apc_trunc']['median']:.2f} (n = {_bg['apc_trunc']['n']}); CTNNB1 hotspot {_bg['ctnnb1_hotspot_only']['median']:.2f} (n = {_bg['ctnnb1_hotspot_only']['n']}); neither {_bg['neither']['median']:.2f} (n = {_bg['neither']['n']}) vs other lineages {_bg['other_lineages_without_apc_trunc_or_ctnnb1_hotspot']['median']:.2f}, p = {_bg['neither_vs_other_lineages_p_one_sided']:.2f}"})
    for _k, _x in _r6j["lineage_given_genotype"].items():
        _add.append({"Analysis": f"Lineage coefficient adjusted for genotype, {_k.replace('_', ' ')}", "Result": f"ratio {_x['ratio']:.2f} (p = {_x['p_adjusted']:.2g}); interaction model: lineage effect in non-carriers {_x['interaction_lineage_effect_in_noncarriers']:.2f} (p = {_x['interaction_lineage_effect_in_noncarriers_p']:.2g}), lineage x genotype {_x['interaction_term']:.2f} (p = {_x['interaction_term_p']:.2g})"})
    _q = _r6j["cross_release_quantile"]; _add.append({"Analysis": "22Q1 with quantile-matched thresholds (shared lines)", "Result": f"selectivity {_q['matched_threshold_22Q1']:.2f}, median {_q['matched_median_threshold_22Q1']:.2f}: {_q['meet_rule_both_matched']} of {_q['testable']} pairs meet the rule; selectivity only: {_q['meet_rule_matched_threshold']}; original thresholds: {_q['meet_rule_minus0_5']}"})
_rvj = _j("rnai_validation_summary.json")
if _rvj:
    _bt = _rvj["by_tier"]; _add.append({"Analysis": "RNAi screens (DEMETER2), prioritised pairs", "Result": f"{_rvj['confirmed']} of {_rvj['testable']} testable pairs with q < 0.05 and negative selectivity; reference {_rvj['reference_confirmed']}/{_rvj['reference_testable']}; high/medium/low robustness tier {_bt['high']['confirmed']}/{_bt['high']['testable']}, {_bt['medium']['confirmed']}/{_bt['medium']['testable']}, {_bt['low']['confirmed']}/{_bt['low']['testable']}"})
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
