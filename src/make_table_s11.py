"""S11 Table: the pre-specified checks of the SNAP23 dependency in bowel, the 50 strongest
co-dependencies and the lineage selectivity of the SNARE partners.
Environment: RESULTS_DIR (reads snap23_validation_summary.json and snap23_codependency_all_genes.csv)."""
import os, json, pandas as pd
O = os.environ.get("RESULTS_DIR", "results")
R = json.load(open(O + "/snap23_validation_summary.json")); P = R["post_hoc_additions"]
V1, V2, V3, V4 = R["V1_primary"], R["V2_genotype"], R["V3_culture_format"], R["V4_subtype"]
V5, V7, V8 = R["V5_paralogue_buffering"], R["V7_orthogonal"], R["V8_expression"]
rows = [
 ("V1", "Primary recomputation, bowel vs all other cancer lines",
  f"median {V1['median_bowel']:.2f} vs {V1['median_other']:.2f}; selectivity {V1['selectivity']:.2f}; "
  f"p = {V1['p_one_sided']:.1e}; bootstrap 95% CI {V1['bootstrap_CI'][0]:.2f} to {V1['bootstrap_CI'][1]:.2f}",
  "q < 0.05, median < -1, selectivity < -0.5, bootstrap CI entirely < -0.5", "pass"),
 ("V2", "Independence from genotype (KRAS, APC, TP53, BRAF, PIK3CA)",
  f"coefficient {V2['coef_unadjusted']:.3f} -> {V2['coef_adjusted']:.3f} "
  f"({100*V2['ratio_adj_over_unadj']:.0f}%); p = {V2['p_adjusted']:.1e}",
  "ratio >= 0.70 and p < 0.05", "pass"),
 ("V3", "Culture format: bowel vs adherent non-bowel lines only",
  f"selectivity {V3['selectivity_vs_adherent']:.2f}; p = {V3['p_one_sided']:.1e}; n = {V3['n_adherent_other']}",
  "selectivity < -0.5 and p < 0.05", "pass"),
 ("V4", "Disease subtypes with at least five lines",
  "; ".join(f"{k}: n = {v['n']}, selectivity {v['selectivity']:.2f}, p = {v['p_one_sided']:.1e}"
            for k, v in V4["subtypes"].items()),
  "selectivity < -0.5 in every subtype", "pass"),
 ("V5", "Paralogue buffering (direction fixed in advance: positive)",
  f"SNAP25 rho = {V5['SNAP25']['spearman_rho']:.2f} (p = {V5['SNAP25']['p_two_sided']:.1e}); "
  f"SNAP29 rho = {V5['SNAP29']['spearman_rho']:.2f} (p = {V5['SNAP29']['p_two_sided']:.1e})",
  "positive rho with p < 0.05 for at least one paralogue", "pass"),
 ("V6", "Co-dependency with all other genes after removing lineage means",
  f"STXBP3 r = 0.481 (rank 1) and STX4 r = 0.477 (rank 2) of {P['V6_enrichment_posthoc']['n_genes']} genes; "
  f"hypergeometric P = {P['V6_enrichment_posthoc']['hypergeometric_P_ge_k']:.1e} (post hoc)",
  "at least one pre-specified SNARE partner in the top 20", "pass"),
 ("V7a", "Sanger Project Score (cell lines not in the Broad screens)",
  f"n = {V7['sanger_project_score']['n_bowel']} bowel vs {V7['sanger_project_score']['n_other']}; "
  f"selectivity {V7['sanger_project_score']['selectivity']:.2f}; "
  f"q = {P['V7_with_q_criterion']['sanger_project_score']['q_value']:.1e}",
  "p < 0.05 with negative selectivity", "pass"),
 ("V7b", "DEMETER2 RNAi screens",
  f"n = {V7['rnai_demeter2']['n_bowel']} bowel; selectivity {V7['rnai_demeter2']['selectivity']:.3f}; "
  f"p = {V7['rnai_demeter2']['p_one_sided']:.3f} but q = {P['V7_with_q_criterion']['rnai_demeter2']['q_value']:.2f}",
  "p < 0.05 with negative selectivity (the manuscript uses q < 0.05)", "pass on p, FAIL on q"),
 ("V8", "Expression of SNAP23 in bowel lines",
  f"median log2(TPM + 1) = {V8['median_log2TPM1_bowel']:.2f} (other lineages {V8['median_log2TPM1_other']:.2f})",
  "median > 1", "pass"),
]
checks = pd.DataFrame(rows, columns=["Check", "Description", "Result", "Pre-specified criterion", "Outcome"])
r = pd.read_csv(O + "/snap23_codependency_all_genes.csv", index_col=0).iloc[:, 0].sort_values(ascending=False).head(50)
cod = pd.DataFrame({"Gene": [g.split(" ")[0] for g in r.index], "Entrez_entry": r.index,
                    "Pearson_r_after_lineage_centering": r.values.round(4), "Rank": range(1, 51)})
mod = pd.DataFrame(P["module_partners_bowel_selectivity_posthoc"]).T.reset_index().rename(columns={"index": "Gene"})
ph = pd.DataFrame([{"Analysis": k, **v} for k, v in P["V5_within_lineage_posthoc"].items()])
with pd.ExcelWriter(O + "/S11_Table.xlsx") as w:
    checks.to_excel(w, sheet_name="Pre-specified checks", index=False)
    cod.to_excel(w, sheet_name="Top 50 co-dependencies", index=False)
    mod.to_excel(w, sheet_name="SNARE module selectivity", index=False)
    ph.to_excel(w, sheet_name="Paralogue corr by lineage", index=False)
print("S11 written")
