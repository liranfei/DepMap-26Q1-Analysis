"""S14 Table: pre-specified checks of FERMT2 (fermt2_validation.py; plan in prespecified_plan_fermt2.md) and post hoc follow-up (fermt2_followup.py).
Environment: RESULTS_DIR."""
import os, json
import pandas as pd
O = os.environ.get("RESULTS_DIR", "results")
R = json.load(open(os.path.join(O, "fermt2_validation_summary.json")))
rows = []
for lin in ["CNS/Brain", "Skin", "Kidney"]:
    x = R[lin]; role = "primary" if lin == "CNS/Brain" else "secondary"
    rows += [
        (lin, role, "0 (reference)", "Primary 2D comparison with all other cancer lines", f"selectivity {x['selectivity_all']:.3f}; p = {x['p_all']:.2g}; n = {x['n_target']}", "-", "-"),
        (lin, role, "1", "Culture format: other lineages, adherent lines only", f"selectivity {x['c1_adherent']['selectivity']:.3f}; p = {x['c1_adherent']['p']:.2g}; n other = {x['c1_adherent']['n_other']}",
         "selectivity <= half of primary and p < 0.05", "pass" if x["c1_pass"] else "fail"),
        (lin, role, "2", "Genotype: TP53, PTEN, NF1 (high/moderate), IDH1, EGFR (hotspot) as covariates",
         f"coefficient {x['c2_genotype']['coef_unadj']:.3f} -> {x['c2_genotype']['coef_adj']:.3f} ({100*x['c2_genotype']['ratio']:.0f}%); p = {x['c2_genotype']['p_adj']:.2g}",
         "ratio >= 0.70 and p < 0.05", "pass" if x["c2_pass"] else "fail"),
        (lin, role, "3", "Disease subtypes with >= 5 lines", "; ".join(f"{k}: n = {v['n']}, selectivity {v['selectivity']:.3f}, p = {v['p']:.2g}" for k, v in x["c3_subtypes"].items()),
         "selectivity < 0 and p < 0.05 in >= 2/3 of subtypes", "pass" if x["c3_pass"] else "fail"),
        (lin, role, "4", "Carcinoma lineages only (the contrast available in the NextGen screens)",
         f"selectivity {x['c4_carcinoma']['selectivity']:.3f}; p = {x['c4_carcinoma']['p']:.2g}; n other = {x['c4_carcinoma']['n_other']}",
         "selectivity < 0 and p < 0.05", "pass" if x["c4_pass"] else "fail")]
P = R["c5_paralogue"]
rows.append(("all lines", "mechanism", "5", "Paralogue buffering (direction fixed in advance): FERMT1 expression vs FERMT2 gene effect",
             f"FERMT1 rho = {P['FERMT1']['rho']:.3f} (p = {P['FERMT1']['p']:.2g}); within lineages rho = {P['FERMT1']['rho_within_lineage']:.3f} (p = {P['FERMT1']['p_within_lineage']:.2g}); FERMT3 rho = {P['FERMT3']['rho']:.3f}, within lineages {P['FERMT3']['rho_within_lineage']:.3f} (p = {P['FERMT3']['p_within_lineage']:.2g}); FERMT1 median log2(TPM+1): " + ", ".join(f"{k} {v:.2f}" for k, v in R["c5_FERMT1_median"].items()),
             "rho > 0 and p < 0.05, across and within lineages", "pass" if R["c5_pass_FERMT1"] else "fail"))
rows.append(("all lines", "mechanism", "6", "Co-dependency with the pre-specified integrin-adhesion set (13 genes), lineage-centred",
             "partner ranks: " + ", ".join(f"{k} {v}" for k, v in sorted(R["c6_partner_ranks"].items(), key=lambda kv: kv[1])),
             ">= 1 partner in the top 20", "pass" if R["c6_pass"] else "fail"))
checks = pd.DataFrame(rows, columns=["Lineage", "Role", "Check", "Description", "Result", "Pre-specified criterion", "Outcome"])
dec = pd.DataFrame([("Decision rule", "Use FERMT2 (CNS/brain) as the example only if checks 1, 2 and 4 pass for CNS/brain"),
                    ("Decision", "met" if R["decision_use_as_example"] else "not met"),
                    ("Plan SHA-256", R["plan_sha256"]),
                    ("Known before the plan (not pre-specified)", "26Q1 2D selectivity and robustness tiers; Project Score result; RNAi (CNS/brain q = 0.10, skin q = 7e-5, kidney q = 3e-4); NextGen CNS/brain vs other NextGen models (mean difference -0.263, q = 4e-4 over 17 pairs); PubMed 0-1 records per pair")],
                   columns=["Item", "Value"])
top = pd.DataFrame(list(R["c6_top20"].items()), columns=["Gene", "Pearson_r_after_lineage_centring"]); top.insert(0, "Rank", range(1, len(top) + 1))
F = json.load(open(os.path.join(O, "fermt2_followup_summary.json"))); fm = F["nextgen_format_matched"]; rn = F["rnai"]; ps = F["project_score_separated_design"]
def ct(c): return f"median {c['median_target']:.3f} (n = {c['n_target']}) vs {c['median_other']:.3f} (n = {c['n_other']}); difference {c['median_difference']:.3f}; one-sided Mann-Whitney p = {c['p_one_sided']:.2g}"
o = fm["organoid_vs_2D_adherent_published"]; c3 = fm["3DN_vs_3DO"]; cal = F["three_d_only_calibration"]
fu = pd.DataFrame([
    ("NextGen, same culture format", "CNS/brain coated-2D (2DN) vs other 2D (2DO) screens", ct(fm["2DN_vs_2DO"])),
    ("NextGen, same culture format", "CNS/brain spheroids (3DN) vs other organoids (3DO)", ct(c3) + f"; bootstrap 95% CI of the median difference {c3['bootstrap_95ci_median_difference'][0]:.2f} to {c3['bootstrap_95ci_median_difference'][1]:.2f}"),
    ("NextGen, culture format", "CNS/brain 3DN vs 2DN (different models)", f"median {fm['CNS_3DN_vs_2DN']['median_3DN']:.3f} vs {fm['CNS_3DN_vs_2DN']['median_2DN']:.3f}; two-sided p = {fm['CNS_3DN_vs_2DN']['p_two_sided']:.2g}; models screened in both formats: {fm['CNS_3DN_vs_2DN']['models_screened_in_both_formats']}"),
    ("NextGen, culture format", "3D organoids vs lineage-matched 2D adherent lines (as published by Neiswender et al.)", f"mean {o['organoid_mean']:.3f} (n = {o['n_organoid']}) vs {o['adherent_2D_mean']:.3f} (n = {o['n_adherent']}); FDR = {o['fdr']:.2g}"),
    ("NextGen, calibration", "Organoid-only contrast (3DO vs 3DO; 3DN for CNS/brain) for the 7 pairs reproduced in S13 Table", f"{cal['reproduced_with_p_lt_0.05_3D_only']} of {cal['reproduced_pairs']} with p < 0.05: " + "; ".join(f"{k} p = {v:.2g}" for k, v in cal["reproduced_pairs_p"].items()) + " (all 17 pairs: nextgen_3d_only_pairs.csv)"),
    ("RNAi (DEMETER2)", "CNS/brain, all lines (as in the RNAi check)", f"n = {rn['cns_brain_all']['n']}; selectivity {rn['cns_brain_all']['selectivity']:.3f}; p = {rn['cns_brain_all']['p_one_sided']:.2g}; q = {rn['cns_brain_all']['q_value_primary_family']:.2g}"),
    ("RNAi (DEMETER2)", "CNS/brain diffuse gliomas only", f"n = {rn['diffuse_glioma']['n']}; selectivity {rn['diffuse_glioma']['selectivity']:.3f}; p = {rn['diffuse_glioma']['p_one_sided']:.2g}"),
    ("RNAi (DEMETER2)", "CNS/brain embryonal and other tumours", f"n = {rn['embryonal_and_other']['n']}; selectivity {rn['embryonal_and_other']['selectivity']:.3f}; p = {rn['embryonal_and_other']['p_one_sided']:.2g}"),
    ("RNAi (DEMETER2)", "CRISPR-RNAi concordance of FERMT2 across shared lines", f"Spearman rho = {rn['crispr_rnai_concordance']['spearman_rho']:.2f}; n = {rn['crispr_rnai_concordance']['n_lines']}"),
    ("Project Score, separated design", "CNS/brain, discovery cohort (not a discovery candidate)", f"selectivity {ps['discovery_selectivity']:.3f}; q = {ps['discovery_q']:.2g}"),
    ("Project Score, separated design", "CNS/brain, Sanger cohort", f"n = {ps['sanger_n']}; selectivity {ps['sanger_selectivity']:.3f}; median {ps['sanger_median']:.3f}; q = {ps['sanger_q']:.2g}"),
    ("Note", "All follow-up analyses", "Post hoc, done after FERMT2 had been chosen; no multiplicity adjustment across them. Per-screen NextGen files: DepMap portal release 'NextGen Model Manuscript 2026' (screen_gene_effect.csv, screen_metadata.csv; SHA-256 in fermt2_followup_summary.json); per-lineage means reproduce the published summaries.")],
    columns=["Data", "Comparison", "Result"])
with pd.ExcelWriter(os.path.join(O, "S14_Table.xlsx")) as w:
    checks.to_excel(w, sheet_name="Checks", index=False)
    dec.to_excel(w, sheet_name="Plan and decision", index=False)
    top.to_excel(w, sheet_name="Top 20 co-dependencies", index=False)
    fu.to_excel(w, sheet_name="Follow-up (post hoc)", index=False)
    pd.read_csv(os.path.join(O, "nextgen_3d_only_pairs.csv")).to_excel(w, sheet_name="NextGen organoid-only, 17 pairs", index=False)
print("S14 written", len(checks))
