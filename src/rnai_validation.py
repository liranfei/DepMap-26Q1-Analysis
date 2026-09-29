"""Orthogonal check of the prioritised pairs in RNAi screens (DEMETER2 combined: Achilles, DRIVE and Marcotte shRNA screens; figshare 10.6084/m9.figshare.6025238).

Rules fixed before the analysis was run:
- cell lines mapped to DepMap ModelIDs through the CCLE name in the 26Q1 Model.csv (fallback: stripped cell-line name, only when unique);
  lineage and non-cancerous status from the 26Q1 model annotation; genes matched on the Entrez identifier;
- the same procedure as the primary analysis (one-sided Mann-Whitney, lineage versus all other cancer lines, at least five non-missing target lines,
  Benjamini-Hochberg over all RNAi tests);
- a prioritised pair is confirmed when q < 0.05 with negative selectivity in the RNAi data; the effect-size rule is reported descriptively only,
  because DEMETER2 scores are not scaled like Chronos scores;
- reference: 26Q1 pairs with q < 0.05 and negative selectivity that were not prioritised; results also by confidence tier and lymphoid/other lineages.
The RNAi screens share many cell lines with the CRISPR data, so this is an orthogonal check of the perturbation method, not a check in independent cell lines.
Input: RNAI_DIR with D2_combined_gene_dep_scores.csv (sample_info.csv is not needed; cell lines are mapped with the 26Q1 Model.csv).  Output: rnai_validation.csv, rnai_all_tests.csv.gz, rnai_validation_summary.json, rnai_cell_line_mapping.csv, rnai_lineage_counts.csv."""
import os, json, re, numpy as np, pandas as pd
from scipy.stats import beta
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; RD = os.environ["RNAI_DIR"]
ci = lambda k, n: [float(beta.ppf(0.025, k, n - k + 1)) if k > 0 else 0.0, float(beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0]
M = pd.read_csv(os.path.join(D, "Model.csv"))
R2 = pd.read_csv(os.path.join(RD, "D2_combined_gene_dep_scores.csv"), index_col=0).T          # lines x genes
by_ccle = M.dropna(subset=["CCLEName"]).drop_duplicates("CCLEName").set_index("CCLEName").ModelID
strip = M.dropna(subset=["StrippedCellLineName"]); strip = strip[~strip.StrippedCellLineName.duplicated(keep=False)].set_index("StrippedCellLineName").ModelID
ids = []
for c in R2.index:
    if c in by_ccle.index: ids.append(by_ccle[c]); continue
    s = c.split("_")[0]; ids.append(strip[s] if s in strip.index else None)
R2["DepMap_ID"] = ids; n_lines = len(R2); unmapped = [c for c, i in zip(R2.index, ids) if i is None]
pd.DataFrame({"CCLE_name": R2.index, "DepMap_ID": ids}).to_csv(O + "rnai_cell_line_mapping.csv", index=False)   # mapping list (unmapped lines have no ID)
R2 = R2[R2.DepMap_ID.notna()].drop_duplicates("DepMap_ID")
meta = M.set_index("ModelID")[["OncotreeLineage", "OncotreePrimaryDisease"]]
R2["lineage"] = R2.DepMap_ID.map(meta.OncotreeLineage).str.strip().str.lower(); R2["primary_disease"] = R2.DepMap_ID.map(meta.OncotreePrimaryDisease)
R2 = R2[R2.lineage.notna() & (R2.primary_disease != "Non-Cancerous")].reset_index(drop=True)
R2.lineage.value_counts().rename_axis("lineage").rename("n_lines").to_csv(O + "rnai_lineage_counts.csv")
genes_all = [c for c in R2.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
def ent(g):
    m = re.fullmatch(r".* \((\d+)\)", g); return m.group(1) if m else None   # entries mapped to several genes ("A&B") have no single Entrez id
genes = [g for g in genes_all if ent(g) is not None]                         # multi-gene entries excluded before testing, so BH runs over single-gene tests only
T = rp.test_all(R2, genes, R2.lineage.values); T.to_csv(O + "rnai_all_tests.csv.gz", index=False)
T["entrez"] = T.Gene.map(ent)
F = pd.read_csv(O + "final_targets.csv"); F["entrez"] = F.Gene.map(ent)
tiers = pd.read_csv(O + "confidence_tiers.csv")[["Gene", "Lineage", "tier"]]
P = F[["Gene", "Lineage", "entrez", "Selectivity"]].merge(T.drop(columns="Gene").rename(columns=lambda c: c + "_RNAi" if c not in ("entrez", "Lineage") else c), on=["entrez", "Lineage"], how="left")
P = P.merge(tiers, on=["Gene", "Lineage"]); P["testable"] = P.q_value_RNAi.notna()
P["confirmed"] = P.testable & (P.q_value_RNAi < 0.05) & (P.Selectivity_RNAi < 0)
P["meets_rule_RNAi"] = P.testable & (P.q_value_RNAi < 0.05) & (P.Chronos_median_RNAi < -1) & (P.Selectivity_RNAi < -0.5)
P.to_csv(O + "rnai_validation.csv", index=False)
# reference: significant, negative, non-prioritised 26Q1 pairs
A = pd.read_csv(O + "all_tests.csv.gz"); A["entrez"] = A.Gene.map(ent); key = set(zip(F.Gene, F.Lineage))
ref = A[(A.q_value < 0.05) & (A.Selectivity < 0) & ~pd.Series([k in key for k in zip(A.Gene, A.Lineage)], index=A.index)]
ref = ref[["entrez", "Lineage"]].merge(T[["entrez", "Lineage", "q_value", "Selectivity"]], on=["entrez", "Lineage"])
ref_conf = int(((ref.q_value < 0.05) & (ref.Selectivity < 0)).sum())
tp = P[P.testable]
out = {"rnai_lines_in_file": n_lines, "unmapped": unmapped, "rnai_cancer_lines_analysed": len(R2), "lines_shared_with_26Q1_crispr": int(R2.DepMap_ID.isin(pd.read_csv(os.path.join(D, "gene_effect.csv"), usecols=[0]).iloc[:, 0]).sum()),
       "genes": len(genes), "multi_gene_entries_excluded": len(genes_all) - len(genes), "eligible_lineages": int((R2.lineage.value_counts() >= rp.MIN_N).sum()), "tests": len(T), "q_lt_0_05": int((T.q_value < 0.05).sum()),
       "testable": int(P.testable.sum()), "confirmed": int(P.confirmed.sum()), "confirmed_ci": ci(int(P.confirmed.sum()), int(P.testable.sum())),
       "negative_selectivity": int((tp.Selectivity_RNAi < 0).sum()), "meets_rule_descriptive": int(P.meets_rule_RNAi.sum()),
       "reference_confirmed": ref_conf, "reference_testable": len(ref),
       "fisher_high_vs_others": float(__import__("scipy.stats", fromlist=["fisher_exact"]).fisher_exact([[int(P[P.tier == "high"].confirmed.sum()), int((P[P.tier == "high"].testable & ~P[P.tier == "high"].confirmed).sum())],
                                                                                                  [int(P[P.tier != "high"].confirmed.sum()), int((P[P.tier != "high"].testable & ~P[P.tier != "high"].confirmed).sum())]])[1]),
       "by_tier": {t: {"confirmed": int(g.confirmed.sum()), "testable": int(g.testable.sum())} for t, g in P.groupby("tier")},
       "by_group": {k: {"confirmed": int(g.confirmed.sum()), "testable": int(g.testable.sum())} for k, g in P.groupby(P.Lineage.map(lambda l: "lymphoid" if l == "lymphoid" else "myeloid" if l == "myeloid" else "other"))},
       "not_testable": [f"{g.split(' (')[0]} {l}" for g, l in P[~P.testable][["Gene", "Lineage"]].values],
       "not_confirmed": [f"{g.split(' (')[0]} {l}" for g, l in P[P.testable & ~P.confirmed][["Gene", "Lineage"]].values],
       "spearman_selectivity_testable": float(tp[["Selectivity", "Selectivity_RNAi"]].corr("spearman").iloc[0, 1])}
json.dump(out, open(O + "rnai_validation_summary.json", "w"), indent=1, default=float); print(json.dumps(out, indent=1, default=float))
