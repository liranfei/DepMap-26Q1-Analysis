"""Supporting tables S6 (26Q1 candidates in 22Q1), S9 (alternative effect-size measures) and S10 (growth pattern and model type) as xlsx.
Run after cross_release.py, benchmark_metrics.py, growth_confound.py and revision3_checks.py (which adds the BH q-values of the growth-adjusted models)."""
import os, pandas as pd, numpy as np
from statsmodels.stats.multitest import multipletests
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; OUT = os.environ.get("FIG_DIR", "figures") + "/supporting_information/"; os.makedirs(OUT, exist_ok=True)
sym = lambda s: s.str.replace(r" \(.*", "", regex=True)

# ---- S6
if os.path.exists(O + "cross_release_candidates.csv"):
    C = pd.read_csv(O + "cross_release_candidates.csv")
    S6 = C[["Gene", "Lineage", "n_target_26", "Chronos_median_26", "Selectivity_26", "q_value_26", "n_target_22", "Chronos_median_22", "Selectivity_22", "q_value_22"]].copy()
    S6.columns = ["Gene", "Lineage", "n target, 26Q1", "Median Chronos, 26Q1", "Selectivity, 26Q1", "q, 26Q1", "n target, 22Q1", "Median Chronos, 22Q1", "Selectivity, 22Q1", "q, 22Q1"]
    with pd.ExcelWriter(OUT + "S6_Table.xlsx") as w: S6.to_excel(w, index=False, sheet_name="26Q1 candidates in 22Q1")

# ---- S9
F = pd.read_csv(O + "final_targets.csv"); BM = pd.read_csv(O + "benchmark_metrics_all_tests.csv.gz"); T = pd.read_csv(O + "all_tests.csv")
S9 = F[["Gene", "Lineage", "n_target", "Chronos_median", "Selectivity"]].merge(BM, on=["Gene", "Lineage"], how="left").merge(T[["Gene", "Lineage", "q_value"]], on=["Gene", "Lineage"], how="left")
S9.insert(0, "Gene_symbol", sym(S9.Gene))
S9 = S9[["Gene_symbol", "Gene", "Lineage", "n_target", "Chronos_median", "Selectivity", "mean_diff", "q_value", "p_t2", "q_t", "cohens_d", "auc_lower"]]
S9.columns = ["Gene_symbol", "Gene", "Lineage", "n_target", "Median Chronos", "Selectivity (diff. medians)", "Difference in means", "q (MWU, BH)", "p (two-tailed t)", "q (BH, t-test, full family)", "Cohen's d", "AUC (target lower)"]
with pd.ExcelWriter(OUT + "S9_Table.xlsx") as w: S9.to_excel(w, index=False, sheet_name="Alternative measures")

# ---- S10 (post-selection diagnostics of the 94 prioritised pairs)
G = pd.read_csv(O + "growth_confound_pairs.csv")
if "q_growth_adjusted" not in G: G["q_growth_adjusted"] = multipletests(G.p_growth_adjusted_two_sided, method="fdr_bh")[1]
ok = G.p_vs_suspension_one_sided.notna(); G["q_vs_suspension"] = np.nan; G.loc[ok, "q_vs_suspension"] = multipletests(G.loc[ok, "p_vs_suspension_one_sided"], method="fdr_bh")[1]
G["q_excl_engineered"] = multipletests(G.p_excl_engineered_one_sided, method="fdr_bh")[1]
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
eng = df.DepMap_ID.map(pd.read_csv(os.path.join(D, "Model.csv")).set_index("ModelID").EngineeredModel).notna().values
G["median_target_excl_engineered"] = [float(np.nanmedian(df.loc[(df.lineage == l).values & ~eng, g].astype(float))) for g, l in zip(G.Gene, G.Lineage)]
cols = ["Gene", "Lineage", "coef_unadjusted", "coef_growth_adjusted", "p_growth_adjusted_two_sided", "q_growth_adjusted", "median_target_excl_engineered", "median_diff_excl_engineered", "p_excl_engineered_one_sided", "q_excl_engineered",
        "n_target_excl_engineered", "n_suspension_comparator", "median_diff_vs_suspension", "p_vs_suspension_one_sided", "q_vs_suspension",
        "coef_unadjusted_excl_other_blood", "coef_growth_adjusted_excl_other_blood", "p_growth_adjusted_excl_other_blood", "q_growth_adjusted_excl_other_blood"]
S10 = G[cols].copy(); S10["Gene_symbol"] = sym(S10.Gene)
note = pd.DataFrame({"Note": ["Post-selection diagnostics of the 94 prioritised pairs; q-values are Benjamini-Hochberg adjusted within the 94 pairs (q_vs_suspension: within the 39 testable lymphoid/myeloid pairs) and are conditional on the selection.",
                              "Growth pattern (Model.csv GrowthPattern: adherent, suspension, mixed, unknown) added as a categorical covariate; engineered or drug-adapted models = EngineeredModel annotated in Model.csv (10 models).",
                              "Suspension comparison: lymphoid/myeloid target lineage versus non-blood cell lines annotated as suspension cultures (41 lines); at least ten comparator lines with data required.",
                              "*_excl_other_blood: lymphoid/myeloid pairs refitted without the other blood lineage (otherwise the suspension term also absorbs dependencies shared by the two blood lineages); q within the 42 blood pairs."]})
with pd.ExcelWriter(OUT + "S10_Table.xlsx") as w: note.to_excel(w, index=False, sheet_name="Read me"); S10.to_excel(w, index=False, sheet_name="Growth pattern and model type")
print("S6, S9, S10 written")
