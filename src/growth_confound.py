"""Culture-format check of the prioritised pairs.  (a) OLS of gene effect on lineage membership with and without growth pattern (adherent / suspension / mixed / unknown) as covariate, for blood pairs also without the other blood lineage;
(b) for lymphoid and myeloid pairs, target lineage compared only with non-blood suspension lines; (c) exclusion of drug-adapted and engineered models.  Growth pattern from DepMap Model.csv."""
import os, json, numpy as np, pandas as pd, statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from scipy.stats import mannwhitneyu
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "primary") + "/"
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
m = pd.read_csv(os.path.join(D, "Model.csv"))[["ModelID", "GrowthPattern", "EngineeredModel"]]; df = df.merge(m, left_on="DepMap_ID", right_on="ModelID", how="left", validate="one_to_one")
F = pd.read_csv(O + "final_targets.csv"); gp = pd.get_dummies(df.GrowthPattern.fillna("Unknown")).astype(float); gp = gp[[c for c in gp.columns if c != "Adherent"]]
susp_nonblood = (df.GrowthPattern == "Suspension") & ~df.lineage.isin(["lymphoid", "myeloid"]); eng = df.EngineeredModel.notna(); rows = []
for _, r in F.iterrows():
    y = pd.to_numeric(df[r.Gene], errors="coerce"); ok = y.notna().values; t = (df.lineage == r.Lineage).astype(float).values
    X0 = sm.add_constant(t[ok]); b0 = sm.OLS(y[ok].values, X0).fit(); X1 = sm.add_constant(np.column_stack([t[ok], gp.values[ok]])); b1 = sm.OLS(y[ok].values, X1).fit()
    row = dict(Gene=r.Gene, Lineage=r.Lineage, coef_unadjusted=b0.params[1], coef_growth_adjusted=b1.params[1], p_growth_adjusted_two_sided=b1.pvalues[1])
    if r.Lineage in ("lymphoid", "myeloid"):
        # same fits without the other blood lineage (otherwise the suspension term also absorbs dependencies shared by the two blood lineages)
        k = ok & (df.lineage != ({"lymphoid": "myeloid", "myeloid": "lymphoid"}[r.Lineage])).values; Xg = gp.values[k]; Xg = Xg[:, Xg.std(0) > 0]
        c0 = sm.OLS(y[k].values, sm.add_constant(t[k])).fit(); c1 = sm.OLS(y[k].values, sm.add_constant(np.column_stack([t[k], Xg]))).fit()
        row.update(coef_unadjusted_excl_other_blood=c0.params[1], coef_growth_adjusted_excl_other_blood=c1.params[1], p_growth_adjusted_excl_other_blood=c1.pvalues[1])
        tg = y[ok & (df.lineage == r.Lineage).values]; cg = y[ok & susp_nonblood.values]; row.update(n_suspension_comparator=len(cg), median_diff_vs_suspension=float(tg.median() - cg.median()) if len(cg) >= 10 else np.nan,
            p_vs_suspension_one_sided=float(mannwhitneyu(tg, cg, alternative="less")[1]) if len(cg) >= 10 else np.nan)
    keep = ok & ~eng.values; tg = y[keep & (df.lineage == r.Lineage).values]; og = y[keep & (df.lineage != r.Lineage).values]
    row.update(median_diff_excl_engineered=float(tg.median() - og.median()), p_excl_engineered_one_sided=float(mannwhitneyu(tg, og, alternative="less")[1]), n_target_excl_engineered=len(tg)); rows.append(row)
G = pd.DataFrame(rows); G["ratio_excl_other_blood"] = G.coef_growth_adjusted_excl_other_blood / G.coef_unadjusted_excl_other_blood
hm = G.Lineage.isin(["lymphoid", "myeloid"]); G.loc[hm, "q_growth_adjusted_excl_other_blood"] = multipletests(G.loc[hm, "p_growth_adjusted_excl_other_blood"], method="fdr_bh")[1]; heme = G[hm]
shared = heme.Gene.map(heme.Gene.value_counts()) == 2
S = dict(n_pairs=len(G), growth_pattern_counts=df.GrowthPattern.fillna("Unknown").value_counts().to_dict(), n_engineered_or_drug_adapted=int(eng.sum()),
    adj_coef_lt_minus0_3=int((G.coef_growth_adjusted < -0.3).sum()), adj_p_lt_0_05=int((G.p_growth_adjusted_two_sided < 0.05).sum()), adj_p_lt_0_05_and_negative=int(((G.p_growth_adjusted_two_sided < 0.05) & (G.coef_growth_adjusted < 0)).sum()),
    median_ratio_adjusted_to_unadjusted=float((G.coef_growth_adjusted / G.coef_unadjusted).median()),
    heme_pairs=len(heme), n_susp_nonblood=int(susp_nonblood.sum()), heme_vs_suspension_sel_lt_minus0_5=int((heme.median_diff_vs_suspension < -0.5).sum()), heme_vs_suspension_sel_lt_0=int((heme.median_diff_vs_suspension < 0).sum()),
    heme_vs_suspension_p_lt_0_05=int((heme.p_vs_suspension_one_sided < 0.05).sum()), heme_vs_suspension_median_diff=float(heme.median_diff_vs_suspension.median()),
    heme_median_ratio=float((heme.coef_growth_adjusted / heme.coef_unadjusted).median()), heme_median_ratio_excl_other_blood=float(heme.ratio_excl_other_blood.median()),
    heme_median_ratio_shared_genes=float((heme.coef_growth_adjusted / heme.coef_unadjusted)[shared].median()), heme_median_ratio_excl_other_blood_shared_genes=float(heme.ratio_excl_other_blood[shared].median()),
    heme_median_ratio_other_genes=float((heme.coef_growth_adjusted / heme.coef_unadjusted)[~shared].median()), heme_median_ratio_excl_other_blood_other_genes=float(heme.ratio_excl_other_blood[~shared].median()),
    heme_n_shared_pairs=int(shared.sum()), heme_ratio_lt_0_5=int(((heme.coef_growth_adjusted / heme.coef_unadjusted) < 0.5).sum()), heme_ratio_lt_0_5_excl_other_blood=int((heme.ratio_excl_other_blood < 0.5).sum()),
    heme_q_lt_0_05_excl_other_blood=int((heme.q_growth_adjusted_excl_other_blood < 0.05).sum()),
    excl_eng_sel_lt_minus0_5=int((G.median_diff_excl_engineered < -0.5).sum()), excl_eng_p_lt_0_05=int((G.p_excl_engineered_one_sided < 0.05).sum()))
G.to_csv(O + "growth_confound_pairs.csv", index=False); json.dump(S, open(O + "growth_confound_summary.json", "w"), indent=1, default=str); print(json.dumps(S, indent=1, default=str))
print(G[G.Lineage.isin(["lymphoid","myeloid"])&(G.median_diff_vs_suspension>-0.5)][["Gene","Lineage","median_diff_vs_suspension","p_vs_suspension_one_sided"]].to_string())
print(G[G.coef_growth_adjusted>-0.3][["Gene","Lineage","coef_unadjusted","coef_growth_adjusted","p_growth_adjusted_two_sided"]].to_string())
