"""Follow-up of FERMT2 in CNS/brain (post hoc; S14 Table, sheet 'Follow-up'), added after the next-generation comparison.
(1) Culture-format-matched comparisons in the next-generation screens, using per-screen gene effects (screen_gene_effect.csv and
    screen_metadata.csv, DepMap portal release "NextGen Model Manuscript 2026"): CNS/brain coated-2D (2DN) vs other 2D (2DO) screens,
    and CNS/brain spheroid (3DN) vs other organoid (3DO) screens; the same organoid-only (3DO vs 3DO) contrast for all 17 pairs of S13 Table
    as a calibration, and the 3D-organoid vs 2D-adherent difference of FERMT2 reported by Neiswender et al.
(2) RNAi (DEMETER2): CNS/brain split into diffuse gliomas and embryonal tumours, and CRISPR-RNAi concordance of FERMT2.
(3) Project Score: the CNS/brain contrast in the Sanger cohort of the separated design (sanger_disjoint_primary94.csv).
Environment: DEPMAP_DIR (26Q1 files), RESULTS_DIR, NEXTGEN_SCREEN_DIR, RNAI_DIR."""
import os, json, hashlib
import numpy as np, pandas as pd
from scipy import stats

DATA = os.environ.get("DEPMAP_DIR", "."); OUT = os.environ.get("RESULTS_DIR", "results")
NGS = os.environ["NEXTGEN_SCREEN_DIR"]; RD = os.environ["RNAI_DIR"]
G = "FERMT2 (10979)"
sha = lambda f: hashlib.sha256(open(f, "rb").read()).hexdigest()
mw = lambda x, y: float(stats.mannwhitneyu(x, y, alternative="less").pvalue)
res = {"gene": G, "status": "post hoc; chosen after the next-generation, Project Score and RNAi results were known"}

# ---- (1) next-generation screens, per screen
fe, fm = os.path.join(NGS, "screen_gene_effect.csv"), os.path.join(NGS, "screen_metadata.csv")
res["nextgen_screen_files_sha256"] = {"screen_gene_effect.csv": sha(fe), "screen_metadata.csv": sha(fm)}
P = pd.read_csv(os.path.join(OUT, "nextgen_check_pairs.csv"))
hdr = pd.read_csv(fe, nrows=0).columns
E = pd.read_csv(fe, usecols=[hdr[0]] + sorted(set(P.Gene))).set_index(hdr[0])
S = pd.read_csv(fm).set_index("ScreenID")
S = S[(S.PassesQC == True) & (S.IsNextGen == True)].join(E, how="inner")
LM = {"bowel": "Bowel", "pancreas": "Pancreas", "esophagus/stomach": "Esophagus/Stomach", "breast": "Breast", "prostate": "Prostate", "cns/brain": "CNS/Brain"}
# consistency with the published lineage summaries (figshare next_gen_dependency_lineage_tests.csv)
dev = [abs(S[S.OncotreeLineage == LM[r.Lineage]][r.Gene].mean() - r.NextGen_target_mean) for _, r in P.iterrows()]
nn = [int((S.OncotreeLineage == LM[r.Lineage]).sum()) == r.NextGen_target_models for _, r in P.iterrows()]
res["consistency_with_published_summaries"] = {"max_abs_mean_difference": float(max(dev)), "all_counts_match": bool(all(nn)), "n_screens": int(len(S))}
cns = S.OncotreeLineage == "CNS/Brain"
def contrast(x, y):
    return {"n_target": int(len(x)), "n_other": int(len(y)), "median_target": float(x.median()), "median_other": float(y.median()),
            "median_difference": float(x.median() - y.median()), "mean_target": float(x.mean()), "mean_other": float(y.mean()), "p_one_sided": mw(x, y)}
fmt = {}
for a, b in [("2DN", "2DO"), ("3DN", "3DO")]:
    fmt[f"{a}_vs_{b}"] = contrast(S[cns & (S.ScreenType == a)][G].dropna(), S[~cns & (S.ScreenType == b)][G].dropna())
x3 = S[cns & (S.ScreenType == "3DN")][G].dropna().values; y3 = S[~cns & (S.ScreenType == "3DO")][G].dropna().values
rng = np.random.default_rng(1); bs = [np.median(rng.choice(x3, len(x3))) - np.median(rng.choice(y3, len(y3))) for _ in range(10000)]
fmt["3DN_vs_3DO"]["bootstrap_95ci_median_difference"] = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
x2, y2 = S[cns & (S.ScreenType == "3DN")][G], S[cns & (S.ScreenType == "2DN")][G]
fmt["CNS_3DN_vs_2DN"] = {"median_3DN": float(x2.median()), "median_2DN": float(y2.median()), "p_two_sided": float(stats.mannwhitneyu(x2, y2).pvalue),
                         "models_screened_in_both_formats": int(S[cns].groupby("ModelID").ScreenType.nunique().gt(1).sum())}
ov = pd.read_csv(os.path.join(OUT, "nextgen_figshare", "organoids_vs_2d_gene_dependency_tests_full.csv")).set_index("Dependency").loc[G]
fmt["organoid_vs_2D_adherent_published"] = {"organoid_mean": float(ov.OrganoidsMean), "n_organoid": int(ov.OrganoidsN), "adherent_2D_mean": float(ov.AdherentMean),
                                            "n_adherent": int(ov.AdherentN), "fdr": float(ov.MannWhitneyFDR)}
res["nextgen_format_matched"] = fmt
cal = []
for _, r in P.iterrows():
    lin = LM[r.Lineage]; t = S.OncotreeLineage == lin
    a_t, b_t = ("3DN", "3DO") if lin == "CNS/Brain" else ("3DO", "3DO")
    x, y = S[t & (S.ScreenType == a_t)][r.Gene].dropna(), S[~t & (S.ScreenType == b_t)][r.Gene].dropna()
    if len(x) < 5: continue
    c = contrast(x, y); c.update(Gene=r.Gene, Lineage=r.Lineage, target_format=a_t, reproduced_in_S13=bool(r["Reproduced_q_lt_0.05"])); cal.append(c)
cal = pd.DataFrame(cal)[["Gene", "Lineage", "target_format", "reproduced_in_S13", "n_target", "n_other", "median_target", "median_other", "median_difference", "p_one_sided"]]
cal.to_csv(os.path.join(OUT, "nextgen_3d_only_pairs.csv"), index=False)
rep = cal[cal.reproduced_in_S13]
res["three_d_only_calibration"] = {"pairs_testable": int(len(cal)), "reproduced_pairs": int(len(rep)), "reproduced_with_p_lt_0.05_3D_only": int((rep.p_one_sided < 0.05).sum()),
                                   "reproduced_pairs_p": {f"{g.split(' ')[0]} {l}": float(p) for g, l, p in zip(rep.Gene, rep.Lineage, rep.p_one_sided)}}

# ---- (2) RNAi
M = pd.read_csv(os.path.join(DATA, "Model.csv")).set_index("ModelID")
R = pd.read_csv(os.path.join(RD, "D2_combined_gene_dep_scores.csv"), index_col=0)
r = R.loc[[i for i in R.index if i.startswith("FERMT2 ")][0]].dropna()
mp = pd.read_csv(os.path.join(OUT, "rnai_cell_line_mapping.csv")).set_index("CCLE_name").DepMap_ID
r.index = r.index.map(mp); r = r[r.index.notna()]; r = r[~r.index.duplicated()]
d = pd.DataFrame({"rnai": r}).join(M[["OncotreeLineage", "OncotreePrimaryDisease"]], how="left")
d = d[d.OncotreeLineage.notna() & (d.OncotreePrimaryDisease != "Non-Cancerous")]
c_ = d.OncotreeLineage.str.strip().str.lower() == "cns/brain"; gl = c_ & d.OncotreePrimaryDisease.str.contains("Glioma")
rest = d[~c_].rnai
rn = {k: {"n": int(m.sum()), "median": float(d[m].rnai.median()), "median_rest": float(rest.median()),
          "selectivity": float(d[m].rnai.median() - rest.median()), "p_one_sided": mw(d[m].rnai, rest)}
      for k, m in [("cns_brain_all", c_), ("diffuse_glioma", gl), ("embryonal_and_other", c_ & ~gl)]}
ref = pd.read_csv(os.path.join(OUT, "rnai_validation.csv")).query("Gene == @G and Lineage == 'cns/brain'").iloc[0]
assert abs(ref.p_value_RNAi - rn["cns_brain_all"]["p_one_sided"]) < 1e-12   # same procedure as rnai_validation.py
rn["cns_brain_all"]["q_value_primary_family"] = float(ref.q_value_RNAi)
hdr = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), nrows=0).columns
cr = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), usecols=[hdr[0], G]).set_index(hdr[0])[G]
j = d.join(cr.rename("crispr"), how="inner").dropna()
rho = stats.spearmanr(j.crispr, j.rnai)
rn["crispr_rnai_concordance"] = {"n_lines": int(len(j)), "spearman_rho": float(rho.statistic), "p": float(rho.pvalue)}
res["rnai"] = rn

# ---- (3) Project Score, separated design
s = pd.read_csv(os.path.join(OUT, "sanger_disjoint_primary94.csv")).query("Gene == @G and Lineage == 'cns/brain'").iloc[0]
res["project_score_separated_design"] = {"discovery_selectivity": float(s.Selectivity_disc), "discovery_q": float(s.q_value_disc), "discovery_candidate": bool(s.candidate_disc),
                                         "sanger_n": int(s.n_target_S), "sanger_selectivity": float(s.Selectivity_S), "sanger_median": float(s.Chronos_median_S),
                                         "sanger_q": float(s.q_value_S)}
json.dump(res, open(os.path.join(OUT, "fermt2_followup_summary.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
