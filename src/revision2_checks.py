"""Second-round checks of downstream analyses (revision 2).

1. Co-dependency after removing lineage structure: candidate genes were selected because they are low in the same lineages, so raw correlations across all
   cell lines are partly induced by the selection.  Gene effects are residualised on lineage (lineage means removed) and correlated again.
2. Over-representation with a background restricted to genes that could have met the median criterion (target-lineage median Chronos < -1 in at least one
   eligible lineage); the primary analysis used all tested genes, which favours broadly essential genes.
3. HNF1B and PAX8 in kidney restricted to renal cell carcinoma lines (the kidney lineage also contains rhabdoid tumour lines; TCGA-KIRC is clear cell RCC only).
4. Subsampling with the significance criterion of the primary analysis: p <= p*, the largest p-value with q < 0.05 in the full test family (instead of nominal p < 0.05).
Outputs: revision2_summary.json, codependency_residual_corr.csv, enrichment_restricted_background.csv, subsampling_power_bh.csv."""
import os, json, numpy as np, pandas as pd, requests
from scipy.stats import mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests
import run_pipeline as rp
D = os.environ.get("DEPMAP_DIR", "."); O = os.environ.get("RESULTS_DIR", "results") + "/"; R = {}
mw = lambda x, y: float(mannwhitneyu(x, y, alternative="less", method="asymptotic")[1])
df = rp.load(D); df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
genes = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
F = pd.read_csv(O + "final_targets.csv"); T = pd.read_csv(O + "all_tests.csv"); tg = sorted(F.Gene.unique()); sym = lambda s: s.split(" (")[0]
col = lambda s: [c for c in genes if c.startswith(s + " (")][0]

# ---- 1. co-dependency with lineage removed
X = df[tg].apply(pd.to_numeric, errors="coerce"); Xr = X - X.groupby(df.lineage).transform("mean")
C0, C1 = X.corr(), Xr.corr(); C1.to_csv(O + "codependency_residual_corr.csv"); iu = np.triu_indices(len(tg), 1)
TF = ["SOX10", "HNF1B", "PAX8", "IRF4", "MYB", "CBFB", "FOXA1"]; MET = ["NAMPT", "NMNAT1", "ADSL", "PAICS", "UMPS", "CTPS1", "DHFR", "TYMS"]
cols = lambda L: [col(g) for g in L if any(c.startswith(g + " (") for c in tg)]
def blk(C, A, B, same):
    v = C.loc[A, B].values; v = v[np.triu_indices(len(A), 1)] if same else v.ravel(); return float(np.median(np.abs(v)))
rng = np.random.default_rng(42); rg = list(rng.choice([g for g in genes if g not in tg], 400, replace=False))
Y = df[rg].apply(pd.to_numeric, errors="coerce"); Yr = Y - Y.groupby(df.lineage).transform("mean")
pick = lambda C, a, b: float(C.loc[col(a), col(b)])
R["codependency_residual"] = {"n_gene_pairs_r_gt_0.5_raw": int((C0.values[iu] > 0.5).sum()), "n_gene_pairs_r_gt_0.5_residual": int((C1.values[iu] > 0.5).sum()),
    "median_abs_r_raw": float(np.median(np.abs(C0.values[iu]))), "median_abs_r_residual": float(np.median(np.abs(C1.values[iu]))),
    "random_median_abs_r_raw": float(np.nanmedian(np.abs(Y.corr().values[np.triu_indices(400, 1)]))), "random_median_abs_r_residual": float(np.nanmedian(np.abs(Yr.corr().values[np.triu_indices(400, 1)]))),
    "metabolic_within_raw": blk(C0, cols(MET), cols(MET), True), "metabolic_within_residual": blk(C1, cols(MET), cols(MET), True),
    "TF_vs_metabolic_raw": blk(C0, cols(TF), cols(MET), False), "TF_vs_metabolic_residual": blk(C1, cols(TF), cols(MET), False),
    **{f"{a}_{b}_{k}": pick(C, a, b) for a, b in [("ADSL", "PAICS"), ("UMPS", "CTPS1"), ("NMNAT1", "NAMPT"), ("SOX10", "HNF1B")] for k, C in [("raw", C0), ("residual", C1)]}}

# ---- 2. enrichment against a background of genes that could meet the median criterion
elig = T[T.Chronos_median < rp.CHRONOS_MAX].Gene.map(sym); universe = set(elig); hits = set(F.Gene.map(sym)); assert hits <= universe
full = set(T.Gene.map(sym)); res = []
for lib in ["GO_Biological_Process_2023", "KEGG_2021_Human"]:
    txt = requests.get("https://maayanlab.cloud/Enrichr/geneSetLibrary", params={"mode": "text", "libraryName": lib}, timeout=120).text; rows = []
    for line in txt.strip().split("\n"):
        p = line.split("\t"); gs = {g.split(",")[0] for g in p[2:] if g}
        if not 5 <= len(gs & full) <= 500: continue          # same set-size rule as the primary analysis
        gs = gs & universe; a = len(gs & hits); b = len(hits) - a; c = len(gs) - a; d = len(universe) - len(gs) - b
        rows.append((lib, p[0], len(gs), a, fisher_exact([[a, b], [c, d]], alternative="greater")[1], ";".join(sorted(gs & hits))))
    E = pd.DataFrame(rows, columns=["library", "term", "set_size_in_background", "n_hits", "p", "genes"]); E["q"] = multipletests(E.p, method="fdr_bh")[1]; res.append(E)
E = pd.concat(res).sort_values("q"); E.to_csv(O + "enrichment_restricted_background.csv", index=False)
R["enrichment_restricted"] = {"background_genes": len(universe), "hit_genes": len(hits), "n_q05": E[E.q < 0.05].groupby("library").size().to_dict(),
    "top": {l: E[E.library == l].head(8)[["term", "n_hits", "set_size_in_background", "q"]].round(5).values.tolist() for l in E.library.unique()}}

# ---- 3. HNF1B / PAX8 in renal cell carcinoma lines only
md = pd.read_csv(D + "/Model.csv")[["ModelID", "OncotreePrimaryDisease", "OncotreeSubtype"]].rename(columns={"ModelID": "DepMap_ID"}); dk = df[["DepMap_ID", "lineage"]].merge(md, on="DepMap_ID")
kid = dk[dk.lineage == "kidney"]; R["kidney_composition"] = kid.OncotreePrimaryDisease.value_counts().to_dict()
rcc = set(kid[kid.OncotreePrimaryDisease == "Renal Cell Carcinoma"].DepMap_ID); oth = df.lineage != "kidney"
for g in ["HNF1B", "PAX8", "CFLAR", "MDM2"]:
    y = pd.to_numeric(df[col(g)], errors="coerce"); a = y[df.DepMap_ID.isin(rcc)].dropna(); nr = y[(df.lineage == "kidney") & ~df.DepMap_ID.isin(rcc)].dropna(); b = y[oth].dropna()
    R[f"{g}_kidney_RCC_only"] = {"n_rcc": len(a), "median_rcc": float(a.median()), "selectivity_rcc": float(a.median() - b.median()), "p_rcc": mw(a, b),
                                 "n_non_rcc_kidney": len(nr), "median_non_rcc_kidney": float(nr.median()) if len(nr) else None}

# ---- 4. subsampling with the BH-equivalent threshold of the primary family
pstar = float(T[T.q_value < rp.Q_MAX].p_value.max()); rows = []; rng = np.random.default_rng(42)
for r in F[F.n_target >= 30].itertuples():
    x = pd.to_numeric(df[r.Gene], errors="coerce"); tv = x[df.lineage == r.Lineage].dropna().values; ov = x[df.lineage != r.Lineage].dropna().values
    for k in [10, 20]:
        hit = 0
        for _ in range(300):
            sm = rng.choice(tv, k, replace=False); hit += (np.median(sm) < rp.CHRONOS_MAX and np.median(sm) - np.median(ov) < rp.SEL_MAX and mw(sm, ov) <= pstar)
        rows.append(dict(Gene=r.Gene, Lineage=r.Lineage, n_full=len(tv), k=k, frac_retained=hit / 300))
S = pd.DataFrame(rows); S.to_csv(O + "subsampling_power_bh.csv", index=False)
R["subsampling_bh"] = {"p_star": pstar, "n_pairs": int(len(S) / 2), "by_k": S.groupby("k").frac_retained.agg(["mean", "median", "min"]).round(3).to_dict("index")}
json.dump(R, open(O + "revision2_summary.json", "w"), indent=1, default=float); print(json.dumps(R, indent=1, default=float))
