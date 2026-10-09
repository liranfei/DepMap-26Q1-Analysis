"""Fig 13: validation of a single candidate of the SNAP23 dependency in the bowel lineage.
Uses the figure conventions of the manuscript (src/figlib.py): Arial 8-12 pt, 300 dpi,
LZW-compressed TIFF, white border, no title inside the file."""
import os, sys, json
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, spearmanr

DATA = os.environ.get("DEPMAP_DIR", "."); DS = os.environ.get("SANGER_DIR", "")
HERE = os.environ.get("RESULTS_DIR", "results")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figlib as fl
import matplotlib.pyplot as plt

GENE, LIN = "SNAP23", "bowel"

# ---------------- data ----------------
hdr = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), nrows=0).columns
need = [hdr[0]] + [c for c in hdr if c.split(" ")[0] in (GENE, "STX4", "STXBP3")]
ge = pd.read_csv(os.path.join(DATA, "gene_effect.csv"), usecols=need).rename(columns={hdr[0]: "DepMap_ID"})
M = pd.read_csv(os.path.join(DATA, "Model.csv"))
meta = M[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "GrowthPattern"]]
meta.columns = ["DepMap_ID", "lineage", "disease", "subtype", "growth"]
meta["lineage"] = meta.lineage.str.strip().str.lower()
d = ge.merge(meta, on="DepMap_ID", how="left", validate="one_to_one")
d = d[d.disease != "Non-Cancerous"].reset_index(drop=True)
gcol = {g: [c for c in need if c.split(" ")[0] == g][0] for g in (GENE, "STX4", "STXBP3")}
for g, c in gcol.items(): d[g] = pd.to_numeric(d[c], errors="coerce")
dd = d[d[GENE].notna()]
bowel = dd[dd.lineage == LIN][GENE].values
oth_all = dd[dd.lineage != LIN][GENE].values
oth_adh = dd[(dd.lineage != LIN) & (dd.growth == "Adherent")][GENE].values

# genotype model
mu = pd.read_csv(os.path.join(DATA, "OmicsSomaticMutations.csv"),
                 usecols=["ModelID", "HugoSymbol", "VepImpact", "Hotspot", "IsDefaultEntryForModel"],
                 low_memory=False)
mu = mu[(mu.IsDefaultEntryForModel == "Yes") & mu.HugoSymbol.isin(["KRAS", "APC", "TP53", "BRAF", "PIK3CA"])]
g = dd.copy(); g["bowel"] = (g.lineage == LIN).astype(int); g["y"] = g[GENE]
for sym, hs in [("KRAS", True), ("BRAF", True), ("PIK3CA", True), ("APC", False), ("TP53", False)]:
    s = mu[mu.HugoSymbol == sym]
    ids = set(s[s.Hotspot == True].ModelID) if hs else set(s[s.VepImpact.isin(["HIGH", "MODERATE"])].ModelID)
    g[sym] = g.DepMap_ID.isin(ids).astype(int)
m0 = smf.ols("y ~ bowel", g).fit(cov_type="HC3", use_t=True)
m1 = smf.ols("y ~ bowel + KRAS + APC + TP53 + BRAF + PIK3CA", g).fit(cov_type="HC3", use_t=True)
coefs = [(float(m.params.bowel), float(m.conf_int().loc["bowel", 0]), float(m.conf_int().loc["bowel", 1])) for m in (m0, m1)]

# Sanger
sg = pd.read_csv(os.path.join(DS, "gene_effect.csv"), low_memory=False)
sg = sg.rename(columns={sg.columns[0]: "DepMap_ID"})
scol = [c for c in sg.columns if c.split(" ")[0] == GENE][0]
s = sg[["DepMap_ID", scol]].copy(); s["y"] = pd.to_numeric(s[scol], errors="coerce")
s = s.merge(meta, on="DepMap_ID", how="left")
s = s[s.y.notna() & s.lineage.notna() & (s.disease != "Non-Cancerous")]
sb, so = s[s.lineage == LIN].y.values, s[s.lineage != LIN].y.values

# co-dependency (stored) and expression
r = pd.read_csv(HERE + "/snap23_codependency_all_genes.csv", index_col=0).iloc[:, 0]
top = r.sort_values(ascending=False).head(20)[::-1]
eh = pd.read_csv(os.path.join(DATA, "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv"), nrows=0).columns
ec = ["ModelID", "IsDefaultEntryForModel"] + [c for c in eh if c.split(" ")[0] == "SNAP25"]
ex = pd.read_csv(os.path.join(DATA, "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv"), usecols=ec)
ex = ex[ex.IsDefaultEntryForModel == "Yes"].rename(columns={"ModelID": "DepMap_ID"})
ex.columns = [c.split(" ")[0] for c in ex.columns]
e = dd[["DepMap_ID", GENE, "lineage"]].merge(ex[["DepMap_ID", "SNAP25"]], on="DepMap_ID")
rho, prho = spearmanr(e.SNAP25, e[GENE])

mw = lambda a, b: float(mannwhitneyu(a, b, alternative="less", method="asymptotic")[1])
def pstr(p):
    m, e = f"{p:.1e}".split("e"); return f"{m} \u00d7 10$^{{{int(e)}}}$"
BOX = dict(facecolor="white", edgecolor="none", alpha=1.0, pad=1.5)

# ---------------- figure ----------------
fig, ax = plt.subplots(2, 3, figsize=(7.2, 6.6), gridspec_kw=dict(height_ratios=[1, 1.45]))   # taller lower row: 20 legible labels in panel E
def box(a, groups, labels, colours, ylab):
    bp = a.boxplot(groups, widths=0.6, showfliers=False, patch_artist=True, medianprops=dict(color="black", lw=1.0))
    for p, c in zip(bp["boxes"], colours): p.set(facecolor=c, alpha=0.35, edgecolor=c, lw=0.8)
    for i, (v, c) in enumerate(zip(groups, colours), start=1):
        x = np.random.default_rng(0).normal(i, 0.07, len(v))
        a.plot(x, v, ".", ms=1.6, color=c, alpha=0.5, zorder=0)
    a.set_xticks(range(1, len(labels) + 1)); a.set_xticklabels(labels)
    a.axhline(-1, ls=":", lw=0.8, color=fl.GREY)
    if ylab: a.set_ylabel(ylab)

# A: bowel vs others
A = ax[0, 0]
box(A, [bowel, oth_adh, oth_all],
    [f"Bowel\n(n = {len(bowel)})", f"Other,\nadherent\n(n = {len(oth_adh)})", f"All other\n(n = {len(oth_all)})"],
    [fl.VERM, fl.BLUE, fl.GREY], "SNAP23 Chronos gene effect")
A.text(0.03, 0.02, f"vs all other: p = {pstr(mw(bowel, oth_all))}\nvs adherent only: p = {pstr(mw(bowel, oth_adh))}",
       transform=A.transAxes, fontsize=8, va="bottom", bbox=BOX, zorder=5)

# B: subtypes
B = ax[0, 1]
sub = dd[(dd.lineage == LIN)].groupby("subtype")[GENE].apply(lambda v: v.dropna().values)
sub = {k: v for k, v in sub.items() if len(v) >= 5}
names = sorted(sub, key=lambda k: -len(sub[k]))
box(B, [sub[n] for n in names] + [oth_all],
    [f"{n.replace(' Adenocarcinoma','')}\nadeno-\ncarcinoma\n(n = {len(sub[n])})" for n in names] + [f"All other\nlineages\n(n = {len(oth_all)})"],
    [fl.VERM, fl.ORANGE, fl.GREY], None)
B.set_ylabel("SNAP23 Chronos gene effect")

# C: genotype adjustment
C = ax[0, 2]
xs = [0, 1]
for i, (c, lo, hi) in enumerate(coefs):
    C.errorbar(i, c, yerr=[[c - lo], [hi - c]], fmt="o", ms=5, lw=1.2, capsize=3,
               color=fl.VERM if i == 0 else fl.BLUE)
C.set_xticks(xs); C.set_xticklabels(["Lineage\nonly", "Adjusted for\nKRAS, APC, TP53,\nBRAF, PIK3CA"])
C.set_xlim(-0.5, 1.5); C.axhline(0, lw=0.8, color="black")
C.set_ylabel("Bowel coefficient (Chronos)")
C.text(0.5, 0.04, f"{100*coefs[1][0]/coefs[0][0]:.0f}% of the\nunadjusted effect retained",
       transform=C.transAxes, fontsize=8, ha="center", bbox=BOX, zorder=5)

# D: Sanger
D = ax[1, 0]
box(D, [sb, so], [f"Bowel\n(n = {len(sb)})", f"All other\n(n = {len(so)})"], [fl.VERM, fl.GREY],
    "SNAP23 Chronos (Project Score)")
D.text(0.03, 0.97, f"selectivity = {np.median(sb)-np.median(so):.2f}\nq = 5.8 × 10$^{{-4}}$",
       transform=D.transAxes, fontsize=8, va="top", bbox=BOX, zorder=5)

# E: co-dependency
E = ax[1, 1]
syms = [t.split(" ")[0] for t in top.index]
cols = [fl.VERM if s_ in ("STX4", "STXBP3") else fl.GREY for s_ in syms]
E.barh(range(len(top)), top.values, color=cols, height=0.7)
E.set_yticks(range(len(top))); E.set_yticklabels(syms, fontsize=8)
E.set_xlabel("Pearson r with SNAP23\n(after removing lineage means)")
E.tick_params(axis="y", length=0)
E.text(0.97, 0.45, f"top 20 of\n{len(r):,} genes", transform=E.transAxes, fontsize=8, ha="right", va="center", bbox=BOX, zorder=5)

# F: paralogue buffering
F = ax[1, 2]
nb = e[e.lineage != LIN]; bw = e[e.lineage == LIN]
F.plot(nb.SNAP25, nb[GENE], ".", ms=2.0, color=fl.GREY, alpha=0.5, label="Other lineages")
F.plot(bw.SNAP25, bw[GENE], ".", ms=3.2, color=fl.VERM, alpha=0.9, label="Bowel")
F.set_xlabel("SNAP25 expression (log$_2$[TPM + 1])"); F.set_ylabel("SNAP23 Chronos gene effect")
F.axhline(-1, ls=":", lw=0.8, color=fl.GREY)
F.legend(loc="lower right", handletextpad=0.2, markerscale=2.5, framealpha=0.85,
         edgecolor="none", facecolor="white")
F.text(0.97, 0.97, f"Spearman $\\rho$ = {rho:.2f}\np = {pstr(prho)}", transform=F.transAxes,
       fontsize=8, va="top", ha="right", bbox=BOX, zorder=5)

for a, lab in zip(ax.ravel(), "ABCDEF"):
    a.text(-0.30, 1.04, lab, transform=a.transAxes, fontsize=10, fontweight="bold", va="bottom")
fig.tight_layout(w_pad=1.6, h_pad=1.4)
fl.save(fig, "Fig13")

json.dump({"panelA_p_vs_all": mw(bowel, oth_all), "panelA_p_vs_adherent": mw(bowel, oth_adh),
           "panelC_coef_unadj": coefs[0], "panelC_coef_adj": coefs[1],
           "panelD_sanger_selectivity": float(np.median(sb) - np.median(so)), "panelD_n": [len(sb), len(so)],
           "panelE_top2": {s_: round(float(v), 3) for s_, v in zip(syms[::-1][:2], top.values[::-1][:2])},
           "panelF_rho": float(rho), "panelF_p": float(prho), "panelF_n": int(len(e))},
          open(HERE + "/Fig13_values.json", "w"), indent=1)
print("figure values written")
