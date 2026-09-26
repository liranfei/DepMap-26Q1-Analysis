"""Fig 12: comparison of the 26Q1 candidate pairs with the independent Sanger Project Score screens (run after independent_sanger.py)."""
import os, json, numpy as np, pandas as pd
from figlib import *
TS = pd.read_csv(O + "all_tests_sanger.csv.gz"); T26 = pd.read_csv(O + ("all_tests.csv" if os.path.exists(O + "all_tests.csv") else "all_tests.csv.gz")); C = pd.read_csv(O + "independent_sanger_candidates.csv")
M = T26.merge(TS, on=["Gene", "Lineage"], suffixes=("_26", "_S"))
cat = np.where(C.q_value_S.isna(), "not testable", np.where((C.q_value_S < 0.05) & (C.Selectivity_S < 0), np.where((C.Chronos_median_S < -1) & (C.Selectivity_S < -0.5), "meets rule", "significant, below rule"), "not significant"))
fig, axs = plt.subplots(1, 2, figsize=(7.3, 3.9), gridspec_kw={"width_ratios": [1.25, 1]}); a, b = axs
a.hexbin(M.Selectivity_26, M.Selectivity_S, gridsize=90, bins="log", cmap="Greys", mincnt=1, extent=(-1.8, 0.8, -1.8, 0.8), linewidths=0)
cc = C[C.q_value_S.notna()]; a.scatter(cc.Selectivity_26, cc.Selectivity_S, s=14, color=VERM, edgecolor="white", linewidth=0.3, zorder=3, label="26Q1 candidates")
a.plot([-1.8, 0.8], [-1.8, 0.8], color="#555555", ls="--", lw=0.8); a.axhline(-0.5, color=BLUE, ls=":", lw=1); a.axvline(-0.5, color=BLUE, ls=":", lw=1)
a.set_xlim(-1.8, 0.8); a.set_ylim(-1.8, 0.8); a.set_xlabel("Selectivity, DepMap 26Q1 (Broad)"); a.set_ylabel("Selectivity, Project Score (Sanger)"); a.legend(frameon=False, loc="lower right")
order = ["meets rule", "significant, below rule", "not significant", "not testable"]; lab = ["Meets the\nsame rule", "Significant,\nbelow rule", "Not\nsignificant", "Not testable\nin Sanger"]; vals = [int((cat == k).sum()) for k in order]
b.bar(range(4), vals, color=[VERM, ORANGE, GREY, "#cccccc"], width=0.7); b.set_xticks(range(4)); b.set_xticklabels(lab); b.set_ylabel("26Q1 candidate pairs")
[b.text(i, v + 1, str(v), ha="center", fontsize=8) for i, v in enumerate(vals)]; b.set_ylim(0, max(vals) * 1.15)
for ax, l, x in ((a, "A", -0.16), (b, "B", -0.2)): ax.text(x, 1.02, l, transform=ax.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout(pad=0.6, w_pad=1.5); save(fig, "Fig12"); print(dict(zip(order, vals)))
