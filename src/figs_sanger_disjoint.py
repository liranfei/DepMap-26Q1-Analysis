"""Fig 12 and S7 Table: candidates from the 26Q1 lines without a Sanger (KY) screen, tested in the Sanger Project Score screens (run after sanger_disjoint.py).
Replaces the earlier comparison (figs_sanger.py), in which the 26Q1 gene effects of the shared cell lines were partly derived from the same Sanger screens."""
import os, json, numpy as np, pandas as pd
from figlib import *
TS = pd.read_csv(O + "all_tests_sanger.csv.gz"); TD = pd.read_csv(O + "all_tests_noKY.csv.gz"); C = pd.read_csv(O + "sanger_disjoint_candidates.csv"); B = pd.read_csv(O + "sanger_disjoint_bands.csv")
M = TD.merge(TS, on=["Gene", "Lineage"], suffixes=("_disc", "_S"))
cat = np.where(C.q_value_S.isna(), "not testable", np.where((C.q_value_S < 0.05) & (C.Selectivity_S < 0), np.where((C.Chronos_median_S < -1) & (C.Selectivity_S < -0.5), "meets rule", "significant, below rule"), "not significant"))
fig = plt.figure(figsize=(7.3, 6.6)); gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1, 0.8]); a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1]); c = fig.add_subplot(gs[1, :])
a.hexbin(M.Selectivity_disc, M.Selectivity_S, gridsize=90, bins="log", cmap="Greys", mincnt=1, extent=(-1.8, 0.8, -1.8, 0.8), linewidths=0)
cc = C[C.q_value_S.notna()]; a.scatter(cc.Selectivity_disc, cc.Selectivity_S, s=14, color=VERM, edgecolor="white", linewidth=0.3, zorder=3, label="Discovery candidates")
a.plot([-1.8, 0.8], [-1.8, 0.8], color="#555555", ls="--", lw=0.8); a.axhline(-0.5, color=BLUE, ls=":", lw=1); a.axvline(-0.5, color=BLUE, ls=":", lw=1)
a.set_xlim(-1.8, 0.8); a.set_ylim(-1.8, 0.8); a.set_xlabel("Selectivity, 26Q1 lines without Sanger screen"); a.set_ylabel("Selectivity, Project Score (Sanger)"); a.legend(frameon=False, loc="lower right")
order = ["meets rule", "significant, below rule", "not significant", "not testable"]; lab = ["Meets the\nsame rule", "Significant,\nbelow rule", "Not\nsignificant", "Not testable\nin Sanger"]; vals = [int((cat == k).sum()) for k in order]
b.bar(range(4), vals, color=[VERM, ORANGE, GREY, "#cccccc"], width=0.7); b.set_xticks(range(4)); b.set_xticklabels(lab); b.set_ylabel("Discovery candidate pairs")
[b.text(i, v + 1, str(v), ha="center", fontsize=8) for i, v in enumerate(vals)]; b.set_ylim(0, max(vals) * 1.15)
labs = ["< −0.7", "−0.7 to −0.5", "−0.5 to −0.4", "−0.4 to −0.3", "−0.3 to −0.2", "−0.2 to 0"]
c.bar(range(len(B)), B.frac * 100, color=[VERM if n > 0 else GREY for n in B.n_candidates], width=0.7); c.set_xticks(range(len(B))); c.set_xticklabels(labs[:len(B)])
c.set_xlabel("Selectivity in 26Q1 lines without Sanger screen"); c.set_ylabel("Replicated in Sanger (%)"); c.set_ylim(0, 118)
[c.text(i, f * 100 + 2, f"n = {int(n):,} ({int(k)} prioritised)", ha="center", fontsize=8) for i, (f, k, n) in enumerate(zip(B.frac, B.n_candidates, B.n))]
for ax, l, x in ((a, "A", -0.16), (b, "B", -0.2), (c, "C", -0.07)): ax.text(x, 1.02, l, transform=ax.transAxes, fontsize=12, fontweight="bold")
fig.tight_layout(pad=0.6, w_pad=1.5); save(fig, "Fig12"); print(dict(zip(order, vals)))

# ---- S7 Table
SI = OUT + "supporting_information/"; os.makedirs(SI, exist_ok=True); sym = lambda s: s.str.replace(r" \(.*", "", regex=True)
C.insert(0, "Gene_symbol", sym(C.Gene)); C["status_in_Sanger"] = cat
ren = {"n_target_disc": "n target (discovery)", "Chronos_median_disc": "Median Chronos (discovery)", "Selectivity_disc": "Selectivity (discovery)", "p_value_disc": "p (discovery)", "q_value_disc": "q (discovery, BH over discovery family)",
       "n_target_S": "n target Sanger", "Chronos_median_S": "Median Chronos Sanger", "Selectivity_S": "Selectivity Sanger", "p_value_S": "p Sanger (one-sided MWU)", "q_value_S": "q Sanger (BH over Sanger family)"}
P = pd.read_csv(O + "sanger_disjoint_primary94.csv"); P.insert(0, "Gene_symbol", sym(P.Gene))
P = P.rename(columns={"Selectivity": "Selectivity (primary)", "Chronos_median": "Median Chronos (primary)", "n_target": "n target (primary)", "q_value": "q (primary)", "candidate_disc": "Candidate in discovery cohort",
                      **{k.replace("_disc", "_disc"): v for k, v in ren.items()}})
O_old = pd.read_csv(O + "independent_sanger_candidates.csv").rename(columns={"confidence": "Low-n flag (fewer than 10 target lines)"}); O_old.insert(0, "Gene_symbol", sym(O_old.Gene))
readme = pd.DataFrame({"Note": [
    "Discovery cohort: the 856 DepMap 26Q1 cancer cell lines without any Sanger (KY-library) screen in the combined 26Q1 data set, after excluding 22 models derived from the same patient as a Project Score model; the complete procedure (one-sided Mann-Whitney, BH over all tests, q < 0.05, median Chronos < -1, selectivity < -0.5) was applied to it.",
    "Replication: the same procedure applied to the Sanger Project Score Chronos data (316 cancer models); discovery and replication share no screens, models or patients.",
    "Sheet 'Primary 94 pairs': the 94 candidate pairs of the primary analysis (all 1,192 lines) with their results in the discovery cohort and in the Sanger data.",
    "Sheet 'Overlapping (not independent)': earlier comparison of the primary candidates with the Sanger data. For the shared cell lines the 26Q1 gene effects include the same Sanger screens, so this comparison is not independent; shown for transparency only."]})
with pd.ExcelWriter(SI + "S7_Table.xlsx") as w:
    readme.to_excel(w, index=False, sheet_name="Read me"); C.rename(columns=ren).to_excel(w, index=False, sheet_name="Discovery candidates vs Sanger")
    P.to_excel(w, index=False, sheet_name="Primary 94 pairs"); O_old.to_excel(w, index=False, sheet_name="Overlapping (not independent)")
print("S7 written")
