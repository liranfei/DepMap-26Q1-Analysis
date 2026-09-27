"""Supporting table S8 (literature and drug-gene annotation) as xlsx.  S7 is written by figs_sanger_disjoint.py."""
import os, pandas as pd
O = os.environ.get("RESULTS_DIR", "results") + "/"; OUT = os.environ.get("FIG_DIR", "figures") + "/supporting_information/"; os.makedirs(OUT, exist_ok=True)
C = pd.read_csv(O + "independent_sanger_candidates.csv"); C["Gene_symbol"] = C.Gene.str.replace(r" \(.*", "", regex=True)
C["status_in_Sanger"] = ["not testable" if pd.isna(q) else ("meets rule" if (q < 0.05 and m < -1 and s < -0.5) else ("q<0.05, negative, below rule" if (q < 0.05 and s < 0) else "not significant")) for q, m, s in zip(C.q_value_S, C.Chronos_median_S, C.Selectivity_S)]
C = C.rename(columns={"Chronos_median_26": "Median Chronos 26Q1", "Selectivity_26": "Selectivity 26Q1", "n_target_26": "n target 26Q1", "q_value_26": "q 26Q1", "Chronos_median_S": "Median Chronos Sanger", "Selectivity_S": "Selectivity Sanger", "n_target_S": "n target Sanger", "p_value_S": "p Sanger (one-sided MWU)", "q_value_S": "q Sanger (BH, full Sanger family)"})
cols = [c for c in ["Gene_symbol", "Gene", "Lineage", "n target 26Q1", "Median Chronos 26Q1", "Selectivity 26Q1", "q 26Q1", "n target Sanger", "Median Chronos Sanger", "Selectivity Sanger", "p Sanger (one-sided MWU)", "q Sanger (BH, full Sanger family)", "status_in_Sanger"] if c in C.columns]
L = pd.read_csv(O + "annotation_pubmed.csv"); G = pd.read_csv(O + "annotation_dgidb.csv")
Z = pd.read_csv(O + "annotation_pubmed_zero_hit_recheck.csv") if os.path.exists(O + "annotation_pubmed_zero_hit_recheck.csv") else None
with pd.ExcelWriter(OUT + "S8_Table.xlsx") as w:
    L.to_excel(w, index=False, sheet_name="PubMed counts per pair"); G.to_excel(w, index=False, sheet_name="DGIdb per gene")
    if Z is not None: Z.to_excel(w, index=False, sheet_name="Recheck of zero-record pairs")
print("S8 written")
