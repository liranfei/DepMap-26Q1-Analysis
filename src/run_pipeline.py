"""Lineage-specific dependency pipeline (DepMap 26Q1), single reproducible script.

Inputs (same folder or --data-dir): gene_effect.csv (CRISPRGeneEffect), Model.csv
Outputs: out dir with lineage_counts.csv, variance_rank.csv, all_tests.csv, final_targets.csv, summary.json
Design: variance top-N -> one-sided Mann-Whitney U (target lineage < all others) for EVERY
gene x lineage (n>=5) -> Benjamini-Hochberg over the FULL family -> then effect-size thresholds.
"""
import argparse, json, os
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

TOP_N, MIN_N = 500, 5
CHRONOS_MAX, SEL_MAX, Q_MAX = -1.0, -0.5, 0.05
EXPECTED_MODELS = 1208


def load(data_dir):
    ge = pd.read_csv(os.path.join(data_dir, "gene_effect.csv"))
    ge = ge.rename(columns={ge.columns[0]: "DepMap_ID"})
    ge["DepMap_ID"] = ge["DepMap_ID"].astype(str).str.strip()
    meta = pd.read_csv(os.path.join(data_dir, "Model.csv"))[["ModelID", "OncotreeLineage", "OncotreePrimaryDisease"]]
    meta.columns = ["DepMap_ID", "lineage", "primary_disease"]
    meta["DepMap_ID"] = meta["DepMap_ID"].astype(str).str.strip()
    df = ge.merge(meta, on="DepMap_ID", how="left", validate="one_to_one")
    n_missing = int(df["lineage"].isna().sum())
    assert n_missing == 0, f"{n_missing} cell lines lack lineage annotation - merge failure"
    assert len(df) == len(ge), "rows lost in merge"
    df["lineage"] = df["lineage"].str.strip().str.lower()
    return df


def test_all(df, genes, lineage):
    rows = []
    for g in genes:
        v = pd.to_numeric(df[g], errors="coerce")
        ok = v.notna().values
        vv, ll = v.values[ok], lineage[ok]
        for c in np.unique(ll):
            m = ll == c
            a, b = vv[m], vv[~m]
            if len(a) < MIN_N or len(b) == 0:
                continue
            _, p = mannwhitneyu(a, b, alternative="less", method="asymptotic")
            rows.append((g, c, float(np.median(a)), float(np.median(a) - np.median(b)), len(a), p))
    t = pd.DataFrame(rows, columns=["Gene", "Lineage", "Chronos_median", "Selectivity", "n_target", "p_value"])
    t["q_value"] = multipletests(t["p_value"], method="fdr_bh")[1]
    return t


def select(t):
    return t[(t.q_value < Q_MAX) & (t.Chronos_median < CHRONOS_MAX) & (t.Selectivity < SEL_MAX)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=".")
    ap.add_argument("--out", default="out")
    ap.add_argument("--top-n", type=int, default=500, help="0 = all genes")
    ap.add_argument("--exclude-noncancerous", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    df = load(a.data_dir)
    assert len(df) == EXPECTED_MODELS, f"expected {EXPECTED_MODELS} cell lines, got {len(df)}"
    n_nc = int((df.primary_disease == "Non-Cancerous").sum())
    if a.exclude_noncancerous:
        df = df[df.primary_disease != "Non-Cancerous"].reset_index(drop=True)
    genes_all = [c for c in df.columns if c not in ("DepMap_ID", "lineage", "primary_disease")]
    var = df[genes_all].apply(pd.to_numeric, errors="coerce").var()
    vr = var.rename("Variance").reset_index().rename(columns={"index": "Gene"})
    vr = vr.sort_values(["Variance", "Gene"], ascending=[False, True], kind="mergesort")
    top_n = len(vr) if a.top_n == 0 else a.top_n
    vr.head(top_n).to_csv(os.path.join(a.out, "variance_rank.csv"), index=False)
    lc = df["lineage"].value_counts()
    lc.to_csv(os.path.join(a.out, "lineage_counts.csv"))
    t = test_all(df, vr.head(top_n)["Gene"].tolist(), df["lineage"].values)
    t.to_csv(os.path.join(a.out, "all_tests.csv"), index=False)
    f = select(t).sort_values("Selectivity")
    f["confidence"] = np.where(f.n_target < 10, "Low", "High")
    f.to_csv(os.path.join(a.out, "final_targets.csv"), index=False)
    s = dict(n_cell_lines_analysed=len(df), n_noncancerous_in_1208=n_nc, excluded_noncancerous=a.exclude_noncancerous, top_n=top_n, n_genes_total=len(genes_all), n_lineages_n_ge_5=int((lc >= MIN_N).sum()),
             n_tests=len(t), n_q_lt_0_05=int((t.q_value < Q_MAX).sum()), n_final_pairs=len(f),
             n_final_genes=int(f.Gene.nunique()), n_final_lineages=int(f.Lineage.nunique()),
             n_low_confidence_pairs=int((f.n_target < 10).sum()))
    json.dump(s, open(os.path.join(a.out, "summary.json"), "w"), indent=2)
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
