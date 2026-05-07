import pandas as pd
import numpy as np
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import os

# === 路径修正：自动获取桌面路径 ===
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
input_file = os.path.join(desktop, "standardized_matrix.csv")

if not os.path.exists(input_file):
    print(f"❌ 错误：在桌面找不到 {input_file}，请先运行 01 脚本生成它。")
else:
    df = pd.read_csv(input_file)
    df = df.dropna(subset=["primary_disease"])
    df["primary_disease"] = df["primary_disease"].str.strip().str.lower()

    # 计算所有基因方差
    all_genes = [c for c in df.columns if c not in ["DepMap_ID", "primary_disease"]]
    print("正在计算全局方差...")
    var_series = df[all_genes].apply(pd.to_numeric, errors="coerce").var().sort_values(ascending=False)

    unique_cancers = df["primary_disease"].unique()
    disease_labels = df["primary_disease"].values

    def run_pipeline(top_n):
        genes = var_series.head(top_n).index.tolist()
        data  = df[genes].values.astype(float)
        idx_dict = {c: np.where(disease_labels == c)[0] for c in unique_cancers}
        p_vals = []
        
        for g_idx in range(len(genes)):
            col = data[:, g_idx]
            for cancer, t_idx in idx_dict.items():
                if len(t_idx) < 5: continue
                target = col[t_idx]
                target = target[~np.isnan(target)]
                if len(target) < 5: continue
                
                o_idx = np.concatenate([v for k, v in idx_dict.items() if k != cancer])
                other = col[o_idx]
                other = other[~np.isnan(other)]
                if len(other) < 5: continue
                
                med_t = np.median(target)
                med_o = np.median(other)
                
                # 使用最终版严格阈值 -0.5
                if med_t < -1.0 and (med_t - med_o) < -0.5:
                    _, p = ranksums(target, other)
                    p_vals.append(p)
        
        if not p_vals: return 0
        _, q_vals, _, _ = multipletests(p_vals, method="fdr_bh")
        return int(np.sum(q_vals < 0.05))

    thresholds = [300, 500, 1000, 2000, 3000]
    results = []
    print("运行敏感性分析（约需5-10分钟，请耐心等待）...")
    for n in thresholds:
        count = run_pipeline(n)
        results.append({"Top_N": n, "Significant_Pairs": count})
        print(f"  ✅ Top {n}: 找到 {count} 对显著关联")

    res_df = pd.DataFrame(results)
    res_df.to_csv(os.path.join(desktop, "sensitivity_results.csv"), index=False)

    # 绘图
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(res_df["Top_N"], res_df["Significant_Pairs"], "o-", color="#2171b5", linewidth=2, markersize=8)
    ax.axvline(500, color="red", linestyle="--", linewidth=1.5, label="Current Selection (Top 500)")

    for _, row in res_df.iterrows():
        ax.annotate(f"{int(row['Significant_Pairs'])}",
                    xy=(row["Top_N"], row["Significant_Pairs"]),
                    xytext=(0, 10), textcoords="offset points", ha="center")

    ax.set_xlabel("Number of Top Variance Genes", fontsize=12)
    ax.set_ylabel("Significant Pairs (q < 0.05, Selectivity < -0.5)", fontsize=11)
    ax.set_title("Sensitivity Analysis: Threshold Robustness Check", fontsize=12, fontweight="bold")
    ax.legend()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    save_path = os.path.join(desktop, "Fig14_Sensitivity.png")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"14 完成，图片保存至: {save_path}")
