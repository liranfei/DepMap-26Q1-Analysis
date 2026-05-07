import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os

# --- 核心修复：自动获取桌面路径 ---
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
input_file = os.path.join(desktop_path, "final_targets.csv")

# 检查文件是否存在
if not os.path.exists(input_file):
    print(f"❌ 错误：在桌面上找不到文件: {input_file}")
    print("提示：请确认 final_targets.csv 是否在桌面上，且没有小云朵图标（已下载）。")
else:
    final = pd.read_csv(input_file)

    # ============ 图3：Top 20 靶点 ============
    top20 = final.nsmallest(20, "Selectivity")

    cancer_list = top20["Cancer"].unique().tolist()
    palette = ["#4472C4","#ED7D31","#70AD47","#C00000","#7030A0",
               "#00B0F0","#FF0000","#FFC000","#00B050","#7F7F7F"]
    color_map = {c: palette[i % len(palette)] for i, c in enumerate(cancer_list)}

    fig, ax = plt.subplots(figsize=(10, 8))
    for _, row in top20.iterrows():
        ax.barh(f"{row['Gene']}", row["Selectivity"],
                color=color_map[row["Cancer"]], edgecolor="white", height=0.7)

    ax.set_xlabel("Selectivity Score (Lower is better)", fontsize=12)
    ax.set_ylabel("Gene", fontsize=12)
    ax.set_title("Top 20 Most Significant Lineage-Specific Dependencies\n(Selectivity < -0.5, q < 0.05)", fontsize=13)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.invert_yaxis()

    patches = [matplotlib.patches.Patch(color=color_map[c], label=c) for c in cancer_list]
    ax.legend(handles=patches, title="Lineage", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.tight_layout()
    
    # 自动保存到桌面
    fig3_path = os.path.join(desktop_path, "Fig3_Top20_v2.png")
    plt.savefig(fig3_path, dpi=300, bbox_inches="tight")
    plt.show()

    # ============ 图4：癌种分布 ============
    cancer_counts = final["Cancer"].value_counts().sort_values()
    n_cancers = len(cancer_counts)
    cmap = matplotlib.colormaps["viridis"]
    colors = [cmap(i / n_cancers) for i in range(n_cancers)]

    fig, ax = plt.subplots(figsize=(10, 8))
    bars = ax.barh(cancer_counts.index, cancer_counts.values, color=colors)
    for bar, val in zip(bars, cancer_counts.values):
        ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                str(val), va="center", fontsize=10)

    ax.set_xlabel("Number of Significant Pairs (q < 0.05, Selectivity < -0.5)", fontsize=11)
    ax.set_ylabel("Cancer Lineage", fontsize=11)
    ax.set_title(f"Distribution of {len(final)} Lineage-Specific Dependencies\nacross {n_cancers} Cancer Lineages", fontsize=13)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    
    # 自动保存到桌面
    fig4_path = os.path.join(desktop_path, "Fig4_CancerDist_v2.png")
    plt.savefig(fig4_path, dpi=300, bbox_inches="tight")
    plt.show()
    
    print("-" * 30)
    print(f"✅ 成功！图3保存至: {fig3_path}")
    print(f"✅ 成功！图4保存至: {fig4_path}")
    print("-" * 30)
