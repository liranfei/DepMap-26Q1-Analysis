import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import os

# === 自动获取桌面路径，避免数字输错 ===
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
final_path = os.path.join(desktop, "final_targets.csv")

# 读取文件
if not os.path.exists(final_path):
    print(f"❌ 错误：在桌面找不到 final_targets.csv，请检查文件名。")
else:
    final = pd.read_csv(final_path)

    # 颜色映射（按癌种）
    cancers = final["Cancer"].unique()
    # 使用正确的 colormap 获取方式 (兼容新旧版本 matplotlib)
    try:
        cmap = plt.get_cmap("tab20")
    except:
        cmap = plt.cm.get_cmap("tab20", len(cancers))
    
    color_map = {c: cmap(i % 20) for i, c in enumerate(cancers)}

    fig, ax = plt.subplots(figsize=(10, 7))

    for cancer in cancers:
        sub = final[final["Cancer"] == cancer]
        ax.scatter(sub["Selectivity"], sub["Chronos_median"],
                   color=color_map[cancer], s=50, alpha=0.8,
                   label=cancer, edgecolors="white", linewidths=0.4, zorder=3)

    # 阈值线
    ax.axvline(-0.5, color="gray", linestyle="--", linewidth=1, alpha=0.7)
    ax.axhline(-1.0, color="gray", linestyle="--", linewidth=1, alpha=0.7)

    # 标注关键靶点 (标注 93 对中的核心基因)
    highlights = ["KRAS (3845)", "SOX10 (6663)", "CTNNB1 (1499)",
                  "MDM2 (4193)", "CBFB (865)", "HNF1B (6928)"]
    for _, row in final[final["Gene"].isin(highlights)].iterrows():
        label = row["Gene"].split(" ")[0]
        ax.annotate(label,
                    xy=(row["Selectivity"], row["Chronos_median"]),
                    xytext=(row["Selectivity"] - 0.05, row["Chronos_median"] - 0.12),
                    fontsize=8, fontweight="bold",
                    path_effects=[pe.withStroke(linewidth=2, foreground="white")],
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.8))

    ax.set_xlabel("Selectivity Score (Lower = More Lineage-Specific)", fontsize=12)
    ax.set_ylabel("Chronos Dependency Score (Lower = More Essential)", fontsize=12)
    ax.set_title("Lineage-Specific Dependencies: Efficacy vs Selectivity\n(93 high-confidence pairs, 20 cancer lineages)",
                 fontsize=13, fontweight="bold")

    # 象限注释
    ax.text(-1.55, -0.2, "High Selectivity\nLow Efficacy",
            fontsize=8, color="gray", alpha=0.7)
    ax.text(-1.55, -2.7, "High Selectivity\nHigh Efficacy ★",
            fontsize=8, color="#C00000", alpha=0.9)

    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left",
              fontsize=7.5, ncol=1, frameon=True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    
    # 保存
    save_path = os.path.join(desktop, "Fig08_Scatter.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"✅ 08 完成，文件已保存至桌面: {save_path}")
