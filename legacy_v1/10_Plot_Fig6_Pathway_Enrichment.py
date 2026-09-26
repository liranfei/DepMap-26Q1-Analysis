import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# === 自动获取桌面路径，避免数字输错 ===
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
final_path = os.path.join(desktop, "final_targets.csv")

# 检查文件是否存在
if not os.path.exists(final_path):
    print(f"❌ 错误：在桌面找不到 final_targets.csv，请检查文件名。")
else:
    final = pd.read_csv(final_path)

    # 手动定义基于已知生物学的通路分组
    pathway_map = {
        "NAD Metabolism":        ["NMNAT1 (64802)", "NAMPT (10135)"],
        "Cell Cycle":            ["CCND1 (595)", "CDK1 (983)", "SKP2 (6502)"],
        "Wnt/β-catenin":         ["CTNNB1 (1499)"],
        "Nucleotide Metabolism": ["ADSL (158)", "UMPS (7372)", "TYMS (7298)",
                                  "PAICS (10606)", "CTPS1 (1503)", "FPGS (2356)"],
        "Folate/1C Metabolism":  ["DHFR (1719)", "YRDC (79693)"],
        "TCA/Mitochondria":      ["SDHB (6390)", "ACO2 (50)"],
        "Transcription Factor":  ["IRF4 (3662)", "SOX10 (6663)", "CBFB (865)",
                                  "HNF1B (6928)"],
        "RAS Signaling":         ["KRAS (3845)", "CRKL (1399)"],
        "p53/MDM2":              ["MDM2 (4193)"],
        "Selenium Metabolism":   ["SEPHS2 (22928)", "TXNRD1 (7296)"],
        "Cytoskeleton/Adhesion": ["FERMT2 (10979)", "TRPM7 (54822)"],
        "Ubiquitin/Proteasome":  ["SNAP23 (8773)", "CHMP4B (128866)"],
    }

    rows = []
    for pathway, genes in pathway_map.items():
        matched = final[final["Gene"].isin(genes)]
        if len(matched) == 0:
            continue
        rows.append({
            "Pathway": pathway,
            "Gene_count": len(matched),
            "Mean_Selectivity": matched["Selectivity"].mean(),
            "Mean_Chronos": matched["Chronos_median"].mean(),
            "Min_qval": matched["q_value"].min(),
        })

    pw_df = pd.DataFrame(rows).sort_values("Mean_Selectivity")
    pw_df["-log10q"] = -np.log10(pw_df["Min_qval"].clip(1e-50))

    fig, ax = plt.subplots(figsize=(10, 7))

    # 绘制气泡
    scatter = ax.scatter(
        pw_df["Mean_Selectivity"],
        pw_df["-log10q"],
        s=pw_df["Gene_count"] * 120,
        c=pw_df["Mean_Chronos"],
        cmap="RdBu",
        vmin=-2.5, vmax=0,
        alpha=0.85,
        edgecolors="black",
        linewidths=0.6,
        zorder=3
    )

    # 添加文字标注
    for _, row in pw_df.iterrows():
        ax.annotate(row["Pathway"],
                    xy=(row["Mean_Selectivity"], row["-log10q"]),
                    xytext=(row["Mean_Selectivity"] + 0.02,
                            row["-log10q"] + 0.3),
                    fontsize=8.5,
                    ha="left")

    # 颜色条
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6)
    cbar.set_label("Mean Chronos Score", fontsize=10)

    # 气泡大小图例
    for size in [1, 3, 5]:
        ax.scatter([], [], s=size*120, color="gray", alpha=0.6,
                   label=f"{size} gene(s)")
    ax.legend(title="Gene count", loc="lower right", fontsize=8)

    ax.axvline(-0.5, color="gray", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xlabel("Mean Selectivity Score", fontsize=12)
    ax.set_ylabel("-log₁₀(minimum q-value)", fontsize=12)
    ax.set_title("Pathway-Level Enrichment of Lineage-Specific Dependencies\n(Bubble size = number of genes per pathway)",
                 fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    
    # 保存
    save_path = os.path.join(desktop, "Fig10_Pathway_Bubble.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"✅ 10 完成，文件已保存至桌面: {save_path}")
