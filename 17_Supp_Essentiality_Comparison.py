import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# === 路径修正：自动获取桌面路径 ===
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
final_path = os.path.join(desktop, "final_targets.csv")

if not os.path.exists(final_path):
    print(f"❌ 错误：在桌面找不到 {final_path}，请确保之前的脚本已成功运行。")
else:
    final = pd.read_csv(final_path)

    # DepMap 公认的 Common Essentials 核心库
    common_essentials = {
        "RPL5","RPL11","RPS14","RPS19","MDM2","CDK4","CDK6",
        "MYC","PCNA","POLA1","RFC1","RRM1","RRM2","TYMS",
        "DHFR","TOP2A","CCNA2","CCNB1","CDC20","BUB1B",
        "MAD2L1","AURKB","PLK1","CENPE","KIF11","KIF2C",
        "NUP98","XPO1","KPNB1","EIF4E","EIF4A1","EIF2S1",
        "ATP1A1","ATP5F1A","VDAC1","TOMM40","HSPA5","HSP90AA1",
        "ACTB","ACTIN","TUBB","VCP","PSMD1","PSMD2","PSMC1"
    }

    # 提取靶点名称
    final["Gene_Symbol"] = final["Gene"].str.extract(r'^(.+?)\s*\(')
    our_genes = set(final["Gene_Symbol"].unique())

    overlap    = our_genes & common_essentials
    unique_our = our_genes - common_essentials

    print(f"--- 必需性校验报告 ---")
    print(f"✅ 总靶点基因: {len(our_genes)}")
    print(f"⚠️ 与核心必需基因重叠: {len(overlap)}")
    if overlap: print(f"   重叠基因包括: {list(overlap)}")
    print(f"🎯 特异性靶点数: {len(unique_our)}")

    # 准备绘图
    final["is_essential"] = final["Gene_Symbol"].isin(common_essentials)
    essential_scores    = final[final["is_essential"]]["Chronos_median"]
    nonessential_scores = final[~final["is_essential"]]["Chronos_median"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 左图：重叠数量 (这里代码做了兼容处理)
    ax1 = axes[0]
    try:
        from matplotlib_venn import venn2
        venn2(subsets=(len(unique_our), len(common_essentials), len(overlap)),
              set_labels=("Our Targets\n(57 genes)", "Common Essentials\n(DepMap)"),
              set_colors=("#C00000", "#2171b5"),
              ax=ax1)
        ax1.set_title("Overlap Analysis", fontsize=12)
    except ImportError:
        # 未安装 venn 库时的条形图替代方案
        ax1.bar(["Unique Targets", "Overlap w/ Essentials"],
                [len(unique_our), len(overlap)],
                color=["#C00000", "#888888"], edgecolor="black")
        ax1.set_ylabel("Number of Genes")
        ax1.set_title("Gene Distribution", fontsize=12)
        for i, v in enumerate([len(unique_our), len(overlap)]):
            ax1.text(i, v + 0.3, str(v), ha="center", fontweight="bold")

    # 右图：Chronos 分布对比
    ax2 = axes[1]
    ax2.hist(nonessential_scores, bins=15, alpha=0.7, color="#C00000", 
             label=f"Specific Targets (n={len(nonessential_scores)})", edgecolor="white")
    
    if not essential_scores.empty:
        ax2.hist(essential_scores, bins=10, alpha=0.7, color="#888888",
                 label=f"Essential Overlap (n={len(essential_scores)})", edgecolor="white")

    ax2.axvline(np.median(nonessential_scores), color="#C00000", linestyle="--", 
                linewidth=1.5, label=f"Median = {np.median(nonessential_scores):.2f}")
    
    ax2.set_xlabel("Chronos Dependency Score", fontsize=11)
    ax2.set_ylabel("Frequency", fontsize=11)
    ax2.set_title("Target Potency Distribution", fontsize=12)
    ax2.legend(fontsize=9)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    save_path = os.path.join(desktop, "Fig17_Essentiality.png")
    plt.suptitle("Target Specificity Validation", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
    print(f"17 完成，图片已保存至桌面: {save_path}")
