import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import os  # 新增：用于自动获取路径

# === 路径修改（自动识别桌面路径） ===
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
# 如果文件在桌面上，下面的代码会自动找到它
sig_path = os.path.join(desktop, "significant_pairs.csv")
final_path = os.path.join(desktop, "final_targets.csv")

# 读取文件
sig = pd.read_csv(sig_path)
final = pd.read_csv(final_path)

# 计算 -log10(q)
sig["-log10q"] = -np.log10(sig["q_value"].clip(lower=1e-50))
sig["final"] = sig["Gene"].isin(final["Gene"])

fig, ax = plt.subplots(figsize=(10, 7))

# 背景点（未入选）
bg = sig[~sig["final"]]
ax.scatter(bg["Selectivity"], bg["-log10q"],
           color="#cccccc", s=15, alpha=0.5, linewidths=0, zorder=1)

# 最终93对
hi = sig[sig["final"]]
ax.scatter(hi["Selectivity"], hi["-log10q"],
           color="#C00000", s=25, alpha=0.8, linewidths=0, zorder=2)

# 标注Top10
top10 = hi.nlargest(10, "-log10q")
for _, row in top10.iterrows():
    label = row["Gene"].split(" ")[0]
    ax.annotate(label,
                xy=(row["Selectivity"], row["-log10q"]),
                xytext=(row["Selectivity"] + 0.03, row["-log10q"] + 0.5),
                fontsize=7.5, color="#333333",
                path_effects=[pe.withStroke(linewidth=2, foreground="white")],
                arrowprops=dict(arrowstyle="-", color="#aaaaaa", lw=0.8))

# 阈值线
ax.axvline(-0.5, color="#2171b5", linestyle="--", linewidth=1,
           label="Selectivity = -0.5")
ax.axhline(-np.log10(0.05), color="orange", linestyle="--", linewidth=1,
           label="FDR = 0.05")

ax.set_xlabel("Selectivity Score (Target vs Others, median difference)",
              fontsize=12)
ax.set_ylabel("-log₁₀(q-value)", fontsize=12)
ax.set_title("Volcano Plot: Lineage-Specific Dependency Screen\n(DepMap, n=1,208 cell lines)",
             fontsize=13, fontweight="bold")

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#C00000',
           markersize=8, label=f'Final targets (n={len(final)})'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#cccccc',
           markersize=8, label='Below threshold'),
    Line2D([0], [0], color='#2171b5', linestyle='--', label='Selectivity < -0.5'),
    Line2D([0], [0], color='orange', linestyle='--', label='FDR = 0.05'),
]
ax.legend(handles=legend_elements, loc="upper left", fontsize=9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
# 自动保存到桌面
save_path = os.path.join(desktop, "Fig07_Volcano.png")
plt.savefig(save_path, dpi=300, bbox_inches="tight")
plt.show()

print(f"07 完成，文件保存在: {save_path}")
