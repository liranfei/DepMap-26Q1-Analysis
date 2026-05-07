import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

# --- 自动处理路径 ---
# 这样无论你的用户 ID 是 3 还是 4，代码都能自动找到桌面
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
save_path = os.path.join(desktop_path, "Fig_Pipeline_v2.png")

# --- 绘图逻辑 ---
fig, ax = plt.subplots(figsize=(10, 6))

steps = ["Genome-wide\n(Genes)", "High Variance\n(Top 500 Genes)", 
         "Significant Pairs\n(Gene-Lineage, q<0.05)", "Specific Targets\n(Unique Genes)"]
counts = [18532, 500, 93, 57]
colors = ["#f0f0f0", "#c6dbef", "#6baed6", "#2171b5"]

bars = ax.bar(steps, counts, color=colors, edgecolor="black", linewidth=0.8, width=0.55)

# 数值标注
for bar, count in zip(bars, counts):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.15,
            f"{count:,}", ha="center", va="bottom", fontsize=13, fontweight="bold")

# 箭头标注
for i in range(len(steps)-1):
    x1 = bars[i].get_x() + bars[i].get_width()
    x2 = bars[i+1].get_x()
    ax.annotate("", xy=(x2 - 0.02, 0.5), xytext=(x1 + 0.02, 0.5),
                xycoords=("data", "axes fraction"),
                textcoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="->", color="gray", lw=1.5))

ax.set_yscale("log")
ax.set_ylabel("Count (Log10 Scale)", fontsize=12)
ax.set_title("Systematic Identification of Lineage-Specific Dependencies", 
             fontsize=14, fontweight="bold", pad=15)
ax.set_ylim(1, 10**5.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis='x', labelsize=11)

plt.tight_layout()

# --- 核心修复点：使用自动生成的路径保存 ---
try:
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    print(f"✅ 图2已成功保存至桌面: {save_path}")
except Exception as e:
    print(f"❌ 保存失败，错误原因: {e}")

plt.show()
