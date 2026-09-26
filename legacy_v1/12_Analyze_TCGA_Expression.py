import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

# --- 核心修改：自动获取当前用户的桌面路径 ---
# 这样代码在任何人的 Mac 或 Windows 上都能直接运行，无需手动修改用户名数字
desktop = os.path.join(os.path.expanduser("~"), "Desktop")

print("正在读取数据，请稍候...")

# --- 1. 读取文件 (路径已设为自动识别) ---
path_tcga = os.path.join(desktop, "TCGA-PAAD.star_counts.tsv")
path_dep = os.path.join(desktop, "standardized_matrix.csv")

# 检查文件是否存在
if not os.path.exists(path_tcga) or not os.path.exists(path_dep):
    print("❌ 错误：桌面上缺少必要的文件！")
    print(f"请确保桌面上包含：\n1. TCGA-PAAD.star_counts.tsv\n2. standardized_matrix.csv")
    exit()

tcga = pd.read_csv(path_tcga, sep="\t", index_col=0)
df_dep = pd.read_csv(path_dep)

# --- 2. 提取 KRAS 数据 ---
# KRAS的Ensembl ID (模糊匹配)
kras_candidates = [i for i in tcga.index if "ENSG00000133703" in str(i)]
if not kras_candidates:
    print("❌ 未在 TCGA 数据中找到 KRAS ID")
    exit()

kras_row = kras_candidates[0]
kras_expr = tcga.loc[kras_row].apply(pd.to_numeric, errors="coerce")
kras_log = np.log2(kras_expr + 1)

# 分肿瘤(01)和正常(11)
tumor_cols = [c for c in kras_log.index if "-01" in c]
normal_cols = [c for c in kras_log.index if "-11" in c]
tumor_vals = kras_log[tumor_cols].dropna()
normal_vals = kras_log[normal_cols].dropna() if normal_cols else None

# --- 3. 准备 DepMap KRAS 依赖性数据 ---
df_dep["primary_disease"] = df_dep["primary_disease"].str.strip().str.lower()
kras_gene = "KRAS (3845)"

# 自动兼容列名
if kras_gene not in df_dep.columns:
    kras_cols = [c for c in df_dep.columns if "KRAS" in c]
    if kras_cols: kras_gene = kras_cols[0]

pancreas = df_dep[df_dep["primary_disease"] == "pancreas"][kras_gene].dropna()
others = df_dep[df_dep["primary_disease"] != "pancreas"][kras_gene].dropna()

# --- 4. 绘图 ---
fig, axes = plt.subplots(1, 2, figsize=(13, 6))

# === 左图：TCGA KRAS表达 ===
ax1 = axes[0]
if normal_vals is not None and len(normal_vals) > 0:
    stat, p = mannwhitneyu(tumor_vals, normal_vals, alternative="greater")
    p_label = f"p = {p:.2e}"
    data_list = [normal_vals.values, tumor_vals.values]
    labels_list = [f"Normal\n(n={len(normal_vals)})", f"Tumor\n(n={len(tumor_cols)})"]
    colors_list = ["#6baed6", "#C00000"]
else:
    p_label = ""
    data_list = [tumor_vals.values]
    labels_list = [f"Tumor\n(n={len(tumor_cols)})"]
    colors_list = ["#C00000"]

bp = ax1.boxplot(data_list, patch_artist=True, widths=0.5, medianprops=dict(color="black", linewidth=2))
for box, color in zip(bp["boxes"], colors_list):
    box.set_facecolor(color)
    box.set_alpha(0.75)

for i, (vals, color) in enumerate(zip(data_list, colors_list), 1):
    jitter = np.random.uniform(-0.08, 0.08, len(vals))
    ax1.scatter([i + j for j in jitter], vals, alpha=0.35, s=12, color=color, zorder=3)

ax1.set_ylabel("KRAS Expression\nlog2(counts + 1)", fontsize=11)
ax1.set_title(f"KRAS Expression in TCGA-PAAD\n{p_label}", fontsize=12, fontweight="bold")
ax1.set_xticks(range(1, len(labels_list) + 1))
ax1.set_xticklabels(labels_list)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# === 右图：DepMap KRAS 依赖性 ===
ax2 = axes[1]
stat2, p2 = mannwhitneyu(pancreas, others, alternative="less")
p2_label = f"p = {p2:.2e}"

bp2 = ax2.boxplot([others.values, pancreas.values], patch_artist=True, widths=0.5, medianprops=dict(color="black", linewidth=2))
bp2["boxes"][0].set_facecolor("#cccccc")
bp2["boxes"][1].set_facecolor("#C00000")

for i, (vals, color) in enumerate([(others, "#999999"), (pancreas, "#C00000")], 1):
    jitter = np.random.uniform(-0.08, 0.08, len(vals))
    ax2.scatter([i + j for j in jitter], vals, alpha=0.25, s=8, color=color, zorder=3)

ax2.set_xticks([1, 2])
ax2.set_xticklabels([f"Other lineages\n(n={len(others)})", f"Pancreas\n(n={len(pancreas)})"])
ax2.set_ylabel("Chronos Dependency Score\n(Lower = More Essential)", fontsize=11)
ax2.set_title(f"KRAS Dependency in DepMap\n{p2_label}", fontsize=12, fontweight="bold")
ax2.axhline(-1.0, color="gray", linestyle="--", linewidth=1)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.suptitle("KRAS Cross-Database Validation", fontsize=14, fontweight="bold")
plt.tight_layout()

# 保存到桌面
save_path = os.path.join(desktop, "Fig12_Universal_Validation.png")
plt.savefig(save_path, dpi=300, bbox_inches="tight")
print(f"✅ 运行成功！图片已保存到桌面：{save_path}")
plt.show()
