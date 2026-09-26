import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

desktop = os.path.expanduser("~/Desktop")

df_dep = pd.read_csv(f"{desktop}/standardized_matrix.csv", low_memory=False)
df_mut = pd.read_csv(f"{desktop}/OmicsSomaticMutations.csv", low_memory=False)

# 找ID列
id_col_dep = "DepMap_ID" if "DepMap_ID" in df_dep.columns else df_dep.columns[0]
id_col_mut = "ModelID" if "ModelID" in df_mut.columns else df_mut.columns[0]
gene_col = "HugoSymbol" if "HugoSymbol" in df_mut.columns else "Gene"

kras_mut = set(df_mut[df_mut[gene_col]=="KRAS"][id_col_mut].unique())

kras_dep_col = "KRAS (3845)"
dep = df_dep[[id_col_dep, kras_dep_col]].dropna()
dep["group"] = dep[id_col_dep].apply(
    lambda x: "Mutant" if x in kras_mut else "Wild-type")

mut_vals = dep[dep["group"]=="Mutant"][kras_dep_col].values
wt_vals  = dep[dep["group"]=="Wild-type"][kras_dep_col].values

stat, p = mannwhitneyu(mut_vals, wt_vals, alternative="less")

# 计算中位数
print(f"突变组 n={len(mut_vals)}, 中位数={np.median(mut_vals):.3f}")
print(f"野生组 n={len(wt_vals)}, 中位数={np.median(wt_vals):.3f}")
print(f"Mann-Whitney U p = {p:.2e}")

# 画箱线图
fig, ax = plt.subplots(figsize=(6,5))
bp = ax.boxplot([mut_vals, wt_vals],
                patch_artist=True,
                medianprops=dict(color="black", linewidth=2),
                widths=0.5)
colors = ["#C00000", "#4472C4"]
for patch, color in zip(bp["boxes"], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.8)

ax.set_xticks([1,2])
ax.set_xticklabels([f"KRAS Mutant\n(n={len(mut_vals)})",
                    f"KRAS Wild-type\n(n={len(wt_vals)})"], fontsize=11)
ax.set_ylabel("KRAS Chronos Dependency Score", fontsize=11)
ax.set_title("KRAS Mutation Status vs. Dependency Score", fontsize=12, fontweight="bold")

# 加显著性标注
y_max = max(mut_vals.max(), wt_vals.max()) + 0.1
ax.plot([1,1,2,2],[y_max+0.05, y_max+0.1, y_max+0.1, y_max+0.05], "k-", lw=1.2)
ax.text(1.5, y_max+0.12, f"p = {p:.2e}", ha="center", fontsize=10)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig(f"{desktop}/Fig_KRAS_mutation_dependency.png", dpi=300)
plt.show()
print("图已保存")

# 同样分析CTNNB1
ctnnb1_mut = set(df_mut[df_mut[gene_col]=="CTNNB1"][id_col_mut].unique())
ctnnb1_dep_col = "CTNNB1 (1499)"
if ctnnb1_dep_col in df_dep.columns:
    dep2 = df_dep[[id_col_dep, ctnnb1_dep_col]].dropna()
    dep2["group"] = dep2[id_col_dep].apply(
        lambda x: "Mutant" if x in ctnnb1_mut else "Wild-type")
    m2 = dep2[dep2["group"]=="Mutant"][ctnnb1_dep_col].values
    w2 = dep2[dep2["group"]=="Wild-type"][ctnnb1_dep_col].values
    _, p2 = mannwhitneyu(m2, w2, alternative="less")
    print(f"\nCTNNB1: 突变组n={len(m2)}, 中位数={np.median(m2):.3f}")
    print(f"CTNNB1: 野生组n={len(w2)}, 中位数={np.median(w2):.3f}")
    print(f"CTNNB1 p = {p2:.2e}")
