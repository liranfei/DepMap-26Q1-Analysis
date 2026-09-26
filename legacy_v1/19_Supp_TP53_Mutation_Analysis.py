import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# 1. 自动定位桌面
home = os.path.expanduser("~")
desktop = os.path.join(home, "Desktop")

# 2. 读取数据
print("正在加载 26Q1 突变与矩阵数据...")
df = pd.read_csv(os.path.join(desktop, "standardized_matrix.csv"))

# 读取突变文件 (智能识别 OmicsSomaticMutations.csv)
mut_path = os.path.join(desktop, "OmicsSomaticMutations.csv")
if not os.path.exists(mut_path):
    print(f" 错误：在桌面找不到 OmicsSomaticMutations.csv")
    exit()

mut = pd.read_csv(mut_path, low_memory=False)

# 3. 筛选 TP53 突变状态
# 26Q1 版本的列名通常是 'ModelID' 和 'HugoSymbol'
id_col = 'ModelID' if 'ModelID' in mut.columns else 'DepMap_ID'
gene_col = 'HugoSymbol' if 'HugoSymbol' in mut.columns else 'hugo_symbol'

# 提取所有携带 TP53 突变的 ModelID
tp53_mut_models = mut[mut[gene_col] == "TP53"][id_col].unique()

# 4. 对接矩阵并打标签
# 矩阵中的 ID 列通常是第一列
dep_id_col = df.columns[0] 
mdm2_col = [c for c in df.columns if c.startswith("MDM2 (")][0]

df["TP53_Status"] = np.where(df[dep_id_col].isin(tp53_mut_models), "Mutant", "Wild-type")

# 5. 统计分析
mdm2_wt = df[df["TP53_Status"] == "Wild-type"][mdm2_col].dropna()
mdm2_mut = df[df["TP53_Status"] == "Mutant"][mdm2_col].dropna()

# 生物学逻辑：MDM2 在 WT 中更重要（更负），所以我们检验 WT 是否 < Mut
stat, p = mannwhitneyu(mdm2_wt, mdm2_mut, alternative="two-sided")
p_label = f"p = {p:.2e}" if p >= 0.0001 else "p < 1e-4"

print(f"✅ 数据匹配完成:")
print(f"TP53 WT 组样本量: {len(mdm2_wt)}, 中位数: {np.median(mdm2_wt):.3f}")
print(f"TP53 Mut 组样本量: {len(mdm2_mut)}, 中位数: {np.median(mdm2_mut):.3f}")

# 6. 绘图 (Figure 19)
plt.figure(figsize=(10, 6))
sns.set_style("white")

# 使用专业的混合图：小提琴图 + 箱线图
palette = {"Wild-type": "#6baed6", "Mutant": "#ef3b2c"}
sns.violinplot(x="TP53_Status", y=mdm2_col, data=df, palette=palette, inner=None, alpha=0.3)
sns.boxplot(x="TP53_Status", y=mdm2_col, data=df, palette=palette, width=0.2, showcaps=True, 
            boxprops={'zorder': 2}, showfliers=False)

plt.axhline(-1.0, color="gray", linestyle="--", alpha=0.5)
plt.title(f"Target Validation: MDM2 Dependency vs TP53 Status\n(Mann-Whitney Test {p_label})", fontsize=14)
plt.ylabel("Chronos Dependency Score", fontsize=12)
plt.xlabel("TP53 Mutation Status", fontsize=12)

# 移除多余边框
sns.despine()

# 保存
save_path = os.path.join(desktop, "Fig19_TP53_MDM2_Validation.png")
plt.savefig(save_path, dpi=300, bbox_inches="tight")
plt.show()

print(f" 脚本 19 运行成功！图片已保存至桌面。")
