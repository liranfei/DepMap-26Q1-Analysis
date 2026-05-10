import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 1. 数据准备
data_26q1 = {
    "Gene": ["KRAS", "SOX10", "NMNAT1", "NAMPT", "CTNNB1", "MDM2", "CBFB", "HNF1B", "IRF4", "ADSL"],
    "Target_Cancer": ["Pancreas", "Skin", "Lymphoid", "Lymphoid", "Bowel", "Eye", "Myeloid", "Kidney", "Lymphoid", "Lymphoid"],
    "Selectivity_Score": [-1.63, -1.38, -0.94, -0.86, -1.18, -1.37, -1.03, -1.15, -0.96, -0.82]
}

df_plot = pd.DataFrame(data_26q1).sort_values("Selectivity_Score")

# 2. 绘图
plt.figure(figsize=(10, 6))
sns.set_theme(style="whitegrid")

# 修正后的部分：指定 hue 并关闭 legend
ax = sns.barplot(
    x="Selectivity_Score", 
    y="Gene", 
    data=df_plot, 
    hue="Gene",          # 新增：将 y 变量赋值给 hue
    palette="Reds_r", 
    legend=False         # 新增：关闭不需要的图例
)

# 3. 添加癌症类型标注
for i, (score, cancer) in enumerate(zip(df_plot["Selectivity_Score"], df_plot["Target_Cancer"])):
    ax.annotate(f' {cancer}', 
                (score, i), 
                va='center', fontsize=10, fontweight='bold')

plt.title("Top 10 Lineage-Specific Dependencies (DepMap 2026Q1 Release)", fontsize=14, fontweight='bold')
plt.xlabel("Selectivity Score (Lower is more essential)", fontsize=12)
plt.ylabel("Gene Symbol", fontsize=12)
plt.xlim(df_plot["Selectivity_Score"].min() - 0.3, 0) 

plt.tight_layout()
plt.savefig("Fig_26Q1_Selectivity_Clean.png", dpi=300)
plt.show()
