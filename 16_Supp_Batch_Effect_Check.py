import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal

home = os.path.expanduser("~")
desktop = os.path.join(home, "Desktop")

# 1. 读取文件
print("正在读取 26Q1 核心数据...")
df = pd.read_csv(os.path.join(desktop, "standardized_matrix.csv"))
info = pd.read_csv(os.path.join(desktop, "sample_info.csv"))
ssm = pd.read_csv(os.path.join(desktop, "ScreenSequenceMap.csv")) 
final = pd.read_csv(os.path.join(desktop, "final_targets.csv")) 

# 2. 统一 ID
id_col = "ModelID" if "ModelID" in info.columns else "DepMap_ID"
info = info.rename(columns={id_col: "DepMap_ID"})

# 3. 合并批次信息
merged_info = info.merge(ssm[['ModelID', 'pDNABatch']], left_on='DepMap_ID', right_on='ModelID', how='inner')

# 4. 【核心修复】智能获取基因列名
# 尝试获取第一行第一列的值作为基因名
top_target = str(final.iloc[0, 0]).strip()
print(f"尝试验证的靶点是: {top_target}")

# 在矩阵列名中寻找包含该基因名的列
full_name_list = [c for c in df.columns if str(top_target) in c]

if not full_name_list:
    print(f"❌ 错误：在矩阵中找不到基因 '{top_target}'")
    print(f"矩阵的前 5 个列名参考: {list(df.columns[1:6])}")
    # 兜底方案：如果找不到，就选矩阵里除 ID 外的第一列
    full_name = df.columns[1]
    print(f"⚠️ 自动切换为兜底基因: {full_name}")
else:
    full_name = full_name_list[0]
    print(f"✅ 成功匹配到矩阵列名: {full_name}")

# 5. 合并矩阵数据与批次数据
plot_df = df[['DepMap_ID', full_name]].merge(merged_info[['DepMap_ID', 'pDNABatch']], on='DepMap_ID')

# 6. 【核心修复】过滤异常值：只取样本量大于 10 的大批次
batch_counts = plot_df['pDNABatch'].value_counts()
top_batches = batch_counts[batch_counts > 10].head(5).index.tolist()

if not top_batches:
    print("❌ 错误：没有一个批次的样本量大于 10，无法进行有效对比。")
    exit()

final_plot = plot_df[plot_df['pDNABatch'].isin(top_batches)]

# 7. 统计检验
groups = [final_plot[final_plot['pDNABatch'] == b][full_name].dropna().values for b in top_batches]
if len(groups) >= 2:
    stat, p = kruskal(*groups)
    p_val = f"p = {p:.4f}" if p > 0.001 else "p < 0.001"
else:
    p_val = "N/A (样本量不足)"

# 8. 绘图
plt.figure(figsize=(10, 6))
sns.set_style("whitegrid")
sns.boxplot(x='pDNABatch', y=full_name, data=final_plot, palette="Set3")
# 添加抖动点，增加数据透明度
sns.stripplot(x='pDNABatch', y=full_name, data=final_plot, color=".3", size=3, alpha=0.3)

plt.title(f"Robustness Check (Batch Effect): {full_name}\nKruskal-Wallis {p_val}", fontsize=12)
plt.ylabel("Chronos Dependency Score")
plt.xlabel("Experimental Batch (pDNABatch)")
plt.xticks(rotation=30)
plt.tight_layout()

save_name = "Fig16_Standard_Batch_Analysis.png"
plt.savefig(os.path.join(desktop, save_name), dpi=300)
plt.show()

print(f"🎉 成功！请在桌面查看：{save_name}")
