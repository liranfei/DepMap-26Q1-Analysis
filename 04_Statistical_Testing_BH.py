import pandas as pd
import numpy as np
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import os

# === 路径锁定 ===
desktop_path = "/Users/lrf15336328932/Desktop"
INPUT_DATA = os.path.join(desktop_path, "standardized_matrix.csv")
INPUT_CAND = os.path.join(desktop_path, "candidate_pairs.csv")
OUTPUT_04 = os.path.join(desktop_path, "significant_pairs.csv")

print("--- 正在执行 04: 统计检验 (Rank-sum + BH校正) ---")

# 1. 加载数据
print("正在读取数据...")
df = pd.read_csv(INPUT_DATA, low_memory=False)
cand_df = pd.read_csv(INPUT_CAND)

if cand_df.empty:
    print("❌ 错误：候选列表为空，请检查步骤 03。")
    exit()

# 2. 执行统计检验
pvals = []
print(f"正在对 {len(cand_df)} 个候选对进行秩和检验...")

for i, row in cand_df.iterrows():
    gene = row["Gene"]
    cancer = row["Cancer"]
    
    # 提取目标癌症组和其它组
    target_vals = df[df["primary_disease"] == cancer][gene].dropna()
    other_vals = df[df["primary_disease"] != cancer][gene].dropna()
    
    # 执行 Wilcoxon 秩和检验
    _, p = ranksums(target_vals, other_vals)
    pvals.append(p)

# 3. 多重假设检验校正 (BH法/FDR)
print("正在进行 BH 校正 (FDR)...")
cand_df["p_value"] = pvals
# alpha=0.05
reject, qvals, _, _ = multipletests(pvals, method='fdr_bh')
cand_df["q_value"] = qvals

# 4. 筛选显著结果 (q < 0.05)
sig_df = cand_df[cand_df["q_value"] < 0.05].copy()
# 按照选择性排序，把最强的靶点排在前面
sig_df = sig_df.sort_values("Selectivity")

# 5. 保存最终结果
sig_df.to_csv(OUTPUT_04, index=False)

print(f"\n✅ 任务全部完成！")
print(f"最终显著基因对数量 (q < 0.05): {len(sig_df)}")
print(f"最终结果已保存至: {OUTPUT_04}")

# 预览前 10 个最显著的靶点
if not sig_df.empty:
    print("\n--- 显著靶点预览 (Top 10) ---")
    print(sig_df[['Gene', 'Cancer', 'Chronos_median', 'Selectivity', 'q_value']].head(10))
