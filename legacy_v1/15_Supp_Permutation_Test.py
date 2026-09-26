import pandas as pd
import numpy as np
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

# === 路径配置 ===
# 使用 "." [span_1](start_span)代表当前文件夹，保护隐私并确保代码在任何电脑上都能运行[span_1](end_span)
desktop_path = "." 

INPUT_DATA  = os.path.join(desktop_path, "standardized_matrix.csv")
INPUT_VAR   = os.path.join(desktop_path, "variance_rank.csv")
OUTPUT_FIG  = os.path.join(desktop_path, "Fig_Supp_Permutation_Final.png")
OUTPUT_DIST = os.path.join(desktop_path, "permutation_distribution_1000.csv")

# === 数据加载 ===
# 检查文件是否存在，防止运行报错
if not os.path.exists(INPUT_DATA) or not os.path.exists(INPUT_VAR):
    print(f"错误：在当前目录下找不到输入文件！")
    print(f"请确保 '{os.path.basename(INPUT_DATA)}' 和 '{os.path.basename(INPUT_VAR)}' 已放入代码所在的文件夹。")
else:
    # 读取数据
    data = pd.read_csv(INPUT_DATA)
    var_rank = pd.read_csv(INPUT_VAR)
    print("数据加载成功，开始执行后续分析...")

# === 数据加载 ===
var_df     = pd.read_csv(INPUT_VAR)
top_genes  = var_df["Gene"].tolist()
df_raw     = pd.read_csv(INPUT_DATA, usecols=["primary_disease"] + top_genes)
df_raw     = df_raw.dropna(subset=["primary_disease"])
df_raw["primary_disease"] = df_raw["primary_disease"].str.strip().str.lower()

data_matrix    = df_raw[top_genes].values.astype(float)  # shape: (n_samples, n_genes)
disease_labels = df_raw["primary_disease"].values
unique_cancers = np.unique(disease_labels)
N_SAMPLES      = len(disease_labels)

# 关键优化：预计算每个癌种的"其他样本"索引，置换时只重用结构
# 这样内层循环不再重复创建 other_mask
observed_count = 274
N_PERM         = 1000
perm_results   = []

print(f"样本数: {N_SAMPLES}, 基因数: {len(top_genes)}, 癌种数: {len(unique_cancers)}")
print(f"开始 {N_PERM} 次置换检验...")

for i in tqdm(range(N_PERM)):
    shuffled = np.random.permutation(disease_labels)
    
    # 一次性建好索引字典，避免内层重复计算
    idx_dict = {c: np.where(shuffled == c)[0] for c in unique_cancers}
    
    p_vals = []
    
    for g_idx in range(len(top_genes)):
        col = data_matrix[:, g_idx]
        
        for cancer, t_idx in idx_dict.items():
            if len(t_idx) < 5:
                continue
            
            target = col[t_idx]
            target = target[~np.isnan(target)]
            if len(target) < 5:
                continue
            
            # 直接用补集索引，不新建布尔数组
            o_idx  = np.concatenate([v for k, v in idx_dict.items() if k != cancer])
            other  = col[o_idx]
            other  = other[~np.isnan(other)]
            if len(other) < 5:
                continue
            
            med_t = np.median(target)
            med_o = np.median(other)
            if med_t < -1.0 and (med_t - med_o) < -0.3:
                _, p = ranksums(target, other)
                p_vals.append(p)
    
    if p_vals:
        _, q_vals, _, _ = multipletests(p_vals, method='fdr_bh')
        perm_results.append(int(np.sum(q_vals < 0.05)))
    else:
        perm_results.append(0)
    
    # 每100次自动存档，防止中途崩溃丢数据
    if (i + 1) % 100 == 0:
        pd.DataFrame({"sig_counts": perm_results}).to_csv(OUTPUT_DIST, index=False)
        print(f"  已完成 {i+1} 次，当前均值: {np.mean(perm_results):.2f}")

# === 保存最终结果 ===
perm_arr = np.array(perm_results)
pd.DataFrame({"sig_counts": perm_arr}).to_csv(OUTPUT_DIST, index=False)

mu    = np.mean(perm_arr)
sigma = np.std(perm_arr)
z     = (observed_count - mu) / sigma
p_emp = np.mean(perm_arr >= observed_count)
p_label = f"p < {1/N_PERM:.3f}" if p_emp == 0 else f"p = {p_emp:.4f}"

print(f"\n结果: Observed={observed_count}, Mean={mu:.2f}, SD={sigma:.2f}")
print(f"Z-Score={z:.2f}, Empirical {p_label}")

# === 双图输出 ===
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# 左图：全局视图
ax1.hist(perm_arr, bins=20, color='gray', alpha=0.6, edgecolor='black')
ax1.axvline(observed_count, color='red', linestyle='--', linewidth=2,
            label=f'Observed ({observed_count})')
ax1.set_xlabel("Number of Significant Pairs (FDR < 0.05)")
ax1.set_ylabel("Frequency")
ax1.set_title("Global View")
ax1.legend()

# 右图：零假设分布放大
ax2.hist(perm_arr, bins=20, color='skyblue', edgecolor='black', alpha=0.8)
ax2.axvline(mu, color='navy', linestyle='-', linewidth=1.5,
            label=f'Mean = {mu:.1f}')
ax2.axvline(mu + 2*sigma, color='orange', linestyle=':', linewidth=1.5,
            label=f'Mean + 2SD = {mu+2*sigma:.1f}')
ax2.set_xlabel("Number of Significant Pairs (FDR < 0.05)")
ax2.set_title(f"Null Distribution (Zoomed)\nMean={mu:.1f}, SD={sigma:.1f}, Z={z:.1f}")
ax2.legend()

plt.suptitle(f"Permutation Test (N={N_PERM}):  {p_label}", fontsize=13)
plt.tight_layout()
plt.savefig(OUTPUT_FIG, dpi=300, bbox_inches='tight')
plt.show()
print(f"图片已保存: {OUTPUT_FIG}")
