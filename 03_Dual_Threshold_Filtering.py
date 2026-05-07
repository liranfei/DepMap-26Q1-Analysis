import pandas as pd
import numpy as np
import os

# === 1. 路径精准定位 ===
desktop_path = "/Users/lrf15336328932/Desktop"
INPUT_DATA = os.path.join(desktop_path, "standardized_matrix.csv")
INPUT_VAR = os.path.join(desktop_path, "variance_rank.csv")
OUTPUT = os.path.join(desktop_path, "candidate_pairs.csv")

print("--- Step 03: Dual Threshold Filtering ---")

# --- 路径自动检查模块 ---
if not os.path.exists(INPUT_DATA):
    print(f"\n❌ 找不到文件: {INPUT_DATA}")
    print("👉 正在为你扫描桌面上的所有文件，请检查文件名是否匹配:")
    try:
        files = os.listdir(desktop_path)
        for f in files:
            if "matrix" in f.lower():
                print(f"   找到相似文件: {f}")
    except Exception as e:
        print(f"   读取桌面失败: {e}")
    # 停止运行，防止后续报错
    exit()

# === 2. 读取数据 ===
print("✅ 文件定位成功，正在读取数据...")
df = pd.read_csv(INPUT_DATA, low_memory=False)
var_df = pd.read_csv(INPUT_VAR)

# === 3. 初始化参数 ===
top_genes = var_df["Gene"].tolist()
grouped = df.groupby("primary_disease")
cancer_types = list(grouped.groups.keys())

results = []
total_tests = 0
skipped_small_n = 0
skipped_missing = 0

print(f"开始分析: {len(top_genes)} 基因 × {len(cancer_types)} 癌种")

# === 4. 主循环 (保留你的核心逻辑) ===
for idx, gene in enumerate(top_genes):
    if idx % 50 == 0:
        print(f"进度: {idx}/{len(top_genes)}")

    if gene not in df.columns:
        skipped_missing += len(cancer_types)
        continue

    for cancer in cancer_types:
        total_tests += 1
        
        try:
            target_group = grouped.get_group(cancer)[gene].dropna()
        except KeyError:
            skipped_missing += 1
            continue
            
        other_group = df[df["primary_disease"] != cancer][gene].dropna()

        # 样本量过滤 (n >= 5)
        if len(target_group) < 5:
            skipped_small_n += 1
            continue

        if len(other_group) == 0:
            skipped_missing += 1
            continue

        # --- 核心计算 ---
        med_target = np.median(target_group)
        med_other = np.median(other_group)
        selectivity = med_target - med_other

        # --- 严格双阈值判断 ---
        if med_target < -1.0 and selectivity < -0.3:
            results.append({
                "Gene": gene,
                "Cancer": cancer,
                "Chronos_median": med_target,
                "Selectivity": selectivity,
                "n_target": len(target_group)
            })

# === 5. 保存结果 ===
res_df = pd.DataFrame(results)
res_df.to_csv(OUTPUT, index=False)

print("\n✅ Step 03 处理完成")
print(f"总测试组合: {total_tests}")
print(f"因样本少跳过: {skipped_small_n}")
print(f"最终候选对数量: {len(res_df)}")
print(f"候选列表已保存至: {OUTPUT}")
