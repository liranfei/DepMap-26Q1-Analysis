import pandas as pd
import os
import time

# === 路径锁定 ===
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
INPUT_01 = os.path.join(desktop_path, "standardized_matrix.csv")
OUTPUT_02 = os.path.join(desktop_path, "variance_rank.csv")

print("--- 正在执行 02: 方差特征筛选 ---")

# 增加一个检查机制
if not os.path.exists(INPUT_01):
    print(f"正在等待系统释放文件...")
    time.sleep(2) # 等待2秒，防止系统文件锁

if os.path.exists(INPUT_01):
    print(f"✅ 已找到文件，开始读取: {INPUT_01}")
    
    # 1. 读取数据
    # 指定 low_memory=False 防止大文件读取警告
    df = pd.read_csv(INPUT_01, low_memory=False)
    
    # 2. 识别基因列
    # 排除 ID 和 癌症分类列
    exclude_cols = ["DepMap_ID", "primary_disease"]
    gene_cols = [c for c in df.columns if c not in exclude_cols]
    print(f"数据读取成功！行数: {len(df)}, 基因列数: {len(gene_cols)}")

    # 3. 计算方差
    print("正在计算方差并筛选前 500 个基因...")
    # 强制将基因列转为浮点型计算，确保不报错
    variance_series = df[gene_cols].apply(pd.to_numeric, errors='coerce').var()

    # 4. 排序并取 Top 500
    variance_df = variance_series.sort_values(ascending=False).reset_index()
    variance_df.columns = ["Gene", "Variance"]
    top500 = variance_df.head(500)

    # 5. 保存
    top500.to_csv(OUTPUT_02, index=False)
    print(f"✅ 步骤 02 完成！Top 500 基因列表已保存。")
    print(f"结果路径: {OUTPUT_02}")
    
    print("\n--- 方差排名前 5 的基因预览 ---")
    print(top500.head())
else:
    print(f"❌ 错误：依然找不到文件。")
    print(f"请检查桌面（Desktop）上是否真的有一个名为 standardized_matrix.csv 的文件？")
