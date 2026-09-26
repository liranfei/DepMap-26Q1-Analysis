import pandas as pd
import os

# === 路径锁定 ===
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
INPUT_GENE = os.path.join(desktop_path, "gene_effect.csv")
INPUT_META = os.path.join(desktop_path, "sample_info.csv")
OUTPUT_01 = os.path.join(desktop_path, "standardized_matrix.csv")

print("--- 01: 深度对齐标准化 ---")

# 1. 加载基因数据
print("读取基因文件...")
# 先不设 index_col，看看第一列到底是什么
df_gene = pd.read_csv(INPUT_GENE)
first_col_name = df_gene.columns[0]
print(f"检测到基因文件第一列名为: {first_col_name}")

# 强制重命名第一列为 DepMap_ID 并去空格
df_gene = df_gene.rename(columns={first_col_name: "DepMap_ID"})
df_gene["DepMap_ID"] = df_gene["DepMap_ID"].astype(str).str.strip()

# 2. 加载样本信息
print("读取样本信息...")
meta = pd.read_csv(INPUT_META)

# 自动找 ID 列：只要列名包含 'ID' 的都试一遍
potential_id_cols = [c for c in meta.columns if 'ID' in c.upper()]
if not potential_id_cols:
    id_col = meta.columns[0]
else:
    id_col = potential_id_cols[0]

# 自动找癌症分类列
potential_dis_cols = [c for c in meta.columns if 'LINEAGE' in c.upper() or 'DISEASE' in c.upper()]
dis_col = potential_dis_cols[0] if potential_dis_cols else meta.columns[1]

print(f"匹配 ID 使用列: {id_col}")
print(f"匹配癌症类型使用列: {dis_col}")

meta_clean = meta[[id_col, dis_col]].copy()
meta_clean.columns = ['DepMap_ID', 'primary_disease']
meta_clean['DepMap_ID'] = meta_clean['DepMap_ID'].astype(str).str.strip()

# 3. 核心合并测试
print("正在尝试合并...")
df_final = df_gene.merge(meta_clean, on="DepMap_ID")

# 4. 诊断输出
if len(df_final) == 0:
    print("❌ 合并依然失败！原因诊断：")
    print(f"基因文件ID样板: {list(df_gene['DepMap_ID'].head(3))}")
    print(f"样本文件ID样板: {list(meta_clean['DepMap_ID'].head(3))}")
    print("请对比上面两个 ID 的长相，是不是一个带 ACH- 而另一个不带？")
else:
    print(f"✅ 成功！合并后有 {len(df_final)} 行数据。")
    # 保存结果
    df_final.to_csv(OUTPUT_01, index=False)
    print(f"文件已保存至: {OUTPUT_01}")
