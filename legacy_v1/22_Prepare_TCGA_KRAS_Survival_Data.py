import pandas as pd
import numpy as np
import os

# 1. 路径处理
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
# 检查文件是否存在
tcga_path = os.path.join(desktop, "TCGA-PAAD.star_counts.tsv")
surv_path = os.path.join(desktop, "PAAD_survival.txt")

# 2. 读取数据
# 表达矩阵通常较大，设置 check.names=False 类似逻辑
tcga = pd.read_csv(tcga_path, sep="\t", index_col=0)
surv = pd.read_csv(surv_path, sep="\t")

# 3. 提取 KRAS (ENSG00000133703) 表达量
# 使用更稳健的索引匹配
kras_rows = [i for i in tcga.index if "ENSG00000133703" in str(i)]
if not kras_rows:
    print("错误：未在矩阵中找到 KRAS 基因 ID。")
    exit()

# 转换数据为数值型并计算 log2(count + 1)
kras_expr = pd.to_numeric(tcga.loc[kras_rows[0]], errors="coerce")
kras_expr = np.log2(kras_expr + 1)

# 4. 筛选肿瘤样本并处理 ID 匹配
# 仅保留末尾为 '01' 的肿瘤样本
tumor_mask = kras_expr.index.str.split('-').str[3].str.startswith('01')
kras_tumor = kras_expr[tumor_mask].dropna()

# 将 'TCGA-XX-XXXX-01A...' 统一截取为 'TCGA-XX-XXXX' 以便和生存表匹配
kras_tumor.index = kras_tumor.index.str[:12]

# 5. 上下四分位数分组
q25 = np.percentile(kras_tumor, 25)
q75 = np.percentile(kras_tumor, 75)

high_ids = kras_tumor[kras_tumor >= q75].index.tolist()
low_ids  = kras_tumor[kras_tumor <= q25].index.tolist()

# 6. 生存数据列名自动识别
id_col = surv.columns[0] # 通常是 sample 或 _PATIENT
# 寻找 OS 状态和时间列
os_time_col = [c for c in surv.columns if "time" in c.lower() and "os" in c.lower()][0]
os_stat_col = [c for c in surv.columns if c.lower() in ["os", "vital_status"]][0]

# 7. 合并分组信息
surv["group"] = "Unknown"
# 注意：surv 里的 ID 可能也需要截取（如果它是 15 位或 16 位的话）
surv[id_col] = surv[id_col].str[:12]

surv.loc[surv[id_col].isin(high_ids), "group"] = "High (top 25%)"
surv.loc[surv[id_col].isin(low_ids),  "group"] = "Low (bottom 25%)"

# 8. 清理最终分析集
sf = surv[surv["group"] != "Unknown"].copy()
sf[os_time_col] = pd.to_numeric(sf[os_time_col], errors="coerce")
sf[os_stat_col] = pd.to_numeric(sf[os_stat_col], errors="coerce")
sf = sf.dropna(subset=[os_time_col, os_stat_col])

# 9. 输出结果预览
print("-" * 30)
print(f"KRAS 表达量分组完成：")
print(f"High 组样本数: {(sf['group']=='High (top 25%)').sum()}")
print(f"Low 组样本数: {(sf['group']=='Low (bottom 25%)').sum()}")
print("\n各组生存时间中位数（天）:")
print(sf.groupby("group")[os_time_col].median())
print("-" * 30)

# 保存结果供脚本 13 画图使用
sf.to_csv(os.path.join(desktop, "KRAS_survival_final.csv"), index=False)
print("已生成 KRAS_survival_final.csv，可直接用于生存曲线绘制。")
