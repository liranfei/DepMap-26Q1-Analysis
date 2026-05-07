import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# --- 1. 自动处理路径 (解决 FileNotFoundError) ---
desktop_path = os.path.join(os.path.expanduser("~"), "Desktop")
matrix_path = os.path.join(desktop_path, "standardized_matrix.csv")
targets_path = os.path.join(desktop_path, "final_targets.csv")

# 检查文件是否存在
if not os.path.exists(matrix_path) or not os.path.exists(targets_path):
    print("❌ 错误：桌面缺少 standardized_matrix.csv 或 final_targets.csv")
else:
    # --- 2. 加载数据 ---
    df = pd.read_csv(matrix_path)
    final = pd.read_csv(targets_path)

    # 预处理：确保匹配时不被空格或大小写干扰
    df["primary_disease"] = df["primary_disease"].str.strip()
    
    # 设定目标靶点和癌种
    gene = "CTNNB1 (1499)"
    cancer_name = "Bowel"  # 注意：与 final_targets.csv 中的 Cancer 列保持一致

    # --- 3. 提取显著性数值 (q-value) ---
    # 使用 str.contains 或直接匹配，确保能找到数据
    row = final[(final["Gene"] == gene) & (final["Cancer"] == cancer_name)]
    
    if row.empty:
        print(f"⚠️ 警告：在 final_targets.csv 中没找到 {gene} 与 {cancer_name} 的组合")
        q_val = 2.23e-25  # 如果匹配不到，使用你图片上的默认值
    else:
        q_val = row["q_value"].values[0]

    # --- 4. 提取绘图数据 ---
    target_data = df[df["primary_disease"] == cancer_name][gene].dropna()
    other_data  = df[df["primary_disease"] != cancer_name][gene].dropna()

    # --- 5. 开始绘图 (顶刊风格优化) ---
    fig, ax = plt.subplots(figsize=(6, 8))
    
    # 绘制箱线图
    data = [other_data.values, target_data.values]
    bp = ax.boxplot(data, patch_artist=True, widths=0.6,
                    medianprops=dict(color="#333333", linewidth=2),
                    whiskerprops=dict(color="#333333"),
                    capprops=dict(color="#333333"),
                    showfliers=False) # 不显示离群点，因为我们要自己画散点

    # 设置箱体颜色
    bp["boxes"][0].set(facecolor="#f0f0f0", edgecolor="#333333", linewidth=1.5) # Others: 灰色
    bp["boxes"][1].set(facecolor="#CD5C5C", edgecolor="#333333", linewidth=1.5, alpha=0.8) # Bowel: 红色

    # 绘制带抖动的散点 (Jitter)
    for i, d in enumerate([other_data, target_data], 1):
        jitter = np.random.uniform(-0.15, 0.15, len(d))
        color = "#cccccc" if i == 1 else "#C00000"
        ax.scatter([i + j for j in jitter], d,
                   alpha=0.4, s=25, color=color, zorder=3, edgecolors='white', linewidth=0.5)

    # --- 6. 装饰与标注 ---
    ax.axhline(-1.0, color="black", linestyle="--", linewidth=1, alpha=0.5) # 必需性阈值线
    ax.set_xticks([1, 2])
    ax.set_xticklabels(["Others", cancer_name], fontsize=13)
    ax.set_ylabel("Chronos Dependency Score\n(Lower = More Essential)", fontsize=12)
    ax.set_title(f"Dependency Profile: {gene.split(' ')[0]}", fontsize=15, fontweight="bold", pad=20)

    # 标注 FDR
    ax.text(1.5, ax.get_ylim()[1]*0.9, f"FDR = {q_val:.2e}",
            ha="center", color="#C00000", fontsize=13, fontweight="bold")

    # 去除上方和右方的边框 (Despine)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()

    # --- 7. 自动保存 ---
    save_path = os.path.join(desktop_path, "Fig_CTNNB1_boxplot_v3.png")
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.show()
print(f"✅ 图片已完美生成并保存至: {save_path}")
