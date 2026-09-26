import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import os

desktop = os.path.join(os.path.expanduser("~"), "Desktop")
input_file = os.path.join(desktop, "standardized_matrix.csv")
final_file = os.path.join(desktop, "final_targets.csv")

if not os.path.exists(input_file) or not os.path.exists(final_file):
    print("❌ 错误：缺少必要的 CSV 文件。")
else:
    df = pd.read_csv(input_file)
    final = pd.read_csv(final_file)

    # 提取靶点基因
    target_genes = final["Gene"].unique().tolist()
    available = [g for g in target_genes if g in df.columns]
    
    # 提取矩阵并处理缺失值
    # 只有至少有一定非空重叠的基因才能计算相关性
    gene_matrix = df[available].apply(pd.to_numeric, errors="coerce")
    
    # 计算 Pearson 相关性，处理缺失数据
    corr_matrix = gene_matrix.corr(method="pearson").fillna(0)

    # 层次聚类排序
    # 距离定义为 1 - 相关系数
    dist = 1 - corr_matrix.values
    dist = np.clip(dist, 0, 2) # 确保距离在合理范围
    np.fill_diagonal(dist, 0)
    
    # 执行聚类
    linkage_matrix = linkage(squareform(dist), method="ward")
    dendro = dendrogram(linkage_matrix, no_plot=True)
    order = dendro["leaves"]

    # 重新排序矩阵
    corr_ordered = corr_matrix.iloc[order, order]
    labels = [g.split(" ")[0] for g in corr_ordered.columns]

    # 绘图
    fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
    im = ax.imshow(corr_ordered.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")

    # 坐标轴设置
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=90, fontsize=8)
    ax.set_yticklabels(labels, fontsize=8)

    # 颜色条
    cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("Pearson Correlation (Chronos Score)", fontsize=10)

    ax.set_title("Figure 18: Co-dependency Network of Target Genes", 
                 fontsize=14, fontweight="bold", pad=20)

    plt.tight_layout()
    save_path = os.path.join(desktop, "Fig18_Heatmap.png")
    plt.savefig(save_path, bbox_inches='tight')
    plt.show()
    print(f"✅ 18 完成，热图保存至: {save_path}")
