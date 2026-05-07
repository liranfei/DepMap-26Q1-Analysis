# DepMap-26Q1-Analysis Repository for CBC Submission

这个仓库包含了针对 DepMap (24Q2/26Q1) 数据集进行谱系特异性依赖 (Lineage-Specific Dependency) 筛选的完整分析工作流。

## 1. 运行环境 (Environment)
* **语言**: Python 3.9.13+
* **核心依赖**: 
    * `pandas`, `numpy` (数据处理)
    * `scipy`, `statsmodels` (统计检验与多重校正)
    * `matplotlib`, `matplotlib-venn` (数据可视化)
    * `openpyxl` (Excel 报表生成)
* **安装方法**: `pip install -r requirements.txt`

## 2. 数据来源 (Data Source)
* **原始数据**: 来源于 [DepMap Portal](https://depmap.org/portal/download/).
* **涉及文件**: 
    1. `standardized_matrix.csv` (标准化后的 Chronos 效应评分)
    2. 相关的癌症谱系标注 (Lineage annotations)

## 3. 脚本执行顺序 (Workflow)
请按以下顺序运行脚本以复现论文中的结果：
* **[01-02]**: 预处理与高方差特征筛选 (Top 500 Genes).
* **[05-09]**: 核心筛选算法与 Wilcoxon 秩和检验，生成 `final_targets.csv`.
* **[07-10]**: 生成火山图、散点图及通路富集气泡图.
* **[14-15]**: 敏感性分析 (Robustness Check) 与实验流程图绘制.
* **[17-18]**: 靶点必需性验证 (Essentiality) 与协同依赖热图分析.
* **[S1-Export]**: 导出符合要求的 Table S1 附表.

## 4. 可复现性声明 (Reproducibility)
所有脚本均已处理路径依赖（采用 `os.path` 自动识别环境）。用户只需将原始数据置于指定目录，依次运行脚本，即可完全复现论文中 Figure 1 及其相关补充材料的所有结果与数值。

---
*如有任何疑问，请通过论文通讯作者邮箱或仓库 Issue 提交。*
