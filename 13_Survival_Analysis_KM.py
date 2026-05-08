import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

# 1. 自动定位路径
desktop = os.path.join(os.path.expanduser("~"), "Desktop")
file_surv = os.path.join(desktop, "PAAD_survival.txt")
file_tcga = os.path.join(desktop, "TCGA-PAAD.star_counts.tsv")

def run_final_km_with_p():
    # 2. 读取数据
    surv = pd.read_csv(file_surv, sep="\t")
    tcga = pd.read_csv(file_tcga, sep="\t", index_col=0)

    # 3. 提取 KRAS 并对齐 ID
    kras_row = [i for i in tcga.index if "ENSG00000133703" in str(i)][0]
    expr = np.log2(tcga.loc[kras_row].apply(pd.to_numeric, errors="coerce") + 1)
    
    # 匹配病人 ID (前 12 位)
    exp_df = pd.DataFrame(expr[expr.index.str.contains("-01")]).reset_index()
    exp_df.columns = ['sample', 'val']
    exp_df['pid'] = exp_df['sample'].str[:12]
    surv['pid'] = surv[surv.columns[0]].str[:12]

    merged = pd.merge(exp_df, surv, on='pid')
    merged['group'] = np.where(merged['val'] >= merged['val'].median(), "High", "Low")

    # 4. 统计 Log-rank P 值
    def logrank_calc(df):
        # 识别列
        t_col = [c for c in df.columns if "time" in c.lower() and "os" in c.lower()][0]
        e_col = [c for c in df.columns if c.lower() in ["os", "vital_status", "_os"]][0]
        df[t_col] = pd.to_numeric(df[t_col], errors='coerce')
        df[e_col] = pd.to_numeric(df[e_col], errors='coerce')
        df = df.dropna(subset=[t_col, e_col])
        
        times = np.sort(df[df[e_col] == 1][t_col].unique())
        O1, E1, V = 0, 0, 0
        for t in times:
            n1 = ((df['group']=="High") & (df[t_col]>=t)).sum()
            n2 = ((df['group']=="Low")  & (df[t_col]>=t)).sum()
            d1 = ((df['group']=="High") & (df[t_col]==t) & (df[e_col]==1)).sum()
            d2 = ((df['group']=="Low")  & (df[t_col]==t) & (df[e_col]==1)).sum()
            n, d = n1 + n2, d1 + d2
            if n > 1:
                O1 += d1
                E1 += n1 * d / n
                V += n1 * n2 * d * (n - d) / (n**2 * (n - 1))
        p = 1 - chi2.cdf((O1 - E1)**2 / V, df=1) if V > 0 else 1
        return p, df, t_col, e_col

    p_val, final_df, os_t, os_e = logrank_calc(merged)

    # 5. 绘图
    plt.figure(figsize=(9, 6))
    for g, c in zip(["High", "Low"], ["#C00000", "#2171b5"]):
        d = final_df[final_df['group'] == g].sort_values(os_t)
        t_plot, s_plot = [0], [1.0]
        n, s = len(d), 1.0
        for _, r in d.iterrows():
            if r[os_e] == 1:
                s *= (1 - 1/n)
                t_plot.append(r[os_t] / 30.44)
                s_plot.append(s)
            n -= 1
        plt.step(t_plot, s_plot, where="post", color=c, linewidth=2.5, label=f"KRAS {g} (n={len(d)})")

    # 在图上强制加上 P 值标注
    plt.text(0.05, 0.15, f"Log-rank P = {p_val:.4f}", transform=plt.gca().transAxes, 
             fontsize=12, fontweight='bold', bbox=dict(facecolor='white', alpha=0.8))

    plt.title("Kaplan-Meier Survival Analysis: KRAS Expression (TCGA-PAAD)", fontweight='bold')
    plt.xlabel("Time (Months)")
    plt.ylabel("Survival Probability")
    plt.legend()
    plt.savefig(os.path.join(desktop, "Fig13_Final_with_P_Value.png"), dpi=300)
    plt.show()

run_final_km_with_p()
