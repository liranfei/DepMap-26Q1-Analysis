import os
import json, numpy as np, pandas as pd, requests
from scipy.stats import mannwhitneyu, fisher_exact
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests
O=os.environ.get("RESULTS_DIR","results")+"/"; R={}
mw=lambda x,y,alt="two-sided": float(mannwhitneyu(x,y,alternative=alt)[1])
def tcga(expr_path,surv_path,ensg,label):
    X=pd.read_csv(expr_path,sep="\t",index_col=0); row=[i for i in X.index if ensg in str(i)][0]
    # Xena GDC "STAR - Counts" are log2(raw count + 1) without between-sample normalisation: convert back to counts and normalise to
    # counts per million (library size = sum of counts over all genes of the sample), then log2(CPM + 1)
    C=np.power(2.0,X.apply(pd.to_numeric,errors="coerce"))-1; lib=C.sum(axis=0); e=np.log2(C.loc[row]/lib*1e6+1)
    pd.DataFrame({"sample":e.index,"log2_cpm":e.values,"log2_count_raw":pd.to_numeric(X.loc[row],errors="coerce").values,"library_size":lib.values}).to_csv(O+f"tcga_{label}_expression_samples.csv",index=False)
    code=lambda c:c.split("-")[3][:2]
    tum=e[[c for c in e.index if code(c)=="01"]].dropna(); nor=e[[c for c in e.index if code(c)=="11"]].dropna()
    out={"n_tumor_samples":len(tum),"n_normal_samples":len(nor),"median_tumor":float(tum.median()),"median_normal":float(nor.median()),"p_tumor_vs_normal_two_sided":mw(tum,nor),
         "p_tumor_lt_normal":mw(tum,nor,"less"),"p_tumor_gt_normal":mw(tum,nor,"greater")}
    d=tum.rename("x").reset_index(); d["pid"]=d["index"].str[:12]; d=d.groupby("pid",as_index=False).x.mean()
    s=pd.read_csv(surv_path,sep="\t"); s["pid"]=s["sample"].str[:12]; s=s[s["sample"].str.split("-").str[3].str[:2]=="01"].drop_duplicates("pid")
    m=d.merge(s[["pid","OS.time","OS"]],on="pid").dropna(); m=m[m["OS.time"]>0]
    m["hi"]=(m.x>=m.x.median()).astype(int); lr=logrank_test(m[m.hi==1]["OS.time"],m[m.hi==0]["OS.time"],m[m.hi==1].OS,m[m.hi==0].OS)
    m["z"]=(m.x-m.x.mean())/m.x.std(); c=CoxPHFitter().fit(m[["OS.time","OS","z"]],"OS.time","OS"); q=c.summary.loc["z"]
    c2=CoxPHFitter().fit(m[["OS.time","OS","hi"]],"OS.time","OS"); q2=c2.summary.loc["hi"]
    out.update({"n_survival":len(m),"n_events":int(m.OS.sum()),"n_high":int(m.hi.sum()),"logrank_p_mediansplit":float(lr.p_value),"cox_binary_HR":float(np.exp(q2["coef"])),
      "cox_binary_CI":[float(np.exp(q2["coef lower 95%"])),float(np.exp(q2["coef upper 95%"]))],"cox_binary_p":float(q2["p"]),
      "cox_continuous_HR_per_SD":float(np.exp(q["coef"])),"cox_continuous_CI":[float(np.exp(q["coef lower 95%"])),float(np.exp(q["coef upper 95%"]))],"cox_continuous_p":float(q["p"])})
    m.to_csv(O+f"tcga_{label}_survival_table.csv",index=False); return out
R["TCGA_KIRC_HNF1B"]=tcga(os.environ.get("TCGA_DIR",".")+"/TCGA-KIRC.star_counts.tsv",os.environ.get("TCGA_DIR",".")+"/TCGA-KIRC.survival.tsv","ENSG00000275410","KIRC_HNF1B")
R["TCGA_PAAD_KRAS"]=tcga(os.environ.get("TCGA_DIR",".")+"/TCGA-PAAD.star_counts.tsv",os.environ.get("TCGA_DIR",".")+"/TCGA-PAAD.survival.tsv","ENSG00000133703","PAAD_KRAS")
if os.environ.get("SKIP_ENRICH"): json.dump(R,open(O+"part3_summary_tcga_only.json","w"),indent=1,default=float); print(json.dumps(R,indent=1,default=float)); raise SystemExit
# ---- enrichment with the tested-gene universe as background (own Fisher test + BH)
T=pd.read_csv(O+"all_tests.csv"); F=pd.read_csv(O+"final_targets.csv"); sym=lambda s:s.split(" (")[0]
universe=set(T.Gene.map(sym)); hits=set(F.Gene.map(sym)); res=[]
for lib in ["GO_Biological_Process_2023","KEGG_2021_Human"]:
    txt=requests.get("https://maayanlab.cloud/Enrichr/geneSetLibrary",params={"mode":"text","libraryName":lib},timeout=120).text
    open(O+f"{lib}.gmt","w").write(txt); rows=[]
    for line in txt.strip().split("\n"):
        p=line.split("\t"); term=p[0]; gs={g.split(",")[0] for g in p[2:] if g}&universe
        if 5<=len(gs)<=500:
            a=len(gs&hits); b=len(hits)-a; c=len(gs)-a; d=len(universe)-len(gs)-b; rows.append((lib,term,len(gs),a,fisher_exact([[a,b],[c,d]],alternative="greater")[1],";".join(sorted(gs&hits))))
    E=pd.DataFrame(rows,columns=["library","term","set_size","n_hits","p","genes"]); E["q"]=multipletests(E.p,method="fdr_bh")[1]; res.append(E)
E=pd.concat(res); E.sort_values("q").to_csv(O+"enrichment_custom_background.csv",index=False)
R["enrichment"]={"universe":len(universe),"n_hit_genes":len(hits),"by_library_n_q05":E[E.q<0.05].groupby("library").size().to_dict(),
   "top":{lib:E[E.library==lib].sort_values("q").head(6)[["term","n_hits","set_size","q","genes"]].round(6).values.tolist() for lib in E.library.unique()}}
json.dump(R,open(O+"part3_summary.json","w"),indent=1,default=float)
print(json.dumps(R,indent=1,default=float)[:6000])
