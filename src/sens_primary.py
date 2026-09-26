import os
import pandas as pd, numpy as np, json
from statsmodels.stats.multitest import multipletests
import run_pipeline as rp
O=os.environ.get("RESULTS_DIR","results")+"/"; T=pd.read_csv(O+"all_tests.csv"); F=pd.read_csv(O+"final_targets.csv"); R={}
R["families"]={"n_tests":len(T),"BH_q<0.05":int((T.q_value<0.05).sum()),"final":len(F)}
# alternative multiplicity procedures on the full family, then the same effect-size rule
eff=(T.Chronos_median<-1)&(T.Selectivity<-0.5)
for name,m in [("BH","fdr_bh"),("BY","fdr_by"),("Holm","holm"),("Bonferroni","bonferroni")]:
    q=multipletests(T.p_value,method=m)[1]; sig=q<0.05; f=T[sig&eff]
    R["multiplicity_"+name]={"n_significant":int(sig.sum()),"n_final":len(f),"n_genes":f.Gene.nunique(),"n_lineages":f.Lineage.nunique()}
# effect-size threshold sensitivity (BH q<0.05)
for sel in [-0.3,-0.5,-0.7,-1.0]:
    f=T[(T.q_value<0.05)&(T.Chronos_median<-1)&(T.Selectivity<sel)]; R[f"selectivity<{sel}"]={"n_final":len(f),"n_genes":f.Gene.nunique(),"n_lineages":f.Lineage.nunique()}
for cm in [-0.5,-1.0,-1.5]:
    f=T[(T.q_value<0.05)&(T.Chronos_median<cm)&(T.Selectivity<-0.5)]; R[f"chronos_median<{cm}"]={"n_final":len(f),"n_genes":f.Gene.nunique()}
for mn in [5,10,20]:
    f=T[(T.q_value<0.05)&eff&(T.n_target>=mn)]; R[f"min_n_target>={mn}"]={"n_final":len(f),"n_genes":f.Gene.nunique(),"n_lineages":f.Lineage.nunique()}
# top-N sensitivity on the cancer-only cohort (same variance ranking rule)
df=rp.load(os.environ.get("DEPMAP_DIR",".")); df=df[df.primary_disease!="Non-Cancerous"].reset_index(drop=True)
genes=[c for c in df.columns if c not in ("DepMap_ID","lineage","primary_disease")]
v=df[genes].apply(pd.to_numeric,errors="coerce").var().rename("Variance").reset_index().rename(columns={"index":"Gene"})
rank=v.sort_values(["Variance","Gene"],ascending=[False,True],kind="mergesort"); rows=[]
for n in [300,500,1000,2000,3000,5000]:
    G=set(rank.head(n).Gene); S=T[T.Gene.isin(G)].copy(); S["q"]=multipletests(S.p_value,method="fdr_bh")[1]
    f=S[(S.q<0.05)&(S.Chronos_median<-1)&(S.Selectivity<-0.5)]; rows.append(dict(top_n=n,n_tests=len(S),n_q05=int((S.q<0.05).sum()),n_final=len(f),n_genes=f.Gene.nunique(),n_lineages=f.Lineage.nunique()))
rows.append(dict(top_n="all",n_tests=len(T),n_q05=int((T.q_value<0.05).sum()),n_final=len(F),n_genes=F.Gene.nunique(),n_lineages=F.Lineage.nunique()))
pd.DataFrame(rows).to_csv(O+"sensitivity_topN.csv",index=False); R["topN"]=rows
json.dump(R,open(O+"sensitivity_summary.json","w"),indent=1,default=float)
for k,v_ in R.items(): print(k,v_)
