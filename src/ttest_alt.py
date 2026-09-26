import os
"""Orthogonal method check: Meyers-2017-style lineage-vs-rest Student's t-test on mean differences (whole genome, cancer-only cohort, same n>=5 rule, BH over full family)."""
import numpy as np, pandas as pd, json
from scipy.stats import t as tdist
import run_pipeline as rp
D=os.environ.get("DEPMAP_DIR","."); RD=os.environ.get("RESULTS_DIR","results"); df=rp.load(D); df=df[df.primary_disease!="Non-Cancerous"].reset_index(drop=True)
genes=[c for c in df.columns if c not in ("DepMap_ID","lineage","primary_disease")]
X=df[genes].apply(pd.to_numeric,errors="coerce").values; ok=~np.isnan(X); Xz=np.where(ok,X,0.0)
lin=df.lineage.values; lc=pd.Series(lin).value_counts(); elig=list(lc[lc>=5].index)
L=np.stack([(lin==e).astype(float) for e in elig]); n1=L@ok; S1=L@Xz; Q1=L@(Xz**2)
N=ok.sum(0)[None,:].astype(float); S=Xz.sum(0)[None,:]; Q=(Xz**2).sum(0)[None,:]; n2=N-n1; S2=S-S1; Q2=Q-Q1
valid=(n1>=5)&(n2>=2)
with np.errstate(all="ignore"):
    m1=S1/n1; m2=S2/n2; v1=(Q1-n1*m1**2)/(n1-1); v2=(Q2-n2*m2**2)/(n2-1)
    sp=((n1-1)*v1+(n2-1)*v2)/(n1+n2-2); se=np.sqrt(sp*(1/n1+1/n2)); tt=(m1-m2)/se; p_two=2*tdist.sf(np.abs(tt),n1+n2-2)
li,gj=np.where(valid); p=p_two[li,gj]; d=(m1-m2)[li,gj]; mt=m1[li,gj]
o=np.argsort(p); m=len(p); q=np.empty(m); r=p[o]*m/(np.arange(m)+1); q[o]=np.minimum.accumulate(r[::-1])[::-1]; q=np.minimum(q,1)
A=pd.DataFrame({"Gene":np.array(genes)[gj],"Lineage":np.array(elig)[li],"mean_target":mt,"mean_diff":d,"p_two_sided":p,"q":q,"n_target":n1[li,gj]})
A.to_csv(RD+"/ttest_alt_all_tests.csv",index=False)
F=pd.read_csv(RD+"/final_targets.csv"); k=F.merge(A,on=["Gene","Lineage"],how="left")
R={"n_tests":int(m),"q_lt_0.05_negative_effect":int(((A.q<0.05)&(A.mean_diff<0)).sum()),
   "meyers_style_calls(q<0.05,mean_diff<0)":int(((A.q<0.05)&(A.mean_diff<0)).sum()),
   "mean_based_candidates(q<0.05,mean_target<-1,mean_diff<-0.5)":int(((A.q<0.05)&(A.mean_target<-1)&(A.mean_diff<-0.5)).sum()),
   "our94_with_ttest_q<0.05_negative":int(((k.q<0.05)&(k.mean_diff<0)).sum()),"our94_meeting_mean_based_rule":int(((k.q<0.05)&(k.mean_target<-1)&(k.mean_diff<-0.5)).sum()),
   "our94_median_ttest_q":float(k.q.median()),"our94_max_ttest_q":float(k.q.max())}
mb=A[(A.q<0.05)&(A.mean_target<-1)&(A.mean_diff<-0.5)]; ks=set(zip(F.Gene,F.Lineage)); R["mean_based_shared_with_our94"]=len(ks&set(zip(mb.Gene,mb.Lineage)))
R["our94_failing_mean_rule"]=k[~((k.q<0.05)&(k.mean_target<-1)&(k.mean_diff<-0.5))][["Gene","Lineage","n_target_x","q","mean_target","mean_diff"]].round(3).values.tolist()
json.dump(R,open(RD+"/ttest_alt_summary.json","w"),indent=1,default=float); print(json.dumps(R,indent=1,default=float))
