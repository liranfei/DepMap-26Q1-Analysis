import os
"""Method 1: second, differently-coded implementation of the downstream analyses, compared with stored results (part2/part3)."""
import json, numpy as np, pandas as pd
from scipy.stats import norm, hypergeom, rankdata
import statsmodels.api as sm
D=os.environ.get("DEPMAP_DIR","."); O=os.environ.get("RESULTS_DIR","results")+"/"; P2=json.load(open(O+"part2_summary.json")); P3=json.load(open(O+"part3_summary.json")); out=[]
def rec(name,mine,other,tol=1e-6,rel=True):
    ok=abs(mine-other)<=tol*(abs(mine) if rel and mine else 1); out.append((name,mine,other,"PASS" if ok else "FAIL")); return ok
# ---- own MWU (one-sided less) written from scratch
def mwu_less(x,y):
    x=np.asarray(x,float); y=np.asarray(y,float); n1,n2=len(x),len(y); a=np.concatenate([x,y]); r=rankdata(a); R1=r[:n1].sum(); U1=R1-n1*(n1+1)/2; U2=n1*n2-U1
    _,c=np.unique(a,return_counts=True); N=n1+n2; s=np.sqrt(n1*n2/12*((N+1)-((c**3-c).sum())/(N*(N-1)))); return float(norm.sf((U2-n1*n2/2-0.5)/s))
# ---- data by different route (numpy arrays, dict mapping)
ge=pd.read_csv(D+"/gene_effect.csv",index_col=0); mo=pd.read_csv(D+"/Model.csv").set_index("ModelID")
keep=[i for i in ge.index if mo.loc[i,"OncotreePrimaryDisease"]!="Non-Cancerous"]; ge=ge.loc[keep]; lin=mo.loc[keep,"OncotreeLineage"].str.strip().str.lower().values; ids=np.array(keep)
def gcol(sym): return ge[[c for c in ge.columns if c.startswith(sym+" (")][0]].values.astype(float)
# ---- DepMap simple comparisons
for sym,l,key in [("KRAS","pancreas","KRAS_pancreas"),("HNF1B","kidney","HNF1B_kidney")]:
    y=gcol(sym); ok=~np.isnan(y); t=y[ok&(lin==l)]; o=y[ok&(lin!=l)]; ref=P2[key]
    rec(f"{key} n_target",len(t),ref["n_pancreas" if sym=="KRAS" else "n_kidney"],0); rec(f"{key} median_target",np.median(t),ref["med_pancreas" if sym=="KRAS" else "med_kidney"],1e-9)
    p=mwu_less(t,o); rec(f"{key} p(log10)",np.log10(p),np.log10(ref["p"]),1e-4)
# ---- mutations from raw file, different route
raw=pd.read_csv(D+"/OmicsSomaticMutations.csv",usecols=["ModelID","HugoSymbol","VepImpact","Hotspot","LikelyLoF","IsDefaultEntryForModel"],low_memory=False)
raw=raw[(raw.IsDefaultEntryForModel=="Yes")&raw.HugoSymbol.isin(["KRAS","CTNNB1","TP53"])]
def frisch(y,x,g):   # OLS coefficient of x controlling for group dummies, via demeaning (Frisch-Waugh)
    d=pd.DataFrame({"y":y,"x":x,"g":g}).dropna(); yy=d.y-d.groupby("g").y.transform("mean"); xx=d.x-d.groupby("g").x.transform("mean"); return float((xx*yy).sum()/(xx*xx).sum())
for sym in ["KRAS","CTNNB1"]:
    m=raw[raw.HugoSymbol==sym]; hot=set(m.loc[m.Hotspot==True,"ModelID"]); alt=set(m.loc[m.VepImpact.isin(["HIGH","MODERATE"]),"ModelID"])
    isH=np.array([i in hot for i in ids]); isA=np.array([(i in alt) and (i not in hot) for i in ids]); y=gcol(sym); ok=~np.isnan(y); ref=P2[f"{sym}_genotype"]
    rec(f"{sym} hotspot n",int((isH&ok).sum()),ref["n"]["hotspot"],0); rec(f"{sym} none n",int((~isH&~isA&ok).sum()),ref["n"]["none"],0); rec(f"{sym} other n",int((isA&ok).sum()),ref["n"]["other_protein_altering"],0)
    rec(f"{sym} hotspot median",np.median(y[isH&ok]),ref["median"]["hotspot"],1e-3); rec(f"{sym} none median",np.median(y[~isH&~isA&ok]),ref["median"]["none"],1e-3)
    p=mwu_less(y[isH&ok],y[~isH&~isA&ok]); rec(f"{sym} hotspot-vs-none p(log10)",np.log10(p),np.log10(ref["hotspot_vs_none_p"]),1e-4)
    sel=ok&~isA; coef=frisch(y[sel],isH[sel].astype(float),lin[sel]); rec(f"{sym} lineage-adjusted coef",coef,ref["hotspot_adj_lineage_coef"],1e-6)
t=raw[raw.HugoSymbol=="TP53"]; pa=set(t.loc[t.VepImpact.isin(["HIGH","MODERATE"]),"ModelID"]); mt=np.array([i in pa for i in ids]); y=gcol("MDM2"); ok=~np.isnan(y); ref=P2["MDM2_TP53_protein_altering"]
rec("TP53 wt n",int((~mt&ok).sum()),ref["n_wt"],0); rec("TP53 mut n",int((mt&ok).sum()),ref["n_mut"],0); rec("MDM2 median wt",np.median(y[~mt&ok]),ref["med_wt"],1e-9); rec("MDM2 adj coef",frisch(y[ok],mt[ok].astype(float),lin[ok]),ref["adj_lineage_coef"],1e-6)
rec("MDM2/TP53 p(log10)",np.log10(mwu_less(y[~mt&ok],y[mt&ok])),np.log10(ref["p_one_sided"]),1e-4)
# ---- batch/library adjustment for a few pairs via least squares on dummies
ssm=pd.read_csv(D+"/ScreenSequenceMap.csv"); ssm=ssm[(ssm.ScreenType=="2DS")&(ssm.PassesQC==True)&(ssm.DrugTreated==False)&(ssm.IsEngineered==False)&(ssm.ExcludeFromCRISPRCombined==False)]
bm=ssm.groupby("ModelID").pDNABatch.agg(lambda s:s.value_counts().index[0]); lm=ssm.groupby("ModelID").Library.agg(lambda s:s.value_counts().index[0])
bat=pd.Series([bm.get(i,"unknown") for i in ids]); lib=pd.Series([lm.get(i,"unknown") for i in ids]); Bt=pd.read_csv(O+"batch_library_adjusted_pairs.csv")
for r in Bt.sort_values("coef_unadj").head(5).itertuples():
    y=gcol(r.Gene.split(" (")[0]); ok=~np.isnan(y); tg=(lin==r.Lineage).astype(float)
    M=np.column_stack([tg[ok],pd.get_dummies(bat[ok],drop_first=True).values.astype(float),pd.get_dummies(lib[ok],drop_first=True).values.astype(float),np.ones(ok.sum())])
    beta=np.linalg.lstsq(M,y[ok],rcond=None)[0][0]; rec(f"batch/library adj coef {r.Gene} {r.Lineage}",float(beta),float(r.coef_adj),1e-6)
# ---- TCGA: hand-written log-rank + statsmodels PHReg vs stored (lifelines)
def tcga(expr,surv,ensg,label):
    # independent route: stream the file, accumulate library sizes (sum of 2^x - 1 over all genes) and normalise to log2(CPM + 1)
    with open(expr) as f:
        hdr=next(f).rstrip("\n").split("\t"); lib=np.zeros(len(hdr)-1); vals=None
        for line in f:
            parts=line.rstrip("\n").split("\t"); cnt=np.exp2(np.array(parts[1:],float))-1; lib+=cnt
            if ensg in parts[0]: vals=cnt
    e=np.log2(vals/lib*1e6+1)
    names=np.array(hdr[1:]); code=np.array([n.split("-")[3][:2] for n in names]); tum=code=="01"
    d=pd.DataFrame({"pid":[n[:12] for n in names[tum]],"x":e[tum]}).groupby("pid",as_index=False).x.mean()
    s=pd.read_csv(surv,sep="\t"); s=s[s["sample"].str.split("-").str[3].str[:2]=="01"].drop_duplicates(subset="sample"); s["pid"]=s["sample"].str[:12]; s=s.drop_duplicates("pid")
    m=d.merge(s[["pid","OS.time","OS"]],on="pid").dropna(); m=m[m["OS.time"]>0]; m["hi"]=(m.x>=m.x.median()).astype(int)
    # log-rank from scratch
    tt=m["OS.time"].values; ev=m.OS.values.astype(int); g=m.hi.values; O1=E1=V=0.0
    for u in np.unique(tt[ev==1]):
        atr=tt>=u; n=atr.sum(); n1=(atr&(g==1)).sum(); dd=((tt==u)&(ev==1)).sum(); d1=((tt==u)&(ev==1)&(g==1)).sum(); O1+=d1; E1+=dd*n1/n
        if n>1: V+=dd*(n1/n)*(1-n1/n)*(n-dd)/(n-1)
    from scipy.stats import chi2
    plr=float(chi2.sf((O1-E1)**2/V,1)); ref=P3[label]
    rec(f"{label} n_survival",len(m),ref["n_survival"],0); rec(f"{label} logrank p (hand-coded)",plr,ref["logrank_p_mediansplit"],1e-6)
    z=((m.x-m.x.mean())/m.x.std()).values; res=sm.PHReg(tt,z[:,None],status=ev,ties="efron").fit(); rec(f"{label} Cox HR/SD (statsmodels)",float(np.exp(res.params[0])),ref["cox_continuous_HR_per_SD"],1e-4); rec(f"{label} Cox p (statsmodels)",float(res.pvalues[0]),ref["cox_continuous_p"],1e-3)
tcga(os.environ.get("TCGA_DIR",".")+"/TCGA-KIRC.star_counts.tsv",os.environ.get("TCGA_DIR",".")+"/TCGA-KIRC.survival.tsv","ENSG00000275410","TCGA_KIRC_HNF1B")
tcga(os.environ.get("TCGA_DIR",".")+"/TCGA-PAAD.star_counts.tsv",os.environ.get("TCGA_DIR",".")+"/TCGA-PAAD.survival.tsv","ENSG00000133703","TCGA_PAAD_KRAS")
# ---- enrichment: hypergeometric survival function vs stored Fisher p
E=pd.read_csv(O+"enrichment_custom_background.csv"); T=pd.read_csv(O+"all_tests.csv"); uni=set(T.Gene.str.split(" (",regex=False).str[0]); F=pd.read_csv(O+"final_targets.csv"); hits=set(F.Gene.str.split(" (",regex=False).str[0])
for r in E.sort_values("q").head(6).itertuples():
    k=r.n_hits; K=r.set_size; p=float(hypergeom.sf(k-1,len(uni),K,len(hits))); rec(f"enrichment {r.term[:28]} p(log10)",np.log10(p),np.log10(r.p),1e-6)
# ---- common essential overlap and named co-dependencies (numpy masked correlation)
ce={g.split(" (")[0] for g in pd.read_csv(D+"/CRISPRInferredCommonEssentials.csv").iloc[:,0]}; rec("common essential overlap (symbols)",len(hits&ce),P2["common_essential_annotation"]["n_overlap"],0)
def pcorr(a,b):
    x,y=gcol(a),gcol(b); ok=~np.isnan(x)&~np.isnan(y); return float(np.corrcoef(x[ok],y[ok])[0,1])
for a,b,key in [("ADSL","PAICS","ADSL_PAICS"),("NMNAT1","NAMPT","NMNAT1_NAMPT"),("UMPS","CTPS1","UMPS_CTPS1")]: rec(f"co-dependency {key}",pcorr(a,b),P2["codependency"][key],1e-6)
df_=pd.DataFrame(out,columns=["check","this_route","stored","status"]); df_.to_csv(O+"cross_check_results.csv",index=False)
print(df_.to_string(index=False)); print("\nPASS:",(df_.status=="PASS").sum(),"FAIL:",(df_.status=="FAIL").sum())
