import os
import json, numpy as np, pandas as pd, statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, spearmanr, percentileofscore
from scipy.cluster.hierarchy import linkage, leaves_list
from statsmodels.stats.multitest import multipletests
import run_pipeline as rp
D=os.environ.get("DEPMAP_DIR","."); O=os.environ.get("RESULTS_DIR","results")+"/"; R={}
mw=lambda x,y,alt="less": float(mannwhitneyu(x,y,alternative=alt,method="asymptotic")[1])
df=rp.load(D); df=df[df.primary_disease!="Non-Cancerous"].reset_index(drop=True)
genes=[c for c in df.columns if c not in ("DepMap_ID","lineage","primary_disease")]
F=pd.read_csv(O+"final_targets.csv"); T=pd.read_csv(O+"all_tests.csv"); tg=sorted(F.Gene.unique())
col=lambda s:[c for c in genes if c.startswith(s+" (")][0]
meta=pd.read_csv(D+"/Model.csv")[["ModelID","OncotreeSubtype"]].rename(columns={"ModelID":"DepMap_ID"}); df=df.merge(meta,on="DepMap_ID")
# ---- batch / library
s=pd.read_csv(D+"/ScreenSequenceMap.csv"); s=s[(s.ScreenType=="2DS")&(s.PassesQC==True)&(s.DrugTreated==False)&(s.IsEngineered==False)&(s.ExcludeFromCRISPRCombined==False)]
mode=lambda x:x.value_counts().index[0]
df=df.merge(s.groupby("ModelID").pDNABatch.agg(mode).rename("batch"),left_on="DepMap_ID",right_index=True,how="left").merge(s.groupby("ModelID").Library.agg(mode).rename("library"),left_on="DepMap_ID",right_index=True,how="left")
R["batch_library_coverage"]={"n":len(df),"with_batch":int(df.batch.notna().sum()),"library_counts":df.library.fillna("unknown").value_counts().to_dict()}
df["batch"]=df.batch.fillna("unknown"); df["library"]=df.library.fillna("unknown")
rows=[]
for r in F.itertuples():
    d=pd.DataFrame({"y":pd.to_numeric(df[r.Gene],errors="coerce"),"tgt":(df.lineage==r.Lineage).astype(int),"batch":df.batch,"library":df.library}).dropna()
    m0=smf.ols("y~tgt",d).fit(); m1=smf.ols("y~tgt+C(batch)+C(library)",d).fit(); mb=smf.ols("y~C(batch)+C(library)",d).fit()
    rows.append(dict(Gene=r.Gene,Lineage=r.Lineage,coef_unadj=m0.params["tgt"],coef_adj=m1.params["tgt"],p_adj=m1.pvalues["tgt"],batch_library_R2=mb.rsquared))
B=pd.DataFrame(rows); B["q_adj"]=multipletests(B.p_adj,method="fdr_bh")[1]; B["ratio"]=B.coef_adj/B.coef_unadj; B["retained"]=(B.q_adj<0.05)&(B.coef_adj<-0.3)
B.to_csv(O+"batch_library_adjusted_pairs.csv",index=False)
R["batch_library"]={"n_pairs":len(B),"median_R2":float(B.batch_library_R2.median()),"max_R2":float(B.batch_library_R2.max()),"n_q05_after_adj":int((B.q_adj<0.05).sum()),
  "n_retained_coef_lt_-0.3":int(B.retained.sum()),"median_ratio":float(B.ratio.median()),"min_ratio":float(B.ratio.min()),"weakened":B[~B.retained][["Gene","Lineage","coef_unadj","coef_adj","q_adj"]].round(3).values.tolist()}
# ---- genotype
mut=pd.read_csv(os.environ.get("RESULTS_DIR","results")+"/mut_3genes.csv"); mut=mut[mut.IsDefaultEntryForModel=="Yes"]
def groups(sym):
    m=mut[mut.HugoSymbol==sym]; hot=set(m[m.Hotspot==True].ModelID); pa=set(m[m.VepImpact.isin(["HIGH","MODERATE"])].ModelID)
    g=np.where(df.DepMap_ID.isin(hot),"hotspot",np.where(df.DepMap_ID.isin(pa),"other_protein_altering","none")); return g
def geno(sym):
    y=pd.to_numeric(df[col(sym)],errors="coerce"); g=groups(sym); d=pd.DataFrame({"y":y,"g":g,"lin":df.lineage}).dropna()
    out={"n":d.g.value_counts().to_dict(),"median":d.groupby("g").y.median().round(3).to_dict()}
    h=d[d.g=="hotspot"].y; n0=d[d.g=="none"].y; pa=d[d.g!="none"].y
    out["hotspot_vs_none_p"]=mw(h,n0); out["any_protein_altering_vs_none_p"]=mw(pa,n0)
    d2=d[d.g!="other_protein_altering"].copy(); d2["mut"]=(d2.g=="hotspot").astype(int); m1=smf.ols("y~mut+C(lin)",d2).fit(); out["hotspot_adj_lineage_coef"]=float(m1.params["mut"]); out["hotspot_adj_lineage_p_two_sided"]=float(m1.pvalues["mut"])
    d3=d.copy(); d3["mut"]=(d3.g!="none").astype(int); m2=smf.ols("y~mut+C(lin)",d3).fit(); out["any_adj_lineage_coef"]=float(m2.params["mut"]); out["any_adj_lineage_p_two_sided"]=float(m2.pvalues["mut"])
    return out,d
for sym,lin in [("KRAS","pancreas"),("CTNNB1","bowel")]:
    o,d=geno(sym); w=d[d.lin==lin]; o["within_lineage"]={lin:{"n":w.g.value_counts().to_dict(),"median":w.groupby("g").y.median().round(3).to_dict()}}
    nonl=d[d.lin!=lin]; o["excluding_lineage_hotspot_vs_none_p"]=mw(nonl[nonl.g=="hotspot"].y,nonl[nonl.g=="none"].y); R[f"{sym}_genotype"]=o
k=df[["lineage",col("KRAS")]].dropna(); kp=k[k.lineage=="pancreas"].iloc[:,1]; ko=k[k.lineage!="pancreas"].iloc[:,1]
R["KRAS_pancreas"]={"n_pancreas":len(kp),"n_other":len(ko),"med_pancreas":float(kp.median()),"med_other":float(ko.median()),"p":mw(kp,ko)}
# TP53 / MDM2
t=mut[mut.HugoSymbol=="TP53"]; pa=set(t[t.VepImpact.isin(["HIGH","MODERATE"])].ModelID); lof=set(t[(t.LikelyLoF==True)|(t.Hotspot==True)].ModelID)
y=pd.to_numeric(df[col("MDM2")],errors="coerce")
for name,S in [("protein_altering",pa),("hotspot_or_LoF",lof)]:
    mm=df.DepMap_ID.isin(S); d=pd.DataFrame({"y":y,"mut":mm.astype(int),"lin":df.lineage}).dropna(); a=d[d.mut==0].y; b=d[d.mut==1].y
    m1=smf.ols("y~mut+C(lin)",d).fit(); R[f"MDM2_TP53_{name}"]={"n_wt":len(a),"n_mut":len(b),"med_wt":float(a.median()),"med_mut":float(b.median()),"p_one_sided":mw(a,b),"adj_lineage_coef":float(m1.params["mut"]),"adj_lineage_p":float(m1.pvalues["mut"])}
# HNF1B kidney
h=df[["lineage",col("HNF1B")]].dropna(); hk=h[h.lineage=="kidney"].iloc[:,1]; ho=h[h.lineage!="kidney"].iloc[:,1]
R["HNF1B_kidney"]={"n_kidney":len(hk),"n_other":len(ho),"med_kidney":float(hk.median()),"med_other":float(ho.median()),"p":mw(hk,ho)}
# ---- common essentials annotation (official list, annotation only)
ce=set(pd.read_csv(D+"/CRISPRInferredCommonEssentials.csv").iloc[:,0]); tested=T.Gene.unique()
ov=[g for g in tg if g in ce]; R["common_essential_annotation"]={"n_list":len(ce),"n_target_genes":len(tg),"n_overlap":len(ov),"frac":len(ov)/len(tg),"pairs_with_CE_gene":int(F.Gene.isin(ce).sum()),"n_pairs":len(F),
  "tested_genes_frac_CE":float(np.mean([g in ce for g in tested])),"overlap_genes":ov}
# ---- co-dependency
X=df[tg].apply(pd.to_numeric,errors="coerce"); C=X.corr(); C.to_csv(O+"codependency_corr.csv"); iu=np.triu_indices(len(tg),1)
pairs=sorted([(C.index[i],C.columns[j],float(C.values[i,j])) for i,j in zip(*iu)],key=lambda z:-z[2])
rng=np.random.default_rng(42); rg=rng.choice(genes,400,replace=False); RC=df[list(rg)].apply(pd.to_numeric,errors="coerce").corr().values[np.triu_indices(400,1)]
gc=lambda a_,b_:(float(C.loc[[g for g in tg if g.startswith(a_+" (")][0],[g for g in tg if g.startswith(b_+" (")][0]]) if any(g.startswith(a_+" (") for g in tg) and any(g.startswith(b_+" (") for g in tg) else None)
R["codependency"]={"n_genes":len(tg),"n_pairs_r_gt_0.5":int((C.values[iu]>0.5).sum()),"top10":[(a_,b_,round(c,3)) for a_,b_,c in pairs[:10]],"random_pair_mean_r":float(np.nanmean(RC)),"random_pair_sd_r":float(np.nanstd(RC)),
  "NMNAT1_NAMPT":gc("NMNAT1","NAMPT"),"ADSL_PAICS":gc("ADSL","PAICS"),"UMPS_CTPS1":gc("UMPS","CTPS1"),"SOX10_HNF1B":gc("SOX10","HNF1B")}
# ---- skewness concordance (descriptive; same data, not independent validation)
sk=df[genes].skew(); best=F.groupby("Gene").Selectivity.min(); rho,p=spearmanr(best[tg],sk[tg]); bg=sk.drop(tg).dropna()
R["skewness_descriptive"]={"median_targets":float(sk[tg].median()),"median_background":float(bg.median()),"p":mw(sk[tg].dropna(),bg),"pct_rank_target_median":float(percentileofscore(sk.dropna(),sk[tg].median())),"spearman_rho":float(rho),"spearman_p":float(p)}
# ---- hematologic background sensitivity (descriptive, candidates only)
rows=[]
for r in F[F.Lineage.isin(["lymphoid","myeloid"])].itertuples():
    y=pd.to_numeric(df[r.Gene],errors="coerce"); ok=y.notna()
    for name,mask in [("exclude_other_heme",~df.lineage.isin(["lymphoid","myeloid"])),("all_others",df.lineage!=r.Lineage)]:
        tv=y[ok&(df.lineage==r.Lineage)]; ov=y[ok&mask&(df.lineage!=r.Lineage)] if name=="all_others" else y[ok&mask]
        rows.append(dict(Gene=r.Gene,Lineage=r.Lineage,background=name,n_target=len(tv),n_other=len(ov),selectivity=float(tv.median()-ov.median()),p=mw(tv,ov)))
H=pd.DataFrame(rows); H.to_csv(O+"heme_background_sensitivity.csv",index=False)
w=H.pivot_table(index=["Gene","Lineage"],columns="background",values="selectivity"); R["heme_background"]={"n_pairs":len(w),"median_selectivity_all_others":float(w.all_others.median()),"median_selectivity_excl_other_heme":float(w.exclude_other_heme.median()),"n_still_lt_-0.5":int((w.exclude_other_heme<-0.5).sum())}
# ---- subsampling power (pairs with n_target>=30)
rows=[]; rng=np.random.default_rng(42)
for r in F[F.n_target>=30].itertuples():
    x=pd.to_numeric(df[r.Gene],errors="coerce"); tv=x[df.lineage==r.Lineage].dropna().values; ov=x[df.lineage!=r.Lineage].dropna().values
    for kk in [10,20]:
        hit=0
        for _ in range(300):
            sm=rng.choice(tv,kk,replace=False); hit+=(np.median(sm)<-1 and np.median(sm)-np.median(ov)<-0.5 and mw(sm,ov)<0.05)
        rows.append(dict(Gene=r.Gene,Lineage=r.Lineage,n_full=len(tv),k=kk,frac_retained=hit/300))
S=pd.DataFrame(rows); S.to_csv(O+"subsampling_power.csv",index=False); R["subsampling"]={"n_pairs":int(S.Gene.size/2),"by_k":S.groupby("k").frac_retained.agg(["mean","median","min"]).round(3).to_dict("index")}
# ---- subtype sensitivity in lymphoid / myeloid
rows=[]
for r in F[F.Lineage.isin(["lymphoid","myeloid"])].itertuples():
    x=pd.to_numeric(df[r.Gene],errors="coerce"); oth=x[df.lineage!=r.Lineage].dropna()
    for st,idx in df[df.lineage==r.Lineage].groupby("OncotreeSubtype").groups.items():
        v=x.loc[idx].dropna()
        if len(v)>=5: rows.append(dict(Gene=r.Gene,Lineage=r.Lineage,Subtype=st,n=len(v),median=float(v.median()),selectivity=float(v.median()-oth.median()),p=mw(v,oth)))
ST=pd.DataFrame(rows); ST.to_csv(O+"subtype_sensitivity.csv",index=False)
R["subtype"]={"n_tests":len(ST),"n_pairs":int(ST[["Gene","Lineage"]].drop_duplicates().shape[0]),"frac_sel_lt_-0.5":float((ST.selectivity<-0.5).mean()),"frac_p_lt_0.05":float((ST.p<0.05).mean()),"n_subtypes":int(ST.Subtype.nunique())}
json.dump(R,open(O+"part2_summary.json","w"),indent=1,default=float)
for k,v in R.items(): print("\n##",k); print(json.dumps(v,default=float)[:700])
