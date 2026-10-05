"""Post hoc quantification added AFTER the pre-specified checks were run (labelled as such)."""
import os, json, numpy as np, pandas as pd
from scipy.stats import hypergeom, spearmanr
DATA=os.environ.get("DEPMAP_DIR","."); REPO=os.environ.get("RESULTS_DIR","results")
R=json.load(open(REPO+"/snap23_validation_summary.json")); P={}

# 1. V6: how unlikely is it that 2 of the pre-specified SNARE partners rank in the top 20?
r=pd.read_csv(REPO+"/snap23_codependency_all_genes.csv",index_col=0).iloc[:,0]
PART={"STX1A","STX2","STX3","STX4","STX7","STX8","VAMP2","VAMP3","VAMP4","VAMP7","VAMP8",
      "SNAP25","SNAP29","NAPA","NAPB","NSF","STXBP1","STXBP2","STXBP3","SEC22B"}
syms=pd.Index([g.split(" ")[0] for g in r.index])
N=len(r); K=int(syms.isin(PART).sum()); n=20
k=int(pd.Index([g.split(" ")[0] for g in r.sort_values(ascending=False).head(n).index]).isin(PART).sum())
P["V6_enrichment_posthoc"]={"n_genes":N,"n_prespecified_partners_present":K,"top_n":n,
  "partners_in_top_n":k,"hypergeometric_P_ge_k":float(hypergeom.sf(k-1,N,K,n)),
  "note":"STXBP3 and STX4 rank 1 and 2 of %d genes" % N}

# 2. V5 within the bowel lineage only
hdr=pd.read_csv(DATA+"/gene_effect.csv",nrows=0).columns
need=[hdr[0]]+[c for c in hdr if c.split(" ")[0]=="SNAP23"]
ge=pd.read_csv(DATA+"/gene_effect.csv",usecols=need).rename(columns={hdr[0]:"DepMap_ID"})
M=pd.read_csv(DATA+"/Model.csv")[["ModelID","OncotreeLineage","OncotreePrimaryDisease"]]
M.columns=["DepMap_ID","lineage","pd_"]; M["lineage"]=M.lineage.str.strip().str.lower()
d=ge.merge(M,on="DepMap_ID").query('pd_ != "Non-Cancerous"')
d["y"]=pd.to_numeric(d[need[1]],errors="coerce")
eh=pd.read_csv(DATA+"/OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv",nrows=0).columns
ec=["ModelID","IsDefaultEntryForModel"]+[c for c in eh if c.split(" ")[0] in ("SNAP25","SNAP29")]
ex=pd.read_csv(DATA+"/OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv",usecols=ec)
ex=ex[ex.IsDefaultEntryForModel=="Yes"].rename(columns={"ModelID":"DepMap_ID"})
ex.columns=[c.split(" ")[0] for c in ex.columns]
e=d[["DepMap_ID","y","lineage"]].merge(ex,on="DepMap_ID")
out={}
for scope in ("bowel","non-bowel"):
    sub=e[e.lineage=="bowel"] if scope=="bowel" else e[e.lineage!="bowel"]
    for p_ in ("SNAP25","SNAP29"):
        ok=sub[[p_,"y"]].dropna(); rho,pv=spearmanr(ok[p_],ok.y)
        out[f"{scope}_{p_}"]={"n":int(len(ok)),"spearman_rho":float(rho),"p_two_sided":float(pv)}
P["V5_within_lineage_posthoc"]=out

# 3. V7 with the manuscript's criterion (BH q over the whole test family)
for name,f in [("sanger_project_score","all_tests_sanger.csv.gz"),("rnai_demeter2","rnai_all_tests.csv.gz")]:
    t=pd.read_csv(REPO+"/"+f)
    row=t[(t.Gene=="SNAP23 (8773)")&(t.Lineage=="bowel")].iloc[0]
    P.setdefault("V7_with_q_criterion",{})[name]={"p_value":float(row.p_value),"q_value":float(row.q_value),
        "selectivity":float(row.Selectivity),"n_target":int(row.n_target),
        "supported_at_q<0.05":bool(row.q_value<0.05 and row.Selectivity<0)}
# 4. are the SNARE partners themselves more required in bowel?
from scipy.stats import mannwhitneyu
import numpy as np
h2=pd.read_csv(os.path.join(DATA,"gene_effect.csv"),nrows=0).columns
nd=[h2[0]]+[c for c in h2 if c.split(" ")[0] in ("STX4","STXBP3","SNAP23")]
g2=pd.read_csv(os.path.join(DATA,"gene_effect.csv"),usecols=nd).rename(columns={h2[0]:"DepMap_ID"})
m2=pd.read_csv(os.path.join(DATA,"Model.csv"))[["ModelID","OncotreeLineage","OncotreePrimaryDisease"]]
m2.columns=["DepMap_ID","lineage","pd_"]; m2["lineage"]=m2.lineage.str.strip().str.lower()
d2=g2.merge(m2,on="DepMap_ID").query('pd_ != "Non-Cancerous"')
rep={}
for c in nd[1:]:
    y=pd.to_numeric(d2[c],errors="coerce"); ok=d2.assign(y=y).dropna(subset=["y"])
    a=ok[ok.lineage=="bowel"].y.values; b=ok[ok.lineage!="bowel"].y.values
    rep[c.split(" ")[0]]={"n_bowel":int(len(a)),"median_bowel":round(float(np.median(a)),3),
      "median_other":round(float(np.median(b)),3),"selectivity":round(float(np.median(a)-np.median(b)),3),
      "p_one_sided":float(mannwhitneyu(a,b,alternative="less",method="asymptotic")[1])}
P["module_partners_bowel_selectivity_posthoc"]=rep

R["post_hoc_additions"]=P
json.dump(R,open(REPO+"/snap23_validation_summary.json","w"),indent=1,ensure_ascii=False)
print(json.dumps(P,indent=1,ensure_ascii=False))
