from figlib import *
import textwrap
from lifelines import KaplanMeierFitter
from matplotlib.patches import Patch
import run_pipeline as rp
F,T,lc=load_all(); SS=json.load(open(O+"sensitivity_summary.json")); TT=json.load(open(O+"ttest_alt_summary.json")); D=os.environ.get("DEPMAP_DIR","."); TD=os.environ.get("TCGA_DIR",".")
def letters(axs,dx=-0.02,dy=1.02):
    for a,l in zip(axs,"ABCDEFGH"): a.text(dx,dy,l,transform=a.transAxes,fontsize=12,fontweight="bold",va="bottom",ha="right")
def bars(ax,labels,vals,colors,ylabel,fmt="{:,}"):
    x=np.arange(len(vals)); ax.bar(x,vals,color=colors,width=0.7); ax.set_xticks(x); ax.set_xticklabels(labels); ax.set_ylabel(ylabel)
    for i,v in enumerate(vals): ax.text(i,v+max(vals)*0.02,fmt.format(v),ha="center",va="bottom",fontsize=8)
    ax.set_ylim(0,max(vals)*1.15)
# ================= Fig 6 sensitivity
fig,axs=plt.subplots(2,2,figsize=(7.3,6.0)); a,b,c,d=axs.ravel()
m=[SS["multiplicity_BH"]["n_final"],SS["multiplicity_BY"]["n_final"],SS["multiplicity_Holm"]["n_final"]]; bars(a,["BH\n(primary)","BY","Holm"],m,[BLUE,SKY,GREY],"Candidate pairs")
sel=[SS["selectivity<-0.3"]["n_final"],SS["selectivity<-0.5"]["n_final"],SS["selectivity<-0.7"]["n_final"],SS["selectivity<-1.0"]["n_final"]]
bars(b,["−0.3","−0.5\n(primary)","−0.7","−1.0"],sel,[SKY,BLUE,SKY,SKY],"Candidate pairs"); b.set_xlabel("Selectivity threshold")
tn=pd.DataFrame(SS["topN"]); x=np.arange(len(tn)); c.plot(x,tn.n_final,"o-",color=BLUE,lw=1.5,ms=5); c.plot([len(tn)-1],[tn.n_final.iloc[-1]],"o",color=VERM,ms=7)
c.set_xticks(x); c.set_xticklabels([f"{v:,}" if v!="all" else "All\n(primary)" for v in tn.top_n]); c.set_xlabel("Genes tested (top-N by variance)"); c.set_ylabel("Candidate pairs"); c.set_ylim(70,110)
for i,v in enumerate(tn.n_final): c.text(i,v+1.8,str(v),ha="center",fontsize=8)
n94=94; both=TT["mean_based_shared_with_our94"]; mean_only=TT["mean_based_candidates(q<0.05,mean_target<-1,mean_diff<-0.5)"]-both
d.barh([2,1,0],[n94-both,both,mean_only],left=[0,n94-both,n94],color=[BLUE,GREEN,ORANGE],height=0.5)
d.barh([0],[both],left=[0],color="none"); d.set_yticks([]); d.set_xlabel("Candidate pairs")
d.text(0,2.45,f"Median rule only: {n94-both}",fontsize=8,va="bottom"); d.text(0,1.45,f"Both rules: {both}",fontsize=8,va="bottom"); d.text(0,0.45,"",fontsize=8)
d.clear(); d.axis("off")
# panel D: stacked single bar with legend
d.axis("on"); d.barh(0,n94-both,color=BLUE,height=0.5); d.barh(0,both,left=n94-both,color=GREEN,height=0.5); d.barh(0,mean_only,left=n94,color=ORANGE,height=0.5)
d.set_yticks([]); d.set_xlim(0,n94+mean_only+8); d.set_ylim(-0.6,0.9); d.set_xlabel("Candidate pairs")
d.legend(handles=[Patch(color=BLUE,label=f"Median rule only (n = {n94-both})"),Patch(color=GREEN,label=f"Both rules (n = {both})"),Patch(color=ORANGE,label=f"Mean rule only (n = {mean_only})")],frameon=False,loc="upper center",bbox_to_anchor=(0.5,1.0),ncol=1)
d.spines["left"].set_visible(False)
letters(axs.ravel(),-0.16); fig.tight_layout(pad=0.6,h_pad=1.6); save(fig,"Fig6")
# ================= Fig 7 enrichment (formal over-representation only)
E=pd.read_csv(O+"enrichment_custom_background.csv"); fig,axs=plt.subplots(1,2,figsize=(7.3,5.4))
for ax,lib,title in zip(axs,["GO_Biological_Process_2023","KEGG_2021_Human"],["GO Biological Process","KEGG pathways"]):
    e=E[(E.library==lib)&(E.q<0.05)].sort_values("q").head(8).iloc[::-1]; y=np.arange(len(e))
    ax.scatter(-np.log10(e.q),y,s=e.n_hits*22,color=BLUE if lib.startswith("GO") else VERM,edgecolor="white",linewidth=0.5,zorder=3); ax.hlines(y,0,-np.log10(e.q),color="#cccccc",lw=1,zorder=1)
    ax.set_yticks(y); ax.set_yticklabels(["\n".join(textwrap.wrap(t.split(" (GO")[0],26)) for t in e.term]); ax.set_xlabel("−log$_{10}$ (BH-adjusted q)"); ax.set_xlim(0,max(-np.log10(e.q))*1.55); ax.margins(y=0.08)
    ax.set_title(title,fontsize=9,loc="center"); [ax.annotate(f"{h}/{s_}",(-np.log10(q),i),xytext=(np.sqrt(h*22)/2+7,0),textcoords="offset points",va="center",ha="left",fontsize=8) for i,(q,h,s_) in enumerate(zip(e.q,e.n_hits,e.set_size))]
letters(axs,-0.02); fig.tight_layout(pad=0.6,w_pad=1.2); save(fig,"Fig7")
# ================= data for Fig 8-9
df=rp.load(D); df=df[df.primary_disease!="Non-Cancerous"].reset_index(drop=True); genes=[c for c in df.columns if c not in ("DepMap_ID","lineage","primary_disease")]
col=lambda s:[c for c in genes if c.startswith(s+" (")][0]; mut=pd.read_csv(O+"mut_3genes.csv"); mut=mut[mut.IsDefaultEntryForModel=="Yes"]
def grp(sym):
    m=mut[mut.HugoSymbol==sym]; hot=set(m[m.Hotspot==True].ModelID); pa=set(m[m.VepImpact.isin(["HIGH","MODERATE"])].ModelID)
    return np.where(df.DepMap_ID.isin(hot),"Hotspot",np.where(df.DepMap_ID.isin(pa),"Other\nprotein-altering","None"))
rng=np.random.default_rng(1)
def box(ax,data,labels,colors,ylabel=None,pos=None):
    pos=pos or list(range(1,len(data)+1)); bp=ax.boxplot(data,positions=pos,widths=0.55,showfliers=False,patch_artist=True,medianprops=dict(color="black",lw=1.5),boxprops=dict(lw=0.8),whiskerprops=dict(lw=0.8),capprops=dict(lw=0.8))
    for p,cc in zip(bp["boxes"],colors): p.set_facecolor(cc); p.set_alpha(0.35)
    for xx,v,cc in zip(pos,data,colors): ax.scatter(xx+rng.uniform(-0.18,0.18,len(v)),v,s=4,color=cc,alpha=0.5,linewidths=0,rasterized=True,zorder=3)
    ax.set_xticks(pos); ax.set_xticklabels(labels)
    if ylabel: ax.set_ylabel(ylabel)
# ================= Fig 8 genotype
fig,axs=plt.subplots(1,3,figsize=(7.3,3.6)); P2=json.load(open(O+"part2_summary.json"))
for ax,sym in zip(axs[:2],["KRAS","CTNNB1"]):
    y=pd.to_numeric(df[col(sym)],errors="coerce"); g=grp(sym); names=["Hotspot","Other\nprotein-altering","None"]; dat=[y[(g==n)&y.notna()].values for n in names]
    box(ax,dat,[f"{n}\n(n = {len(v)})" for n,v in zip(names,dat)],[VERM,ORANGE,GREY],f"{sym} Chronos" if sym=="KRAS" else None); ax.set_title(f"{sym}\nHotspot vs. none: p = {P2[sym+'_genotype']['hotspot_vs_none_p']:.1e}",fontsize=9); ax.tick_params(axis="x",labelsize=8)
axs[0].set_ylabel("Chronos gene effect")
t53=set(mut[(mut.HugoSymbol=="TP53")&mut.VepImpact.isin(["HIGH","MODERATE"])].ModelID); yy=pd.to_numeric(df[col("MDM2")],errors="coerce"); mm=df.DepMap_ID.isin(t53)
dat=[yy[(~mm)&yy.notna()].values,yy[mm&yy.notna()].values]; box(axs[2],dat,[f"No TP53\nvariant\n(n = {len(dat[0])})",f"TP53\nmutant\n(n = {len(dat[1])})"],[BLUE,VERM]); axs[2].set_title(f"MDM2\nNo variant vs. mutant: p = {P2['MDM2_TP53_protein_altering']['p_one_sided']:.1e}",fontsize=9)
letters(axs,-0.05); fig.tight_layout(pad=0.6,w_pad=0.8); save(fig,"Fig8")
# ================= Fig 9 TCGA + DepMap
P3=json.load(open(O+"part3_summary.json")); HAVE_TCGA=all(os.path.exists(O+f"tcga_{c}_expression_samples.csv") for c in ("KIRC_HNF1B","PAAD_KRAS"))
if not HAVE_TCGA: print("normalised TCGA expression tables not found - Fig 9 skipped")
fig,axs=plt.subplots(2,3,figsize=(7.3,5.6))
def expr_groups(label):
    E=pd.read_csv(O+f"tcga_{label}_expression_samples.csv"); code=E["sample"].str.split("-").str[3].str[:2]; return E.log2_cpm[code=="01"].values,E.log2_cpm[code=="11"].values
def km(ax,tab,label,key):
    kmf=KaplanMeierFitter(); tab=tab.copy(); tab["m"]=tab["OS.time"]/30.44
    for hi,cc,nm_ in [(1,VERM,"High"),(0,BLUE,"Low")]:
        s=tab[tab.hi==hi]; kmf.fit(s.m,s.OS,label=f"{label} {nm_} (n = {len(s)})"); kmf.plot_survival_function(ax=ax,ci_show=True,ci_alpha=0.12,color=cc,lw=1.6,show_censors=False)
    ax.set_xlabel("Time (months)"); ax.set_ylabel("Overall survival"); ax.set_ylim(0,1.02)
    r=P3[key]; txt=f"Log-rank p = {r['logrank_p_mediansplit']:.3f}\nCox HR per SD = {r['cox_continuous_HR_per_SD']:.2f}\n(p = {r['cox_continuous_p']:.3f})"
    if key=="TCGA_KIRC_HNF1B": ax.legend(frameon=False,loc="lower left",fontsize=8); ax.text(0.98,0.97,txt,transform=ax.transAxes,ha="right",va="top",fontsize=8)
    else: ax.legend(frameon=False,loc="upper right",fontsize=8,bbox_to_anchor=(1.0,1.0)); ax.text(0.98,0.74,txt,transform=ax.transAxes,ha="right",va="top",fontsize=8)
for row,(sym,lin,ens,path,tabf,key) in enumerate([] if not HAVE_TCGA else [("HNF1B","kidney","ENSG00000275410",TD+"/TCGA-KIRC.star_counts.tsv","tcga_KIRC_HNF1B_survival_table.csv","TCGA_KIRC_HNF1B"),("KRAS","pancreas","ENSG00000133703",TD+"/TCGA-PAAD.star_counts.tsv","tcga_PAAD_KRAS_survival_table.csv","TCGA_PAAD_KRAS")]):
    y=pd.to_numeric(df[col(sym)],errors="coerce"); a=y[(df.lineage==lin)&y.notna()].values; o=y[(df.lineage!=lin)&y.notna()].values
    box(axs[row,0],[o,a],[f"Other\nlineages\n(n = {len(o):,})",f"{nm(lin)}\n(n = {len(a)})"],[GREY,VERM],f"{sym} Chronos"); axs[row,0].axhline(-1,color="#555555",ls="--",lw=0.8)
    t,nrm=expr_groups(tabf.replace("tcga_","").replace("_survival_table.csv","")); r=P3[key]; box(axs[row,1],[nrm,t],[f"Normal\n(n = {len(nrm)})",f"Tumor\n(n = {len(t)})"],[BLUE,VERM],f"{sym} log$_2$(CPM + 1)")
    if row==0: KP=json.load(open(O+"kirc_paired.json")); axs[row,1].set_title(f"TCGA-KIRC\npaired p = {KP['wilcoxon_p_two_sided']:.3f} (MW p = {r['p_tumor_vs_normal_two_sided']:.1e})",fontsize=8)
    else: axs[row,1].set_title("TCGA-PAAD\n4 normal samples (descriptive)",fontsize=9)
    km(axs[row,2],pd.read_csv(O+tabf),sym,key)
axs[0,0].set_title("DepMap",fontsize=9); axs[0,2].set_title("TCGA-KIRC",fontsize=9); axs[1,2].set_title("TCGA-PAAD",fontsize=9); axs[1,0].set_title("DepMap",fontsize=9)
if HAVE_TCGA: letters(axs.ravel(),-0.12); fig.tight_layout(pad=0.6,h_pad=1.0,w_pad=0.8); save(fig,"Fig9")
else: plt.close(fig)
# ================= Fig 10 robustness of biology
B=pd.read_csv(O+"batch_library_adjusted_pairs.csv"); H=pd.read_csv(O+"heme_background_sensitivity.csv"); S=pd.read_csv(O+"subsampling_power_bh.csv"); ST=pd.read_csv(O+"subtype_sensitivity.csv")
fig,axs=plt.subplots(3,2,figsize=(7.3,8.6)); a,b,c,d,e,f=axs.ravel()
a.scatter(B.coef_unadj,B.coef_adj,s=12,color=BLUE,edgecolor="white",linewidth=0.3); lim=[min(B.coef_unadj.min(),B.coef_adj.min())-0.05,-0.2]; a.plot(lim,lim,color="#555555",ls="--",lw=0.8); a.set_xlim(lim); a.set_ylim(lim)
a.set_xlabel("Lineage effect, unadjusted"); a.set_ylabel("Lineage effect, batch-adjusted")
w=H.pivot_table(index=["Gene","Lineage"],columns="background",values="selectivity"); b.scatter(w.all_others,w.exclude_other_heme,s=14,color=VERM,edgecolor="white",linewidth=0.3); l2=[w.min().min()-0.05,-0.3]; b.plot(l2,l2,color="#555555",ls="--",lw=0.8)
b.axhline(-0.5,color=BLUE,ls=":",lw=1); b.set_xlim(l2); b.set_ylim(l2); b.set_xlabel("Selectivity, all other lineages"); b.set_ylabel("Selectivity, excluding the other\nblood lineage")
vs=[S[S.k==k_].frac_retained.values for k_ in (10,20)]; c.hist(vs,bins=np.linspace(0,1,11),color=[SKY,BLUE],label=[f"{k_} lines per subsample (median {np.median(v):.2f})" for k_,v in zip((10,20),vs)],edgecolor="white")
c.set_xlabel("Fraction of subsamples retaining the candidate"); c.set_ylabel("Candidate pairs"); c.legend(frameon=False,loc="upper left"); c.yaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
g=ST.groupby("Subtype").agg(n=("selectivity","size"),med=("selectivity","median")).sort_values("med"); yv=np.arange(len(g))
d.barh(yv,g.med,color=[VERM if s_ in set(ST[ST.Lineage=="lymphoid"].Subtype) else ORANGE for s_ in g.index],height=0.7); d.set_yticks(yv); d.set_yticklabels(["\n".join(textwrap.wrap(str(s_),24)) for s_ in g.index])
d.axvline(-0.5,color=BLUE,ls=":",lw=1); d.set_xlabel("Median selectivity of candidate pairs"); [d.text(m-0.02,i,f"n = {n_}",va="center",ha="right",fontsize=8) for i,(m,n_) in enumerate(zip(g.med,g.n))]; d.set_xlim(min(g.med)-0.35,0); d.legend(handles=[Patch(color=VERM,label="Lymphoid subtypes"),Patch(color=ORANGE,label="Myeloid subtypes")],frameon=False,loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=2)
GC=pd.read_csv(O+"growth_confound_pairs.csv"); GC["heme"]=GC.Lineage.isin(["lymphoid","myeloid"]); FT=pd.read_csv(O+"final_targets.csv")[["Gene","Lineage","Selectivity"]]; GC=GC.merge(FT,on=["Gene","Lineage"])
for flag,col,lab in ((False,SKY,"Non-blood lineages"),(True,VERM,"Lymphoid and myeloid")):
    q=GC[GC.heme==flag]; e.scatter(q.coef_unadjusted,q.coef_growth_adjusted,s=12,color=col,edgecolor="white",linewidth=0.3,label=lab)
lim=[GC[["coef_unadjusted","coef_growth_adjusted"]].min().min()-0.05,0.0]; e.plot(lim,lim,color="#555555",ls="--",lw=0.8); e.set_xlim(lim); e.set_ylim(lim); e.set_xlabel("Lineage effect, unadjusted"); e.set_ylabel("Lineage effect, growth-pattern-adjusted"); e.legend(frameon=False,loc="lower right")
q=GC[GC.heme]; f.scatter(q.Selectivity,q.median_diff_vs_suspension,s=14,color=VERM,edgecolor="white",linewidth=0.3); l3=[min(q.Selectivity.min(),q.median_diff_vs_suspension.min())-0.05,0.0]; f.plot(l3,l3,color="#555555",ls="--",lw=0.8); f.axhline(-0.5,color=BLUE,ls=":",lw=1)
f.set_xlim(l3); f.set_ylim(l3); f.set_xlabel("Selectivity, all other lineages"); f.set_ylabel("Selectivity, non-blood suspension\nlines only")
letters(axs.ravel(),-0.16); fig.tight_layout(pad=0.6,h_pad=1.6,w_pad=1.0); save(fig,"Fig10")
