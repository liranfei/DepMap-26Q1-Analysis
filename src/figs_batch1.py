from figlib import *
from matplotlib.patches import FancyBboxPatch
from adjustText import adjust_text
from matplotlib import cm
F,T,lc=load_all(); S=json.load(open(O+"summary.json")); lc=lc[lc.n>=5]
# ================= Fig 1 workflow
fig,ax=plt.subplots(figsize=(7.3,5.6)); ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")
def box(x,y,w,h,txt,fc="#f2f2f2",ec="#444444",fs=8,bold=False):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.4,rounding_size=1.2",fc=fc,ec=ec,lw=0.9)); ax.text(x+w/2,y+h/2,txt,ha="center",va="center",fontsize=fs,fontweight="bold" if bold else "normal",linespacing=1.25)
def arrow(x1,y1,x2,y2): ax.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="-|>",color="#444444",lw=0.9))
box(12,88,76,9,"DepMap 26Q1 CRISPR gene effect (Chronos)\n1,208 models × 18,531 genes; Model.csv metadata (all IDs matched)",fc="#e8f1f8")
box(12,74,76,9,"Cohort: remove 16 models annotated Non-Cancerous\n1,192 cancer cell lines in 29 lineages; 26 lineages with ≥ 5 lines are targets")
box(12,60,76,9,"Test every gene in every eligible lineage (target vs. all other lines)\none-sided Mann–Whitney U; missing values removed per gene\n467,091 tests")
box(12,46,76,9,"Benjamini–Hochberg correction over the whole test family\n9,805 pairs with q < 0.05")
box(12,32,76,9,"Effect-size filter: target median Chronos < −1 and selectivity < −0.5\n94 candidate gene–lineage pairs (58 genes, 21 lineages)",fc="#fbe9dc",bold=False)
for y in (88,74,60,46): arrow(50,y-0.4,50,y-4.2)
box(1,3,22.5,22,"Permutation tests\n1,000 shuffles each,\nfull pipeline rerun;\nunrestricted and\npatient-block null:\n0 candidates",fc="#e9f5ef",fs=8)
box(26,3,22.5,22,"Sensitivity\nBY / Holm, thresholds,\ntop-N genes, t-tests,\none model per patient,\nbootstrap; batch and\ngrowth checks",fc="#e9f5ef",fs=8)
box(51,3,22.5,22,"Replication and\nstability\nSanger Project Score\n(separate models\nand patients);\nDepMap 22Q1;\nRNAi (DEMETER2)",fc="#e9f5ef",fs=8)
box(76,3,22.5,22,"Exploratory\nannotation\ngenotype, TCGA,\nenrichment,\nco-dependency,\nliterature",fc="#f2f2f2",fs=8)
for x in (12.25,37.25,62.25,87.25): arrow(50,31.6,x,25.6)
save(fig,"Fig1")
# ================= Fig 2 lineages
D=lc.set_index("lineage").join(F.groupby("Lineage").size().rename("pairs")).fillna(0).astype(int).sort_values(["pairs","n"])
fig,(a1,a2)=plt.subplots(1,2,figsize=(7.3,6.4),sharey=True,gridspec_kw={"width_ratios":[1,1.1]}); y=np.arange(len(D)); low=(D.n<10).values
a1.barh(y,D.n,color=GREY,height=0.75); [a1.text(n+2,i,str(n),va="center",fontsize=8) for i,n in enumerate(D.n)]
a1.set_yticks(y); a1.set_yticklabels([nm(i)+(" *" if l else "") for i,l in zip(D.index,low)]); a1.set_xlabel("Cancer cell lines per lineage"); a1.set_xlim(0,145)
cols=cm.cividis(np.linspace(0.1,0.95,len(D))); a2.barh(y,D.pairs,color=cols,height=0.75); [a2.text(p+0.3,i,str(p),va="center",fontsize=8) for i,p in enumerate(D.pairs)]
a2.set_xlabel("Candidate gene–lineage pairs"); a2.set_xlim(0,26)
fig.tight_layout(pad=0.6,rect=(0,0,1,0.97))
for a,l,x in ((a1,"A",0.01),(a2,"B",a2.get_position().x0-0.02)): fig.text(x,0.975,l,fontsize=12,fontweight="bold",va="bottom")   # panel letters in figure coordinates (inside the canvas)
save(fig,"Fig2")
# ================= Fig 3 volcano
T["nlq"]=-np.log10(T.q_value.clip(lower=1e-300)); fin=set(zip(F.Gene,F.Lineage)); T["final"]=[(g,l) in fin for g,l in zip(T.Gene,T.Lineage)]
fig,ax=plt.subplots(figsize=(7.3,5.3)); g=T[~T.final]; r=T[T.final]
ax.scatter(g.Selectivity,g.nlq,s=3,c="#c8c8c8",alpha=0.6,linewidths=0,rasterized=True,label=f"Other tested pairs (n = {len(g):,})")
ax.scatter(r.Selectivity,r.nlq,s=16,c=VERM,edgecolor="white",linewidth=0.3,label=f"Candidate pairs (n = {len(r)})")
ax.axvline(-0.5,color=BLUE,ls="--",lw=1,label="Selectivity = −0.5"); ax.axhline(-np.log10(0.05),color="#555555",ls=":",lw=1,label="q = 0.05")
key={"KRAS","SOX10","CTNNB1","NMNAT1","NAMPT","HNF1B","CBFB","IRF4","MDM2","PAX8","BRAF","MYB"}
lab=r[(r.sym.isin(key))|(r.nlq.rank(ascending=False)<=6)].sort_values("nlq",ascending=False).drop_duplicates("sym")
tx=[ax.text(x,yv,s + (f" ({l})" if s=="MDM2" else ""),fontsize=8) for x,yv,s,l in zip(lab.Selectivity,lab.nlq,lab.sym,lab.Lineage)]   # MDM2 is a candidate in four lineages: name the one labelled
adjust_text(tx,x=lab.Selectivity.values,y=lab.nlq.values,ax=ax,arrowprops=dict(arrowstyle="-",color="0.4",lw=0.5),expand=(1.4,1.6),iter_lim=1000)   # fixed iteration count: the default (1 s time limit) makes label placement machine-dependent
ax.set_xlim(-1.75,0.7); ax.set_xlabel("Selectivity (median Chronos, target lineage − all other cancer cell lines)"); ax.set_ylabel("−log$_{10}$ (BH-adjusted q)"); ax.legend(loc="upper right",frameon=False)
fig.tight_layout(pad=0.6); save(fig,"Fig3")
# ================= Fig 4 efficacy vs selectivity
fig,ax=plt.subplots(figsize=(7.3,5.3)); h=F[F.Lineage.isin(HEME)]; s=F[~F.Lineage.isin(HEME)]
ax.scatter(s.Selectivity,s.Chronos_median,s=22,c=BLUE,edgecolor="white",linewidth=0.3,label="Solid-tumor lineages")
ax.scatter(h.Selectivity,h.Chronos_median,s=22,c=VERM,edgecolor="white",linewidth=0.3,label="Lymphoid / myeloid")
lo=F[F.n_target<10]; ax.scatter(lo.Selectivity,lo.Chronos_median,s=60,facecolors="none",edgecolors="black",linewidth=0.8,label="Pair with < 10 target observations")
ax.axvline(-0.5,color="#555555",ls="--",lw=0.9); ax.axhline(-1,color="#555555",ls="--",lw=0.9)
lab=F[F.sym.isin({"KRAS","SOX10","HNF1B","MDM2","CTNNB1","CFLAR","MYB","CBFB","PAX8","NMNAT1","MYC","FOXA1","DHFR","TUBB4B","DCAF7"})].drop_duplicates("sym")
tx=[ax.text(x,yv,sy + (f" ({l})" if sy=="MDM2" else ""),fontsize=8) for x,yv,sy,l in zip(lab.Selectivity,lab.Chronos_median,lab.sym,lab.Lineage)]
adjust_text(tx,x=lab.Selectivity.values,y=lab.Chronos_median.values,ax=ax,arrowprops=dict(arrowstyle="-",color="0.4",lw=0.5),iter_lim=1000)
ax.set_xlabel("Selectivity (target lineage − all other cancer cell lines, median Chronos)"); ax.set_ylabel("Median Chronos, target lineage"); ax.legend(frameon=False,loc="lower center",bbox_to_anchor=(0.42,0.0))
fig.tight_layout(pad=0.6); save(fig,"Fig4")
# ================= Fig 5 top 20
t=F.sort_values("Selectivity").head(20).iloc[::-1]; fig,ax=plt.subplots(figsize=(7.3,5.6)); y=np.arange(len(t))
ax.barh(y,-t.Selectivity,color=[VERM if l in HEME else BLUE for l in t.Lineage],height=0.72)
ax.set_yticks(y); ax.set_yticklabels([f"{a} ({nm(b)})"+(" *" if n<10 else "") for a,b,n in zip(t.sym,t.Lineage,t.n_target)])
for i,(s_,q_) in enumerate(zip(t.Selectivity,t.q_value)): ax.text(-s_+0.02,i,f"{s_:.2f}".replace("-","−"),va="center",fontsize=8)
ax.set_xlabel("−Selectivity (larger = more lineage-selective)"); ax.set_xlim(0,1.9)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=BLUE,label="Solid-tumor lineage"),Patch(color=VERM,label="Lymphoid / myeloid")],frameon=False,loc="lower right")
fig.tight_layout(pad=0.6); save(fig,"Fig5")
