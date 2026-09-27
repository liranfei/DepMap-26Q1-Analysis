import os
from figlib import *
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform
from matplotlib.patches import Patch
F,T,lc=load_all(); P2=json.load(open(O+"part2_summary.json")); PP=json.load(open(O+"permutation_summary_genome.json")); Pm=pd.read_csv(O+"permutation_counts_genome.csv")
def letters(axs,dx=-0.02,dy=1.02):
    for a,l in zip(axs,"ABCDEFGH"): a.text(dx,dy,l,transform=a.transAxes,fontsize=12,fontweight="bold",va="bottom",ha="right")
# ---------- S1 Fig permutation
fig,axs=plt.subplots(1,2,figsize=(7.3,3.4)); a,b=axs
for ax,col_,obs,lab,cc in [(a,"n_q_lt_0_05",PP["observed_q_lt_0_05"],"Pairs with q < 0.05 (permuted labels)",BLUE),(b,"n_final_pairs",PP["observed_final"],"Candidate pairs (permuted labels)",VERM)]:
    v=Pm[col_].values; ax.hist(v,bins=np.arange(-0.5,6.5,1),color=GREY,edgecolor="white"); ax.set_xlim(-0.5,6.5); ax.set_xlabel(lab); ax.set_ylabel("Permutations (n = 1,000)"); ax.set_ylim(0,1150)
    ax.annotate(f"Observed = {obs:,}\n(off scale)",xy=(6.4,600),xytext=(3.0,600),ha="center",va="center",fontsize=8,color=cc,arrowprops=dict(arrowstyle="-|>",color=cc,lw=1.4))
    ax.text(0.97,0.95,f"Null maximum = {int(v.max())}",transform=ax.transAxes,ha="right",va="top",fontsize=8)
letters(axs,-0.16); fig.tight_layout(pad=0.6,w_pad=1.5); save(fig,"S1_Fig")
# ---------- S2 Fig co-dependency
C=pd.read_csv(O+"codependency_corr.csv",index_col=0); sym=[g.split(" (")[0] for g in C.index]; C.index=sym; C.columns=sym
d=1-C.values; np.fill_diagonal(d,0); d=(d+d.T)/2; order=leaves_list(linkage(squareform(np.clip(d,0,None),checks=False),"average")); Cc=C.iloc[order,order]  # average linkage on correlation distance 1-r (Ward requires Euclidean distances)
fig,ax=plt.subplots(figsize=(7.3,7.0)); im=ax.imshow(Cc.values,cmap="RdBu_r",vmin=-1,vmax=1); ax.set_xticks(range(len(Cc))); ax.set_yticks(range(len(Cc)))
ax.set_xticklabels(Cc.columns,rotation=90,fontsize=8); ax.set_yticklabels(Cc.index,fontsize=8); cb=fig.colorbar(im,ax=ax,fraction=0.035,pad=0.02); cb.set_label("Pearson r"); cb.ax.tick_params(labelsize=8)
for s_ in ax.spines.values(): s_.set_visible(True)
fig.tight_layout(pad=0.6); save(fig,"S4_Fig")
# ---------- S3 Fig common essential annotation
tg=sorted(F.Sym if False else F.sym.unique()); ce={g.split(" (")[0] for g in pd.read_csv(os.environ.get("DEPMAP_DIR",".")+"/CRISPRInferredCommonEssentials.csv").iloc[:,0]}
uni=T.sym.unique(); ce_t=np.mean([g in ce for g in tg]); ce_u=np.mean([g in ce for g in uni]); frac_lines=None
fig,axs=plt.subplots(1,2,figsize=(7.3,3.6)); a,b=axs
a.bar([0,1],[ce_u*100,ce_t*100],color=[GREY,VERM],width=0.6); a.set_xticks([0,1]); a.set_xticklabels([f"All tested genes\n(n = {len(uni):,})",f"Candidate genes\n(n = {len(tg)})"]); a.set_ylabel("Genes in official common-essential list (%)")
[a.text(i,v*100+1.5,f"{v*100:.1f}%",ha="center",fontsize=8) for i,v in enumerate([ce_u,ce_t])]; a.set_ylim(0,70)
import run_pipeline as rp
df=rp.load(os.environ.get("DEPMAP_DIR",".")); df=df[df.primary_disease!="Non-Cancerous"].reset_index(drop=True); genes=[c for c in df.columns if c not in ("DepMap_ID","lineage","primary_disease")]
X=df[[g for g in genes if g.split(" (")[0] in set(tg)]].apply(pd.to_numeric,errors="coerce"); fr=(X<-1).mean(); isce=[c.split(" (")[0] in ce for c in X.columns]
b.hist([fr[[not i for i in isce]]*100,fr[isce]*100],bins=np.linspace(0,100,11),color=[BLUE,VERM],label=["Not in list","In list"],edgecolor="white"); b.set_xlabel("Cancer cell lines with Chronos < −1 (%)"); b.set_ylabel("Candidate genes"); b.legend(frameon=False)
letters(axs,-0.16); fig.tight_layout(pad=0.6,w_pad=1.5); save(fig,"S3_Fig")
# ---------- S4 Fig Mann-Whitney vs t-test on the 94
A=pd.read_csv(O+"ttest_alt_all_tests.csv.gz"); K=F.merge(A,on=["Gene","Lineage"],how="left")
fig,axs=plt.subplots(1,2,figsize=(7.3,3.6)); a,b=axs
a.scatter(-np.log10(K.q_value.clip(lower=1e-300)),-np.log10(K.q.clip(lower=1e-300)),s=14,color=BLUE,edgecolor="white",linewidth=0.3); a.plot([0,40],[0,40],color="#555555",ls="--",lw=0.8); a.set_xlim(0,40); a.set_ylim(0,160)
a.axhline(-np.log10(0.05),color=VERM,ls=":",lw=1); a.set_xlabel("−log$_{10}$ q, Mann–Whitney (primary)"); a.set_ylabel("−log$_{10}$ q, Student t-test")
b.scatter(K.Selectivity,K.mean_diff,s=14,color=VERM,edgecolor="white",linewidth=0.3); lim=[-1.7,-0.2]; b.plot(lim,lim,color="#555555",ls="--",lw=0.8); b.axhline(-0.5,color=BLUE,ls=":",lw=1); b.axvline(-0.5,color=BLUE,ls=":",lw=1); b.set_xlim(lim); b.set_ylim(lim)
b.set_xlabel("Difference in medians (primary)"); b.set_ylabel("Difference in means")
letters(axs,-0.16); fig.tight_layout(pad=0.6,w_pad=1.5); save(fig,"S2_Fig")
# ---------- S5 Fig skewness (descriptive)
sk=pd.read_csv(O+"selectivity_vs_skewness.csv") if False else None
gen_all=[c for c in genes]; skew=df[gen_all].skew(); skew.index=[c.split(" (")[0] for c in skew.index]; best=F.groupby("sym").Selectivity.min()
fig,axs=plt.subplots(1,2,figsize=(7.3,3.6)); a,b=axs; bg=skew.drop(list(best.index),errors="ignore").dropna(); tv=skew[best.index].dropna()
vp=a.violinplot([bg.values,tv.values],showmedians=True,widths=0.8); [ (p.set_facecolor(cc),p.set_alpha(0.6)) for p,cc in zip(vp["bodies"],[GREY,VERM]) ]
a.set_xticks([1,2]); a.set_xticklabels([f"Other genes\n(n = {len(bg):,})",f"Candidate genes\n(n = {len(tv)})"]); a.set_ylabel("Skewness of gene effect across cell lines")
b.scatter(best[tv.index],tv,s=16,color=VERM,edgecolor="white",linewidth=0.3); b.set_xlabel("Selectivity (most negative per gene)"); b.set_ylabel("Skewness"); b.text(0.03,0.05,f"Spearman ρ = {P2['skewness_descriptive']['spearman_rho']:.2f}",transform=b.transAxes,fontsize=8)
letters(axs,-0.16); fig.tight_layout(pad=0.6,w_pad=1.5); save(fig,"S5_Fig")
