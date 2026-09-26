import os
import json, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from PIL import Image
mpl.rcParams.update({"font.family":"Arial","font.size":8,"axes.labelsize":9,"xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,"axes.linewidth":0.8,
  "pdf.fonttype":42,"mathtext.fontset":"custom","mathtext.rm":"Arial","mathtext.it":"Arial:italic","mathtext.default":"regular","axes.spines.top":False,"axes.spines.right":False})
BLUE,VERM,ORANGE,GREEN,SKY,PINK,GREY="#0072B2","#D55E00","#E69F00","#009E73","#56B4E9","#CC79A7","#8c8c8c"
O=os.environ.get("RESULTS_DIR","results")+"/"; OUT=os.environ.get("FIG_DIR","figures")+"/"
HEME={"lymphoid","myeloid"}
def save(fig,name):
    import matplotlib.text as mt
    bad=[(t.get_text()[:25],t.get_fontsize()) for t in fig.findobj(mt.Text) if t.get_text().strip() and t.get_visible() and not (8<=t.get_fontsize()<=12)]
    assert not bad, (name,"font size outside 8-12 pt:",bad[:5])
    p=OUT+name+".tif"; fig.savefig(p,dpi=300,pil_kwargs={"compression":"tiff_lzw"},facecolor="white"); plt.close(fig)
    im=Image.open(p).convert("RGB"); _a=np.asarray(im); _nz=(_a<245).any(axis=2); _ys,_xs=np.where(_nz)
    im=im.crop((_xs.min(),_ys.min(),_xs.max()+1,_ys.max()+1)); _new=Image.new("RGB",(im.width+24,im.height+24),"white"); _new.paste(im,(12,12)); im=_new
    im.save(p,compression="tiff_lzw",dpi=(300,300)); w,h=im.size
    assert 789<=w<=2250 and h<=2625, (name,w,h); print(name,w,h,"OK")
def load_all():
    F=pd.read_csv(O+"final_targets.csv"); T=pd.read_csv(O+"all_tests.csv"); lc=pd.read_csv(O+"lineage_counts.csv"); lc.columns=["lineage","n"]
    F["sym"]=F.Gene.str.split(" (",regex=False).str[0]; T["sym"]=T.Gene.str.split(" (",regex=False).str[0]; return F,T,lc
nm=lambda s:{"cns/brain":"CNS/Brain"}.get(s,s.title().replace("And","and").replace("Of","of"))
