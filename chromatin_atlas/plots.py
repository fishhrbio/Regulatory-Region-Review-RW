"""Publication-exportable original charts with embedded authorship marks."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from . import AUTHOR

COLORS = ["#c34154","#df7b88","#9766b0","#e49330","#d5a43a","#d5bc57","#dcd488",
          "#5688b0","#468b78","#286b5b","#9bc6a4","#7f8091","#d3d7db","#e9e9eb","#f2f2f4"]


def save(fig,path):
    fig.text(.5,.024,AUTHOR,ha="center",fontsize=9,color="#495a62")
    fig.savefig(path.with_suffix(".png"),dpi=180,facecolor="white")
    fig.savefig(path.with_suffix(".svg"),metadata={"Date":None,"Creator":AUTHOR})
    plt.close(fig)


def make_figures(result,folder):
    folder.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"axes.spines.top":False,
                         "axes.spines.right":False,"svg.hashsalt":"waite-chromatin-atlas-1"})
    primary = result["scenarios"]["primary"]
    fig,ax = plt.subplots(figsize=(12,6.3))
    fig.subplots_adjust(left=.12,right=.70,top=.80,bottom=.22)
    left = [0.,0.]
    for r,color,name in zip(primary["coverage"],COLORS,result["state_names"]):
        values = [r["gm_fraction"]*100,r["k_fraction"]*100]
        ax.barh([1,0],values,left=left,color=color,label=f'{r["state"]}  {name}',height=.42)
        left = [a+b for a,b in zip(left,values)]
    ax.set(yticks=[1,0],yticklabels=["GM12878","K562"],xlim=(0,100),ylim=(-.6,1.6),xlabel="Coverage (% of common eligible territory)")
    ax.legend(bbox_to_anchor=(1.02,1.05),loc="upper left",frameon=False,fontsize=8.5)
    fig.suptitle("Chromatin composition on chromosome 22",x=.12,ha="left",y=.96,fontsize=18,weight="bold")
    fig.text(.12,.87,f'hg19 | {primary["eligible_bp"]:,} eligible bp | states 14/15 excluded in either line',fontsize=10)
    fig.text(.12,.10,"Published ChIP-seq-derived ChromHMM states; base-pair weighted. Low-signal state 13 is retained.",fontsize=9)
    save(fig,folder/"01_state_coverage")
    fig,ax = plt.subplots(figsize=(12,5.7))
    fig.subplots_adjust(left=.20,right=.95,top=.77,bottom=.31)
    labels = ["Active promoter (1)","Strong enhancer (4+5)"]
    for y,name in enumerate(["strong_enhancer","active_promoter"]):
        r = primary["overlap"][name]
        start = 0
        for key,color in zip(["gm_only_bp","shared_bp","k_only_bp"],["#276f83","#373d59","#da9b3a"]):
            v = r[key]/r["union_bp"]*100
            ax.barh(y,v,left=start,color=color,height=.48)
            if v > 5: ax.text(start+v/2,y,f"{v:.1f}%",ha="center",va="center",color="white",weight="bold")
            start += v
        ax.text(100,y+.36,f'J = {r["jaccard"]:.3f}  |  union {r["union_bp"]:,} bp',ha="right",fontsize=10)
    ax.set(yticks=[1,0],yticklabels=labels,xlim=(0,100),ylim=(-.6,1.7),xlabel="Fraction of each regulatory class union (%)")
    handles = [plt.Rectangle((0,0),1,1,color=c) for c in ["#276f83","#373d59","#da9b3a"]]
    ax.legend(handles,["GM12878 only","Shared","K562 only"],loc="upper center",bbox_to_anchor=(.5,-.23),ncol=3,frameon=False)
    fig.suptitle("Sharing of regulatory annotations",x=.08,ha="left",y=.96,fontsize=18,weight="bold")
    fig.text(.08,.86,"Exact overlap on common eligible hg19 chr22 territory; J = shared bp / union bp",fontsize=10)
    fig.text(.08,.10,"Cross-sectional annotation specificity does not demonstrate enhancer function or temporal switching.",fontsize=9)
    save(fig,folder/"02_regulatory_sharing")
    fig,ax = plt.subplots(figsize=(9,8))
    fig.subplots_adjust(left=.12,right=.90,top=.84,bottom=.16)
    matrix = np.array(primary["state_matrix_bp"],dtype=float)
    totals = matrix.sum(axis=1,keepdims=True)
    fractions = np.divide(matrix,totals,out=np.zeros_like(matrix),where=totals>0)
    image = ax.imshow(fractions,vmin=0,vmax=1,cmap="Blues")
    for i in range(15):
        if totals[i,0] == 0:
            ax.axhspan(i-.5,i+.5,color="#ddd")
    ax.set(xticks=range(15),xticklabels=range(1,16),yticks=range(15),yticklabels=range(1,16),xlabel="K562 state",ylabel="GM12878 state")
    fig.colorbar(image,ax=ax,label="Fraction of each GM12878 state row")
    fig.suptitle("Paired state map",x=.12,ha="left",y=.96,fontsize=18,weight="bold")
    fig.text(.12,.89,"Same genomic bases across cell lines; this is not a temporal transition matrix.",fontsize=10)
    fig.text(.12,.08,"Grey rows: no eligible bases. States 14/15 excluded. See README for state labels.",fontsize=9)
    save(fig,folder/"03_paired_states")
