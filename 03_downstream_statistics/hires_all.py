import os, re, numpy as np, pandas as pd, pickle
from scipy import stats
from statsmodels.stats.multitest import multipletests
BASE="/sessions/modest-sleepy-pasteur/mnt/NanoString GeoMx data analysis"
HD=f"{BASE}/data/2026-01-28_melanoma_GeoMx_hires"
OUT="/sessions/modest-sleepy-pasteur/mnt/outputs/hires_out"
meta=pd.read_csv(f"{BASE}/data/2025/2025-01-27_GeoMxmetadata.txt",sep="\t")
meta["dcc"]=meta["ID"].astype(str); meta["AOI"]=meta["Segment_Tags"].replace({"neg/neg":"neg_neg"})
meta=meta.set_index("dcc")
SUBS=["LMM","SITU","SSM","NOD","AMEL","MET"]
PAIR={"Melanoma":"S100","CAF":"neg_neg","Endothelial":"neg_neg","Mural_cell":"neg_neg",
      "T-CD8":"CD45","T-CD4":"CD45","NK":"CD45","B_cell":"CD45","Macrophage":"CD45"}
# epidermal / keratinocyte contamination gene set (spillover, not cell-intrinsic)
EPI_PAT=re.compile(r"^(KRT\d|KRT[0-9]|KRTAP|KRTDAP|LCE\d|SPRR|FLG|LOR|IVL|TGM[135]|CDSN|CASP14|"
 r"DSG[1-4]|DSC[1-3]|DSP$|PKP[1-3]|PPL|EVPL|PI3$|SERPINB[3-9]|CSTA|CSTB|SBSN|DMKN|CALML[35]|"
 r"S100A7|S100A8|S100A9|GJB[26]|CLDN1$|PERP|CNFN|ASPRV1|ABCA12|ALOX12B|HRNR|LGALS7|SLURP1|A2ML1|KLK[5-9]|DSG|KRT)")
EPI_EXTRA={"KRT1","KRT2","KRT5","KRT6A","KRT6B","KRT9","KRT10","KRT14","KRT15","KRT16","KRT17","DSG1","DSC1","SFN","S100A7","DMKN","GJB2","CALML5","CALML3","SPRR1B","SPRR2A","LCE3D","FLG2","KRTDAP","CSTA","SBSN","PI3"}
def is_epi(g): return bool(EPI_PAT.match(g)) or g in EPI_EXTRA
summary={}
allcand=[]
for cell,comp in PAIR.items():
    fn=f"CIBERSORTxHiRes_NA_{cell}_Window36.txt"
    df=pd.read_csv(f"{HD}/{fn}",sep="\t",index_col=0); df.columns=df.columns.astype(str)
    keep=[c for c in df.columns if c in meta.index and meta.loc[c,"Pathophysiology"]=="Melanoma"
          and meta.loc[c,"AOI"]==comp and meta.loc[c,"Path_sub_type"] in SUBS]
    sub=meta.loc[keep,"Path_sub_type"]
    X=df[keep].astype(float)
    X=X[(X==1).mean(1)<0.4]; X=np.log2(X+1); X=X[X.std(1)>0]
    grp={s:[c for c in keep if sub[c]==s] for s in SUBS if (sub==s).sum()>=3}
    labs=list(grp.keys())
    rows=[]
    for g in X.index:
        vals=[X.loc[g,grp[s]].values for s in labs]
        try: H,p=stats.kruskal(*vals)
        except Exception: H,p=np.nan,1.0
        meds={s:np.median(X.loc[g,grp[s]]) for s in labs}
        hi=max(meds,key=meds.get); lo=min(meds,key=meds.get)
        rows.append([g,p,hi,lo,meds[hi]-meds[lo],is_epi(g)])
    R=pd.DataFrame(rows,columns=["gene","p","peak_subtype","low_subtype","range_log2","epidermal"])
    R["q"]=multipletests(R["p"].fillna(1),method="fdr_bh")[1]
    R["cell"]=cell; R["compartment"]=comp
    R=R.sort_values("p")
    R.to_csv(f"{OUT}/scan_{cell}_{comp}.csv",index=False)
    Rci=R[~R["epidermal"]]      # cell-intrinsic only
    nsig=(R["q"]<0.05).sum(); nci=(Rci["q"]<0.05).sum()
    print(f"{cell:11s} {comp:8s} AOIs={len(keep):3d} genes={X.shape[0]:4d}  q<.05 all={nsig:3d} cell-intrinsic={nci:3d}  top: "
          + ", ".join(f"{r.gene}({r.peak_subtype})" for r in Rci.head(5).itertuples()))
    summary[cell]=(comp,labs,grp,X,R)
    for r in Rci[Rci["q"]<0.1].head(8).itertuples():
        allcand.append([cell,comp,r.gene,r.peak_subtype,r.low_subtype,round(r.range_log2,3),float(r.p),float(r.q)])
pickle.dump(summary,open(f"{OUT}/summary_all.pkl","wb"))
cand=pd.DataFrame(allcand,columns=["cell","compartment","gene","peak_subtype","low_subtype","range_log2","p","q"])
cand.to_csv(f"{OUT}/candidate_genes_all.csv",index=False)
print("\ncandidates (cell-intrinsic, q<0.1):", cand.shape[0])
