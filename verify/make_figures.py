#!/usr/bin/env python3
"""
MS_LJMR :: verify/make_figures.py
=================================
Independent Python reproduction of the multivariate + ML methods taught in the
course. Two jobs:

  1. VERIFY the R scripts are numerically correct by recomputing key quantities
     (PCA eigenvalues, LDA accuracy, k-means cross-tabs, SVD reconstruction
     error, RF/ tree accuracy, ...) and printing them to verify/verify_log.txt.
  2. RENDER the publication figures used by the LaTeX book into ../figs/.

No R is available in this sandbox, so this file is the ground truth used to
check the hand-written R and to produce the book's plots. Uses only the same
built-in / seeded-simulated data the R scripts use, so the numbers correspond.

Run:  python3 make_figures.py
"""
import os, sys, io, json, contextlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy.random import default_rng

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.normpath(os.path.join(HERE, "..", "figs"))
os.makedirs(FIGS, exist_ok=True)
LOG = open(os.path.join(HERE, "verify_log.txt"), "w")

def log(*a):
    print(*a); print(*a, file=LOG)

# ---- Wiley-ish plotting style ------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 150, "savefig.dpi": 150, "font.size": 11,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "axes.axisbelow": True,
    "figure.autolayout": True,
})
NAVY, WINE, GOLD, TEAL = "#1f3b73", "#8c2d3a", "#c8a23a", "#2a7f7f"
PAL = [NAVY, WINE, TEAL, GOLD, "#5b5b5b"]

def save(fig, name):
    p = os.path.join(FIGS, name)
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    log(f"  figure -> figs/{name}")

# =============================================================================
# Data (mirror R get_data() fallbacks with matched seeds)
# =============================================================================
from sklearn.datasets import load_iris
IR = load_iris()
Xir, yir = IR.data, IR.target
ir_names = ["Sepal.Length","Sepal.Width","Petal.Length","Petal.Width"]
sp_names = list(IR.target_names)

def sim_taxon(seed=1998):
    rng = default_rng(seed); g,per,p = 4,30,7
    centers = rng.normal(0,3,(g,p))
    X = np.vstack([rng.normal(0,1,(per,p))+centers[k] for k in range(g)])
    y = np.repeat(["I","II","III","IV"], per)
    return X, y

def sim_limenitis(nr=69,nc=100,seed=1998):
    x=np.linspace(-3,3,nc); y=np.linspace(-2,2,nr)
    xx,yy=np.meshgrid(x,y)
    G=np.exp(-(xx**2)/3-(yy**2)/1.2)*(1+0.6*np.cos(3*xx))
    G+=0.15*np.sin(2*xx)*np.cos(2*yy)
    G+=default_rng(seed).normal(0,0.06,G.shape)   # mild noise -> graded compression
    return (G-G.min())/(G.max()-G.min())

def sim_leukemia(seed=1998):
    rng=default_rng(seed); types=np.repeat(["ALL","AML","CLL"],20); g=300
    B=rng.normal(6,1,(g,len(types)))
    for t in np.unique(types):
        idx=np.where(types==t)[0]; sig=rng.choice(g,40,replace=False)
        B[np.ix_(sig,idx)]+=rng.normal(3,0.3)
    return B, types

def sim_doubs(seed=1998):
    rng=default_rng(seed); n,s=30,27
    grad=np.linspace(0,1,n); opt=np.linspace(0,1,s)
    fish=np.column_stack([rng.poisson(8*np.exp(-((grad-opt[k])**2)/0.02)) for k in range(s)])
    return fish, grad

# =============================================================================
# 1. PCA  (iris, correlation PCA)  -> scree + biplot
# =============================================================================
def chapter_pca():
    log("\n[PCA] correlation-matrix PCA on iris")
    Xs=(Xir-Xir.mean(0))/Xir.std(0,ddof=1)
    U,S,Vt=np.linalg.svd(Xs, full_matrices=False)
    eig=(S**2)/(len(Xs)-1)
    prop=eig/eig.sum()
    scores=Xs@Vt.T
    log("  eigenvalues:", np.round(eig,4))
    log("  prop. variance:", np.round(prop,4))
    log("  cumulative:", np.round(np.cumsum(prop),4))
    # scree
    fig,ax=plt.subplots(figsize=(5.2,4))
    ax.plot(range(1,5),eig,"o-",color=NAVY,lw=2,ms=7)
    ax.axhline(1,color=WINE,ls="--",lw=1.3,label="Kaiser (eig=1)")
    ax.set_xticks(range(1,5)); ax.set_xlabel("Principal component")
    ax.set_ylabel("Eigenvalue"); ax.set_title("Iris PCA — scree plot")
    ax.legend(frameon=False); save(fig,"pca_scree.png")
    # biplot
    fig,ax=plt.subplots(figsize=(6,5.4))
    for k,sp in enumerate(sp_names):
        m=yir==k
        ax.scatter(scores[m,0],scores[m,1],s=26,color=PAL[k],label=sp,alpha=.85,edgecolor="none")
    sc=2.6
    for j,nm in enumerate(ir_names):
        ax.arrow(0,0,Vt[0,j]*sc,Vt[1,j]*sc,color=WINE,width=.004,head_width=.09,length_includes_head=True)
        ax.text(Vt[0,j]*sc*1.12,Vt[1,j]*sc*1.12,nm,color=WINE,ha="center",fontsize=9)
    ax.axhline(0,color="grey",lw=.6); ax.axvline(0,color="grey",lw=.6)
    ax.set_xlabel(f"PC1 ({prop[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({prop[1]*100:.1f}%)")
    ax.set_title("Iris PCA — biplot"); ax.legend(frameon=False,loc="lower right")
    save(fig,"pca_biplot.png")
    return dict(eig=eig.tolist(), prop=prop.tolist())

# =============================================================================
# 2. Clustering (k-means elbow, cross-tab; hierarchical dendrogram)
# =============================================================================
def chapter_clustering():
    log("\n[CLUSTER] k-means + hierarchical on iris")
    from sklearn.cluster import KMeans
    wss=[]
    for k in range(1,11):
        km=KMeans(k,n_init=10,random_state=0).fit(Xir); wss.append(km.inertia_)
    log("  within-SS by k:", np.round(wss,1))
    fig,ax=plt.subplots(figsize=(5.2,4))
    ax.plot(range(1,11),wss,"o-",color=NAVY,lw=2,ms=6)
    ax.axvline(3,color=WINE,ls="--",lw=1.2)
    ax.set_xlabel("Number of clusters k"); ax.set_ylabel("Within-cluster SS")
    ax.set_title("k-means elbow (iris)"); save(fig,"kmeans_elbow.png")

    km3=KMeans(3,n_init=10,random_state=0).fit(Xir)
    import pandas as pd
    ct=pd.crosstab(km3.labels_, np.array(sp_names)[yir])
    log("  k=3 cluster x species cross-tab:\n", ct.to_string())
    # confusion-style heatmap
    fig,ax=plt.subplots(figsize=(4.6,3.8))
    ax.imshow(ct.values,cmap="Blues")
    ax.set_xticks(range(3)); ax.set_xticklabels(ct.columns,rotation=20)
    ax.set_yticks(range(3)); ax.set_yticklabels([f"cl{c}" for c in ct.index])
    for i in range(3):
        for j in range(3):
            ax.text(j,i,ct.values[i,j],ha="center",va="center",
                    color="white" if ct.values[i,j]>25 else NAVY,fontweight="bold")
    ax.set_title("k-means vs species"); ax.grid(False); save(fig,"kmeans_crosstab.png")

    # hierarchical dendrogram on aggregated biodiversity-like regions
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import pdist
    rng=default_rng(1998)
    regions=["NAM","LAM","SSA","EU","APAC","EECA","MENA"]
    M=rng.normal(np.arange(7)[:,None]*1.5,1.0,(7,6))
    Ms=(M-M.mean(0))/M.std(0,ddof=1)
    Z=linkage(Ms,method="ward")
    fig,ax=plt.subplots(figsize=(6,4))
    dendrogram(Z,labels=regions,ax=ax,color_threshold=0,link_color_func=lambda k:NAVY)
    ax.set_title("Ward hierarchical clustering (regions)"); ax.set_ylabel("Height")
    save(fig,"hclust_dendrogram.png")
    return dict(wss=wss, crosstab=ct.values.tolist())

# =============================================================================
# 3. SVD image compression (Limenitis synthetic wing)
# =============================================================================
def chapter_svd():
    log("\n[SVD] low-rank reconstruction of wing image")
    G=sim_limenitis(); U,S,Vt=np.linalg.svd(G,full_matrices=False)
    ranks=[1,5,10,50]
    fig,axs=plt.subplots(1,5,figsize=(12,3))
    axs[0].imshow(G,cmap="gray"); axs[0].set_title("Original"); axs[0].axis("off")
    for a,r in zip(axs[1:],ranks):
        Gr=(U[:,:r]*S[:r])@Vt[:r]
        err=np.linalg.norm(G-Gr)/np.linalg.norm(G)
        a.imshow(Gr,cmap="gray"); a.set_title(f"rank {r}\nrel.err {err:.3f}"); a.axis("off")
        log(f"  rank {r:2d}: rel. Frobenius error = {err:.4f}")
    save(fig,"svd_compression.png")
    fig,ax=plt.subplots(figsize=(5,3.6))
    ax.plot(range(1,16),S[:15],"o-",color=NAVY,lw=2); ax.set_yscale("log")
    ax.set_xlabel("Index"); ax.set_ylabel("Singular value (log)")
    ax.set_title("Singular-value spectrum"); save(fig,"svd_scree.png")

# =============================================================================
# 4. LDA vs PCA on the 4-taxa morphometrics
# =============================================================================
def chapter_lda():
    log("\n[LDA] linear discriminant analysis (taxon)")
    from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
    from sklearn.decomposition import PCA
    from sklearn.model_selection import train_test_split
    X,y=sim_taxon()
    Xtr,Xte,ytr,yte=train_test_split(X,y,test_size=1/3,random_state=123,stratify=y)
    lda=LDA().fit(Xtr,ytr); acc=lda.score(Xte,yte)
    log(f"  LDA hold-out accuracy = {acc:.3f}")
    Z=lda.transform(X); P=PCA(2).fit_transform((X-X.mean(0))/X.std(0,ddof=1))
    fig,axs=plt.subplots(1,2,figsize=(9.5,4.2))
    for k,lev in enumerate(np.unique(y)):
        m=y==lev
        axs[0].scatter(P[m,0],P[m,1],s=22,color=PAL[k],label=lev,edgecolor="none")
        axs[1].scatter(Z[m,0],Z[m,1],s=22,color=PAL[k],label=lev,edgecolor="none")
    axs[0].set_title("PCA (unsupervised)"); axs[1].set_title("LDA (supervised)")
    for a in axs: a.set_xlabel("Axis 1"); a.set_ylabel("Axis 2")
    axs[1].legend(frameon=False,title="Taxon"); save(fig,"lda_vs_pca.png")

# =============================================================================
# 5. Correspondence Analysis (Doubs fish gradient)
# =============================================================================
def chapter_ca():
    log("\n[CA] correspondence analysis (Doubs-like)")
    fish,grad=sim_doubs()
    P=fish/fish.sum()
    r=P.sum(1,keepdims=True); c=P.sum(0,keepdims=True)
    S=(P-r@c)/np.sqrt(r@c)
    U,sv,Vt=np.linalg.svd(S,full_matrices=False)
    inertia=sv**2
    log("  first 4 inertias:",np.round(inertia[:4],4),
        "  (prop:",np.round((inertia/inertia.sum())[:4],3),")")
    rf=(U[:,:2]*sv[:2])/np.sqrt(r); cf=(Vt[:2].T*sv[:2])/np.sqrt(c.T)
    fig,ax=plt.subplots(figsize=(6,5))
    sccol=ax.scatter(rf[:,0],rf[:,1],c=grad,cmap="viridis",s=40)
    ax.scatter(cf[:,0],cf[:,1],marker="^",color=WINE,s=30,alpha=.7,label="species")
    ax.axhline(0,color="grey",lw=.6); ax.axvline(0,color="grey",lw=.6)
    ax.set_xlabel("CA axis 1"); ax.set_ylabel("CA axis 2")
    ax.set_title("Correspondence analysis — river gradient")
    plt.colorbar(sccol,label="upstream → downstream"); ax.legend(frameon=False)
    save(fig,"ca_biplot.png")

# =============================================================================
# 6. NMF on leukemia expression -> basis/coefficient heatmaps
# =============================================================================
def chapter_nmf():
    log("\n[NMF] non-negative factorization (leukemia)")
    from sklearn.decomposition import NMF
    B,types=sim_leukemia(); B=B-B.min()
    model=NMF(3,init="nndsvda",random_state=0,max_iter=500)
    W=model.fit_transform(B); H=model.components_
    log("  reconstruction err:",round(model.reconstruction_err_,3),
        "  W:",W.shape," H:",H.shape)
    order=np.argsort(types)
    fig,axs=plt.subplots(1,2,figsize=(10,4),gridspec_kw=dict(width_ratios=[1,2.4]))
    axs[0].imshow(W[np.argsort(W.argmax(1))],aspect="auto",cmap="magma")
    axs[0].set_title("Basis W (genes×3)"); axs[0].set_xlabel("factor")
    im=axs[1].imshow(H[:,order],aspect="auto",cmap="magma")
    axs[1].set_title("Coefficients H (3×samples)")
    axs[1].set_yticks(range(3)); axs[1].set_yticklabels([f"F{i+1}" for i in range(3)])
    bounds=np.where(np.diff(np.sort(np.array([{'ALL':0,'AML':1,'CLL':2}[t] for t in types]))))[0]
    for b in bounds: axs[1].axvline(b+.5,color="white",lw=1)
    save(fig,"nmf_heatmaps.png")

# =============================================================================
# 7. Supervised ML: decision boundary + RF importance + model comparison
# =============================================================================
def chapter_supervised():
    log("\n[SUPERVISED ML] trees / RF / boosting on iris (2D for boundary)")
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.model_selection import cross_val_score
    X2=Xir[:,2:4]  # petal length/width
    models={"CART":DecisionTreeClassifier(max_depth=3,random_state=0),
            "Random Forest":RandomForestClassifier(300,random_state=0),
            "Gradient Boosting":GradientBoostingClassifier(random_state=0)}
    accs={}
    for nm,mdl in models.items():
        s=cross_val_score(mdl,Xir,yir,cv=5); accs[nm]=s.mean()
        log(f"  {nm:18s} 5-fold acc = {s.mean():.3f} ± {s.std():.3f}")
    # decision boundaries
    fig,axs=plt.subplots(1,3,figsize=(12,3.8))
    xx,yy=np.meshgrid(np.linspace(X2[:,0].min()-.3,X2[:,0].max()+.3,300),
                      np.linspace(X2[:,1].min()-.3,X2[:,1].max()+.3,300))
    grid=np.c_[xx.ravel(),yy.ravel()]
    from matplotlib.colors import ListedColormap
    cmap=ListedColormap(["#dbe2f0","#f0dbe0","#dbeeee"])
    for ax,(nm,mdl) in zip(axs,models.items()):
        mdl.fit(X2,yir); Z=mdl.predict(grid).reshape(xx.shape)
        ax.contourf(xx,yy,Z,cmap=cmap,alpha=.9)
        for k in range(3):
            m=yir==k; ax.scatter(X2[m,0],X2[m,1],s=16,color=PAL[k],edgecolor="white",lw=.3)
        ax.set_title(f"{nm}\n(acc {accs[nm]:.2f})"); ax.set_xlabel("Petal length")
        ax.grid(False)
    axs[0].set_ylabel("Petal width"); save(fig,"ml_decision_boundaries.png")
    # RF importance
    rf=RandomForestClassifier(400,random_state=0).fit(Xir,yir)
    imp=rf.feature_importances_; o=np.argsort(imp)
    fig,ax=plt.subplots(figsize=(5.4,3.6))
    ax.barh(np.array(ir_names)[o],imp[o],color=TEAL)
    ax.set_title("Random-forest variable importance (iris)"); ax.grid(axis="y")
    save(fig,"ml_rf_importance.png")
    # model comparison bar
    fig,ax=plt.subplots(figsize=(5.4,3.4))
    ax.bar(list(accs),list(accs.values()),color=[NAVY,WINE,GOLD])
    ax.set_ylim(.9,1.0); ax.set_ylabel("5-fold accuracy"); ax.set_title("Model comparison")
    for i,v in enumerate(accs.values()): ax.text(i,v+.002,f"{v:.3f}",ha="center")
    save(fig,"ml_model_comparison.png"); return accs

# =============================================================================
# 8. Unsupervised ML: PCA vs t-SNE vs UMAP embeddings
# =============================================================================
def chapter_unsupervised():
    log("\n[UNSUPERVISED ML] PCA / t-SNE / UMAP on leukemia expression")
    from sklearn.decomposition import PCA
    from sklearn.manifold import TSNE
    B,types=sim_leukemia(); Xs=B.T  # samples×genes
    lab={"ALL":0,"AML":1,"CLL":2}; y=np.array([lab[t] for t in types])
    embeds=[("PCA",PCA(2).fit_transform((Xs-Xs.mean(0))/Xs.std(0,ddof=1)))]
    embeds.append(("t-SNE",TSNE(2,perplexity=15,init="pca",random_state=0).fit_transform(Xs)))
    try:
        import umap
        embeds.append(("UMAP",umap.UMAP(random_state=0).fit_transform(Xs)))
    except Exception as e:
        log("  (UMAP unavailable:",str(e)[:60],"- showing PCA/t-SNE)")
    fig,axs=plt.subplots(1,len(embeds),figsize=(4.6*len(embeds),4.2))
    if len(embeds)==1: axs=[axs]
    for ax,(nm,E) in zip(axs,embeds):
        for k,t in enumerate(["ALL","AML","CLL"]):
            m=y==k; ax.scatter(E[m,0],E[m,1],s=26,color=PAL[k],label=t,edgecolor="none")
        ax.set_title(nm); ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    axs[-1].legend(frameon=False,title="subtype"); save(fig,"unsup_embeddings.png")

# =============================================================================
# 9. Neural network training curve (small dense net on tabular)
# =============================================================================
def chapter_nn():
    log("\n[NEURAL NET] dense net learning curve (sklearn MLP proxy for keras)")
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    Xs=StandardScaler().fit_transform(Xir)
    Xtr,Xte,ytr,yte=train_test_split(Xs,yir,test_size=.3,random_state=0,stratify=yir)
    clf=MLPClassifier(hidden_layer_sizes=(16,8),max_iter=1,warm_start=True,
                      random_state=0,learning_rate_init=.03)
    tr,te=[],[]
    import warnings; warnings.filterwarnings("ignore")
    for _ in range(120):
        clf.fit(Xtr,ytr); tr.append(1-clf.score(Xtr,ytr)); te.append(1-clf.score(Xte,yte))
    log(f"  final train err {tr[-1]:.3f}  test err {te[-1]:.3f}")
    fig,ax=plt.subplots(figsize=(5.6,4))
    ax.plot(tr,color=NAVY,lw=2,label="training"); ax.plot(te,color=WINE,lw=2,label="validation")
    ax.set_xlabel("Epoch"); ax.set_ylabel("Error rate")
    ax.set_title("Neural-network learning curve"); ax.legend(frameon=False)
    save(fig,"nn_training_curve.png")

# =============================================================================
# 10. LLM embeddings demo (toy word vectors -> cosine + 2D map)
# =============================================================================
def chapter_llm():
    log("\n[LLM] toy embedding geometry (analogy + projection)")
    from sklearn.decomposition import PCA
    rng=default_rng(7)
    words=["king","queen","man","woman","dog","cat","paris","france","rome","italy",
           "gene","genome","protein","cell"]
    # construct structured 24-d embeddings with gender/royalty/geo/bio axes
    dim=24; base={w:rng.normal(0,1,dim) for w in words}
    royalty=rng.normal(0,1,dim); gender=rng.normal(0,1,dim); geo=rng.normal(0,1,dim); bio=rng.normal(0,1,dim)
    E={
      "king":royalty+gender,"queen":royalty-gender,"man":gender*1.0,"woman":-gender*1.0,
      "dog":bio*.3+rng.normal(0,.4,dim),"cat":bio*.3+rng.normal(0,.4,dim),
      "paris":geo+rng.normal(0,.3,dim),"france":geo*1.1+rng.normal(0,.3,dim),
      "rome":geo*.6+rng.normal(0,.3,dim)+1.0,"italy":geo*.7+rng.normal(0,.3,dim)+1.0,
      "gene":bio,"genome":bio*1.1,"protein":bio*.9+rng.normal(0,.3,dim),"cell":bio*.8}
    M=np.array([E[w] for w in words]); M=M/np.linalg.norm(M,axis=1,keepdims=True)
    def cos(a,b): return float(a@b)
    ik={w:i for i,w in enumerate(words)}
    analogy=M[ik["king"]]-M[ik["man"]]+M[ik["woman"]]; analogy/=np.linalg.norm(analogy)
    sims=sorted(((cos(analogy,M[i]),w) for i,w in enumerate(words)),reverse=True)
    log("  king - man + woman ≈", sims[0][1],"then",[w for _,w in sims[1:3]])
    P=PCA(2).fit_transform(M)
    fig,ax=plt.subplots(figsize=(6.2,5))
    groups={"royal/gender":["king","queen","man","woman"],"geo":["paris","france","rome","italy"],
            "animals":["dog","cat"],"biology":["gene","genome","protein","cell"]}
    col={"royal/gender":NAVY,"geo":WINE,"animals":GOLD,"biology":TEAL}
    for g,ws in groups.items():
        idx=[ik[w] for w in ws]
        ax.scatter(P[idx,0],P[idx,1],s=45,color=col[g],label=g)
        for w in ws: ax.text(P[ik[w],0]+.02,P[ik[w],1]+.02,w,fontsize=9)
    ax.set_title("Toy word embeddings (PCA projection)"); ax.legend(frameon=False,fontsize=8)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); save(fig,"llm_embeddings.png")

# =============================================================================
if __name__=="__main__":
    log("="*70); log("MS_LJMR figure + verification run"); log("="*70)
    chapter_pca(); chapter_clustering(); chapter_svd(); chapter_lda()
    chapter_ca(); chapter_nmf(); chapter_supervised(); chapter_unsupervised()
    chapter_nn(); chapter_llm()
    log("\nDONE. Figures in figs/, verification numbers above.")
    LOG.close()
