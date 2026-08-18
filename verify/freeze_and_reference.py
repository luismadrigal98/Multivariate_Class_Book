#!/usr/bin/env python3
"""
MS_LJMR :: verify/freeze_and_reference.py
=========================================
Single source of truth for the DEMO data and the reference OUTPUTS quoted in the
book. It (1) writes frozen demo datasets to ../data/demo/ so the R scripts and
this verifier read *identical bytes* (deterministic methods then match exactly),
and (2) recomputes every numeric result the book prints, in R-like formatting,
into reference_outputs.txt.  It also (re)renders the figures that depend on the
simulated data so figures and quoted numbers agree.

Run:  python3 freeze_and_reference.py
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from numpy.random import default_rng

HERE=os.path.dirname(os.path.abspath(__file__))
DEMO=os.path.normpath(os.path.join(HERE,"..","data","demo")); os.makedirs(DEMO,exist_ok=True)
FIGS=os.path.normpath(os.path.join(HERE,"..","figs")); os.makedirs(FIGS,exist_ok=True)
REF=open(os.path.join(HERE,"reference_outputs.txt"),"w")
def out(*a):
    print(*a); print(*a,file=REF)
def head(t): out("\n"+"="*68); out(t); out("="*68)

NAVY,WINE,GOLD,TEAL="#1f3b73","#8c2d3a","#c8a23a","#2a7f7f"; PAL=[NAVY,WINE,TEAL,GOLD]
plt.rcParams.update({"figure.dpi":150,"savefig.dpi":150,"font.size":11,
    "axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,"grid.alpha":.25})
def save(fig,n): fig.savefig(os.path.join(FIGS,n),bbox_inches="tight"); plt.close(fig); out("  fig ->",n)

# ---------------------------------------------------------------- freeze data
def freeze():
    head("FREEZING DEMO DATA -> data/demo/  (R and Python read these same files)")
    # taxon: 4 groups x 30, 7 vars
    rng=default_rng(1998); g,per,p=4,30,7
    centers=rng.normal(0,3,(g,p))
    X=np.vstack([rng.normal(0,1,(per,p))+centers[k] for k in range(g)])
    vn=["Petals","Internode","Sepal","Bract","Petiole","Leaf","Fruit"]
    tx=pd.DataFrame(X,columns=vn); tx.insert(0,"Taxon",np.repeat(["I","II","III","IV"],per))
    tx.to_csv(os.path.join(DEMO,"taxon.csv"),index=False)
    # doubs fish + env
    rng=default_rng(1998); n,s=30,27; grad=np.linspace(0,1,n); opt=np.linspace(0,1,s)
    fish=np.column_stack([rng.poisson(8*np.exp(-((grad-opt[k])**2)/0.02)) for k in range(s)])
    pd.DataFrame(fish,columns=[f"Sp{i+1}" for i in range(s)]).to_csv(os.path.join(DEMO,"doubs_fish.csv"),index=False)
    env=pd.DataFrame(dict(alt=900-850*grad+rng.normal(0,10,n),slo=np.exp(4-4*grad+rng.normal(0,.3,n)),
        flo=10+400*grad+rng.normal(0,20,n),pH=8-grad+rng.normal(0,.1,n),
        oxy=12-4*grad+rng.normal(0,.5,n),nit=1+3*grad+rng.normal(0,.2,n)))
    env.to_csv(os.path.join(DEMO,"doubs_env.csv"),index=False)
    # leukemia expression genes x samples
    rng=default_rng(1998); types=np.repeat(["ALL","AML","CLL"],20); gg=300
    B=rng.normal(6,1,(gg,len(types)))
    for t in np.unique(types):
        idx=np.where(types==t)[0]; sig=rng.choice(gg,40,replace=False); B[np.ix_(sig,idx)]+=rng.normal(3,0.3)
    cols=[f"{t}.{i}" for t,i in zip(types,[list(range(1,21))*3][0])]
    pd.DataFrame(B,columns=cols,index=[f"g{i}" for i in range(gg)]).to_csv(os.path.join(DEMO,"leukemia.csv"))
    # limenitis image
    x=np.linspace(-3,3,100); y=np.linspace(-2,2,69); xx,yy=np.meshgrid(x,y)
    G=np.exp(-(xx**2)/3-(yy**2)/1.2)*(1+0.6*np.cos(3*xx))+0.15*np.sin(2*xx)*np.cos(2*yy)
    G=G+default_rng(1998).normal(0,0.06,G.shape); G=(G-G.min())/(G.max()-G.min())
    pd.DataFrame(G).to_csv(os.path.join(DEMO,"limenitis.csv"),index=False,header=False)
    # crawley mixed (Biomass right-skewed positive -> Gamma; Species counts -> Poisson)
    rng=default_rng(1998); nc=90
    Tmp=rng.uniform(-2,12,nc); Precip=rng.gamma(2,60,nc); pH=rng.uniform(3.5,8.5,nc)
    Soil=rng.choice(["clay","loam","sand"],nc)
    mu_b=np.exp(1.6+0.04*Tmp-0.06*(pH-6)**2)       # positive mean surface
    Biomass=rng.gamma(shape=3.0,scale=mu_b/3.0)     # Gamma -> right-skewed, non-normal
    Species=rng.poisson(np.exp(1.2+.06*Tmp+.15*(pH-4)))
    pd.DataFrame(dict(Species=Species,Biomass=Biomass,Tmp=Tmp,Precip=Precip,pH=pH,Soil=Soil)).to_csv(
        os.path.join(DEMO,"crawley.csv"),index=False)
    # hanta (ENM): presence/absence with 2 bioclim predictors
    rng=default_rng(1998); nh=400
    bio1=rng.uniform(-5,30,nh); bio12=rng.uniform(100,2500,nh)
    eta=-2.0+0.15*bio1-0.0015*bio12+0.00000 -0.004*(bio1-15)**2
    p=1/(1+np.exp(-eta)); Sp=rng.binomial(1,p)
    pd.DataFrame(dict(Sp=Sp,bio_1=bio1,bio_12=bio12)).to_csv(os.path.join(DEMO,"hanta.csv"),index=False)
    # PAM: sites x species presence/absence on a latitudinal gradient
    rng=default_rng(1998); nsite,nsp=80,120
    lat=np.linspace(0,1,nsite); optS=rng.uniform(0,1,nsp); breadth=rng.uniform(.05,.25,nsp)
    prob=np.exp(-((lat[:,None]-optS[None])**2)/(2*breadth[None]**2))
    pam=(rng.random((nsite,nsp))<prob).astype(int)
    dfp=pd.DataFrame(pam,columns=[f"sp{i+1}" for i in range(nsp)]); dfp.insert(0,"lat",lat)
    dfp.to_csv(os.path.join(DEMO,"pam.csv"),index=False)
    out("  wrote: taxon, doubs_fish, doubs_env, leukemia, limenitis, crawley, hanta, pam")
    return dict(taxon=tx,fish=fish,env=env,leuk=(B,types),lim=G,craw=None,
                hanta=pd.read_csv(os.path.join(DEMO,"hanta.csv")),pam=(pam,lat))

D=freeze()

# =============================================================== iris built-ins
from sklearn.datasets import load_iris
IR=load_iris(); Xir=IR.data; yir=IR.target
irn=["Sepal.Length","Sepal.Width","Petal.Length","Petal.Width"]; spn=list(IR.target_names)

head("CH.1  iris R-mode correlation matrix (deterministic; built-in data)")
R=np.corrcoef(((Xir-Xir.mean(0))/Xir.std(0,ddof=1)).T)
out("             "+ "".join(f"{n:>13}" for n in irn))
for i,n in enumerate(irn): out(f"{n:<13}"+"".join(f"{R[i,j]:>13.2f}" for j in range(4)))

head("CH.2  SVD low-rank reconstruction error (limenitis frozen image)")
G=D["lim"]; U,S,Vt=np.linalg.svd(G,full_matrices=False)
errs=[np.linalg.norm(G-(U[:,:k]*S[:k])@Vt[:k])/np.linalg.norm(G) for k in (1,5,10,50)]
out("  rank 1,5,10,50 rel.Frobenius error:", " ".join(f"{e:.3f}" for e in errs))

head("CH.3  iris PCA eigenvalues / variance (built-in)")
Xs=(Xir-Xir.mean(0))/Xir.std(0,ddof=1); ev=np.linalg.svd(Xs,full_matrices=False)[1]**2/(len(Xs)-1)
out("  eigenvalue :", " ".join(f"{v:.6f}" for v in ev))
out("  prop.var   :", " ".join(f"{v:.6f}" for v in ev/ev.sum()))
out("  cumulative :", " ".join(f"{v:.6f}" for v in np.cumsum(ev/ev.sum())))

head("CH.4  Correspondence analysis inertia (doubs_fish frozen)")
fish=D["fish"]; P=fish/fish.sum(); r=P.sum(1,keepdims=True); c=P.sum(0,keepdims=True)
Sca=(P-r@c)/np.sqrt(r@c); sv=np.linalg.svd(Sca,full_matrices=False)[1]; inr=sv**2
out("  first 4 inertias :", " ".join(f"{v:.4f}" for v in inr[:4]))
out("  prop. of total   :", " ".join(f"{v:.3f}" for v in (inr/inr.sum())[:4]))
grad=np.linspace(0,1,fish.shape[0])
U2,s2,V2t=np.linalg.svd(Sca,full_matrices=False)
rf=(U2[:,:2]*s2[:2])/np.sqrt(r); cf=(V2t[:2].T*s2[:2])/np.sqrt(c.T)
fig,ax=plt.subplots(figsize=(6,5))
sc=ax.scatter(rf[:,0],rf[:,1],c=grad,cmap="viridis",s=40)
ax.scatter(cf[:,0],cf[:,1],marker="^",color=WINE,s=30,alpha=.7,label="species")
ax.axhline(0,color="grey",lw=.6); ax.axvline(0,color="grey",lw=.6)
ax.set_xlabel("CA axis 1"); ax.set_ylabel("CA axis 2"); ax.set_title("Correspondence analysis — river gradient")
plt.colorbar(sc,label="upstream → downstream"); ax.legend(frameon=False); save(fig,"ca_biplot.png")

head("CH.6  NMF reconstruction (leukemia frozen)")
from sklearn.decomposition import NMF
B,types=D["leuk"]; Bn=B-B.min()
nmf=NMF(3,init="nndsvda",random_state=0,max_iter=500); W=nmf.fit_transform(Bn); Hh=nmf.components_
out(f"  rank-3 reconstruction error = {nmf.reconstruction_err_:.3f}   W {W.shape}  H {Hh.shape}")
order=np.argsort([{'ALL':0,'AML':1,'CLL':2}[t] for t in types])
fig,axs=plt.subplots(1,2,figsize=(10,4),gridspec_kw=dict(width_ratios=[1,2.4]))
axs[0].imshow(W[np.argsort(W.argmax(1))],aspect="auto",cmap="magma"); axs[0].set_title("Basis W (genes×3)")
axs[1].imshow(Hh[:,order],aspect="auto",cmap="magma"); axs[1].set_title("Coefficients H (3×samples)")
for b in [20,40]: axs[1].axvline(b-.5,color="white",lw=1)
save(fig,"nmf_heatmaps.png")

head("CH.7  iris k-means (k=3) cross-tab + agreement (built-in)")
from sklearn.cluster import KMeans
km=KMeans(3,n_init=25,random_state=0).fit(Xir); ct=pd.crosstab(km.labels_,np.array(spn)[yir])
out(ct.to_string()); out("  agreement (max-assignment) =",round(ct.values.max(1).sum()/150,3))

head("CH.8  LDA hold-out accuracy on frozen taxon")
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.model_selection import train_test_split
tx=D["taxon"]; Xt=tx.iloc[:,1:].values; yt=tx["Taxon"].values
Xtr,Xte,ytr,yte=train_test_split(Xt,yt,test_size=1/3,random_state=123,stratify=yt)
lda=LDA().fit(Xtr,ytr); out(f"  LDA hold-out accuracy = {lda.score(Xte,yte):.3f}")
from sklearn.decomposition import PCA
Z=LDA().fit(Xt,yt).transform(Xt); Pp=PCA(2).fit_transform((Xt-Xt.mean(0))/Xt.std(0,ddof=1))
fig,axs=plt.subplots(1,2,figsize=(9.5,4.2))
for k,lev in enumerate(np.unique(yt)):
    m=yt==lev; axs[0].scatter(Pp[m,0],Pp[m,1],s=22,color=PAL[k],label=lev,edgecolor="none")
    axs[1].scatter(Z[m,0],Z[m,1],s=22,color=PAL[k],label=lev,edgecolor="none")
axs[0].set_title("PCA (unsupervised)"); axs[1].set_title("LDA (supervised)")
for a in axs: a.set_xlabel("Axis 1"); a.set_ylabel("Axis 2")
axs[1].legend(frameon=False,title="Taxon"); save(fig,"lda_vs_pca.png")

head("CH.9  supervised ML 5-fold accuracy (iris, built-in)")
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.model_selection import cross_val_score
for nm,mdl in [("CART (rpart)",DecisionTreeClassifier(max_depth=3,random_state=0)),
               ("Random Forest",RandomForestClassifier(500,random_state=0)),
               ("Gradient Boosting",GradientBoostingClassifier(random_state=0))]:
    s=cross_val_score(mdl,Xir,yir,cv=5); out(f"  {nm:18s} {s.mean():.3f} ± {s.std():.3f}")

head("CH.10  embedding neighbourhood purity (leukemia frozen)")
from sklearn.manifold import TSNE
Xsz=(B.T-B.T.mean(0))/B.T.std(0,ddof=1); ylab=np.array([{'ALL':0,'AML':1,'CLL':2}[t] for t in types])
def purity(Y,lab,k=10):
    Dm=np.sqrt(((Y[:,None]-Y[None])**2).sum(-1)); s=0
    for i in range(len(Y)): s+=np.mean(lab[np.argsort(Dm[i])[1:k+1]]==lab[i])
    return s/len(Y)
Ppca=PCA(2).fit_transform(Xsz); Pts=TSNE(2,perplexity=15,init="pca",random_state=0).fit_transform(B.T)
try:
    import umap; Pum=umap.UMAP(random_state=0).fit_transform(B.T); um_ok=True
except Exception: um_ok=False
out(f"  PCA purity {purity(Ppca,ylab):.3f}   t-SNE {purity(Pts,ylab):.3f}"+(f"   UMAP {purity(Pum,ylab):.3f}" if um_ok else ""))

head("CH.11  nnet-style test accuracy (iris, built-in)")
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
import warnings; warnings.filterwarnings("ignore")
Xz=StandardScaler().fit_transform(Xir)
Xtr,Xte,ytr,yte=train_test_split(Xz,yir,test_size=.3,random_state=0,stratify=yir)
clf=MLPClassifier((8,),alpha=5e-3,max_iter=1500,random_state=0).fit(Xtr,ytr)
out(f"  test accuracy = {clf.score(Xte,yte):.3f}")

head("CH.ENM  logistic ENM coefficients + AUC (hanta frozen)")
import statsmodels.api as sm
from sklearn.metrics import roc_auc_score
h=D["hanta"]; Xd=sm.add_constant(h[["bio_1","bio_12"]]); m=sm.Logit(h["Sp"],Xd).fit(disp=0)
out("  coefficients:");
for nm,co in m.params.items(): out(f"    {nm:8s} {co:+.5f}")
out(f"  training AUC = {roc_auc_score(h['Sp'],m.predict(Xd)):.3f}")
# ENM response curve + suitability map
b1=np.linspace(-5,30,100); base=np.median(h["bio_12"])
resp=m.predict(sm.add_constant(pd.DataFrame({"bio_1":b1,"bio_12":base}),has_constant="add"))
fig,axs=plt.subplots(1,2,figsize=(10,4))
axs[0].plot(b1,resp,color=NAVY,lw=2); axs[0].set_xlabel("bio_1 (temperature)")
axs[0].set_ylabel("P(presence)"); axs[0].set_title("ENM response curve")
gx,gy=np.meshgrid(np.linspace(-5,30,120),np.linspace(100,2500,120))
gg=m.predict(sm.add_constant(pd.DataFrame({"bio_1":gx.ravel(),"bio_12":gy.ravel()}),has_constant="add")).values.reshape(gx.shape)
im=axs[1].contourf(gx,gy,gg,levels=12,cmap="YlGnBu"); axs[1].scatter(h.loc[h.Sp==1,"bio_1"],h.loc[h.Sp==1,"bio_12"],s=6,color=WINE)
axs[1].set_xlabel("bio_1"); axs[1].set_ylabel("bio_12"); axs[1].set_title("Suitability surface"); axs[1].grid(False)
plt.colorbar(im,ax=axs[1],label="suitability"); save(fig,"enm.png")

head("CH.PAM  presence-absence matrix summaries (pam frozen)")
pam,lat=D["pam"]
alpha=pam.sum(1); omega=pam.sum(0)
out(f"  sites={pam.shape[0]}  species={pam.shape[1]}  fill={pam.mean():.3f}")
out(f"  mean alpha (richness/site) = {alpha.mean():.2f}   mean omega (range/species) = {omega.mean():.2f}")
fig,axs=plt.subplots(1,3,figsize=(13,3.8))
axs[0].imshow(pam[np.argsort(lat)],aspect="auto",cmap="Greys"); axs[0].set_title("PAM (sites×species)"); axs[0].grid(False)
axs[1].plot(lat,alpha,color=NAVY); axs[1].set_xlabel("latitude"); axs[1].set_ylabel("richness α"); axs[1].set_title("Richness gradient")
axs[2].hist(omega,bins=20,color=TEAL); axs[2].set_xlabel("range size ω"); axs[2].set_title("Range-size frequency")
save(fig,"pam.png")

head("CH.FA  factor analysis (frozen taxon, ML, varimax) — uniquenesses/loadings")
from sklearn.decomposition import FactorAnalysis
Xz2=(Xt-Xt.mean(0))/Xt.std(0,ddof=1)
fa=FactorAnalysis(n_components=2,rotation="varimax",random_state=0).fit(Xz2)
L=fa.components_.T
out("  loadings (varimax, 2 factors):")
for i,nm in enumerate(tx.columns[1:]): out(f"    {nm:10s} {L[i,0]:+.2f} {L[i,1]:+.2f}")
ev2=np.linalg.eigvalsh(np.corrcoef(Xz2.T))[::-1]
out("  correlation-matrix eigenvalues:", " ".join(f"{v:.3f}" for v in ev2))
fig,ax=plt.subplots(figsize=(5.2,4))
ax.plot(range(1,len(ev2)+1),ev2,"o-",color=NAVY,lw=2); ax.axhline(1,color=WINE,ls="--")
ax.set_xlabel("factor"); ax.set_ylabel("eigenvalue"); ax.set_title("Factor analysis — scree (taxon)")
save(fig,"factor_scree.png")

out("\nDONE. Frozen data in data/demo/, reference outputs above, figures in figs/.")
REF.close()
