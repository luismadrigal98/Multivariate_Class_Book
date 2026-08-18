#!/usr/bin/env python3
"""Extra figures: dissimilarity heatmap (Ch.1) and MDS + Shepard (Ch. MDS)."""
import os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from numpy.random import default_rng
from sklearn.datasets import load_iris
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import MDS

HERE=os.path.dirname(os.path.abspath(__file__)); FIGS=os.path.normpath(os.path.join(HERE,"..","figs"))
NAVY,WINE,GOLD,TEAL="#1f3b73","#8c2d3a","#c8a23a","#2a7f7f"; PAL=[NAVY,WINE,TEAL,GOLD]
plt.rcParams.update({"figure.dpi":150,"savefig.dpi":150,"font.size":11,
    "axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,"grid.alpha":.25})
def save(fig,n): fig.savefig(os.path.join(FIGS,n),bbox_inches="tight"); plt.close(fig); print("figs/"+n)

IR=load_iris(); X=IR.data; y=IR.target

# ---- 1. Q-mode dissimilarity heatmap (Euclidean, standardized iris) ----------
Xs=(X-X.mean(0))/X.std(0,ddof=1)
D=squareform(pdist(Xs,"euclidean"))
o=np.argsort(y)                                   # order by species -> block structure
fig,ax=plt.subplots(figsize=(5.4,4.6))
im=ax.imshow(D[np.ix_(o,o)],cmap="magma")
for b in [50,100]: ax.axhline(b-.5,color="w",lw=.8); ax.axvline(b-.5,color="w",lw=.8)
ax.set_title("Q-mode Euclidean dissimilarity (iris, ordered by species)")
ax.set_xlabel("observation"); ax.set_ylabel("observation"); ax.grid(False)
plt.colorbar(im,label="distance"); save(fig,"diss_heatmap.png")

# ---- 2. R-mode correlation heatmap ------------------------------------------
R=np.corrcoef(Xs.T); names=["Sep.L","Sep.W","Pet.L","Pet.W"]
fig,ax=plt.subplots(figsize=(4.2,3.8))
im=ax.imshow(R,cmap="RdBu_r",vmin=-1,vmax=1)
ax.set_xticks(range(4)); ax.set_xticklabels(names); ax.set_yticks(range(4)); ax.set_yticklabels(names)
for i in range(4):
    for j in range(4): ax.text(j,i,f"{R[i,j]:.2f}",ha="center",va="center",
        color="white" if abs(R[i,j])>.6 else "black",fontsize=9)
ax.set_title("R-mode correlations"); ax.grid(False); plt.colorbar(im,shrink=.8); save(fig,"corr_heatmap.png")

# ---- 3. Non-metric MDS of a seeded multivariate sample + Shepard diagram -----
# (offline stand-in for the mtcars example used in the R script: 32 objects,
#  a few correlated latent gradients + noise, standardized)
rng=default_rng(1998); n=32
g1=rng.normal(0,1,n); g2=rng.normal(0,1,n)
Z=np.column_stack([g1, g1+0.3*rng.normal(0,1,n), -g1+0.4*rng.normal(0,1,n),
                   g2, g2+0.3*rng.normal(0,1,n), 0.6*g1-0.6*g2+0.3*rng.normal(0,1,n)])
Z=(Z-Z.mean(0))/Z.std(0,ddof=1)
labels=[f"obj{i+1:02d}" for i in range(n)]
Dm=squareform(pdist(Z,"euclidean"))
mds=MDS(n_components=2,dissimilarity="precomputed",random_state=1,n_init=4,
        normalized_stress="auto")
Y=mds.fit_transform(Dm); stress=mds.stress_
print("MDS raw stress:",round(stress,2))
fig,axs=plt.subplots(1,2,figsize=(10,4.4))
axs[0].scatter(Y[:,0],Y[:,1],s=28,color=NAVY)
for i,nm in enumerate(labels): axs[0].text(Y[i,0],Y[i,1],nm,fontsize=6,alpha=.8)
axs[0].set_title("nMDS configuration"); axs[0].set_xlabel("dim 1"); axs[0].set_ylabel("dim 2")
# Shepard: original distance vs configuration distance
dcfg=squareform(pdist(Y)); iu=np.triu_indices_from(Dm,1)
axs[1].scatter(Dm[iu],dcfg[iu],s=8,color=TEAL,alpha=.5)
lim=[0,Dm[iu].max()*1.05]; axs[1].plot(lim,lim,color=WINE,lw=1.5)
axs[1].set_xlabel("original distance"); axs[1].set_ylabel("configuration distance")
axs[1].set_title("Shepard diagram"); save(fig,"mds_shepard.png")
print("done")
