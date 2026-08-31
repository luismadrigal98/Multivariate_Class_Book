## ============================================================================
##  MS_LJMR :: 18_LLM_Embeddings.R — Large Language Models & Embeddings
##
##  Author: Luis J. Madrigal-Roca  (new session for the MS_LJMR edition)
##
##  NEW session (revised Fall-2026 schedule).
##  Reading: Vaswani et al. 2017 ("Attention is all you need").
##
##  Bridge to the course: an *embedding* is exactly the multivariate object we
##  have studied all semester — a high-dimensional vector per item. An LLM turns
##  each token/word/document into such a vector so that geometry encodes meaning
##  (cosine similarity ~ semantic similarity). We can then apply every tool from
##  this book to embeddings: distances, PCA/UMAP, clustering, classification.
##
##  This session is self-contained: instead of calling an external LLM API we
##  build small STRUCTURED embeddings (gender / royalty / geography / biology
##  axes) so the classic analogies work and reproduce offline. The exact same
##  code applies to real embeddings — just replace `E` with a matrix returned by
##  a model (e.g. text-embedding-3, or a local sentence-transformer).
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/18_LLM_Embeddings.R by scripts/make_class_src.py --
##  edit R/18_LLM_Embeddings.R and regenerate; changes made here will be overwritten.
## ============================================================================

set.seed(7)

words <- c("king","queen","man","woman","dog","cat",
           "paris","france","rome","italy","gene","genome","protein","cell")
d <- 24
royalty <- rnorm(d); gender <- rnorm(d); geo <- rnorm(d); bio <- rnorm(d)
mk <- function(v) v / sqrt(sum(v^2))
E <- rbind(
  king = royalty + gender, queen = royalty - gender,
  man = gender, woman = -gender,
  dog = 0.3*bio + rnorm(d,0,.4), cat = 0.3*bio + rnorm(d,0,.4),
  paris = geo + rnorm(d,0,.3), france = 1.1*geo + rnorm(d,0,.3),
  rome = 0.6*geo + rnorm(d,0,.3) + 1, italy = 0.7*geo + rnorm(d,0,.3) + 1,
  gene = bio, genome = 1.1*bio, protein = 0.9*bio + rnorm(d,0,.3), cell = 0.8*bio)
E <- t(apply(E, 1, mk))                          # unit-normalise each row

## ---- 1. Cosine similarity is the natural metric for embeddings -------------
cosine <- function(a, b) sum(a * b)              # rows are unit vectors
cat("cos(king, queen)  =", round(cosine(E["king",],  E["queen",]),  3), "\n")
cat("cos(king, cell)   =", round(cosine(E["king",],  E["cell",]),   3), "\n")
cat("cos(gene, genome) =", round(cosine(E["gene",],  E["genome",]), 3), "\n")

## ---- 2. Vector analogy: king - man + woman ~ queen -------------------------
q <- E["king",] - E["man",] + E["woman",]; q <- q / sqrt(sum(q^2))
sims <- sort(E %*% q, decreasing = TRUE)
cat("\nking - man + woman  ->  nearest words:\n")
print(round(head(sims[!rownames(E) %in% c("king","man","woman")], 3), 3))

## ---- 3. Every multivariate tool applies to embeddings ----------------------
## (a) PCA projection to 2-D for visualisation
pc <- prcomp(E)
grp <- c(king="royal",queen="royal",man="royal",woman="royal",
         dog="animal",cat="animal",paris="geo",france="geo",rome="geo",
         italy="geo",gene="bio",genome="bio",protein="bio",cell="bio")
pal <- c(royal="#1f3b73", geo="#8c2d3a", animal="#c8a23a", bio="#2a7f7f")
par(mfrow = c(1, 1))
plot(pc$x[,1:2], pch = 19, col = pal[grp[rownames(E)]], cex = 1.4,
     xlab = "PC1", ylab = "PC2", main = "Word embeddings — PCA projection")
text(pc$x[,1], pc$x[,2], rownames(E), pos = 3, cex = 0.75)
legend("topright", names(pal), col = pal, pch = 19, bty = "n")

## (b) hierarchical clustering on cosine distance recovers the themes
Dcos <- as.dist(1 - E %*% t(E))
cat("\nHierarchical clusters of words (cosine distance):\n")
print(cutree(hclust(Dcos, "average"), k = 4))

cat("\n[18_LLM_Embeddings] done. Embeddings = multivariate vectors; reuse the",
    "whole toolbox.\n")
