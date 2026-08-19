## ============================================================================
##  MS_LJMR :: 09_PCA.R  — Principal Component Analysis (I–III, PCoA folded in)
##
##  Original authors: Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from Soberon & Jimenez 09_PCAI / 10_PCAII / 11_PCAIII_Biplots and
##  12_PCOAI, per the Fall-2026 schedule (PCoA now folded into the PCA arc).
##
##  What changed vs. the originals
##  ------------------------------
##  * No setwd(); data via get_data(); no x11() (uses open_dev()).
##  * princomp / prcomp / SVD are shown to compute the SAME object; the
##    equivalence "PCA = spectral decomposition of the correlation matrix" is
##    demonstrated numerically.
##  * A clean biplot, a scree plot with the Kaiser rule, and a PCoA section that
##    shows metric MDS on a Euclidean distance reproduces PCA scores.
##
##  Run:  source("R/00_utils.R"); source("R/09_PCA.R")
## ============================================================================
## ---------------------------------------------------------------------------
##  STANDALONE USE (no repository needed)
##  ---------------------------------------------------------------------------
##  As shipped, the source() line below borrows three helpers from the course
##  repository: get_data() (loads a data set), need() (loads packages) and
##  with_fig() (opens a plot device). To run this script entirely on its own,
##  delete that line and uncomment the block below. Nothing else changes.
##
##  The standalone with_fig() just draws each figure to the screen, one after
##  the other. To save them as files instead, replace its body with
##      png(paste0(name, ".png")); on.exit(dev.off()); force(expr)
##
##  Files to keep next to this script: none — the data sets used here ship with R itself
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) { utils::data("iris");   iris }
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("MASS")                                   # base + MASS are enough here

## ---- 1. Data ---------------------------------------------------------------
iris <- get_data("iris")
X    <- as.matrix(iris[, 1:4])
sp   <- iris$Species

## ---- 2. Standardize and build the correlation matrix -----------------------
Xs <- scale(X)                                 # z-scores: mean 0, sd 1
R  <- cor(X)                                   # == crossprod(Xs)/(n-1)
cat("Correlation matrix:\n"); print(round(R, 3))

## ---- 3. PCA three equivalent ways ------------------------------------------
pc_princomp <- princomp(Xs)                    # spectral decomposition of cov(Xs)
pc_prcomp   <- prcomp(X, center = TRUE, scale. = TRUE)   # SVD of the data
eig         <- eigen(R)                         # do it "by hand"

## eigenvalues == variances of the components == sdev^2
lambda <- eig$values
cat("\nEigenvalues (variance per PC):\n"); print(round(lambda, 4))
cat("Proportion of variance:\n");          print(round(lambda / sum(lambda), 4))
cat("prcomp sdev^2 (identical):\n");       print(round(pc_prcomp$sdev^2, 4))

## the loadings are the eigenvectors; scores are Xs %*% loadings
loadings <- eig$vectors
scores   <- Xs %*% loadings
## sign is arbitrary; align with prcomp for reproducibility
flip     <- sign(diag(crossprod(scores, pc_prcomp$x)))
scores   <- sweep(scores, 2, flip, `*`)
loadings <- sweep(loadings, 2, flip, `*`)
cat("\nLoadings (eigenvectors):\n")
dimnames(loadings) <- list(colnames(X), paste0("PC", 1:4)); print(round(loadings, 3))

## numerical proof that hand-PCA == prcomp scores
cat("\nmax|hand scores - prcomp scores| =",
    format(max(abs(scores - pc_prcomp$x)), digits = 3), "\n")

## ---- 4. Scree plot (Kaiser rule) -------------------------------------------
with_fig("09_pca_scree", {
  plot(lambda, type = "b", pch = 19, col = "#1f3b73", lwd = 2,
       xlab = "Principal component", ylab = "Eigenvalue",
       main = "Iris PCA — scree plot")
  abline(h = 1, col = "#8c2d3a", lty = 2, lwd = 1.5)   # keep PCs with eig > 1
  legend("topright", "Kaiser (eig = 1)", lty = 2, col = "#8c2d3a", bty = "n")
})

## ---- 5. Biplot -------------------------------------------------------------
cols <- c("#1f3b73", "#8c2d3a", "#2a7f7f")
with_fig("09_pca_biplot", {
  plot(scores[, 1:2], col = cols[sp], pch = 19,
       xlab = sprintf("PC1 (%.1f%%)", 100 * lambda[1] / sum(lambda)),
       ylab = sprintf("PC2 (%.1f%%)", 100 * lambda[2] / sum(lambda)),
       main = "Iris PCA — biplot", asp = 1)
  abline(h = 0, v = 0, col = "grey70")
  s <- 2.6
  arrows(0, 0, loadings[, 1] * s, loadings[, 2] * s,
         col = "#8c2d3a", length = 0.08, lwd = 2)
  text(loadings[, 1] * s * 1.1, loadings[, 2] * s * 1.1,
       rownames(loadings), col = "#8c2d3a", cex = 0.8)
  legend("bottomright", levels(sp), col = cols, pch = 19, bty = "n")
})

## ---- 6. PCoA folded in: metric MDS of Euclidean distance == PCA ------------
## Principal Coordinates Analysis on a Euclidean distance recovers PCA scores
## (Gower 1966). This is why the standalone PCoA session is folded in here.
D      <- dist(Xs)                              # Euclidean distance
pcoa   <- cmdscale(D, k = 2, eig = TRUE)        # classical (metric) MDS
## align signs with PCA for comparison
al     <- sign(diag(crossprod(pcoa$points, scores[, 1:2])))
pcoa_pts <- sweep(pcoa$points, 2, al, `*`)
cat("\nPCoA vs PCA — correlation of axis 1:",
    round(cor(pcoa_pts[, 1], scores[, 1]), 4),
    " axis 2:", round(cor(pcoa_pts[, 2], scores[, 2]), 4), "\n")

with_fig("09_pcoa_vs_pca", {
  plot(pcoa_pts, col = cols[sp], pch = 19, asp = 1,
       xlab = "PCoA axis 1", ylab = "PCoA axis 2",
       main = "PCoA on Euclidean distance = PCA")
  abline(h = 0, v = 0, col = "grey70")
})

cat("\n[09_PCA] done. Figures written to figs/ when run in batch.\n")
