## ============================================================================
##  MS_LJMR :: 26_UnsupervisedML.R — Unsupervised Machine Learning
##
##  Author: Luis J. Madrigal-Roca  (new session for the MS_LJMR edition)
##
##  NEW session (revised Fall-2026 schedule).
##  Readings: ISLR Ch. 12; McInnes et al. 2018 (UMAP); van der Maaten & Hinton
##  2008 (t-SNE).
##
##  Idea: PCA is a *linear* map that preserves global variance. Modern manifold
##  learners (t-SNE, UMAP) are *non-linear* and preserve local neighbourhoods,
##  so they separate clusters that PCA leaves entangled. We show all three on a
##  gene-expression matrix (leukemia subtypes).
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
##  Files to keep next to this script: leukemiaExpressionSubset.rds
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) readRDS("leukemiaExpressionSubset.rds")
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("Rtsne", "uwot")                            # uwot = UMAP in pure R
set.seed(1998)

## ---- Data: genes x samples -> samples x genes ------------------------------
expr  <- get_data("leukemia")                    # matrix, rows = genes
X     <- t(as.matrix(expr))                      # rows = samples
type  <- factor(sub("\\..*$", "", rownames(X)))  # ALL / AML / CLL
cols  <- c("#1f3b73", "#8c2d3a", "#2a7f7f")[type]

## ---- 1. PCA (linear baseline) ----------------------------------------------
pc <- prcomp(X, center = TRUE, scale. = TRUE)
ve <- pc$sdev^2 / sum(pc$sdev^2)
cat("PCA variance explained (PC1-3):", round(ve[1:3], 3), "\n")

## ---- 2. t-SNE (local, probabilistic) ---------------------------------------
ts <- Rtsne::Rtsne(X, dims = 2, perplexity = 15, pca = TRUE, verbose = FALSE)

## ---- 3. UMAP (local, graph-based; usually faster & keeps more structure) ---
um <- uwot::umap(X, n_neighbors = 15, min_dist = 0.1, n_components = 2)

## ---- Compare the three embeddings ------------------------------------------
with_fig("26_embeddings", {
  op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1)); on.exit(par(op))
  plot(pc$x[, 1:2], col = cols, pch = 19, main = "PCA (linear)",
       xlab = sprintf("PC1 (%.0f%%)", 100*ve[1]), ylab = "PC2")
  plot(ts$Y, col = cols, pch = 19, main = "t-SNE", xlab = "dim 1", ylab = "dim 2")
  plot(um,   col = cols, pch = 19, main = "UMAP",  xlab = "dim 1", ylab = "dim 2")
  legend("topright", levels(type), col = c("#1f3b73","#8c2d3a","#2a7f7f"),
         pch = 19, bty = "n")
})

## ---- Quick quantitative check: neighbourhood purity ------------------------
purity <- function(Y, lab, k = 10) {
  D <- as.matrix(dist(Y)); n <- nrow(Y); s <- 0
  for (i in 1:n) {
    nn <- order(D[i, ])[2:(k + 1)]
    s <- s + mean(lab[nn] == lab[i])
  }
  s / n
}
cat("k-NN label purity  PCA:", round(purity(pc$x[,1:2], type), 3),
    " t-SNE:", round(purity(ts$Y, type), 3),
    " UMAP:", round(purity(um, type), 3), "\n")

cat("\n[26_UnsupervisedML] done. Non-linear maps usually raise local purity.\n")
