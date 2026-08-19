## ============================================================================
##  MS_LJMR :: 11_FactorAnalysis.R — Exploratory Factor Analysis
##
##  Original authors: Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from FactorAnalysis.R. Uses frozen taxon data (portable).
##
##  FA vs PCA: PCA summarises TOTAL variance with components that are exact
##  linear combinations of the variables. FA models the COMMON variance with a
##  small number of latent factors plus variable-specific 'uniquenesses':
##      x = Lambda f + e,   Cov(x) = Lambda Lambda' + Psi.
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
##  Files to keep next to this script: taxon.csv
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) read.csv("taxon.csv", stringsAsFactors = TRUE)
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

taxon <- get_data("taxon")
X     <- scale(taxon[, -1])

## ---- how many factors? eigenvalues of the correlation matrix ---------------
ev <- eigen(cor(X))$values
cat("Correlation-matrix eigenvalues:\n"); print(round(ev, 3))
with_fig("11_fa_scree", {
  plot(ev, type = "b", pch = 19, col = "#1f3b73",
       xlab = "factor", ylab = "eigenvalue", main = "Factor analysis scree")
  abline(h = 1, col = "#8c2d3a", lty = 2)          # Kaiser criterion
})

## ---- maximum-likelihood factor analysis with varimax rotation --------------
fa <- factanal(X, factors = 2, rotation = "varimax", scores = "regression")
cat("\nFactor loadings (varimax):\n")
print(round(unclass(fa$loadings), 2))
cat("\nUniquenesses (variable-specific variance):\n")
print(round(fa$uniquenesses, 2))
cat("\nLikelihood-ratio test that 2 factors suffice: p =",
    round(fa$PVAL, 4), "\n")

with_fig("11_fa_loadings", {
  L <- fa$loadings
  plot(L[, 1], L[, 2], pch = 19, col = "#1f3b73", xlim = c(-1, 1), ylim = c(-1, 1),
       xlab = "Factor 1", ylab = "Factor 2", main = "Rotated loadings")
  abline(h = 0, v = 0, col = "grey70")
  text(L[, 1], L[, 2], rownames(L), pos = 3, cex = .8)
})

cat("\n[11_FactorAnalysis] FA models shared variance via latent factors +",
    "uniquenesses; rotation (varimax) aids interpretation.\n")
