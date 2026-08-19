## ============================================================================
##  MS_LJMR :: 17_NMF.R — Non-negative Matrix Factorization
##
##  Original author: Jorge Soberón
##                   (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 17_NNF.R. Uses frozen leukemia expression matrix (no GitHub
##  download, no absolute paths).
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
need("NMF")

expr <- get_data("leukemia")             # genes x samples (>= 0)
A    <- expr - min(expr)                  # ensure non-negativity
type <- sub("\\..*$", "", colnames(expr)) # ALL / AML / CLL

## A ~ W H  with W (genes x rank) and H (rank x samples), both >= 0
set.seed(1998)
res <- NMF::nmf(A, rank = 3, seed = "nndsvd")
W   <- NMF::basis(res); H <- NMF::coef(res)
cat("NMF rank 3:  W", dim(W)[1], "x", dim(W)[2],
    "  H", dim(H)[1], "x", dim(H)[2],
    "  residual", round(sqrt(sum((A - W %*% H)^2)), 1), "\n")

with_fig("17_nmf", width = 10, height = 4, {
  op <- par(mfrow = c(1, 2)); on.exit(par(op))
  image(t(W[order(max.col(W)), ]), col = hcl.colors(30, "Inferno"),
        main = "Basis W (gene signatures)", axes = FALSE)
  image(H[, order(type)], col = hcl.colors(30, "Inferno"),
        main = "Coefficients H (samples)", axes = FALSE)
})

## the three factors should line up with the three clinical subtypes
assign3 <- max.col(t(H))
cat("\nFactor vs subtype (columns = dominant NMF factor):\n")
print(table(subtype = type, factor = assign3))

cat("\n[17_NMF] non-negativity makes the factors additive 'parts' (metagenes)",
    "that recover subtype structure without labels.\n")
