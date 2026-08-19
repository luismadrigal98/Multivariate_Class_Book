## ============================================================================
##  MS_LJMR :: 15_MDS.R — (non-metric) Multidimensional Scaling
##
##  Original authors: Marlon E. Cobos & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 15_MultidimensionalScaling.R. Uses mtcars (built-in).
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
# get_data <- function(name) { utils::data("mtcars"); mtcars }
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("vegan")

cars <- get_data("mtcars")

## nMDS preserves the RANK ORDER of dissimilarities by minimizing 'stress'
set.seed(1998)
md <- vegan::metaMDS(cars, distance = "euclidean", k = 2, trace = FALSE,
                     autotransform = FALSE)
cat("nMDS stress:", round(md$stress, 4),
    "  (<0.05 excellent, <0.1 good, <0.2 usable)\n")

with_fig("15_mds", width = 10, height = 4.4, {
  op <- par(mfrow = c(1, 2)); on.exit(par(op))
  plot(md$points, pch = 19, col = "#1f3b73", main = "nMDS configuration",
       xlab = "dim 1", ylab = "dim 2")
  text(md$points, labels = rownames(cars), cex = .5, pos = 3)
  vegan::stressplot(md, main = "Shepard diagram")   # fit of distances
})

## compare two distance choices with a Procrustes test
mdE <- vegan::metaMDS(cars, distance = "euclidean", k = 2, trace = FALSE,
                      autotransform = FALSE)
mdG <- vegan::metaMDS(cars, distance = "gower", k = 2, trace = FALSE,
                      autotransform = FALSE)
pr <- vegan::protest(mdE, mdG, permutations = 999)
cat("\nProcrustes correlation (Euclidean vs Gower):", round(pr$t0, 3),
    "  p =", pr$signif, "\n")

cat("\n[15_MDS] stress is an optimisation score, not a probability; a tight,",
    "monotone Shepard diagram means the low-D map is faithful.\n")
