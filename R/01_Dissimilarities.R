## ============================================================================
##  MS_LJMR :: 01_Dissimilarities.R — the Q and R views of a data matrix
##
##  Original authors: Jorge Soberón & Laura Jiménez
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 01_QandR_dissimilarities.R.
##  Fixes: the original referenced an undefined object `d_euc2`; here the
##  standardized distance is computed and compared correctly. No setwd()/x11().
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
need("vegan")

iris <- get_data("iris")
X    <- iris[, 1:4]

## ---- R view: covariance / correlation among VARIABLES -----------------------
m.cov <- cov(X); m.cor <- cor(X)
cat("Correlation matrix (R view):\n"); print(round(m.cor, 2))

## ---- standardization --------------------------------------------------------
Xs <- scale(X)                                   # mean 0, sd 1 per column
stopifnot(all(abs(colMeans(Xs)) < 1e-12))
stopifnot(all(abs(apply(Xs, 2, sd) - 1) < 1e-12))

## ---- Q view: dissimilarity among OBJECTS ------------------------------------
d_euc  <- vegan::vegdist(X,  method = "euclidean")   # raw units
d_eucS <- vegan::vegdist(Xs, method = "euclidean")   # standardized (fixed!)
d_bray <- vegan::vegdist(X,  method = "bray")
d_gow  <- vegan::vegdist(X,  method = "gower")

## raw vs standardized Euclidean: standardization changes the geometry
with_fig("01_euc_vs_eucS", {
  plot(as.vector(d_euc), as.vector(d_eucS), pch = 19, col = "#1f3b7355",
       xlab = "Euclidean (raw)", ylab = "Euclidean (standardized)",
       main = "Standardization changes distances")
  abline(0, 1, col = "#8c2d3a", lwd = 2)
})

## a distance heatmap (image of the Q-mode matrix)
with_fig("01_diss_image", {
  image(as.matrix(d_euc), col = hcl.colors(20, "YlOrRd", rev = TRUE),
        main = "Q-mode Euclidean dissimilarity (iris)")
})

cat("\n[01_Dissimilarities] Q view = objects (distances); R view = variables",
    "(correlations). Choice of distance is the first modelling decision.\n")
