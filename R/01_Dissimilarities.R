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
##  As shipped, the source() line below borrows a few helpers from the course
##  repository: get_data() (loads a data set), need() (loads packages),
##  with_fig() (opens a plot device) and has_pkg()/skip_note() (let an optional
##  section be skipped rather than crash). To run this script entirely on its
##  own, delete that line and uncomment the block below. Nothing else changes.
##
##  The standalone with_fig() draws each figure to the current device -- inside
##  RStudio that is the Plot pane, so figures accumulate in the plot history.
##  To save them as files instead, replace its body with
##      png(paste0(name, ".png")); on.exit(dev.off()); force(expr)
##
##  Files to keep next to this script: none -- the data sets used here ship with R itself
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) {          # draw on the current device;
#   op <- par(no.readonly = TRUE); on.exit(par(op))  # in RStudio that is the Plot pane
#   invisible(force(expr))
# }
# has_pkg <- function(...) all(vapply(c(...), requireNamespace, logical(1), quietly = TRUE))
# skip_note <- function(what, pkgs) {
#   message("  [skipped] ", what, " -- needs ", paste(pkgs, collapse = ", "),
#           ":  install.packages(c(",
#           paste(sprintf('"%s"', pkgs), collapse = ", "), "))")
#   invisible(FALSE)
# }
# get_data <- function(name) { utils::data("iris");   iris }
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("vegan")

## ============================================================================
##  1. THE DATA
## ============================================================================
iris <- get_data("iris")
head(iris)
unique(iris[, 5])                      # which species are in the table
table(iris[, 5])                       # ... and how many of each
cat("\nThe table is", nrow(iris), "objects x", ncol(iris) - 1, "variables\n")

X <- iris[, 1:4]

## ============================================================================
##  2. R MODE -- comparing VARIABLES
## ============================================================================
##  In the R view the units of each variable do not matter to the *question*,
##  but they matter enormously to the *answer*: covariance is in squared units,
##  so a variable measured in millimetres dominates one measured in metres.
m.cov <- cov(X)
cat("\nCovariance matrix (R view):\n"); print(round(m.cov, 3))

m.cor <- cor(X)
cat("\nCorrelation matrix (R view):\n"); print(round(m.cor, 3))

## ---- standardization --------------------------------------------------------
Xs <- scale(X)                                   # z-transform: mean 0, sd 1

cat("\nRaw column means / sds:\n")
print(round(rbind(mean = colMeans(X),  sd = apply(X,  2, sd)), 3))
cat("\nStandardized column means / sds:\n")
print(round(rbind(mean = colMeans(Xs), sd = apply(Xs, 2, sd)), 3))

stopifnot(all(abs(colMeans(Xs)) < 1e-12))
stopifnot(all(abs(apply(Xs, 2, sd) - 1) < 1e-12))

## The point of standardizing, stated as an identity: the COVARIANCE of the
## standardized data IS the CORRELATION of the original data.
ms.cov <- cov(Xs)
cat("\ncov(scale(X)) equals cor(X):", isTRUE(all.equal(ms.cov, m.cor,
                                                       check.attributes = FALSE)), "\n")

## ============================================================================
##  3. Q MODE -- comparing OBJECTS
## ============================================================================
##  There are 150 objects, so 150*149/2 = 11175 pairwise distances.
##  vegan::vegdist() offers dozens of coefficients:
##
##    "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
##    "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
##    "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger",
##    "aitchison", "robust.aitchison"
##
##  A metric satisfies three conditions:  d(x,y) >= 0;  d(x,x) = 0;
##  and the triangle inequality  d(x,y) + d(y,z) >= d(x,z).
##  Several of the coefficients above are *not* metrics -- worth knowing before
##  handing one to a method that assumes Euclidean geometry.
cat("\nNumber of pairwise distances:", choose(nrow(X), 2), "\n")

d_euc  <- vegan::vegdist(X,  method = "euclidean", diag = TRUE)   # raw units
d_eucS <- vegan::vegdist(Xs, method = "euclidean", diag = TRUE)   # standardized
d_bray <- vegan::vegdist(X,  method = "bray",      diag = TRUE)
d_gow  <- vegan::vegdist(X,  method = "gower",     diag = TRUE)

cat("\nRange of each coefficient:\n")
print(round(sapply(list(euclidean = d_euc, euclidean_std = d_eucS,
                        bray = d_bray, gower = d_gow), range), 3))

## ---- raw vs standardized Euclidean ------------------------------------------
## The red line is y = x. Points depart from it because standardizing gives
## sepal width -- the smallest-variance column -- the same say as petal length.
with_fig("01_euc_vs_eucS", {
  plot(as.vector(d_euc), as.vector(d_eucS), pch = 19, col = "#1f3b7355",
       xlab = "Euclidean (raw)", ylab = "Euclidean (standardized)",
       main = "Standardization changes distances")
  abline(0, 1, col = "#8c2d3a", lwd = 2)
})

## ---- one distance against another -------------------------------------------
with_fig("01_dist_comparisons", {
  par(mfrow = c(1, 2))
  plot(as.vector(d_euc), as.vector(d_bray), pch = 19, col = "#1f3b7355",
       xlab = "Euclidean", ylab = "Bray-Curtis", main = "Euclidean vs Bray-Curtis")
  abline(0, 1, col = "#8c2d3a", lwd = 2)
  plot(as.vector(d_euc), as.vector(d_gow), pch = 19, col = "#1f3b7355",
       xlab = "Euclidean", ylab = "Gower", main = "Euclidean vs Gower")
  abline(0, 1, col = "#8c2d3a", lwd = 2)
}, width = 10, height = 5)

## ============================================================================
##  4. LOOKING AT A DISTANCE MATRIX
## ============================================================================
m_euc   <- as.matrix(d_euc)
m_eucS  <- as.matrix(d_eucS)
m_bray  <- as.matrix(d_bray)

## The same matrix under two palettes. The palette is not cosmetic: a sequential
## ramp shows the gradient, a rainbow invents boundaries that are not there.
with_fig("01_diss_palettes", {
  par(mfrow = c(1, 2))
  image(m_euc, col = heat.colors(20),       main = "heat.colors")
  image(m_euc, col = rainbow(20),           main = "rainbow (avoid)")
}, width = 10, height = 5)

## the figure used in the book
with_fig("01_diss_image", {
  image(m_euc, col = hcl.colors(20, "YlOrRd", rev = TRUE),
        main = "Q-mode Euclidean dissimilarity (iris)")
})

## ---- ordering the objects reveals the block structure -----------------------
##  Rows arrive in whatever order the file had. Sorting them by mean distance to
##  everything else turns a speckled square into a gradient: this is the whole
##  idea behind seriation, and the reason a clustered heatmap is readable.
##
##  (The original script computed iiE/jjE and then indexed with ii/jj, and drew
##  setosaE in both Bray-Curtis panels. Both are fixed here.)
order_block <- function(M) M[order(rowMeans(M)), order(colMeans(M))]

setosaE <- m_euc[1:50, 1:50]
setosaB <- m_bray[1:50, 1:50]

with_fig("01_setosa_ordered", {
  par(mfrow = c(2, 2))
  image(setosaE,              col = rainbow(20), main = "setosa Euclidean, as given")
  image(order_block(setosaE), col = rainbow(20), main = "setosa Euclidean, ordered")
  image(setosaB,              col = rainbow(20), main = "setosa Bray-Curtis, as given")
  image(order_block(setosaB), col = rainbow(20), main = "setosa Bray-Curtis, ordered")
}, width = 9, height = 8)

## ---- raw vs standardized, side by side --------------------------------------
with_fig("01_raw_vs_std_image", {
  par(mfrow = c(1, 2))
  image(m_euc,  col = rainbow(20), main = "Euclidean, raw")
  image(m_eucS, col = rainbow(20), main = "Euclidean, standardized")
}, width = 10, height = 5)

## ---- and the R-mode matrices ------------------------------------------------
with_fig("01_cov_cor_images", {
  par(mfrow = c(1, 2))
  image(m.cov, col = heat.colors(20), main = "Covariance")
  image(m.cor, col = heat.colors(20), main = "Correlation")
}, width = 10, height = 5)

cat("\n[01_Dissimilarities] Q view = objects (distances); R view = variables",
    "(correlations). Choice of distance is the first modelling decision.\n")
