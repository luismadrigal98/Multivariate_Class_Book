## ============================================================================
##  MS_LJMR :: 07_LinearAlgebra.R — transformations, eigenvectors, ellipses
##
##  Original authors: Jorge Soberón & Laura Jiménez
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 07_Transformations.R. No x11(); portable device via with_fig().
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
##  Files to keep next to this script: none — this session builds its own data
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("ellipse")

## ---- linear maps are matrix multiplications --------------------------------
X <- matrix(c(1,1, 10,1, 1,10, 10,10), nrow = 2)   # four corner points
colnames(X) <- c("a","b","c","d")

trans <- function(m1, m2, title) {
  plot(t(m1), xlim = c(-20, 30), ylim = c(-20, 30), asp = 1, main = title,
       xlab = "x", ylab = "y"); text(t(m1), colnames(m1), pos = 4, cex = .8)
  points(t(m2), pch = 19, col = "#8c2d3a")
  text(t(m2), colnames(m2), pos = 4, cex = .8, col = "#8c2d3a")
}
U_stretch <- matrix(c(2, 0, 0, 1.2), 2)            # stretch
rot <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2)  # rotate
U_rot <- rot(pi/4)

with_fig("07_transformations", {
  op <- par(mfrow = c(1, 2)); on.exit(par(op))
  trans(X, U_stretch %*% X, "Stretch")
  trans(X, U_rot %*% X,     "Rotate (45 deg)")
})

## ---- eigen-decomposition recovers the axes of a data ellipsoid -------------
set.seed(1998)
cloud <- t(rot(pi/4) %*% rbind(rnorm(2000, 0, 3), rnorm(2000, 0, 1)))
S     <- cov(cloud)
eg    <- eigen(S)                                  # eigenvectors = axes
cat("Eigenvalues (variances along principal axes):\n"); print(round(eg$values, 3))
len   <- 2 * sqrt(eg$values)                        # axis half-lengths (2 sd)

with_fig("07_ellipse_axes", {
  plot(cloud, pch = ".", col = "#c8a23a", asp = 1,
       main = "Eigenvectors are the axes of dispersion")
  arrows(0, 0, eg$vectors[1,] * len, eg$vectors[2,] * len,
         col = "#1f3b73", lwd = 3, length = .1)
  lines(ellipse::ellipse(S, centre = colMeans(cloud), level = .95), lwd = 2)
})

cat("\n[07_LinearAlgebra] eigenvectors of the covariance matrix are the",
    "principal axes; eigenvalues are the variances along them (this IS PCA).\n")
