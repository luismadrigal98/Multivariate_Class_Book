## ============================================================================
##  MS_LJMR :: 05_LinearAlgebra.R — transformations, eigenvectors, ellipses
##
##  Original authors: Jorge Soberón & Laura Jiménez
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 07_Transformations.R. No x11(); portable device via with_fig().
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/05_LinearAlgebra.R by scripts/make_class_src.py --
##  edit R/05_LinearAlgebra.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("ellipse")
library(ellipse)

## ============================================================================
##  1. A LINEAR MAP IS A MATRIX MULTIPLICATION
## ============================================================================
##  X is a swarm of four points in two dimensions -- the corners of a square.
##  Columns are points, rows are coordinates, so a 2x2 matrix acting on the left
##  transforms every point at once.
X <- matrix(c(1,1, 10,1, 1,10, 10,10), nrow = 2)
colnames(X) <- c("a", "b", "c", "d")
print(X)

par(mfrow = c(1, 1))
plot(t(X), xlim = c(-20, 30), ylim = c(-20, 30), xlab = "x", ylab = "y",
     pch = 19, main = "Four points before anything happens")
text(t(X), colnames(X), pos = 4, cex = 1.2)

##  A helper that draws the original swarm joined into a square (black) and the
##  transformed one on top (red). Drawing the EDGES is what makes the difference
##  between a rotation and a shear visible -- four loose points do not show it.
trans <- function(m1, m2, title) {
  edges <- function(m, col) {
    p <- t(m)
    for (e in list(c(1,2), c(1,3), c(3,4), c(4,2)))
      lines(p[e, 1], p[e, 2], col = col)
  }
  plot(t(m1), xlim = c(-20, 30), ylim = c(-20, 30), asp = 1, main = title,
       xlab = "x", ylab = "y")
  text(t(m1), colnames(m1), pos = 4, cex = .8)
  edges(m1, "black")
  points(t(m2), pch = 19, col = "#8c2d3a")
  text(t(m2), colnames(m2), pos = 4, cex = .8, col = "#8c2d3a")
  edges(m2, "#8c2d3a")
}

## four maps: identity, pure stretch, pure rotation, and both at once
U1 <- matrix(c(1, 0, 0, 1),      nrow = 2)      # identity: does nothing
U2 <- matrix(c(2, 0, 0, 1.2),    nrow = 2)      # stretch along the axes
U3 <- matrix(c(.8, .6, -.6, .8), nrow = 2)      # rotation (columns orthonormal)
U4 <- matrix(c(1.5, .1, .9, 1),  nrow = 2)      # stretch + rotation + shear

Y1 <- U1 %*% X; Y2 <- U2 %*% X; Y3 <- U3 %*% X; Y4 <- U4 %*% X

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 10, height = 9)
op <- par(mfrow = c(2, 2))
trans(X, Y1, "Nothing (identity)")
trans(X, Y2, "Stretches")
trans(X, Y3, "Rotates")
trans(X, Y4, "Stretches and rotates")
par(op)

## Why is U3 a rotation? Because its columns are orthonormal, so it preserves
## lengths and angles. That is checkable, not a matter of opinion:
cat("\nU3'U3 = I ?\n"); print(round(t(U3) %*% U3, 10))
cat("determinants -- identity:", det(U1), " stretch:", det(U2),
    " rotation:", round(det(U3), 3), " general:", det(U4), "\n")
cat("(a determinant is the factor by which AREA changes; 1 means area-preserving)\n")

## ---- every 2-D rotation, in one function ------------------------------------
##  CAUTION: R's cos() and sin() take RADIANS. The original script commented
##  "angle in degrees" and then called rot(35) and rot(45), which are 35 and 45
##  *radians* -- a rotation of about 2005 and 2578 degrees. The version below
##  takes degrees and converts, so the name and the behaviour agree.
rot_rad <- function(a)
  matrix(c(cos(a), sin(a), -sin(a), cos(a)), nrow = 2)
rot <- function(degrees) round(rot_rad(degrees * pi / 180), 3)

cat("\nrot(45):\n");  print(rot(45))
cat("rot(90) applied to (1,0) should be (0,1):\n")
print(round(rot(90) %*% c(1, 0), 10))

## ============================================================================
##  2. A CLOUD, ROTATED
## ============================================================================
##  Draw 20000 points from independent normals with different sds: an
##  axis-aligned ellipse. Then rotate it 45 degrees, and the two variables
##  become correlated -- correlation is what a rotation looks like in a table.
set.seed(1998)
xx  <- rnorm(20000, mean = 0, sd = 3)
yy  <- rnorm(20000, mean = 0, sd = 1)
mat <- cbind(xx, yy)

mat2 <- t(rot(45) %*% t(mat))

par(mfrow = c(1, 1))
plot(mat, xlim = c(-10, 10), ylim = c(-10, 10), pch = ".", asp = 1,
     xlab = "x", ylab = "y", main = "Before (grey) and after a 45 deg rotation")
points(mat2, pch = ".", col = "orange")

cat("\nCorrelation before rotation:", round(cor(mat)[1, 2], 3),
    "  after:", round(cor(mat2)[1, 2], 3), "\n")

## ============================================================================
##  3. UNDOING IT: EIGENVECTORS OF THE COVARIANCE MATRIX
## ============================================================================
cov2 <- cov(mat2)
cat("\nCovariance of the rotated cloud:\n"); print(round(cov2, 3))

##  The axes of the cloud are the eigenvectors of the covariance matrix, and the
##  variances along them are the eigenvalues.
eg <- eigen(cov2)
cat("\nEigenvalues (variances along the principal axes):\n")
print(round(eg$values, 3))
cat("Eigenvectors (the axes, as columns):\n"); print(round(eg$vectors, 3))
cat("sqrt of eigenvalues:", round(sqrt(eg$values), 3),
    " -- compare with the sds we simulated, 3 and 1\n")
cat("angle of the first axis:",
    round(atan2(eg$vectors[2, 1], eg$vectors[1, 1]) * 180 / pi, 1), "degrees\n")

##  The original worked with eigen(solve(cov2)) -- the inverse covariance, or
##  precision matrix. It has the SAME eigenvectors with reciprocal eigenvalues,
##  which is why its "stretchers" were 2/sqrt(values) rather than 2*sqrt(values).
##  Both routes draw the same axes; going through cov() directly is easier to
##  connect to PCA in the next session.
##  eigen() always sorts by DECREASING eigenvalue, so inverting the matrix
##  reverses the column order -- the axes are the same, listed the other way up.
eg_inv <- eigen(solve(cov2))
cat("\nsame axes from the inverse covariance:",
    isTRUE(all.equal(abs(eg$vectors),
                     abs(eg_inv$vectors[, 2:1]), check.attributes = FALSE)), "\n")
cat("its eigenvalues are the reciprocals:",
    round(rev(1 / eg_inv$values), 3), "\n")

len <- 2 * sqrt(eg$values)                       # 2 sd along each axis

par(mfrow = c(1, 1))
plot(mat2, pch = ".", col = "#c8a23a", asp = 1,
     xlab = "x", ylab = "y",
     main = "Eigenvectors are the axes of dispersion")
arrows(0, 0, eg$vectors[1, ] * len, eg$vectors[2, ] * len,
       col = "#1f3b73", lwd = 3, length = .1)
## statistical ellipses: contours containing 25% and 95% of a Gaussian with
## this covariance. The 95% ellipse is the picture of "two standard deviations"
## generalised to two dimensions.
lines(ellipse::ellipse(cov2, centre = c(0, 0), level = 0.25), lwd = 2)
lines(ellipse::ellipse(cov2, centre = c(0, 0), level = 0.95), lwd = 2)

cat("\n[05_LinearAlgebra] eigenvectors of the covariance matrix are the",
    "principal axes; eigenvalues are the variances along them (this IS PCA).\n")
