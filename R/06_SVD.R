## ============================================================================
##  MS_LJMR :: 06_SVD.R — Singular Value Decomposition & low-rank approximation
##
##  Original authors: Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 08_SingularValueDecomposition.R. Uses frozen limenitis image.
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
##  Files to keep next to this script: Limenitis_archippus.csv
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
# get_data <- function(name) as.matrix(read.csv("Limenitis_archippus.csv", header = FALSE))
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

## ============================================================================
##  1. AN IMAGE IS A MATRIX
## ============================================================================
##  Limenitis_archippus.csv is a headerless numeric grid: a greyscale photograph
##  of a viceroy butterfly wing, one intensity per pixel. Nothing about the SVD
##  cares that it is a picture -- but a picture makes the rank-k approximation
##  visible in a way a table of numbers never does.
G  <- get_data("limenitis")
nr <- nrow(G); nc <- ncol(G)
cat("Image matrix:", nr, "x", nc, " (", nr * nc, "cells )\n")
cat("intensity range:", round(range(G), 3), "\n")

with_fig("06_image_raw", {
  par(mfrow = c(1, 2), mar = c(2, 2, 3, 1))
  image(G, main = "default palette")
  image(G, col = gray(1:50 / 50), main = "greyscale")
}, width = 10, height = 5)

## ============================================================================
##  2. THE DECOMPOSITION
## ============================================================================
##      G = U D V'
##  where D is diagonal with the singular values, U is an orthonormal basis for
##  the columns and V an orthonormal basis for the rows.
s <- svd(G)
str(s)

## The number of non-zero singular values IS the rank of the matrix. qr() will
## tell you the rank independently, which is a useful cross-check.
rnk <- qr(G)$rank
cat("\nRank from qr():", rnk, "  non-negligible singular values:",
    sum(s$d > max(dim(G)) * .Machine$double.eps * max(s$d)), "\n")

d <- s$d[seq_len(min(nr, nc))]
u <- as.matrix(s$u)
v <- as.matrix(s$v)
options(digits = 3)
cat("first ten singular values:\n"); print(d[1:10])

## ---- scree: how fast do they decay? -----------------------------------------
with_fig("06_svd_scree", {
  par(mfrow = c(1, 2))
  plot(d, type = "b", pch = 19, col = "#1f3b73",
       xlab = "index", ylab = "singular value", main = "All singular values")
  plot(d[1:10], type = "b", pch = 19, col = "#1f3b73",
       xlab = "index", ylab = "singular value", main = "First ten")
}, width = 10, height = 5)

cat("\nshare of squared 'energy' in the first k:\n")
cum <- cumsum(d^2) / sum(d^2)
print(round(cum[c(1, 5, 10, 50, 100)], 4))

## ============================================================================
##  3. REBUILDING THE MATRIX FROM ITS PARTS
## ============================================================================
##  G = sum over i of  d[i] * u[,i] %*% t(v[,i])
##  Each term is a full-size matrix of rank one -- an outer product of a column
##  pattern and a row pattern -- weighted by its singular value. Truncating the
##  sum after n terms gives the best rank-n approximation there is, in the
##  least-squares sense (Eckart & Young 1936).
compact <- function(u, d, v, n) {
  result <- matrix(0, nrow = nrow(u), ncol = nrow(v))
  for (i in seq_len(n))
    result <- result + d[i] * (u[, i] %*% t(v[, i]))   # outer product, weighted
  result
}

## the same thing without the loop, for comparison
compact_fast <- function(u, d, v, n)
  u[, 1:n, drop = FALSE] %*% diag(d[1:n], n, n) %*% t(v[, 1:n, drop = FALSE])

c5   <- compact(u, d, v, 5)
c10  <- compact(u, d, v, 10)
c50  <- compact(u, d, v, 50)
c100 <- compact(u, d, v, 100)
cat("\nloop and matrix forms agree:",
    isTRUE(all.equal(c10, compact_fast(u, d, v, 10))), "\n")

relerr <- function(A) norm(G - A, "F") / norm(G, "F")
cat("\nrelative error of the rank-k approximation:\n")
print(round(c(k5 = relerr(c5), k10 = relerr(c10),
              k50 = relerr(c50), k100 = relerr(c100)), 4))

with_fig("06_svd_compression", {
  par(mfrow = c(1, 5), mar = c(1, 1, 2, 1))
  image(c5,   col = gray(1:50 / 50), asp = 1, axes = FALSE, main = "5 terms")
  image(c10,  col = gray(1:50 / 50), asp = 1, axes = FALSE, main = "10 terms")
  image(c50,  col = gray(1:50 / 50), asp = 1, axes = FALSE, main = "50 terms")
  image(c100, col = gray(1:50 / 50), asp = 1, axes = FALSE, main = "100 terms")
  image(G,    col = gray(1:50 / 50), asp = 1, axes = FALSE,
        main = paste("rank =", rnk))
}, width = 13, height = 3.4)

## ---- what did that buy? ------------------------------------------------------
##  A rank-k approximation stores k*(nr + nc + 1) numbers instead of nr*nc.
store <- function(k) k * (nr + nc + 1)
cat("\nnumbers stored, full vs rank-k:\n")
print(data.frame(k = c(5, 10, 50, 100),
                 stored = store(c(5, 10, 50, 100)),
                 full   = nr * nc,
                 ratio  = round(store(c(5, 10, 50, 100)) / (nr * nc), 3)))

cat("\n[06_SVD] every matrix is a weighted sum of rank-one outer products;",
    "keeping the largest few is the best low-rank approximation possible.\n")
