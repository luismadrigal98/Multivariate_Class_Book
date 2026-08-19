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
##  As shipped, the source() line below borrows three helpers from the course
##  repository: get_data() (loads a data set), need() (loads packages) and
##  with_fig() (opens a plot device). To run this script entirely on its own,
##  delete that line and uncomment the block below. Nothing else changes.
##
##  The standalone with_fig() just draws each figure to the screen, one after
##  the other. To save them as files instead, replace its body with
##      png(paste0(name, ".png")); on.exit(dev.off()); force(expr)
##
##  Files to keep next to this script: Limenitis_archippus.csv
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) as.matrix(read.csv("Limenitis_archippus.csv", header = FALSE))
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

G <- get_data("limenitis")               # 69 x 100 grayscale intensity matrix
s <- svd(G)                              # G = U D V'
d <- s$d

## rank-k reconstruction and its relative error
approx <- function(k) s$u[, 1:k] %*% diag(d[1:k], k, k) %*% t(s$v[, 1:k])
relerr <- function(k) norm(G - approx(k), "F") / norm(G, "F")
ks  <- c(1, 5, 10, 50)
err <- sapply(ks, relerr)
cat("Relative Frobenius error at ranks", paste(ks, collapse = ", "), ":\n")
print(round(err, 3))

with_fig("06_svd_compression", width = 12, height = 3, {
  op <- par(mfrow = c(1, 5), mar = c(1,1,2,1)); on.exit(par(op))
  image(t(G[nrow(G):1, ]), col = gray(0:50/50), axes = FALSE, main = "Original")
  for (k in ks) {
    Gk <- approx(k)
    image(t(Gk[nrow(Gk):1, ]), col = gray(0:50/50), axes = FALSE,
          main = sprintf("rank %d (err %.3f)", k, relerr(k)))
  }
})

with_fig("06_svd_scree", {
  plot(d[1:15], type = "b", pch = 19, col = "#1f3b73", log = "y",
       xlab = "Index", ylab = "Singular value (log)",
       main = "Singular-value spectrum")
})

cat("\n[06_SVD] the first few singular triplets carry most of the matrix;",
    "truncating them is the best rank-k approximation (Eckart-Young).\n")
