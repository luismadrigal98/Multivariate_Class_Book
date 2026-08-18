## ============================================================================
##  MS_LJMR :: 08_SVD.R — Singular Value Decomposition & low-rank approximation
##
##  Original authors: Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 08_SingularValueDecomposition.R. Uses frozen limenitis image.
## ============================================================================
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

with_fig("08_svd_compression", width = 12, height = 3, {
  op <- par(mfrow = c(1, 5), mar = c(1,1,2,1)); on.exit(par(op))
  image(t(G[nrow(G):1, ]), col = gray(0:50/50), axes = FALSE, main = "Original")
  for (k in ks) {
    Gk <- approx(k)
    image(t(Gk[nrow(Gk):1, ]), col = gray(0:50/50), axes = FALSE,
          main = sprintf("rank %d (err %.3f)", k, relerr(k)))
  }
})

with_fig("08_svd_scree", {
  plot(d[1:15], type = "b", pch = 19, col = "#1f3b73", log = "y",
       xlab = "Index", ylab = "Singular value (log)",
       main = "Singular-value spectrum")
})

cat("\n[08_SVD] the first few singular triplets carry most of the matrix;",
    "truncating them is the best rank-k approximation (Eckart-Young).\n")
