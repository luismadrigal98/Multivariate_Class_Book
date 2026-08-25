## ============================================================================
##  MS_LJMR :: 10_NMF.R — Non-negative Matrix Factorization
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
##  Files to keep next to this script: CountriesToLive.csv; leukemiaExpressionSubset.rds
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
# get_data <- function(name) switch(name,
#   countries_live = read.csv("CountriesToLive.csv", stringsAsFactors = TRUE),
#   doubs          = { utils::data("doubs", package = "ade4")
#                     list(fish = doubs$fish, env = doubs$env, xy = doubs$xy) },
#   leukemia       = readRDS("leukemiaExpressionSubset.rds"),
#   stop("unknown data set: ", name))
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("NMF")

## ============================================================================
##  1. THE IDEA
## ============================================================================
##  Let A be the n x v data matrix. We want
##
##        A  ~  W H       with W (n x rank) and H (rank x v), both >= 0
##
##  The columns of W are the BASIS, the columns of H are the WEIGHTS.
##  A = WH means every column of A is a weighted sum of the columns of W, the
##  weights being the corresponding column of H.
##
##  Two properties do the interpretive work:
##    * non-negativity -- parts add, they never cancel, so a basis vector is a
##      "part of the whole" rather than a contrast the way a PC is;
##    * sparsity of W -- each part involves few variables.
##
##  Unlike the SVD of session 06, this problem has no closed-form solution and
##  is not convex in W and H jointly, so nmf() runs an iterative algorithm from
##  a seed. Fix the seed or the answer moves.
expr <- get_data("leukemia")             # genes x samples
A    <- expr - min(expr)                 # ensure non-negativity
type <- sub("\\..*$", "", colnames(expr))
cat("expression matrix:", nrow(A), "genes x", ncol(A), "samples\n")
cat("leukemia subtypes:\n"); print(table(type))

## ============================================================================
##  2. THE FACTORIZATION
## ============================================================================
set.seed(1998)
res <- NMF::nmf(A, rank = 3, seed = "nndsvd")
w   <- NMF::basis(res)                   # W: genes x 3
h   <- NMF::coef(res)                    # H: 3 x samples
cat("\nW:", paste(dim(w), collapse = " x "),
    "  H:", paste(dim(h), collapse = " x "), "\n")
cat("reconstruction error (Frobenius):",
    round(norm(A - w %*% h, "F") / norm(A, "F"), 4), "\n")
cat("fraction of exact zeros in W:", round(mean(w == 0), 3), "\n")

## ---- the weights separate the subtypes --------------------------------------
##  Each sample is a point in "amount of each part" space. If the parts are
##  biologically real, the subtypes should separate without ever being told.
pal <- c("#8c2d3a", "#1f3b73", "#2a7f7f")
with_fig("10_nmf", {
  par(mfrow = c(1, 2))
  image(t(w[order(w[, 1]), ]), col = hcl.colors(30, "YlGnBu", rev = TRUE),
        main = "W: basis (genes x 3)", axes = FALSE)
  plot(h[2, ], h[3, ], col = pal[factor(type)], pch = 19,
       xlab = "factor 2", ylab = "factor 3",
       main = "H: samples in factor space")
  legend("topright", legend = levels(factor(type)), fill = pal, bty = "n")
}, width = 11, height = 5)

## all three pairs of factors at once
with_fig("10_nmf_factor_pairs", {
  pairs(t(h), col = pal[factor(type)], pch = 19, cex = .8,
        labels = paste("factor", 1:3),
        main = "Samples in the space of the three parts")
})

## ---- basismap / coefmap ------------------------------------------------------
##  The NMF package's own displays: basismap() clusters the rows of W, coefmap()
##  the columns of H, both with subtype annotation. They need Biobase.
if (has_pkg("Biobase")) {
  with_fig("10_nmf_maps", {
    layout(cbind(1, 2))
    NMF::basismap(res)
    NMF::coefmap(res)
  }, width = 12, height = 6)
} else skip_note("basismap / coefmap displays", "Biobase")

## ============================================================================
##  3. A SMALL, READABLE EXAMPLE: WHERE WOULD YOU LIVE?
## ============================================================================
##  Thirteen countries ranked on six criteria. Small enough that you can read
##  the parts straight off the heatmaps.
ctl2 <- get_data("countries_live")
paises <- as.character(ctl2[[1]])
ctl    <- as.matrix(t(ctl2[, -1]))
colnames(ctl) <- paises
cat("\nliveability matrix:", nrow(ctl), "criteria x", ncol(ctl), "countries\n")
print(ctl[, 1:min(6, ncol(ctl))])

set.seed(1998)
nmfCtl <- NMF::nmf(ctl, 3)               # rank = 3 means three basis vectors
with_fig("10_nmf_countries", {
  par(mfrow = c(1, 2))
  ## dark = this criterion loads heavily on that part
  heatmap(NMF::basis(nmfCtl), main = "W, basis", Colv = NA, Rowv = NA,
          scale = "none")
  ## dark = this country carries a lot of that part
  heatmap(NMF::coef(nmfCtl),  main = "H, weights", Colv = NA, Rowv = NA,
          scale = "none")
}, width = 12, height = 6)

if (has_pkg("Biobase")) {
  with_fig("10_nmf_countries_maps", {
    layout(cbind(1, 2))
    NMF::basismap(nmfCtl)
    NMF::coefmap(nmfCtl)
  }, width = 12, height = 6)
}

## ============================================================================
##  4. COUNT DATA: BIRDS AND FISH
## ============================================================================
##  NMF is a natural fit for species-by-site count tables: counts are already
##  non-negative, and "parts" read as assemblages.
if (has_pkg("ade4")) {
  e <- new.env(); utils::data("atlas", package = "ade4", envir = e)
  birds <- as.matrix(t(e$atlas$birds))       # 19 species x 23 sites
  cat("\natlas birds:", nrow(birds), "species x", ncol(birds), "sites\n")
  set.seed(1998)
  nmfAtlas <- NMF::nmf(birds, 3)
  with_fig("10_nmf_birds", {
    par(mfrow = c(1, 2))
    heatmap(NMF::basis(nmfAtlas), main = "W, basis",   Colv = NA, Rowv = NA,
            scale = "none")
    heatmap(NMF::coef(nmfAtlas),  main = "H, weights", Colv = NA, Rowv = NA,
            scale = "none")
  }, width = 12, height = 6)
} else skip_note("the ade4 'atlas' bird example", "ade4")

## ---- the Doubs fish, and the parts mapped back onto the river ----------------
doubs <- get_data("doubs")
keep  <- rowSums(doubs$fish) > 0
fish  <- t(doubs$fish[keep, ])
fish  <- fish[rowSums(fish) > 0, ]
cat("\ndoubs fish (transposed):", nrow(fish), "species x", ncol(fish), "sites\n")

set.seed(1998)
nmfF <- NMF::nmf(fish, 3)
hF   <- NMF::coef(nmfF)

with_fig("10_nmf_fish", {
  par(mfrow = c(1, 2))
  heatmap(NMF::basis(nmfF), main = "W, basis (species)",  Colv = NA, Rowv = NA,
          scale = "none")
  heatmap(hF,               main = "H, weights (sites)",  Colv = NA, Rowv = NA,
          scale = "none")
}, width = 12, height = 6)

##  The payoff: colour each site on the real river map by how much of part 1 it
##  carries. A latent factor found with no spatial information at all turns out
##  to be the upstream-downstream gradient.
if (!is.null(doubs$xy)) {
  xy <- doubs$xy[keep, ]
  with_fig("10_nmf_fish_map", {
    par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    for (k in 1:3)
      plot(xy, pch = 19, cex = 1.6, asp = 1,
           col = hsv(hF[k, ] / max(hF[k, ]), 1, 1),
           main = paste("part", k, "along the river"))
  }, width = 13, height = 5)
}

## ============================================================================
##  5. CHOOSING THE RANK
## ============================================================================
##  There is no eigenvalue to look at. The usual approach is to fit several
##  ranks and watch the residual and the cophenetic correlation (how stable the
##  sample clustering is across runs). A drop in cophenetic correlation is the
##  conventional signal that you have gone one rank too far.
set.seed(1998)
ranks <- 2:5
resid <- sapply(ranks, function(k) {
  f <- NMF::nmf(A, k, seed = "nndsvd")
  norm(A - NMF::basis(f) %*% NMF::coef(f), "F") / norm(A, "F")
})
cat("\nrelative residual by rank:\n")
print(round(setNames(resid, paste0("rank", ranks)), 4))

with_fig("10_nmf_rank", {
  plot(ranks, resid, type = "b", pch = 19, col = "#1f3b73",
       xlab = "rank", ylab = "relative residual",
       main = "Residual falls with rank -- always")
})

cat("\n[10_NMF] A ~ WH with everything non-negative: parts that add rather than",
    "contrasts that cancel, which is why the factors are readable as",
    "assemblages, gene programmes or criteria.\n")
