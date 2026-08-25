## ============================================================================
##  MS_LJMR :: 21_PAMs.R — Presence-Absence Matrices (biodiversity structure)
##
##  Original author: Jorge Soberón
##                   (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 24_PAMs.R / 25_AnalysisOfPAMs.R. Uses frozen 'pam' matrix.
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
##  Files to keep next to this script: pam.csv
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
# get_data <- function(name) read.csv("pam.csv")
# data_dir <- function() "."                       # look for data files here
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

## ============================================================================
##  1. A PAM IS A BINARY MATRIX, AND ITS MARGINS ARE BIODIVERSITY
## ============================================================================
##  A presence-absence matrix records which species occur at which sites. If
##  rows are sites and columns are species, then
##      row sums  = how many species live at each site      -> alpha (richness)
##      col sums  = at how many sites each species is found -> omega (range size)
pamdf <- get_data("pam")
lat   <- pamdf$lat
pam   <- as.matrix(pamdf[, -1])            # sites x species, 0/1
n     <- nrow(pam); s <- ncol(pam)
cat("PAM:", n, "sites x", s, "species;  fill =", round(mean(pam), 3), "\n")

## ---- the marginal sums, three ways ------------------------------------------
##  rowSums/colSums are the obvious route. The matrix-product route is worth
##  seeing once, because it is the same algebra as everything else in this book:
##  multiplying by a vector of ones IS summing.
uc <- matrix(1, nrow = s, ncol = 1)        # column of ones, length = species
ur <- matrix(1, nrow = 1, ncol = n)        # row of ones,    length = sites
alphas <- pam %*% uc                        # n x 1
omegas <- ur  %*% pam                       # 1 x s
stopifnot(all.equal(as.vector(alphas), unname(rowSums(pam))),
          all.equal(as.vector(omegas), unname(colSums(pam))))
alpha <- as.vector(alphas); omega <- as.vector(omegas)
cat("mean alpha (richness/site) =", round(mean(alpha), 2),
    "   mean omega (range/species) =", round(mean(omega), 2), "\n")

## normalised versions: proportions rather than counts, so the two margins
## become comparable and dimensionless
alphast <- alpha / s                        # proportion of the species pool
omegast <- omega / n                        # proportion of the sites

with_fig("21_pam_margins", {
  par(mfrow = c(2, 2))
  hist(alpha,   col = "grey85", main = "richness per site (alpha)",  xlab = "")
  hist(omega,   col = "grey85", main = "range size per species (omega)", xlab = "")
  plot(density(alphast), main = "alpha, normalised", xlab = "")
  plot(density(omegast), main = "omega, normalised", xlab = "")
}, width = 10, height = 8)

## the book's figure
with_fig("21_pam", {
  par(mfrow = c(1, 3))
  image(t(pam[order(lat), ]), col = c("white", "#1f3b73"),
        main = "PAM (sites x species)", axes = FALSE)
  plot(lat, alpha, type = "l", col = "#1f3b73", lwd = 2,
       xlab = "latitude", ylab = "richness (alpha)", main = "Richness gradient")
  hist(omega, breaks = 20, col = "#2a7f7f", border = "white",
       xlab = "range size (omega)", main = "Range-size frequency")
}, width = 13, height = 3.8)

## ============================================================================
##  2. THE DISPERSION FIELD
## ============================================================================
##  alpha says how many species a site has. It says nothing about WHICH. Two
##  sites with ten species each are very different if one holds ten widespread
##  species and the other ten narrow endemics.
##
##  The dispersion field of a site is the mean normalised range size of the
##  species living there:   fi = (P omega*) / alpha .  A low value means the
##  site is full of restricted species -- an inverse measure of endemism.
fist     <- pam %*% omegast                 # summed range sizes per site
fistprom <- fist / alpha                    # averaged
cat("\ndispersion field: min", round(min(fistprom), 3),
    " max", round(max(fistprom), 3), "\n")

with_fig("21_dispersion_field", {
  par(mfrow = c(1, 2))
  plot(lat, fistprom, pch = 19, cex = .6, col = "#8c2d3a",
       xlab = "latitude", ylab = "mean range size of resident species",
       main = "Dispersion field")
  plot(alphast, fistprom, pch = 19, cex = .6, col = "#1f3b73",
       xlab = "proportional richness", ylab = "dispersion field",
       main = "Richness vs. endemism")
}, width = 11, height = 5)

##  With real coordinates this is a map -- richness in one panel, inverse
##  endemism in the other. The frozen demo PAM carries latitude only; drop
##  PAM_NA_WOGreenland.csv or pamMamms5.csv (site, long, lat, species...) into
##  data/ and the geographic version below runs instead.
pam_geo <- file.path(data_dir(), "PAM_NA_WOGreenland.csv")
if (file.exists(pam_geo)) {
  g    <- read.csv(pam_geo)
  crds <- as.matrix(g[, 2:3])
  M    <- as.matrix(g[, -(1:3)])
  ok   <- rowSums(M) > 0
  M    <- M[ok, colSums(M[ok, ]) > 0]; crds <- crds[ok, ]
  a    <- rowSums(M); o <- colSums(M) / nrow(M)
  fp   <- (M %*% o) / a
  with_fig("21_pam_maps", {
    par(mfrow = c(1, 2))
    plot(crds, col = hsv(a / max(a), 1, 1), pch = 19, asp = 1,
         main = "Richness")
    if (has_pkg("maps")) maps::map("world", add = TRUE, col = "grey40")
    plot(crds, col = hsv(fp / max(fp), 1, 1), pch = 19, asp = 1,
         main = "Inverse endemism (dispersion field)")
    if (has_pkg("maps")) maps::map("world", add = TRUE, col = "grey40")
  }, width = 13, height = 6)
} else {
  message("  [note] no georeferenced PAM in data/ -- plotting against latitude ",
          "only.\n         Add PAM_NA_WOGreenland.csv for the mapped version.")
}

## ============================================================================
##  3. RANGE-DIVERSITY PLOTS
## ============================================================================
##  rdp() is Jorge Soberon's range-diversity plot. Every site (or species) is
##  one point: proportional richness against dispersion field. The vertical line
##  sits at 1/beta, where beta is Whittaker's beta diversity, and the two
##  hyperbolas bound the region the points can occupy -- they are not fitted,
##  they are algebraic limits.
##
##  The original wrote all its intermediate quantities to the global environment
##  with `<<-`. This version returns them in a list instead, so nothing leaks
##  and the function can be called twice without the second call overwriting the
##  first one's results. The mathematics is unchanged.
rdp <- function(mat, view = 1, limits = 2) {
  mat <- as.matrix(mat)
  mat <- mat[, colSums(mat) > 0, drop = FALSE]     # drop empty species
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]     # drop empty sites
  if (view != 1) mat <- t(mat)                     # per-species view

  n <- nrow(mat); s <- ncol(mat)
  a  <- rowSums(mat); o <- colSums(mat)
  at <- a / s                                      # proportional richness
  ot <- o / n                                      # proportional range
  fi <- as.vector((mat %*% ot) / a)                # dispersion field

  beta <- 1 / mean(at)                             # Whittaker's beta
  rho  <- at * (fi - 1 / beta)                     # covariance term

  xl <- if (limits == 1) 1 else max(fi) * 1.1
  yl <- if (limits == 1) 1 else max(at) * 1.1
  xM <- seq(1.01 / beta, 1, length = 100); yM <- max(rho) / (xM - 1 / beta)
  xm <- seq(0, 0.99 / beta, length = 50);  ym <- min(rho) / (xm - 1 / beta)

  plot(fi, at, xlim = c(0, xl), ylim = c(0, yl), pch = 19, cex = .6,
       col = "#1f3b7399",
       xlab = if (view == 1) "Dispersion field" else "Richness field",
       ylab = if (view == 1) "Proportional richness" else "Proportional range")
  abline(v = 1 / beta, col = "grey40")
  lines(xM, yM); lines(xm, ym)
  legend("topright", legend = sprintf("Beta = %.3f", beta), bty = "n")
  invisible(list(alpha = a, omega = o, alphast = at, omegast = ot,
                 field = fi, beta = beta, rho = rho))
}

with_fig("21_range_diversity", {
  par(mfrow = c(1, 2))
  rs <- rdp(pam, view = 1, limits = 2)     # per sites
  title("Per site")
  rp <- rdp(pam, view = 2, limits = 2)     # per species
  title("Per species")
}, width = 12, height = 6)

## the same two views with axes fixed to the unit square, which is how they are
## usually published -- comparable across data sets
with_fig("21_range_diversity_unit", {
  par(mfrow = c(1, 2))
  rdp(pam, view = 1, limits = 1); title("Per site (unit axes)")
  rdp(pam, view = 2, limits = 1); title("Per species (unit axes)")
}, width = 12, height = 6)

rs <- rdp(pam, 1, 2)
cat("\nWhittaker's beta (sites view):", round(rs$beta, 3), "\n")
cat("1/beta, the vertical line:", round(1 / rs$beta, 3), "\n")

## ============================================================================
##  4. THE TWO MATRIX PRODUCTS
## ============================================================================
##  P P' counts SHARED SPECIES between every pair of sites.
##  P' P counts CO-OCCURRENCES between every pair of species.
##  Everything above falls out of these two, which closes the loop back to the
##  Q and R views of session 01.
PPt <- pam %*% t(pam)                       # sites x sites
PtP <- t(pam) %*% pam                       # species x species
cat("\nP P' is", paste(dim(PPt), collapse = " x "),
    " and its diagonal IS alpha:", isTRUE(all.equal(unname(diag(PPt)), alpha)), "\n")
cat("P'P is", paste(dim(PtP), collapse = " x "),
    " and its diagonal IS omega:", isTRUE(all.equal(unname(diag(PtP)), omega)), "\n")

with_fig("21_pam_products", {
  par(mfrow = c(1, 2))
  o1 <- order(lat)
  image(PPt[o1, o1], col = hcl.colors(24, "YlGnBu", rev = TRUE), axes = FALSE,
        main = "P P' : species shared between sites")
  o2 <- order(omega)
  image(PtP[o2, o2], col = hcl.colors(24, "YlOrRd", rev = TRUE), axes = FALSE,
        main = "P'P : species co-occurrence")
}, width = 11, height = 5.5)

##  And because they are just Q-mode and R-mode matrices, every method in this
##  book applies to them: a similarity from P P' can be clustered (session 04),
##  ordinated (07), or scaled (09).
sim <- PPt / outer(alpha, alpha, pmin)      # proportion of the smaller list shared
d   <- as.dist(1 - sim)
tree <- hclust(d, method = "ward.D2")
with_fig("21_pam_clustering", {
  plot(tree, labels = FALSE, main = "Sites clustered by shared species",
       xlab = "", sub = "")
  rect.hclust(tree, k = 3, border = c("#1f3b73", "#8c2d3a", "#2a7f7f"))
})
cat("\nclusters of sites vs latitude:\n")
print(table(cluster = cutree(tree, 3), band = cut(lat, 3)))

cat("\n[21_PAMs] a PAM's row and column sums (alpha, omega) drive richness and",
    "range-size patterns; the dispersion field adds WHICH species, not just how",
    "many; and range-diversity plots put both on one pair of axes.\n")
