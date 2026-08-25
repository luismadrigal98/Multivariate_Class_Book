## ============================================================================
##  MS_LJMR :: 08_CorrespondenceAnalysis.R — Correspondence Analysis (CA)
##
##  Original authors: Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 14_CorrespondaceAnalysis.R. Uses frozen Doubs fish counts.
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
##  Files to keep next to this script: ButterfliesQRoo2.csv
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
#   butterflies = read.csv("ButterfliesQRoo2.csv", stringsAsFactors = TRUE),
#   doubs       = { utils::data("doubs", package = "ade4")
#                     list(fish = doubs$fish, env = doubs$env, xy = doubs$xy) },
#   stop("unknown data set: ", name))
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("vegan")

## ============================================================================
##  1. WHY A DIFFERENT DISTANCE
## ============================================================================
##  Correspondence analysis is PCA adapted to COUNTS -- non-negative integers,
##  which are not normally distributed and have no meaningful zero-centred
##  ellipse. It has been called reciprocal averaging and dual scaling; the
##  mathematics is the same in every case.
##
##  The appropriate distance is chi-square, not Euclidean. CA eigen-decomposes
##  the chi-square distances between ROWS (sites), and the same eigenvalues
##  scale to the chi-square distances between COLUMNS (species). So CA is an
##  R-mode and a Q-mode analysis at once -- which follows from the fact, from
##  session 01, that XX' and X'X share their non-zero eigenvalues.
##
##  Data: Verneaux (1973) -- 27 fish species, 11 environmental variables and
##  coordinates for 30 sites along the Doubs river in the Jura.
doubs <- get_data("doubs")
cat("tables in doubs:", paste(names(doubs), collapse = ", "), "\n")
print(head(doubs$fish[, 1:8]))
print(head(doubs$env))

## ---- the sites are a river, so they are spatially ordered -------------------
if (!is.null(doubs$xy)) {
  with_fig("08_doubs_geography", {
    plot(doubs$xy, pch = 19, col = "seagreen", asp = 1,
         main = "The 30 sites, in real geography")
    text(doubs$xy, labels = seq_len(nrow(doubs$xy)), pos = 3, cex = .7)
  })
  cat("\nThe sites trace the river downstream -- so spatial covariance is not a",
      "hypothesis, it is a certainty.\n")
}

## CA cannot handle an all-zero row or column
cat("\nempty sites:",   paste(which(rowSums(doubs$fish) == 0), collapse = ", "), "\n")
cat("empty species:", paste(which(colSums(doubs$fish) == 0), collapse = ", "), "\n")
keep <- rowSums(doubs$fish) > 0
fish <- doubs$fish[keep, colSums(doubs$fish[keep, ]) > 0]
env  <- doubs$env[keep, ]
cat("after removing empty rows:", nrow(fish), "sites x", ncol(fish), "species\n")

## ============================================================================
##  2. RUNNING IT, AND THE TWO SCALINGS
## ============================================================================
##  vegan::cca() with a single table is a CA. Results can be presented two ways:
##
##  scaling = 1 (by sites)   chi-square distance is preserved among SITES.
##      Sites near each other have similar relative species composition; a site
##      near a species probably has a high proportion of it.
##  scaling = 2 (by species) chi-square distance is preserved among SPECIES.
##      Species near each other occur at similar relative frequencies; a species
##      near a site is probably frequent there.
##
##  The scaling changes the eigenVECTORS, never the eigenVALUES.
f.ca <- vegan::cca(fish)
print(summary(f.ca))                 # default scaling 2
print(summary(f.ca, scaling = 1))

ev <- f.ca$CA$eig
cat("\neigenvalues (identical under either scaling):\n"); print(round(ev, 4))
cat("total inertia:", round(sum(ev), 4), "\n")
cat("In CA an eigenvalue above 0.6 signals a very strong gradient.\n")

## ---- how many axes? Kaiser and broken stick --------------------------------
##  evplot() is Francois Gillet's function from Numerical Ecology in R. The
##  original session source()d it from a local RCode folder; it is defined here
##  so the script has no external file to find.
evplot <- function(ev, main = "") {
  n   <- length(ev)
  bsm <- data.frame(j = seq_len(n), p = 0)
  bsm$p[1] <- 1 / n
  for (i in 2:n) bsm$p[i] <- bsm$p[i - 1] + (1 / (n + 1 - i))
  bsm$p <- 100 * bsm$p / n
  par(mfrow = c(2, 1))
  barplot(ev, main = paste("Eigenvalues", main), col = "bisque", las = 2)
  abline(h = mean(ev), col = "red")
  legend("topright", "Average eigenvalue", lwd = 1, col = 2, bty = "n")
  barplot(t(cbind(100 * ev / sum(ev), bsm$p[n:1])), beside = TRUE,
          main = "% variation", col = c("bisque", 2), las = 2)
  legend("topright", c("% eigenvalue", "Broken stick model"),
         pch = 15, col = c("bisque", 2), bty = "n")
}
with_fig("08_ca_evplot", evplot(ev, "(Doubs fish)"), width = 8, height = 7)

## ============================================================================
##  3. THE BIPLOTS
## ============================================================================
with_fig("08_ca_biplot", {
  par(mfrow = c(1, 2))
  plot(f.ca, scaling = 1, main = "Scaling 1: sites are centroids of species")
  plot(f.ca, scaling = 2, main = "Scaling 2: species are centroids of sites")
}, width = 12, height = 6)

##  Axis 1 separates sites 19-30, the downstream ones, and most species sit near
##  them -- they are frequent there. Axis 2 picks out intermediate sites and the
##  species characteristic of them. Interpretation needs care: by construction
##  the site scores and species scores maximise their correlation, so proximity
##  on the plot is a statement about relative frequency, not abundance.

## ---- overlaying the environment --------------------------------------------
##  envfit() regresses each environmental variable on the ordination axes and
##  tests it by permutation, which turns "the axis looks like a gradient" into
##  a number.
set.seed(1998)
spe.ca.env <- vegan::envfit(f.ca, env, permutations = 999)
print(spe.ca.env)

with_fig("08_ca_envfit", {
  plot(f.ca, scaling = 2, main = "CA with environmental vectors")
  plot(spe.ca.env)
})

## ============================================================================
##  4. THE ARCH EFFECT, AND DETRENDING
## ============================================================================
##  A single strong gradient makes CA fold the second axis into a horseshoe --
##  the arch effect. It is an artefact of the method, not structure in the data.
##  Detrended correspondence analysis (decorana) removes it by flattening the
##  second axis in segments.
f.dca <- vegan::decorana(fish)
print(f.dca)

with_fig("08_ca_vs_dca", {
  par(mfrow = c(1, 2))
  plot(f.ca,  main = "CA: note the arch")
  plot(f.dca, main = "Detrended (decorana)")
}, width = 12, height = 6)

## ============================================================================
##  5. THE SAME ANALYSIS ON THE BUTTERFLIES
## ============================================================================
##  Sites x wing patterns, from session 04. Here the gradient is succession, and
##  the arch is almost entirely a consequence of it.
qroo3 <- get_data("butterflies")
sites <- c("HD", "SD", "GA", "YA", "MA", "OA", "PF")
qroo2 <- aggregate(qroo3[, sites], by = list(Pattern = qroo3$Pattern), FUN = sum)
qroo  <- as.matrix(t(qroo2[, -1]))
colnames(qroo) <- as.character(qroo2$Pattern)
qroo  <- qroo[rowSums(qroo) > 0, colSums(qroo) > 0]
cat("\nbutterfly table:", nrow(qroo), "sites x", ncol(qroo), "patterns\n")

qr.ca  <- vegan::cca(qroo)
qr.dca <- vegan::decorana(qroo)
with_fig("08_ca_butterflies_evplot", evplot(qr.ca$CA$eig, "(Quintana Roo)"),
         width = 8, height = 7)

with_fig("08_ca_butterflies", {
  par(mfrow = c(1, 2))
  plot(qr.ca,  main = "CA: the arch is the successional gradient")
  plot(qr.dca, main = "Detrended")
}, width = 12, height = 6)

cat("\n[08_CorrespondenceAnalysis] CA is PCA under chi-square distance, it",
    "ordinates rows and columns simultaneously, and its arch is an artefact",
    "of a dominant gradient rather than a second gradient.\n")
