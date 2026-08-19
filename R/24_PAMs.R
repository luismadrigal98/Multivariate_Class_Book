## ============================================================================
##  MS_LJMR :: 24_PAMs.R — Presence-Absence Matrices (biodiversity structure)
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
##  As shipped, the source() line below borrows three helpers from the course
##  repository: get_data() (loads a data set), need() (loads packages) and
##  with_fig() (opens a plot device). To run this script entirely on its own,
##  delete that line and uncomment the block below. Nothing else changes.
##
##  The standalone with_fig() just draws each figure to the screen, one after
##  the other. To save them as files instead, replace its body with
##      png(paste0(name, ".png")); on.exit(dev.off()); force(expr)
##
##  Files to keep next to this script: pam.csv
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) read.csv("pam.csv")
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

pamdf <- get_data("pam")
lat   <- pamdf$lat
pam   <- as.matrix(pamdf[, -1])            # sites x species, 0/1

## ---- marginal sums are the fundamental biodiversity quantities -------------
alpha <- rowSums(pam)                       # richness per site
omega <- colSums(pam)                       # range size per species
cat("PAM:", nrow(pam), "sites x", ncol(pam), "species;",
    "fill =", round(mean(pam), 3), "\n")
cat("mean alpha (richness/site) =", round(mean(alpha), 2),
    "   mean omega (range/species) =", round(mean(omega), 2), "\n")

with_fig("24_pam", width = 13, height = 3.8, {
  op <- par(mfrow = c(1, 3)); on.exit(par(op))
  image(t(pam[order(lat), ]), col = c("white", "#1f3b73"),
        main = "PAM (sites x species)", axes = FALSE)
  plot(lat, alpha, type = "l", col = "#1f3b73", lwd = 2,
       xlab = "latitude", ylab = "richness (alpha)", main = "Richness gradient")
  hist(omega, breaks = 20, col = "#2a7f7f", border = "white",
       xlab = "range size (omega)", main = "Range-size frequency")
})

## ---- dispersion field: the beta-diversity relationship ---------------------
## fill-corrected richness vs. mean co-occurrence (Soberon's "range-diversity")
betty  <- 1 / mean(alpha / ncol(pam))
disp   <- (pam %*% omega) / (alpha * nrow(pam))
with_fig("24_range_diversity", {
  plot(disp, alpha / ncol(pam), pch = 19, col = "#8c2d3a55",
       xlab = "dispersion field (mean range of co-occurring species)",
       ylab = "proportional richness",
       main = "Range-diversity plot")
})

cat("\n[24_PAMs] a PAM's row/column sums (alpha, omega) drive richness and",
    "range-size patterns; range-diversity plots relate the two.\n")
