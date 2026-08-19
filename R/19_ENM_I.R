## ============================================================================
##  MS_LJMR :: 19_ENM_I.R — Ecological Niche Modeling I: fitting the model
##
##  Original authors: Marlon E. Cobos, Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 21_ENMI.R. The heavy spatial pipeline (terra, dismo, worldclim
##  downloads) is replaced by a self-contained binomial ENM on the real 'hanta'
##  occurrence table with two bioclim predictors, so the core idea runs without
##  GIS rasters or network access. Session II (20_ENM_II.R) takes the model
##  fitted here and asks whether it is any good.
##
##  The idea
##  --------
##  An ENM is a function from ENVIRONMENTAL space to suitability. Each site is
##  a point in the multivariate space of climate variables — exactly the data
##  matrix of Chapter 1 — and the model is a response surface over it. What
##  makes the problem ecological rather than merely statistical is that the
##  surface is then read back into GEOGRAPHIC space, where a suitable climate
##  may still be unoccupied because the species never got there (dispersal) or
##  cannot persist there (biotic interactions).
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
##  Files to keep next to this script: hanta_virtual.csv
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) read.csv("hanta_virtual.csv")
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

occ <- get_data("hanta")                  # Sp (0/1), bio_1, bio_12
cat("Detections:", sum(occ$Sp), "of", nrow(occ), "sites",
    "(prevalence =", round(mean(occ$Sp), 3), ")\n")

## ---- 1. The niche in environmental space -----------------------------------
## Before fitting anything: do the presences occupy a restricted region of the
## two-variable environmental space? If they scatter like the background, no
## model will help.
with_fig("19_enm_espace", {
  plot(occ$bio_1, occ$bio_12, pch = 19, cex = .35, col = "grey75",
       xlab = "bio_1 (mean annual temperature)",
       ylab = "bio_12 (annual precipitation)",
       main = "Environmental space: background vs. detections")
  points(occ$bio_1[occ$Sp == 1], occ$bio_12[occ$Sp == 1],
         pch = 19, cex = .5, col = "#8c2d3a")
  legend("topright", c("all sites", "detections"), pch = 19,
         col = c("grey75", "#8c2d3a"), bty = "n")
})

cat("\nEnvironment at detections vs. everywhere:\n")
print(round(rbind(
  detections = colMeans(occ[occ$Sp == 1, c("bio_1", "bio_12")]),
  all        = colMeans(occ[, c("bio_1", "bio_12")])), 2))

## ---- 2. Suitability as a binomial response surface -------------------------
## The response is 0/1, so the model is binomial with a logit link — the same
## GLM machinery as session 19, now with a spatial reading.
enm1 <- glm(Sp ~ bio_1,                    data = occ, family = binomial)
enm2 <- glm(Sp ~ bio_1 + bio_12,           data = occ, family = binomial)
enmq <- glm(Sp ~ poly(bio_1, 2) + bio_12,  data = occ, family = binomial)

cat("\nAIC (linear1, linear2, quadratic):\n")
print(AIC(enm1, enm2, enmq))
cat("\nCoefficients of the additive two-predictor model:\n")
print(round(coef(enm2), 5))

## A quadratic term in temperature is the difference between a monotone
## preference ("warmer is always better") and a genuine OPTIMUM — the shape
## Hutchinson's niche concept actually predicts. Compare the two curves below
## before deciding which model to carry into session II.
cat("\nLikelihood-ratio test, additive vs. quadratic in bio_1:\n")
print(anova(enm2, enmq, test = "Chisq"))

## ---- 3. Response curves ----------------------------------------------------
## A response curve varies ONE predictor with the others held at their median.
## It is a slice through the surface, not the surface.
with_fig("19_enm_response", {
  op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  b1   <- seq(min(occ$bio_1), max(occ$bio_1), length = 200)
  mods <- list("additive" = enm2, "quadratic in bio_1" = enmq)
  for (nm in names(mods)) {
    pr <- predict(mods[[nm]],
                  data.frame(bio_1 = b1, bio_12 = median(occ$bio_12)),
                  type = "response")
    plot(b1, pr, type = "l", lwd = 2, col = "#1f3b73", ylim = c(0, 1),
         xlab = "bio_1 (temperature)", ylab = "P(presence)", main = nm)
    rug(occ$bio_1[occ$Sp == 1], col = "#8c2d3a")
  }
  par(op)
}, width = 10, height = 4.6)

## ---- 4. The suitability surface --------------------------------------------
grid_pred <- function(mod, n = 120) {
  gx <- seq(min(occ$bio_1),  max(occ$bio_1),  length = n)
  gy <- seq(min(occ$bio_12), max(occ$bio_12), length = n)
  z  <- outer(gx, gy, function(a, b)
    predict(mod, data.frame(bio_1 = a, bio_12 = b), type = "response"))
  list(x = gx, y = gy, z = z)
}
g <- grid_pred(enmq)

with_fig("19_enm_surface", {
  image(g$x, g$y, g$z, col = hcl.colors(24, "YlGnBu", rev = TRUE),
        xlab = "bio_1 (temperature)", ylab = "bio_12 (precipitation)",
        main = "Modelled suitability over environmental space")
  contour(g$x, g$y, g$z, add = TRUE, levels = c(.25, .5, .75), col = "grey30")
  points(occ$bio_1[occ$Sp == 1], occ$bio_12[occ$Sp == 1],
         pch = 19, cex = .4, col = "#8c2d3a")
})

## Carry the fitted objects forward: 20_ENM_II.R re-fits from scratch so it can
## be run standalone, but the model formulas are the ones chosen here.
cat("\n[19_ENM_I] a niche model is a binomial response surface over",
    "environmental predictors. Session II asks whether it generalizes.\n")
