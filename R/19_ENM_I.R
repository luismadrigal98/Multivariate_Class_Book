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
##  Files to keep next to this script: hanta_virtual.csv; vars.tif
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
# get_data <- function(name) read.csv("hanta_virtual.csv")
# data_dir <- function() "."                       # look for data files here
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

## ============================================================================
##  1. THE DATA
## ============================================================================
##  Detection / non-detection records from tests for Hantavirus, with two
##  bioclimatic predictors: bio_1 (mean annual temperature) and bio_12 (annual
##  precipitation). The raster vars.tif carries the same two variables over the
##  continental USA, so a model fitted in environmental space can be projected
##  back into geography.
occ <- get_data("hanta")
cat("Detections:", sum(occ$Sp), "of", nrow(occ), "sites",
    "(prevalence =", round(mean(occ$Sp), 3), ")\n")
print(head(occ))
cat("\ncolumns:", paste(colnames(occ), collapse = ", "),
    "\n  Sp is the response; the rest are predictors\n")

env_path <- file.path(data_dir(), "vars.tif")
have_rast <- file.exists(env_path) && has_pkg("terra")
if (have_rast) {
  env_vars <- terra::rast(env_path)
  print(env_vars)
  with_fig("19_env_layers", terra::plot(env_vars), width = 10, height = 5)
} else {
  skip_note("geographic projection onto vars.tif", "terra")
}

## ============================================================================
##  2. EXPLORATION IN ENVIRONMENTAL SPACE
## ============================================================================
with_fig("19_enm_espace", {
  plot(occ[, 2:3], pch = 16, cex = .35, col = "grey75",
       xlab = "bio_1 (mean annual temperature)",
       ylab = "bio_12 (annual precipitation)",
       main = "Environmental space: background vs. detections")
  points(occ[occ$Sp == 1, 2:3], pch = 19, cex = .5, col = "red")
  legend("topright", c("all sites", "detections"), pch = 19,
         col = c("grey75", "red"), bty = "n")
})

##  The records alone do not say whether the species is picky -- that depends on
##  what was AVAILABLE. Plotting every cell of the continental USA behind them
##  supplies the missing denominator.
if (have_rast) {
  bg_cells <- na.omit(as.data.frame(env_vars))
  with_fig("19_enm_espace_background", {
    plot(bg_cells, col = "gray", pch = ".",
         xlab = "bio_1", ylab = "bio_12",
         main = "All continental-US conditions, with the records on top")
    points(occ[, 2:3], pch = 16, cex = .3)
    points(occ[occ$Sp == 1, 2:3], pch = 19, cex = .45, col = "red")
    legend("topright", c("available (USA)", "sampled", "detections"),
           pch = c(20, 16, 19), col = c("gray", "black", "red"), bty = "n")
  })
}

cat("\nEnvironment at detections vs. everywhere sampled:\n")
print(round(rbind(detections = colMeans(occ[occ$Sp == 1, c("bio_1", "bio_12")]),
                  all        = colMeans(occ[, c("bio_1", "bio_12")])), 2))

## ============================================================================
##  3. A LADDER OF MODELS
## ============================================================================
##  An ENM is a statistical representation of "A" in the BAM diagram -- the
##  abiotically suitable region. It is fitted with a binomial or binomial-like
##  model in which 1 = presence or detection and 0 = absence, non-detection, or
##  background. Start simple and add structure one step at a time.
enmb1  <- glm(Sp ~ bio_1,  data = occ, family = binomial(link = "logit"))
enmb12 <- glm(Sp ~ bio_12, data = occ, family = binomial(link = "logit"))
print(summary(enmb1)); print(summary(enmb12))
cat("\none variable at a time:\n"); print(AIC(enmb1, enmb12))

enm2 <- glm(Sp ~ bio_1 + bio_12, data = occ, family = binomial(link = "logit"))
print(summary(enm2))
cat("\nadding the second variable:\n"); print(AIC(enmb1, enmb12, enm2))

## quadratic terms let the response be hump-shaped -- an OPTIMUM rather than a
## monotone preference, which is what Hutchinson's niche concept predicts
enm2q <- glm(Sp ~ bio_1 + bio_12 + I(bio_1^2) + I(bio_12^2),
             data = occ, family = binomial(link = "logit"))
print(summary(enm2q))

## an interaction lets the response to one variable depend on the other
enm2qp <- glm(Sp ~ bio_1 + bio_12 + I(bio_1^2) + I(bio_12^2) + bio_1:bio_12,
              data = occ, family = binomial(link = "logit"))
print(summary(enm2qp))
cat("\nthe whole ladder:\n"); print(AIC(enmb1, enmb12, enm2, enm2q, enm2qp))

## ============================================================================
##  4. THE SURPRISE, AND WHY IT MATTERS
## ============================================================================
##  The most complex model wins. Does that mean complexity always wins? Add a
##  variable that is pure noise -- a reshuffle of the real predictors, so it has
##  the same distribution and no relationship to Sp whatsoever -- and see.
set.seed(3)
occ1 <- cbind(occ, rvar = sample(c(occ$bio_1, occ$bio_12), nrow(occ)))
cat("\ncorrelation of the junk variable with the response:",
    round(cor(occ1$rvar, occ1$Sp), 4), "\n")

enm3l   <- glm(Sp ~ bio_1 + bio_12 + rvar + I(bio_1^2) + I(bio_12^2) +
                 bio_1:bio_12, data = occ1, family = binomial(link = "logit"))
enm3lq  <- glm(Sp ~ bio_1 + bio_12 + rvar + I(bio_1^2) + I(bio_12^2) +
                 I(rvar^2) + bio_1:bio_12, data = occ1,
               family = binomial(link = "logit"))
enm3lqp <- glm(Sp ~ bio_1 + bio_12 + rvar + I(bio_1^2) + I(bio_12^2) +
                 I(rvar^2) + bio_1:bio_12 + bio_1:rvar + rvar:bio_12,
               data = occ1, family = binomial(link = "logit"))
cat("\nwith a junk variable added:\n")
print(AIC(enmb1, enmb12, enm2, enm2q, enm2qp, enm3l, enm3lq, enm3lqp))

##  AIC penalises parameters, but only mildly: the junk models land close to the
##  honest one and can even edge past it. AIC ranks models, it does not tell you
##  whether a term MEANS anything. A sequential deviance table does:
cat("\ndeviance contributions, honest model:\n")
print(anova(enm2qp, test = "Chi"))
cat("\ndeviance contributions, with junk:\n")
print(anova(enm3lqp, test = "Chi"))
cat("\nThe junk terms contribute essentially no deviance, whatever AIC says.\n")

## ============================================================================
##  5. RESPONSE CURVES AND THE SUITABILITY SURFACE
## ============================================================================
##  A response curve varies ONE predictor with the others held at their median:
##  a slice through the surface, not the surface.
with_fig("19_enm_response", {
  par(mfrow = c(1, 2))
  b1   <- seq(min(occ$bio_1), max(occ$bio_1), length = 200)
  mods <- list("additive" = enm2, "quadratic + interaction" = enm2qp)
  for (nm in names(mods)) {
    pr <- predict(mods[[nm]],
                  data.frame(bio_1 = b1, bio_12 = median(occ$bio_12)),
                  type = "response")
    plot(b1, pr, type = "l", lwd = 2, col = "#1f3b73", ylim = c(0, 1),
         xlab = "bio_1 (temperature)", ylab = "P(presence)", main = nm)
    rug(occ$bio_1[occ$Sp == 1], col = "#8c2d3a")
  }
}, width = 10, height = 4.6)

grid_pred <- function(mod, n = 120) {
  gx <- seq(min(occ$bio_1),  max(occ$bio_1),  length = n)
  gy <- seq(min(occ$bio_12), max(occ$bio_12), length = n)
  z  <- outer(gx, gy, function(a, b)
    predict(mod, data.frame(bio_1 = a, bio_12 = b), type = "response"))
  list(x = gx, y = gy, z = z)
}
g <- grid_pred(enm2qp)
with_fig("19_enm_surface", {
  image(g$x, g$y, g$z, col = hcl.colors(24, "YlGnBu", rev = TRUE),
        xlab = "bio_1 (temperature)", ylab = "bio_12 (precipitation)",
        main = "Modelled suitability over environmental space")
  contour(g$x, g$y, g$z, add = TRUE, levels = c(.25, .5, .75), col = "grey30")
  points(occ$bio_1[occ$Sp == 1], occ$bio_12[occ$Sp == 1],
         pch = 19, cex = .4, col = "#8c2d3a")
})

## ============================================================================
##  6. BACK INTO GEOGRAPHY
## ============================================================================
##  Which predictors you use depends on the goal. Principal components often
##  predict better; RAW variables are what you need if you want to interpret the
##  response. This model uses raw variables, so the curves above mean something.
if (have_rast) {
  enmpre <- terra::predict(env_vars, enm2qp, type = "response")
  with_fig("19_enm_prediction_layers", {
    terra::plot(c(env_vars, enmpre))
  }, width = 10, height = 7)

  ##  The same prediction seen twice: as a map, and as a colouring of
  ##  environmental space. Points of the two figures correspond one to one.
  bgdf  <- as.data.frame(env_vars, na.rm = TRUE)
  suit  <- as.data.frame(enmpre, na.rm = TRUE)[, 1]
  cols  <- rev(terrain.colors(64))[cut(suit, 64, labels = FALSE)]

  with_fig("19_enm_geo_vs_env", {
    par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))
    terra::plot(enmpre, main = "Suitability for Hantavirus, continental USA")
    plot(bgdf, pch = 19, cex = .3, col = cols,
         xlab = "bio_1", ylab = "bio_12", main = "The same, in E-space")
    plot(bgdf, pch = 19, cex = .3, col = cols,
         xlab = "bio_1", ylab = "bio_12", main = "... with NEGATIVE records")
    points(occ[occ$Sp == 0, 2:3], col = "gray30", pch = "x", cex = .5)
    plot(bgdf, pch = 19, cex = .3, col = cols,
         xlab = "bio_1", ylab = "bio_12", main = "... with POSITIVE records")
    points(occ[occ$Sp == 1, 2:3], col = "purple", pch = "+", cex = 1)
  }, width = 12, height = 10)
}

cat("\n[19_ENM_I] a niche model is a binomial response surface over",
    "environmental predictors. AIC ranks models but does not detect nonsense;",
    "check each term's deviance. Session II asks whether it generalizes.\n")
