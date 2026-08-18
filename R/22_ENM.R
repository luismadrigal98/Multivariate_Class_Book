## ============================================================================
##  MS_LJMR :: 22_ENM.R — Ecological Niche Modeling (bases)
##
##  Original authors: Marlon E. Cobos, Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 21_ENMI.R / 22_ENMII.R. The heavy spatial pipeline (terra,
##  dismo, worldclim downloads) is replaced by a self-contained binomial ENM on
##  the frozen 'hanta' presence/absence data with two bioclim predictors, so the
##  core idea runs without GIS rasters or network access.
## ============================================================================
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

occ <- get_data("hanta")                  # Sp (0/1), bio_1, bio_12
cat("Detections:", sum(occ$Sp), "of", nrow(occ), "sites\n")

## ---- ENMs are binomial(-like) models of suitability ------------------------
enm1 <- glm(Sp ~ bio_1,            data = occ, family = binomial)
enm2 <- glm(Sp ~ bio_1 + bio_12,   data = occ, family = binomial)
enmq <- glm(Sp ~ poly(bio_1, 2) + bio_12, data = occ, family = binomial)
cat("\nAIC (linear1, linear2, quadratic):\n")
print(AIC(enm1, enm2, enmq))
cat("\nBest-model coefficients:\n"); print(round(coef(enm2), 5))

## ---- discrimination: AUC on the training data ------------------------------
auc <- function(p, y) {                    # rank-based AUC, no extra packages
  r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
cat("Training AUC (bio_1 + bio_12):", round(auc(fitted(enm2), occ$Sp), 3), "\n")

## ---- response curve + suitability surface ----------------------------------
with_fig("22_enm", width = 10, height = 4.4, {
  op <- par(mfrow = c(1, 2)); on.exit(par(op))
  b1 <- seq(-5, 30, length = 100)
  pr <- predict(enm2, data.frame(bio_1 = b1, bio_12 = median(occ$bio_12)),
                type = "response")
  plot(b1, pr, type = "l", lwd = 2, col = "#1f3b73",
       xlab = "bio_1 (temperature)", ylab = "P(presence)",
       main = "ENM response curve")
  gx <- seq(-5, 30, length = 80); gy <- seq(100, 2500, length = 80)
  gg <- outer(gx, gy, function(a, b)
    predict(enm2, data.frame(bio_1 = a, bio_12 = b), type = "response"))
  image(gx, gy, gg, col = hcl.colors(20, "YlGnBu", rev = TRUE),
        xlab = "bio_1", ylab = "bio_12", main = "Suitability surface")
  points(occ$bio_1[occ$Sp == 1], occ$bio_12[occ$Sp == 1], pch = 19,
         cex = .4, col = "#8c2d3a")
})

cat("\n[22_ENM] a niche model is a binomial response surface over environmental",
    "predictors; evaluate it with AIC (fit) and AUC (discrimination).\n")
