## ============================================================================
##  MS_LJMR :: 20_ENM_II.R — Ecological Niche Modeling II: evaluation & transfer
##
##  Original authors: Marlon E. Cobos, Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 22_ENMII.R. The original session was about downloading and
##  cleaning rasters with dismo/geodata/rgl; that work is real but it is data
##  engineering, and it cannot run without network access or a GIS stack. The
##  statistically substantive half of niche modelling — how you know a model is
##  any good, and when you are allowed to project it somewhere else — is what
##  is kept and expanded here, on the same 'hanta' data as session I.
##
##  Three questions
##  ---------------
##  1. Does the model DISCRIMINATE on data it has not seen?  (cross-validated AUC)
##  2. Where do you cut a continuous suitability into presence/absence?  (thresholds)
##  3. Is the projection an interpolation or an extrapolation?  (novel environments)
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/20_ENM_II.R by scripts/make_class_src.py --
##  edit R/20_ENM_II.R and regenerate; changes made here will be overwritten.
## ============================================================================

# Working directory
#   Point this at the folder that holds this session's data. Every file
#   name below is resolved relative to it, so the script and its data have
#   to travel together -- or at least stay in step.
#
#   This session reads: hanta_virtual.csv
setwd("YOUR/DIRECTORY")

# Area under the ROC curve, from the rank-sum identity (no extra package)
auc <- function(p, y) {
  r <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ############################################################################
##  PART A -- MUCH OF NICHE MODELLING IS DATA PREPARATION
## ############################################################################
##  The original session downloaded WorldClim layers and pulled occurrences from
##  GBIF. Both need a live network, so the demonstrations below use the example
##  data bundled with dismo -- the same workflow, reproducible offline. The
##  network versions are kept as commented recipes where they belong.
have_dismo <- requireNamespace("dismo", quietly = TRUE) && requireNamespace("raster", quietly = TRUE)

if (have_dismo) {
  ## dismo defines S4 plot/predict methods, so these two must be ATTACHED for
  ## dispatch to find them -- dismo::plot(model) alone falls through to
  ## plot.default and errors on the model-specific arguments.
  library(dismo)
  library(raster)

  ## ---- a clean file, to start with ------------------------------------------
  ##  Bradypus variegatus, the brown-throated sloth: longitude/latitude only,
  ##  already tidied. Real data rarely arrives like this.
  bv2 <- read.csv(system.file("ex/bradypus.csv", package = "dismo"))
  bv  <- bv2[, 2:3]
  cat("Bradypus records:", nrow(bv), "\n"); print(head(bv2))

  par(mfrow = c(1, 1))
  if (requireNamespace("maps", quietly = TRUE)) maps::map("world", xlim = c(-100, -35), ylim = c(-40, 25),
                                 col = "grey70", fill = TRUE, bg = "white") else plot(bv, type = "n")
  points(bv, pch = 19, col = "orange", cex = .75)
  points(bv, cex = 1.2)
  title("Bradypus variegatus")


  ## ---- and a messy one -------------------------------------------------------
  ##  Solanum acaule as it comes out of a herbarium aggregator: duplicates,
  ##  missing coordinates, and points in impossible places. This is what a raw
  ##  GBIF download looks like.
  ##
  ##  The live version:
  ##      ep <- dismo::gbif("Echinacea", "purpurea*", geo = TRUE)   # * = wildcard
  ##      vm <- dismo::gbif("Vespa", "mandarinia", geo = TRUE)      # killer hornet
  ##  Always record the DOI from the GBIF portal for a citable, repeatable query.
  ac <- read.csv(system.file("ex/acaule.csv", package = "dismo"))
  cat("\nraw records:", nrow(ac), "\n")
  cat("countries:", paste(head(sort(unique(na.omit(ac$country))), 8),
                          collapse = ", "), "...\n")

  xy3 <- cbind(ac$lon, ac$lat)
  xy2 <- na.omit(xy3)
  xy  <- unique(xy2)
  cat("after dropping NA coordinates:", nrow(xy2),
      "  after removing duplicates:", nrow(xy), "\n")

  ## Two separate figures rather than one two-panel one: maps::map() fixes the
  ## aspect ratio to the projection, so a world map squeezed into half a device
  ## demands a plot region larger than the device and errors out.
  par(mfrow = c(1, 1))
  # Wide figure: widen the Plot pane, or open a sized device first --
  #   dev.new(width = 11, height = 6)
  if (requireNamespace("maps", quietly = TRUE)) maps::map("world", col = "grey75") else plot(xy, type = "n", xlab = "lon", ylab = "lat")
  points(xy, pch = 19, col = "red", cex = .6)
  title("Everything, worldwide")


  par(mfrow = c(1, 1))
  if (requireNamespace("maps", quietly = TRUE)) maps::map("world", xlim = c(-85, -55), ylim = c(-40, 5),
                                 col = "grey75") else plot(xy, type = "n", xlab = "lon", ylab = "lat")
  points(xy, pch = 19, col = "orange", cex = .8); points(xy, cex = 1.1)
  title("The intended range")


  cat("\nRecords far outside the known range are the point of the exercise:\n",
      "  they may be cultivated material, transposed coordinates, a country\n",
      "  centroid, or simply lon/lat swapped. See Hijmans & Elith (2017) sec 2.3\n",
      "  and Chapman (2005) on cleaning primary occurrence data.\n")

  ## ############################################################################
  ##  PART B -- PREDICTORS
  ## ############################################################################
  ##  dismo ships a small set of bioclimatic layers, which stand in for the
  ##  WorldClim download:
  ##      bios <- geodata::worldclim_global(var = "bio", res = 10, path = ".")
  ##  BIO5 max temperature of the warmest month | BIO6 min temp of the coldest
  ##  BIO13 precipitation of the wettest month  | BIO14 of the driest
  fnames <- list.files(system.file("ex", package = "dismo"), "grd$",
                       full.names = TRUE)
  bios   <- stack(fnames)
  names(bios) <- sub("\\.grd$", "", basename(fnames))
  print(bios)

  par(mfrow = c(1, 1))
  plot(bios)

  ## ---- correlated predictors can still define a niche ------------------------
  ##  Minimum and maximum temperature are strongly correlated, as you would
  ##  expect. That does not make them useless: the SCATTER around the line is
  ##  where the niche lives.
  pts <- rasterToPoints(bios[["bio5"]])
  set.seed(1998)
  iran <- sample(seq_len(nrow(pts)), min(10000, nrow(pts)))
  b5P  <- extract(bios[["bio5"]], pts[iran, 1:2])
  b6P  <- extract(bios[["bio6"]], pts[iran, 1:2])
  ok   <- complete.cases(b5P, b6P)
  mod  <- lm(b6P[ok] ~ b5P[ok])
  print(summary(mod))

  par(mfrow = c(1, 1))
  plot(b5P[ok], b6P[ok], pch = ".", col = "grey40",
       xlab = "BIO5 (max temp, warmest month)",
       ylab = "BIO6 (min temp, coldest month)",
       main = sprintf("r = %.2f, R2 = %.2f -- and still room for a niche",
                      cor(b5P[ok], b6P[ok]), summary(mod)$r.squared))
  abline(mod, col = "red", lwd = 2)


  ## ---- the species in E-space --------------------------------------------------
  st  <- bios[[c("bio5", "bio6", "bio12", "bio16")]]
  vrs <- na.omit(extract(st, bv))
  cat("\npredictors at the occurrence points:", nrow(vrs), "usable records\n")
  print(head(vrs))

  par(mfrow = c(1, 1))
  pairs(vrs, pch = 19, cex = .4, col = "#1f3b73",
        main = "Bradypus in environmental space")


  if (requireNamespace("car", quietly = TRUE)) {
    par(mfrow = c(1, 1))
    car::scatterplotMatrix(vrs, pch = 19, col = "blue", regLine = FALSE,
                           smooth = FALSE,
                           ellipse = list(levels = c(0.95, 0.999),
                                          robust = FALSE, fill = TRUE))

  }

  ## ############################################################################
  ##  PART C -- BIOCLIM: THE VENERABLE ENVELOPE
  ## ############################################################################
  ##  bioclim scores a site by where each of its variables falls in the
  ##  percentile distribution of the occurrence points -- a box, not a surface.
  ##  It is presence-only and it is old, which makes it a good baseline.
  bio_vm <- bioclim(st, bv)
  # Wide figure: widen the Plot pane, or open a sized device first --
  #   dev.new(width = 11, height = 5.5)
  op <- par(mfrow = c(1, 2))
  plot(bio_vm, p = 0.85, main = "85% envelope")
  plot(bio_vm, p = 0.95, main = "95% envelope")
  par(op)

  pred_vm <- predict(st, bio_vm, progress = "")
  par(mfrow = c(1, 1))
  plot(pred_vm, main = "bioclim suitability (raw)")
  if (requireNamespace("maps", quietly = TRUE)) maps::map("world", add = TRUE, col = "grey40")


  ## ---- from continuous to binary: pick a threshold ---------------------------
  qq <- quantile(pred_vm, probs = (1:20) / 20, na.rm = TRUE)
  cat("\nquantiles of the prediction:\n"); print(round(qq, 4))

  # Wide figure: widen the Plot pane, or open a sized device first --
  #   dev.new(width = 13, height = 6)
  op <- par(mfrow = c(1, 2))
  plot(pred_vm > qq[18], main = "threshold = top 90%")
  if (requireNamespace("maps", quietly = TRUE)) maps::map("world", add = TRUE, col = "grey40")
  plot(pred_vm > qq[19], main = "threshold = top 95%")
  if (requireNamespace("maps", quietly = TRUE)) maps::map("world", add = TRUE, col = "grey40")
  par(op)

  ## ############################################################################
  ##  PART D -- ELLIPSOIDS: THE MAHALANOBIS MODEL
  ## ############################################################################
  ##  mahal() measures distance to the centroid of the occurrence cloud in units
  ##  of its own covariance -- the ellipsoid of session 07, used as a niche
  ##  model. It respects correlations among predictors, which a bioclim box
  ##  cannot.
  maha_vm <- mahal(st, bv)
  ext     <- extent(-90, -32, -35, 15)          # restrict: it is slow
  pred_mh <- predict(crop(st, ext), maha_vm, progress = "")

  # Wide figure: widen the Plot pane, or open a sized device first --
  #   dev.new(width = 13, height = 6)
  op <- par(mfrow = c(1, 2))
  plot(pred_vm, ext = ext, main = "bioclim (box)")
  if (requireNamespace("maps", quietly = TRUE)) maps::map("world", add = TRUE, col = "grey40")
  plot(pred_mh, main = "Mahalanobis (ellipsoid)")
  if (requireNamespace("maps", quietly = TRUE)) maps::map("world", add = TRUE, col = "grey40")
  par(op)
} else {
  message("Skipped: the dismo workflow (bioclim, Mahalanobis, raster projection).  Install with: install.packages(c(\"dismo\", \"raster\"))")
}

## ############################################################################
##  PART E -- EVALUATION AND TRANSFER
## ############################################################################
##  Back to the Hantavirus data of session I. Three questions decide whether a
##  model is worth anything:
##    1. does it DISCRIMINATE on data it has not seen?  (cross-validated AUC)
##    2. where do you cut suitability into presence/absence?  (thresholds)
##    3. is a projection interpolation or extrapolation?  (novel environments)
occ  <- read.csv("hanta_virtual.csv")
form <- Sp ~ poly(bio_1, 2) + bio_12
full <- glm(form, data = occ, family = binomial)

cat("\nSites:", nrow(occ), " detections:", sum(occ$Sp), "\n")
cat("Training AUC (optimistic):", round(auc(fitted(full), occ$Sp), 3), "\n")

## ---- 1. cross-validated discrimination --------------------------------------
set.seed(1998)
k    <- 5
fold <- sample(rep_len(seq_len(k), nrow(occ)))
cv   <- vapply(seq_len(k), function(i) {
  m <- glm(form, data = occ[fold != i, ], family = binomial)
  auc(predict(m, occ[fold == i, ], type = "response"), occ$Sp[fold == i])
}, numeric(1))
cat("\n", k, "-fold CV AUC per fold: ", paste(round(cv, 3), collapse = ", "),
    "\n  mean = ", round(mean(cv), 3), "  sd = ", round(sd(cv), 3), "\n", sep = "")

## These folds are RANDOM, so a held-out site usually has a training site a few
## kilometres away. With spatially autocorrelated environments that makes the
## test easy and the AUC optimistic even after cross-validation. Spatially
## blocked partitions give the honest number for a model meant to be projected.

## ---- 2. from continuous suitability to a binary prediction -------------------
p    <- fitted(full)
pres <- p[occ$Sp == 1]
thr  <- c(
  "minimum training presence"   = min(pres),
  "10th percentile presence"    = unname(quantile(pres, 0.10)),
  "max sensitivity+specificity" = {
    cand <- sort(unique(round(p, 3)))
    ss <- vapply(cand, function(t)
      mean(p[occ$Sp == 1] >= t) + mean(p[occ$Sp == 0] < t), numeric(1))
    cand[which.max(ss)]
  })

confusion <- function(t) {
  pred <- p >= t
  c(sensitivity = mean(pred[occ$Sp == 1]),
    specificity = mean(!pred[occ$Sp == 0]),
    omission    = mean(!pred[occ$Sp == 1]),
    area        = mean(pred))
}
tab <- t(vapply(thr, confusion, numeric(4)))
cat("\nThreshold rules (threshold, then consequences):\n")
print(round(cbind(threshold = thr, tab), 3))

## Minimum-training-presence admits every known detection (omission = 0) at the
## cost of calling a large area suitable; the 10th-percentile rule discards the
## most marginal tenth of records -- often the misidentified or mis-georeferenced
## ones -- and halves the predicted area. Neither is correct in the abstract;
## failing to say which you used is not defensible.

op <- par(mfrow = c(1, 2))
ord <- order(-p)
tpr <- cumsum(occ$Sp[ord] == 1) / sum(occ$Sp == 1)
fpr <- cumsum(occ$Sp[ord] == 0) / sum(occ$Sp == 0)
plot(fpr, tpr, type = "l", lwd = 2, col = "#1f3b73",
     xlab = "false positive rate", ylab = "true positive rate",
     main = paste0("ROC (AUC = ", round(auc(p, occ$Sp), 3), ")"))
abline(0, 1, lty = 2, col = "grey60")
ts <- seq(0, 1, length = 200)
plot(ts, vapply(ts, function(t) mean(p >= t), numeric(1)), type = "l", lwd = 2,
     col = "#1f3b73", xlab = "threshold",
     ylab = "fraction of sites called suitable",
     main = "Predicted area vs. threshold")
abline(v = thr, col = c("#8c2d3a", "#2a7f7f", "grey30"), lty = 2)
legend("topright", names(thr), lty = 2, bty = "n", cex = .7,
       col = c("#8c2d3a", "#2a7f7f", "grey30"))
par(op)

## ---- 3. transfer, and how much of it is extrapolation ------------------------
future <- transform(occ, bio_1 = bio_1 + 3)
p_now  <- predict(full, occ,    type = "response")
p_fut  <- predict(full, future, type = "response")
tcut   <- thr[["10th percentile presence"]]
cat("\nSuitable fraction at the 10th-percentile threshold:",
    "\n  current  =", round(mean(p_now >= tcut), 3),
    "\n  +3 C     =", round(mean(p_fut >= tcut), 3), "\n")

in_range <- function(v, ref) v >= min(ref) & v <= max(ref)
novel    <- !(in_range(future$bio_1, occ$bio_1) & in_range(future$bio_12, occ$bio_12))
cat("Projection points in novel environments:", round(100 * mean(novel), 1), "%\n")

par(mfrow = c(1, 1))
plot(occ$bio_1, p_now, pch = 19, cex = .3, col = "grey70",
     xlab = "bio_1 (temperature)", ylab = "predicted suitability",
     main = "Current vs +3 C, novel environments flagged")
points(future$bio_1[!novel], p_fut[!novel], pch = 19, cex = .3, col = "#1f3b73")
points(future$bio_1[novel],  p_fut[novel],  pch = 4,  cex = .5, col = "#8c2d3a")
abline(h = tcut, lty = 2, col = "grey30")
legend("topleft", c("current", "+3 C, within training range",
                    "+3 C, EXTRAPOLATED"),
       pch = c(19, 19, 4), col = c("grey70", "#1f3b73", "#8c2d3a"),
       bty = "n", cex = .75)

cat("\n[20_ENM_II] most of niche modelling is cleaning data; cross-validate",
    "before believing an AUC; state the threshold rule you used; and never",
    "report a projection without saying how much of it is extrapolation.\n")
