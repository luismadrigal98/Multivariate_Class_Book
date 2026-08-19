## ============================================================================
##  MS_LJMR :: 23_ENM_II.R — Ecological Niche Modeling II: evaluation & transfer
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
# auc <- function(p, y) {                 # rank-based AUC, no extra package
#   r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0)
#   (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
# }
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

occ <- get_data("hanta")
form <- Sp ~ poly(bio_1, 2) + bio_12       # the model selected in session I
full <- glm(form, data = occ, family = binomial)

cat("Sites:", nrow(occ), " detections:", sum(occ$Sp), "\n")
cat("Training AUC (optimistic):", round(auc(fitted(full), occ$Sp), 3), "\n")

## ---- 1. Cross-validated discrimination -------------------------------------
## Training AUC is measured on the points that chose the coefficients, so it is
## biased upward. k-fold cross-validation refits the model k times, each time
## scoring the fold that was held out.
set.seed(1998)
k     <- 5
fold  <- sample(rep_len(seq_len(k), nrow(occ)))
cv    <- vapply(seq_len(k), function(i) {
  tr <- occ[fold != i, ]; te <- occ[fold == i, ]
  m  <- glm(form, data = tr, family = binomial)
  auc(predict(m, te, type = "response"), te$Sp)
}, numeric(1))
cat("\n", k, "-fold CV AUC per fold: ", paste(round(cv, 3), collapse = ", "),
    "\n  mean = ", round(mean(cv), 3), "  sd = ", round(sd(cv), 3), "\n", sep = "")

## A caution worth stating in class: these folds are RANDOM, so a held-out site
## usually has a training site a few kilometres away. With spatially
## autocorrelated environments that makes the test easy and the AUC optimistic
## even after cross-validation. Spatially blocked partitions (contiguous
## regions held out together) give the honest number for a model meant to be
## projected to new areas.

## ---- 2. From continuous suitability to a binary prediction -----------------
## Nothing in the model says where "suitable" begins. Three standard rules:
p    <- fitted(full)
pres <- p[occ$Sp == 1]
thr  <- c(
  "minimum training presence" = min(pres),
  "10th percentile presence"   = unname(quantile(pres, 0.10)),
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

## The trade-off is the whole point: the minimum-training-presence rule admits
## every known detection (omission = 0) at the cost of calling a large area
## suitable, while the 10th-percentile rule discards the 10% most marginal
## records — usually the ones most likely to be identification or georeference
## errors — and shrinks the predicted area sharply.

with_fig("23_enm_thresholds", {
  op <- par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  ## ROC curve
  ord <- order(-p)
  tpr <- cumsum(occ$Sp[ord] == 1) / sum(occ$Sp == 1)
  fpr <- cumsum(occ$Sp[ord] == 0) / sum(occ$Sp == 0)
  plot(fpr, tpr, type = "l", lwd = 2, col = "#1f3b73",
       xlab = "false positive rate", ylab = "true positive rate",
       main = paste0("ROC (AUC = ", round(auc(p, occ$Sp), 3), ")"))
  abline(0, 1, lty = 2, col = "grey60")
  ## predicted area as a function of the threshold
  ts <- seq(0, 1, length = 200)
  plot(ts, vapply(ts, function(t) mean(p >= t), numeric(1)), type = "l", lwd = 2,
       col = "#1f3b73", xlab = "threshold", ylab = "fraction of sites called suitable",
       main = "Predicted area vs. threshold")
  abline(v = thr, col = c("#8c2d3a", "#2a7f7f", "grey30"), lty = 2)
  legend("topright", names(thr), lty = 2, bty = "n", cex = .7,
         col = c("#8c2d3a", "#2a7f7f", "grey30"))
  par(op)
}, width = 10, height = 4.6)

## ---- 3. Transfer: projecting to a scenario the model never saw -------------
## Projection is where niche models earn their keep and where they fail most
## quietly. Shift temperature by +3 C and re-score the same sites.
future <- transform(occ, bio_1 = bio_1 + 3)
p_now  <- predict(full, occ,    type = "response")
p_fut  <- predict(full, future, type = "response")
tcut   <- thr[["10th percentile presence"]]
cat("\nSuitable fraction at the 10th-percentile threshold:",
    "\n  current  =", round(mean(p_now >= tcut), 3),
    "\n  +3 C     =", round(mean(p_fut >= tcut), 3), "\n")

## ---- 4. Are we extrapolating? ----------------------------------------------
## A prediction outside the range of the training data is an extrapolation, and
## a quadratic term extrapolates violently. This is the idea behind MESS
## (multivariate environmental similarity surface): flag every projection point
## that falls outside the training envelope of ANY predictor.
in_range <- function(v, ref) v >= min(ref) & v <= max(ref)
novel    <- !(in_range(future$bio_1, occ$bio_1) & in_range(future$bio_12, occ$bio_12))
cat("Projection points in novel environments:",
    round(100 * mean(novel), 1), "%\n")

with_fig("23_enm_transfer", {
  plot(occ$bio_1, p_now, pch = 19, cex = .3, col = "grey70",
       xlab = "bio_1 (temperature)", ylab = "predicted suitability",
       main = "Current vs. +3 C, novel environments flagged")
  points(future$bio_1[!novel], p_fut[!novel], pch = 19, cex = .3, col = "#1f3b73")
  points(future$bio_1[novel],  p_fut[novel],  pch = 4,  cex = .5, col = "#8c2d3a")
  abline(h = tcut, lty = 2, col = "grey30")
  legend("topleft", c("current", "+3 C, within training range",
                      "+3 C, EXTRAPOLATED"),
         pch = c(19, 19, 4), col = c("grey70", "#1f3b73", "#8c2d3a"),
         bty = "n", cex = .75)
})

cat("\n[23_ENM_II] cross-validate before believing an AUC; state the threshold",
    "rule you used; and never report a projection without saying how much of",
    "it is extrapolation.\n")
