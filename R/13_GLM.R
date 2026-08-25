## ============================================================================
##  MS_LJMR :: 13_GLM.R — Generalized Linear Models (multivariate predictors)
##
##  Original authors: Marlon E. Cobos & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 19_GLModels.R. Uses frozen 'crawley' mixed-type table.
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
##  Files to keep next to this script: speciesCrawley3.csv
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
# get_data <- function(name) read.csv("speciesCrawley3.csv", stringsAsFactors = TRUE)
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

## ============================================================================
##  1. WHY GENERALIZED
## ============================================================================
##  A linear model says   Y ~ mu + X + error, and ordinary regression assumes
##  the error is (1) normal with mean 0 and (2) of constant variance across all
##  observations. A GLM relaxes both: the error may follow another distribution,
##  and its variance may depend on the mean. So you must state the FAMILY --
##  gaussian, Gamma, Poisson, binomial -- and the LINK that connects the linear
##  predictor to the mean.
m <- get_data("crawley")
m$Soil <- factor(m$Soil)
str(m)

with_fig("13_glm_pairs", {
  pairs(m[, 2:5], col = factor(m$Species), pch = 19, cex = .6,
        main = "speciesCrawley3: biomass, temperature, precipitation, pH")
})

## ---- is the response even close to normal? ----------------------------------
with_fig("13_glm_qq_raw", {
  par(mfrow = c(1, 2))
  qqnorm(m$Biomass, main = "Variable: Biomass"); qqline(m$Biomass)
  hist(m$Biomass, col = "grey85", main = "Biomass", xlab = "")
}, width = 10, height = 5)
cat("\nShapiro-Wilk on Biomass: p =", signif(shapiro.test(m$Biomass)$p.value, 4),
    " (small p => not normal)\n")

## ---- and are the RESIDUALS of an OLS fit normal? ----------------------------
##  This is the assumption that actually matters; the marginal distribution of
##  the response does not have to be normal for OLS to be valid.
mLM <- lm(Biomass ~ Tmp, data = m)
with_fig("13_glm_qq_resid", {
  qqnorm(residuals(mLM), main = "Residuals, ordinary least squares")
  qqline(residuals(mLM))
})
cat("Shapiro-Wilk on OLS residuals: p =",
    signif(shapiro.test(residuals(mLM))$p.value, 4), "\n")

## ============================================================================
##  2. CHOOSING A FAMILY
## ============================================================================
##  continuous and unbounded      -> gaussian (ordinary regression)
##  continuous and non-negative   -> Gamma, or inverse gaussian
##      (biomass, precipitation, concentrations -- anything with a hard floor)
##  counts in a fixed interval    -> poisson
##      poisson assumes mean == variance; when that fails use quasipoisson or
##      a negative binomial
##  binary / proportions          -> binomial (logistic regression)
##  more than two categories      -> multinomial
mGauss <- glm(Biomass ~ Tmp, data = m, family = gaussian)
mGamma <- glm(Biomass ~ Tmp, data = m, family = Gamma)

with_fig("13_glm_families", {
  par(mfrow = c(1, 2))
  o <- order(m$Tmp)
  plot(m$Tmp, m$Biomass, pch = 19, cex = .6, main = "Gaussian error",
       xlab = "Tmp", ylab = "Biomass")
  lines(m$Tmp[o], fitted(mGauss)[o], col = "#2a7f7f", lwd = 2)
  plot(m$Tmp, m$Biomass, pch = 19, cex = .6, main = "Gamma error",
       xlab = "Tmp", ylab = "Biomass")
  lines(m$Tmp[o], fitted(mGamma)[o], col = "#8c2d3a", lwd = 2)
}, width = 10, height = 5)

##  predict(type = "response") returns the fitted MEAN; without it you get the
##  linear predictor, which for Gamma is on the inverse scale.
cat("\nfirst fitted values, response vs link scale:\n")
print(round(head(cbind(response = predict(mGamma, type = "response"),
                       link     = predict(mGamma)), 4), 4))

with_fig("13_glm", {
  par(mfrow = c(1, 2))
  plot(mGamma, which = 1)      # residuals vs fitted: outliers and structure
  plot(mGamma, which = 2)      # normal Q-Q of deviance residuals
}, width = 10, height = 5)

## ---- adding a second predictor, and an interaction ---------------------------
mG_TP  <- glm(Biomass ~ Tmp + Precip, data = m, family = Gamma)
mG_TxP <- glm(Biomass ~ Tmp * Precip, data = m, family = Gamma)
print(summary(mG_TP))
print(summary(mG_TxP))
cat("\nAIC:\n"); print(AIC(mGamma, mG_TP, mG_TxP))
## (the original script compared an object `mBT` that was never created; the
##  three models above are the intended comparison)

## ============================================================================
##  3. A FACTOR PREDICTOR IS AN ANOVA
## ============================================================================
mBS <- glm(Biomass ~ Soil, data = m, family = gaussian)
print(summary(mBS))
with_fig("13_glm_soil_boxplot", {
  boxplot(Biomass ~ Soil, data = m, col = "grey85",
          main = "Biomass by soil type")
})

##  Exactly the same model, in the vocabulary of experimental design:
aovBS <- aov(Biomass ~ Soil, data = m)
print(summary(aovBS))
cat("\nsame residual deviance:",
    isTRUE(all.equal(deviance(mBS), sum(residuals(aovBS)^2))), "\n")
with_fig("13_aov_diagnostics", {
  par(mfrow = c(2, 2)); plot(aovBS)
}, width = 9, height = 8)

## ============================================================================
##  4. COUNTS: THE POISSON FAMILY
## ============================================================================
##  Species is a count, so the natural family is Poisson with a log link. The
##  coefficients are therefore on the LOG scale -- a unit change in Tmp
##  multiplies the expected count by exp(beta).
mST <- glm(Species ~ Tmp, data = m, family = poisson)
print(summary(mST))
cat("\nexp(coefficients) -- multiplicative effect on the count:\n")
print(round(exp(coef(mST)), 4))

with_fig("13_glm_poisson", {
  par(mfrow = c(1, 2))
  plot(m$Tmp, log(m$Species), pch = 19, col = "#2a7f7f",
       xlab = "Tmp", ylab = "log(Species)", main = "log counts vs temperature")
  abline(coef(mST), lwd = 2)
  plot(m$Species, fitted(mST), pch = 19, col = "#2a7f7f",
       xlab = "observed", ylab = "fitted", main = "Observed vs fitted")
  abline(0, 1, lty = 2)
}, width = 10, height = 5)

##  Poisson assumes mean = variance. Check it before trusting the standard
##  errors: a dispersion far above 1 means the SEs are too small.
disp <- sum(residuals(mST, type = "pearson")^2) / df.residual(mST)
cat("\ndispersion:", round(disp, 3),
    if (disp > 1.5) " -- overdispersed, consider quasipoisson\n" else " -- acceptable\n")

## ---- interactions -----------------------------------------------------------
mSTp  <- glm(Species ~ Tmp + pH, data = m, family = poisson)
mST_p <- glm(Species ~ Tmp * pH, data = m, family = poisson)
print(summary(mSTp))
print(summary(mST_p))
cat("\nAIC, additive vs interaction:\n"); print(AIC(mSTp, mST_p))
with_fig("13_glm_interaction", {
  par(mfrow = c(2, 2)); plot(mST_p)
}, width = 9, height = 8)

## ============================================================================
##  5. ONE MODEL PER GROUP, OR ONE MODEL WITH GROUPS?
## ============================================================================
##  Fitting the three soils separately gives three lines but no test of whether
##  they differ.
mBTl <- glm(Biomass ~ Tmp, data = m, family = gaussian, subset = Soil == "loam")
mBTc <- glm(Biomass ~ Tmp, data = m, family = gaussian, subset = Soil == "clay")
mBTs <- glm(Biomass ~ Tmp, data = m, family = gaussian, subset = Soil == "sand")

with_fig("13_glm_by_soil", {
  cols <- c(clay = "black", loam = "red", sand = "green3")
  plot(m$Tmp, m$Biomass, col = cols[as.character(m$Soil)], pch = 19,
       xlab = "Tmp", ylab = "Biomass", main = "One regression per soil type")
  abline(mBTc, col = cols["clay"], lwd = 2)
  abline(mBTl, col = cols["loam"], lwd = 2)
  abline(mBTs, col = cols["sand"], lwd = 2)
  legend("bottomleft", legend = names(cols), col = cols, lwd = 2, bty = "n")
})

##  Analysis of covariance does it properly, in one model. The interaction term
##  IS the hypothesis "the slopes differ"; dropping it assumes common slope and
##  different intercepts.
mBT_S  <- glm(Biomass ~ Tmp * Soil, data = m, family = gaussian)   # sep. slopes
mBT_Si <- glm(Biomass ~ Tmp + Soil, data = m, family = gaussian)   # common slope
print(summary(mBT_S))
print(summary(mBT_Si))
cat("\nAIC, separate vs common slopes:\n"); print(AIC(mBT_S, mBT_Si))
cat("\nlikelihood-ratio test of the interaction:\n")
print(anova(mBT_Si, mBT_S, test = "F"))

cat("\n[13_GLM] state the family and the link, check the dispersion, and let a",
    "single model with an interaction answer 'do the groups differ' rather",
    "than fitting each group on its own.\n")
