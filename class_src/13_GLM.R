## ============================================================================
##  MS_LJMR :: 13_GLM.R — Generalized Linear Models (multivariate predictors)
##
##  Original authors: Marlon E. Cobos & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 19_GLModels.R. Uses frozen 'crawley' mixed-type table.
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/13_GLM.R by scripts/make_class_src.py --
##  edit R/13_GLM.R and regenerate; changes made here will be overwritten.
## ============================================================================

# Working directory
#   Point this at the folder that holds this session's data. Every file
#   name below is resolved relative to it, so the script and its data have
#   to travel together -- or at least stay in step.
#
#   This session reads: speciesCrawley3.csv
setwd("YOUR/DIRECTORY")

## ============================================================================
##  1. WHY GENERALIZED
## ============================================================================
##  A linear model says   Y ~ mu + X + error, and ordinary regression assumes
##  the error is (1) normal with mean 0 and (2) of constant variance across all
##  observations. A GLM relaxes both: the error may follow another distribution,
##  and its variance may depend on the mean. So you must state the FAMILY --
##  gaussian, Gamma, Poisson, binomial -- and the LINK that connects the linear
##  predictor to the mean.
m <- read.csv("speciesCrawley3.csv", stringsAsFactors = TRUE)
m$Soil <- factor(m$Soil)
str(m)

par(mfrow = c(1, 1))
pairs(m[, 2:5], col = factor(m$Species), pch = 19, cex = .6,
      main = "speciesCrawley3: biomass, temperature, precipitation, pH")

## ---- is the response even close to normal? ----------------------------------
op <- par(mfrow = c(1, 2))
qqnorm(m$Biomass, main = "Variable: Biomass"); qqline(m$Biomass)
hist(m$Biomass, col = "grey85", main = "Biomass", xlab = "")
par(op)
cat("\nShapiro-Wilk on Biomass: p =", signif(shapiro.test(m$Biomass)$p.value, 4),
    " (small p => not normal)\n")

## ---- and are the RESIDUALS of an OLS fit normal? ----------------------------
##  This is the assumption that actually matters; the marginal distribution of
##  the response does not have to be normal for OLS to be valid.
mLM <- lm(Biomass ~ Tmp, data = m)
par(mfrow = c(1, 1))
qqnorm(residuals(mLM), main = "Residuals, ordinary least squares")
qqline(residuals(mLM))
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

op <- par(mfrow = c(1, 2))
o <- order(m$Tmp)
plot(m$Tmp, m$Biomass, pch = 19, cex = .6, main = "Gaussian error",
     xlab = "Tmp", ylab = "Biomass")
lines(m$Tmp[o], fitted(mGauss)[o], col = "#2a7f7f", lwd = 2)
plot(m$Tmp, m$Biomass, pch = 19, cex = .6, main = "Gamma error",
     xlab = "Tmp", ylab = "Biomass")
lines(m$Tmp[o], fitted(mGamma)[o], col = "#8c2d3a", lwd = 2)
par(op)

##  predict(type = "response") returns the fitted MEAN; without it you get the
##  linear predictor, which for Gamma is on the inverse scale.
cat("\nfirst fitted values, response vs link scale:\n")
print(round(head(cbind(response = predict(mGamma, type = "response"),
                       link     = predict(mGamma)), 4), 4))

op <- par(mfrow = c(1, 2))
plot(mGamma, which = 1)      # residuals vs fitted: outliers and structure
plot(mGamma, which = 2)      # normal Q-Q of deviance residuals
par(op)

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
par(mfrow = c(1, 1))
boxplot(Biomass ~ Soil, data = m, col = "grey85",
        main = "Biomass by soil type")

##  Exactly the same model, in the vocabulary of experimental design:
aovBS <- aov(Biomass ~ Soil, data = m)
print(summary(aovBS))
cat("\nsame residual deviance:",
    isTRUE(all.equal(deviance(mBS), sum(residuals(aovBS)^2))), "\n")
op <- par(mfrow = c(2, 2)); plot(aovBS)
par(op)

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

op <- par(mfrow = c(1, 2))
plot(m$Tmp, log(m$Species), pch = 19, col = "#2a7f7f",
     xlab = "Tmp", ylab = "log(Species)", main = "log counts vs temperature")
abline(coef(mST), lwd = 2)
plot(m$Species, fitted(mST), pch = 19, col = "#2a7f7f",
     xlab = "observed", ylab = "fitted", main = "Observed vs fitted")
abline(0, 1, lty = 2)
par(op)

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
op <- par(mfrow = c(2, 2)); plot(mST_p)
par(op)

## ============================================================================
##  5. ONE MODEL PER GROUP, OR ONE MODEL WITH GROUPS?
## ============================================================================
##  Fitting the three soils separately gives three lines but no test of whether
##  they differ.
mBTl <- glm(Biomass ~ Tmp, data = m, family = gaussian, subset = Soil == "loam")
mBTc <- glm(Biomass ~ Tmp, data = m, family = gaussian, subset = Soil == "clay")
mBTs <- glm(Biomass ~ Tmp, data = m, family = gaussian, subset = Soil == "sand")

par(mfrow = c(1, 1))
cols <- c(clay = "black", loam = "red", sand = "green3")
plot(m$Tmp, m$Biomass, col = cols[as.character(m$Soil)], pch = 19,
     xlab = "Tmp", ylab = "Biomass", main = "One regression per soil type")
abline(mBTc, col = cols["clay"], lwd = 2)
abline(mBTl, col = cols["loam"], lwd = 2)
abline(mBTs, col = cols["sand"], lwd = 2)
legend("bottomleft", legend = names(cols), col = cols, lwd = 2, bty = "n")

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
