## ============================================================================
##  MS_LJMR :: 19_GLM.R — Generalized Linear Models (multivariate predictors)
##
##  Original authors: Marlon E. Cobos & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 19_GLModels.R. Uses frozen 'crawley' mixed-type table.
## ============================================================================
if (!exists("get_data")) source(file.path("R", "00_utils.R"))

m <- get_data("crawley")
m$Soil <- factor(m$Soil)
str(m)

## ---- Why GLM? the response is bounded/positive, errors non-Gaussian --------
cat("\nShapiro-Wilk on Biomass:", round(shapiro.test(m$Biomass)$p.value, 4),
    " (small p => non-normal)\n")

## ---- Gamma GLM for a positive continuous response --------------------------
g_gamma <- glm(Biomass ~ Tmp + pH, family = Gamma(link = "log"), data = m)
cat("\nGamma GLM (Biomass ~ Tmp + pH):\n"); print(round(coef(summary(g_gamma)), 4))

## ---- Poisson GLM for counts ------------------------------------------------
g_pois  <- glm(Species ~ Tmp + pH, family = poisson, data = m)
cat("\nPoisson GLM (Species ~ Tmp + pH) coefficients:\n")
print(round(coef(g_pois), 4))

## interaction vs additive: compare with AIC
g_int <- glm(Species ~ Tmp * pH, family = poisson, data = m)
cat("\nAIC additive vs interaction:\n"); print(AIC(g_pois, g_int))

with_fig("19_glm", width = 10, height = 4.4, {
  op <- par(mfrow = c(1, 2)); on.exit(par(op))
  plot(m$Tmp, m$Biomass, pch = 19, col = "#1f3b7355",
       xlab = "Temperature", ylab = "Biomass", main = "Gamma GLM fit")
  o <- order(m$Tmp)
  lines(m$Tmp[o], predict(g_gamma, type = "response")[o], col = "#8c2d3a", lwd = 2)
  plot(m$pH, m$Species, pch = 19, col = "#2a7f7f55",
       xlab = "pH", ylab = "Species (count)", main = "Poisson GLM fit")
})

cat("\n[19_GLM] choose the error family by the response: Gamma for positive",
    "continuous, Poisson for counts, binomial for 0/1. AIC compares models.\n")
