## ============================================================================
##  MS_LJMR :: 11_FactorAnalysis.R — Exploratory Factor Analysis
##
##  Original authors: Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from FactorAnalysis.R. Uses frozen taxon data (portable).
##
##  FA vs PCA: PCA summarises TOTAL variance with components that are exact
##  linear combinations of the variables. FA models the COMMON variance with a
##  small number of latent factors plus variable-specific 'uniquenesses':
##      x = Lambda f + e,   Cov(x) = Lambda Lambda' + Psi.
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/11_FactorAnalysis.R by scripts/make_class_src.py --
##  edit R/11_FactorAnalysis.R and regenerate; changes made here will be overwritten.
## ============================================================================

# Working directory
#   Point this at the folder that holds this session's data. Every file
#   name below is resolved relative to it, so the script and its data have
#   to travel together -- or at least stay in step.
#
#   This session reads: BiodiversityCountriesSSAFactanal2.csv, GraceSEM.csv  (optional), Insatisf2.csv
setwd("YOUR/DIRECTORY")

## ============================================================================
##  1. FA IS NOT PCA
## ============================================================================
##  PCA re-expresses ALL the variance with components that are exact linear
##  combinations of the variables. Factor analysis posits a small number of
##  unobserved COMMON factors that generate the correlations, plus a
##  variable-specific term:
##
##        x = Lambda f + e ,      Cov(x) = Lambda Lambda' + Psi
##
##  Lambda holds the loadings, f are the common factors (uncorrelated, unit
##  variance) and Psi is a diagonal matrix of UNIQUENESSES -- the part of each
##  variable no factor explains. FA does not try to fit that part; PCA has no
##  way of not fitting it. Estimation here is by maximum likelihood, which buys
##  a likelihood-ratio test of "are k factors enough?".

## ============================================================================
##  2. FIRST EXAMPLE: capacity, biodiversity and governance in Africa
## ============================================================================
m  <- read.csv("BiodiversityCountriesSSAFactanal2.csv", stringsAsFactors = TRUE)
m2 <- m[, -1]                       # drop the country label column
cat("table:", nrow(m2), "countries x", ncol(m2), "indicators\n")
print(names(m2))

## ---- how many factors? -------------------------------------------------------
ev <- eigen(cor(m2))
cat("\neigenvalues of the correlation matrix:\n"); print(round(ev$values, 3))
cat("above 1 (Kaiser):", sum(ev$values > 1), "\n")

##  Kaiser's rule is crude. Parallel analysis simulates correlation matrices
##  from RANDOM data of the same shape and keeps only factors whose eigenvalue
##  beats the simulated distribution -- a null model rather than a threshold.
if (requireNamespace("nFactors", quietly = TRUE)) {
  set.seed(1998)
  ap <- nFactors::parallel(subject = nrow(m2), var = ncol(m2),
                           rep = 100, cent = .05)
  ns <- nFactors::nScree(x = ev$values, aparallel = ap$eigen$qevpea)
  print(ns)
  par(mfrow = c(1, 1))
  nFactors::plotnScree(ns)
} else message("Skipped: parallel analysis / scree criteria.  Install with: install.packages(c(\"nFactors\"))")
## the plain scree, for the book
par(mfrow = c(1, 1))
plot(ev$values, type = "b", pch = 19, col = "#1f3b73",
     xlab = "factor", ylab = "eigenvalue", main = "Factor analysis scree")
abline(h = 1, col = "#8c2d3a", lty = 2)          # Kaiser criterion
legend("topright", "Kaiser (eig = 1)", lty = 2, col = "#8c2d3a", bty = "n")

## ---- fitting, and the effect of rotation -------------------------------------
##  The number of factors is a MODEL choice: factanal() reports a chi-square
##  test of the hypothesis that k factors are sufficient.
mf3 <- factanal(m2, factors = 3, rotation = "none", scores = "regression")
print(mf3, digits = 2, cutoff = 0.2)

mf4 <- factanal(m2, factors = 4, rotation = "none")
print(mf4, digits = 2, cutoff = 0.2)
cat("\n3 factors: chi-sq p =", signif(mf3$PVAL, 3),
    " | 4 factors: p =", signif(mf4$PVAL, 3),
    "  (large p = the model is adequate)\n")

##  Rotation does not change the fit at all -- it changes which basis of the
##  same factor space you look at. varimax keeps factors orthogonal and pushes
##  loadings toward 0 or 1; promax allows them to correlate.
mf3a <- factanal(m2, factors = 3, rotation = "varimax")
mf3b <- factanal(m2, factors = 3, rotation = "promax")
print(mf3a, digits = 2, cutoff = 0.2)
print(mf3b, digits = 2, cutoff = 0.2)
cat("\nsame uniquenesses under all three rotations:",
    isTRUE(all.equal(mf3$uniquenesses, mf3a$uniquenesses)), "\n")

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 13, height = 4.8)
op <- par(mfrow = c(1, 3))
for (fit in list(none = mf3, varimax = mf3a, promax = mf3b)) {
  L <- unclass(loadings(fit))
  plot(L[, 1:2], pch = 19, col = "#1f3b73", xlim = c(-1, 1), ylim = c(-1, 1),
       xlab = "Factor 1", ylab = "Factor 2")
  abline(h = 0, v = 0, col = "grey70")
  text(L[, 1:2], labels = rownames(L), cex = .55, pos = 3)
}
par(op)
cat("(panels: no rotation, varimax, promax -- same fit, different basis)\n")

## ============================================================================
##  3. SECOND EXAMPLE: social indicators of the US states
## ============================================================================
d  <- read.csv("Insatisf2.csv", stringsAsFactors = TRUE)
d2 <- na.omit(d[, -1])
rownames(d2) <- d[[1]][as.integer(rownames(d2))]
cat("\nstates table:", nrow(d2), "x", ncol(d2), "\n")

evd <- eigen(cor(d2))
cat("eigenvalues:\n"); print(round(evd$values, 3))
if (requireNamespace("nFactors", quietly = TRUE)) {
  set.seed(1998)
  apd <- nFactors::parallel(subject = nrow(d2), var = ncol(d2),
                            rep = 100, cent = .05)
  nsd <- nFactors::nScree(x = evd$values, aparallel = apd$eigen$qevpea)
  par(mfrow = c(1, 1))
  nFactors::plotnScree(nsd)
}

faI  <- factanal(x = d2, factors = 2, rotation = "none")
print(faI, digits = 2, cutoff = 0.2)

par(mfrow = c(1, 1))
L <- unclass(loadings(faI))
plot(L, pch = 19, col = "#1f3b73", xlim = c(-1, 1), ylim = c(-1, 1),
     main = "Two factors behind nine social indicators")
abline(h = 0, v = 0, col = "grey70")
text(L, labels = rownames(L), pos = 3, cex = .8)

faI2 <- factanal(x = d2, factors = 3, rotation = "none")
print(faI2, digits = 2, cutoff = 0.2, sort = TRUE)

## ============================================================================
##  4. FROM EXPLORATORY TO CONFIRMATORY
## ============================================================================
##  Everything above is EXPLORATORY: the data chose the structure. In a
##  confirmatory factor analysis you write the structure down first and the data
##  gets to reject it. Structural equation modelling generalises that to a whole
##  system of directed relationships.
##
##  The example is from Grace, Scheiner & Schoolmaster, "Structural equation
##  modeling: building and evaluating causal models", ch. 8 of Ecological
##  Statistics (OUP). Wetland variables:
##    landuse   ordinal human development in the watershed
##    buffer    whether development comes within 50 m of the wetland edge
##    hydro     degree of alteration of the natural hydrology
##    flooding  average water depth
##    soil      degree of soil disturbance
##    cond      water conductivity (a quality measure)
##    cattails  cover of Typha, an invader, as a percentage
##    richness  number of native plant species
grace_path <- "GraceSEM.csv"
if (file.exists(grace_path) && requireNamespace("lavaan", quietly = TRUE)) {
  ml <- read.csv(grace_path)

  par(mfrow = c(1, 1))
  plot(richness ~ landuse, data = ml, pch = 16, ylim = c(0, 10),
       xlab = "Land-use intensity", ylab = "Native species richness",
       cex = 1.1, cex.lab = 1.4)
  box(lwd = 1.8)
  text(0, 2, adj = 0, paste("r =", round(cor(ml$landuse, ml$richness), 2)),
       cex = 1.2)


  # Wide figure: widen the Plot pane, or open a sized device first --
  #   dev.new(width = 12, height = 6)
  op <- par(mfrow = c(2, 4), cex.lab = 1.2)
  for (v in c("landuse","buffer","hydro","soil",
              "flooding","richness","cond","cattails"))
    hist(ml[[v]], xlab = v, main = "", col = "light grey")
  par(op)

  ## the model of figure 6 in Grace et al.
  Mod1 <- 'buffer   ~ landuse
           hydro    ~ buffer + landuse
           flooding ~ hydro + landuse
           soil     ~ buffer
           cond     ~ soil + buffer
           cattails ~ cond
           richness ~ flooding + cattails'
  fit1 <- lavaan::sem(Mod1, data = ml)
  print(lavaan::summary(fit1))

  ## modification indices suggest paths the model is missing -- to be used with
  ## judgement, not automatically
  mi <- lavaan::modindices(fit1)
  print(head(mi[order(mi$mi, decreasing = TRUE), ], 8))

  Mod3 <- 'buffer   ~ landuse
           hydro    ~ buffer + landuse
           flooding ~ hydro + landuse
           soil     ~ buffer
           cond     ~ soil + buffer + landuse
           cattails ~ cond
           richness ~ flooding + cattails + landuse'
  print(lavaan::summary(lavaan::sem(Mod3, data = ml)))

  ## drop the non-significant paths for the final model
  Mod5 <- 'buffer   ~ landuse
           hydro    ~ landuse
           flooding ~ hydro
           soil     ~ buffer
           cond     ~ soil + landuse
           cattails ~ cond
           richness ~ flooding + landuse'
  print(lavaan::summary(lavaan::sem(Mod5, data = ml)))
} else {
  message("Skipped: ", paste("the confirmatory SEM example (needs data/GraceSEM.csv",
                  if (file.exists(grace_path)) "" else "-- file not present",
                  ")"), ".  Install with: install.packages(c(\"lavaan\"))")
}

cat("\n[11_FactorAnalysis] FA models the COMMON variance and leaves",
    "uniquenesses alone; rotation changes the basis, never the fit; and the",
    "number of factors is a hypothesis you can test.\n")
