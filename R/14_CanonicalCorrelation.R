## ============================================================================
##  MS_LJMR :: 14_CanonicalCorrelation.R — Canonical Correlation & Redundancy
##
##  Original authors: Marlon E. Cobos, Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 20_Canonical_correlation.R / 20a_...RedundancyAnalysis.R.
##  Two tables from the frozen Doubs data: environment (X) vs fish (Y).
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
##  Files to keep next to this script: BiodiversityCountriesBiGv.csv
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
# get_data <- function(name) switch(name,
#   biodiv_gv = read.csv("BiodiversityCountriesBiGv.csv", stringsAsFactors = TRUE),
#   doubs     = { utils::data("doubs", package = "ade4")
#                     list(fish = doubs$fish, env = doubs$env, xy = doubs$xy) },
#   stop("unknown data set: ", name))
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("vegan")

## ############################################################################
##  PART A -- CANONICAL CORRELATION
## ############################################################################
##  Simple correlation relates two VARIABLES. Canonical correlation relates two
##  TABLES measured on the same objects: it finds linear combinations of the
##  columns of X and of Y whose correlation is as large as possible, then the
##  next such pair uncorrelated with the first, and so on. It is an ordination
##  that maximises CORRELATION rather than variance, and it is symmetric --
##  neither table is the response.
##
##  Redundancy analysis (Part B) is the asymmetric counterpart: it generalises
##  REGRESSION instead of correlation.

## ============================================================================
##  A1. Biodiversity and governance
## ============================================================================
##  BiodiversityCountriesBiGv.csv: complete cases only, four richness variables
##  and four World Bank governance indicators (www.govindicators.org):
##    VA voice and accountability | PS political stability and absence of
##    violence | GE government effectiveness | RL rule of law
b <- get_data("biodiv_gv")
rownames(b) <- make.unique(as.character(b[[2]]))
cat("table:", nrow(b), "countries x", ncol(b), "\n")
print(head(b))
cat("\nregions represented:\n"); print(table(b$RegionCode))

bb <- as.matrix(b[, 3:6])           # biodiversity block
bg <- as.matrix(b[, 7:10])          # governance block
colnames(bg) <- c("Voice&Account", "Stability", "GvtEfect", "RuleLaw")

##  vegan prefers standardized input; some CCA routines need the biodiversity
##  block left positive, which is why it is scaled separately here.
bbs <- scale(bb)
options(digits = 3)

## ---- look at the correlations first -----------------------------------------
##  With two tables the correlation matrix splits into four blocks:
##      {cor(X), cor(XY); cor(YX), cor(Y)}
##  CCA::img.matcor() draws exactly that. type = 1 shows the assembled matrix,
##  type = 2 the three blocks separately.
if (has_pkg("CCA")) {
  corrbg <- CCA::matcor(bbs, bg)
  with_fig("14_matcor", CCA::img.matcor(corrbg, type = 2), width = 10, height = 7)
} else skip_note("the img.matcor correlation display", "CCA")

cb  <- cor(bbs)
cg  <- cor(bg)
cgb <- cor(cbind(bbs, bg))
cat("\nwithin biodiversity:\n"); print(round(cb, 2))
cat("\nwithin governance:\n");   print(round(cg, 2))
cat("\ncross-block correlations:\n"); print(round(cgb[1:4, 5:8], 2))

##  gplots::heatmap.2 gives more control than heatmap(): no trace, no density
##  strip, readable label sizes.
if (has_pkg("gplots")) {
  with_fig("14_cor_heatmaps", {
    gplots::heatmap.2(cgb, trace = "none", col = hcl.colors(32, "RdBu", rev = TRUE),
                      cexRow = .8, cexCol = .8, density.info = "none",
                      main = "Biodiversity + governance")
  }, width = 8, height = 7)
} else skip_note("heatmap.2 correlation maps", "gplots")

## ---- the canonical correlation ----------------------------------------------
set.seed(1998)
ccabg <- vegan::CCorA(bbs, bg, stand.Y = TRUE, stand.X = TRUE, permutations = 999)
print(ccabg)

ne <- ccabg$Eigenvalues / sum(ccabg$Eigenvalues)
with_fig("14_ccora_scree", {
  plot(ne, type = "b", pch = 19, col = "#1f3b73",
       xlab = "canonical axis", ylab = "proportion",
       main = "Canonical correlation scree")
  abline(h = mean(ne), col = "red")
})

##  Two tests of H0 "the two tables are independent": a parametric one based on
##  Pillai's trace, and a permutation test.
cat("\nPillai trace p =", ccabg$p.Pillai, "   permutation p =", ccabg$p.perm, "\n")

##  In the biplot the arrow length is the correlation between a variable and the
##  axis, taken from $corr.Y.Cy and $corr.X.Cx. The two circles have radius 1
##  and 0.5, which is what makes those lengths readable.
with_fig("14_ccora_biplot", {
  biplot(ccabg)
}, width = 11, height = 6)
with_fig("14_ccora_biplot2", {
  biplot(ccabg, plot.type = "biplot")
}, width = 11, height = 6)

## ============================================================================
##  A2. The same question, region by region
## ============================================================================
##  Pooling all countries can hide opposite relationships in different regions.
for (rg in c("LAM", "SSA")) {
  ii <- which(b$RegionCode == rg)
  if (length(ii) < 8) { cat("\n", rg, ": too few countries\n"); next }
  set.seed(1998)
  cc <- vegan::CCorA(bbs[ii, ], bg[ii, ], stand.Y = TRUE, stand.X = TRUE,
                     permutations = 999)
  cat("\n", rg, " (n =", length(ii), ")  permutation p =", cc$p.perm,
      "  Pillai p =", cc$p.Pillai, "\n")
  nex <- cc$Eigenvalues / sum(cc$Eigenvalues)
  with_fig(paste0("14_ccora_", tolower(rg)), {
    par(mfrow = c(1, 2))
    plot(nex, type = "b", pch = 19, main = rg, xlab = "axis", ylab = "proportion")
    abline(h = mean(nex), col = "red")
    biplot(cc)
  }, width = 12, height = 6)
}

## ############################################################################
##  PART B -- REDUNDANCY ANALYSIS
## ############################################################################
##  RDA is multivariate regression followed by a PCA:
##    1. regress every variable in Y on all the variables in X;
##    2. collect the fitted values in Yhat;
##    3. run a PCA on Yhat;
##    4. use those rotations to compute fitted site scores.
##  So it is an ordination of Y CONSTRAINED by X.
doubs <- get_data("doubs")
keep  <- rowSums(doubs$fish) > 0
env   <- doubs$env[keep, ]
spe2  <- doubs$fish[keep, colSums(doubs$fish[keep, ]) > 0]

##  The Hellinger transformation -- square root of row proportions -- makes
##  species counts safe for the Euclidean geometry that RDA assumes.
Y <- vegan::decostand(spe2, "hellinger")

##  Build a predictor table that deliberately includes a FACTOR, so the machinery
##  has to cope with mixed types: slope cut at its quartiles.
q     <- quantile(env$slo)
slop2 <- cut(env$slo, breaks = c(-Inf, q[2], q[3], q[4], Inf),
             labels = c("low", "moderate", "steep", "v_steep"))
cat("\nslope classes:\n"); print(table(slop2))

env.phys  <- env[, c("alt", "slo", "flo")]
env.chem  <- env[, setdiff(names(env), c("alt", "slo", "flo"))]
env.physS <- scale(env.phys)
env.chemS <- scale(env.chem)
X <- data.frame(env.chemS, slope = slop2, env.physS[, c("alt", "flo")])

## ---- the correlation structure between the two tables ------------------------
if (has_pkg("CCA")) {
  with_fig("14_doubs_matcor", {
    CCA::img.matcor(CCA::matcor(as.matrix(Y), as.matrix(env.chemS)), type = 1)
  }, width = 9, height = 8)
}

## ---- symmetric first: chemistry vs physiography ------------------------------
set.seed(1998)
ccrPC <- vegan::CCorA(env.chemS, env.physS, permutations = 999)
print(ccrPC)
with_fig("14_ccora_doubs", {
  par(mfrow = c(1, 2))
  plot(ccrPC$Eigenvalues, type = "b", pch = 19, main = "scree",
       xlab = "axis", ylab = "eigenvalue")
  biplot(ccrPC, "b")
}, width = 12, height = 6)

## ---- then the constrained ordination -----------------------------------------
rda.spe.env <- vegan::rda(Y ~ ., data = X)      # ~. means every column of X
print(rda.spe.env)

##  Read the "Partitioning of variance" block. Inertia is variance; the
##  CONSTRAINED part is what X explains, the UNCONSTRAINED part is what is left.
##  The raw proportion is the multivariate R^2 and it is BIASED upward, so use
##  the adjusted version (Borcard et al. 2011).
cat("\nR-squared, raw and adjusted:\n"); print(vegan::RsquareAdj(rda.spe.env))

with_fig("14_rda_ordiplots", {
  par(mfrow = c(1, 2))
  ## scaling 1: distances among SITES are meaningful
  vegan::ordiplot(rda.spe.env, scaling = 1, type = "text",
                  main = "scaling 1: sites")
  ## scaling 2: angles between arrows and species are meaningful
  vegan::ordiplot(rda.spe.env, scaling = 2, type = "text",
                  main = "scaling 2: variables")
}, width = 13, height = 6.5)

##  The fitted object carries BOTH ordinations: $CCA is the constrained part
##  (the regression), $CA the unconstrained residual PCA.
cat("\nconstrained eigenvalues:\n"); print(round(rda.spe.env$CCA$eig, 4))

evplot <- function(ev, main = "") {
  n <- length(ev); bsm <- data.frame(j = seq_len(n), p = 0)
  bsm$p[1] <- 1 / n
  for (i in 2:n) bsm$p[i] <- bsm$p[i - 1] + (1 / (n + 1 - i))
  bsm$p <- 100 * bsm$p / n
  par(mfrow = c(2, 1))
  barplot(ev, main = paste("Eigenvalues", main), col = "bisque", las = 2)
  abline(h = mean(ev), col = "red")
  legend("topright", "Average eigenvalue", lwd = 1, col = 2, bty = "n")
  barplot(t(cbind(100 * ev / sum(ev), bsm$p[n:1])), beside = TRUE,
          main = "% variation", col = c("bisque", 2), las = 2)
  legend("topright", c("% eigenvalue", "Broken stick model"),
         pch = 15, col = c("bisque", 2), bty = "n")
}
with_fig("14_rda_evplot", evplot(rda.spe.env$CCA$eig, "(constrained axes)"),
         width = 8, height = 7)

## ============================================================================
##  B2. Which predictors are worth keeping?
## ============================================================================
##  ordiR2step() adds variables one at a time, maximising the ADJUSTED R^2, and
##  stops when it starts to fall or the full model's R^2 is reached. Starting
##  from Y ~ 1 means "intercept only", so the search is forward.
set.seed(1998)
rda.step <- vegan::ordiR2step(vegan::rda(Y ~ 1, data = X),
                              scope = formula(rda.spe.env),
                              direction = "forward", R2scope = TRUE,
                              pstep = 1000, trace = FALSE)
cat("\nselected model:\n"); print(rda.step$call)
cat("\nR-squared -- full vs selected:\n")
print(rbind(full = unlist(vegan::RsquareAdj(rda.spe.env)),
            step = unlist(vegan::RsquareAdj(rda.step))))

## two small models, for comparison
mod1 <- vegan::rda(Y ~ alt + oxy, data = X)
mod2 <- vegan::rda(Y ~ pH + flo,  data = X)
with_fig("14_rda_small_models", {
  par(mfrow = c(1, 2))
  vegan::ordiplot(mod1, scaling = 3, type = "text", main = "Y ~ alt + oxy")
  grid()
  vegan::ordiplot(mod2, scaling = 3, type = "text", main = "Y ~ pH + flo")
  grid()
}, width = 13, height = 6.5)

## ---- is any of this significant? ---------------------------------------------
##  anova.cca() is badly named -- it is a PERMUTATION test, not an ANOVA. H0:
##  no linear relationship between Y and X, i.e. the association between fish
##  and chemistry is indistinguishable from a random pattern.
set.seed(1998)
cat("\nfull model:\n");     print(vegan::anova.cca(rda.spe.env))
cat("\nselected model:\n"); print(vegan::anova.cca(rda.step))
cat("\nby axis (selected model):\n")
print(vegan::anova.cca(rda.step, by = "axis"))

## ---- collinearity -------------------------------------------------------------
##  vif.cca() reports how much each coefficient's variance is inflated by
##  correlation with the other predictors. Above ~10 is the usual alarm.
cat("\nvariance inflation, full model:\n");     print(round(vegan::vif.cca(rda.spe.env), 2))
cat("\nvariance inflation, selected model:\n"); print(round(vegan::vif.cca(rda.step), 2))

## the book's figure
with_fig("14_rda", {
  vegan::ordiplot(rda.step, scaling = 2, type = "text",
                  main = "Redundancy analysis: fish constrained by environment")
})

cat("\n[14_CanonicalCorrelation] canonical correlation is symmetric and",
    "generalises correlation; RDA is asymmetric and generalises regression.",
    "Report the ADJUSTED R-squared and a permutation test, never the raw one.\n")
