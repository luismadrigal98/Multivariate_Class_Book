## ============================================================================
##  MS_LJMR :: 20_CanonicalCorrelation.R — Canonical Correlation & Redundancy
##
##  Original authors: Marlon E. Cobos, Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 20_Canonical_correlation.R / 20a_...RedundancyAnalysis.R.
##  Two tables from the frozen Doubs data: environment (X) vs fish (Y).
## ============================================================================
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("vegan")

doubs <- get_data("doubs")
keep  <- rowSums(doubs$fish) > 0
Y     <- doubs$fish[keep, colSums(doubs$fish[keep, ]) > 0]
X     <- scale(doubs$env[keep, ])

## ---- Canonical correlation: symmetric association between two tables -------
cca_sym <- vegan::CCorA(scale(Y), X, permutations = 999)
cat("Canonical correlation (CCorA):\n")
cat("  Pillai trace p (permutation) =", cca_sym$p.perm, "\n")
cat("  canonical correlations       =",
    paste(round(head(cca_sym$CanCorr, 3), 3), collapse = ", "), "\n")

## ---- Redundancy analysis: DIRECTIONAL (Y explained by X) -------------------
Yh  <- vegan::decostand(Y, "hellinger")            # makes RDA appropriate
rda <- vegan::rda(Yh ~ ., data = as.data.frame(X))
cat("\nRDA variance partition:\n")
print(round(c(
  constrained   = rda$CCA$tot.chi / rda$tot.chi,
  unconstrained = rda$CA$tot.chi  / rda$tot.chi), 3))
cat("Adjusted R^2:", round(vegan::RsquareAdj(rda)$adj.r.squared, 3), "\n")
cat("Global test (anova.cca) p =",
    anova.cca(rda, permutations = 199)$`Pr(>F)`[1], "\n")

with_fig("20_rda", { ordiplot(rda, scaling = 2, type = "text",
                              main = "RDA: fish constrained by environment") })

cat("\n[20_CanonicalCorrelation] canonical correlation is symmetric (two-way",
    "association); RDA is directional (Y ~ X), the multivariate regression.\n")
