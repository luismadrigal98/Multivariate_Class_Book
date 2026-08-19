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
##  As shipped, the source() line below borrows three helpers from the course
##  repository: get_data() (loads a data set), need() (loads packages) and
##  with_fig() (opens a plot device). To run this script entirely on its own,
##  delete that line and uncomment the block below. Nothing else changes.
##
##  The standalone with_fig() just draws each figure to the screen, one after
##  the other. To save them as files instead, replace its body with
##      png(paste0(name, ".png")); on.exit(dev.off()); force(expr)
##
##  Files to keep next to this script: none — doubs comes from the ade4 package (install.packages("ade4"))
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) { utils::data("doubs", package = "ade4")
#                  list(fish = doubs$fish, env = doubs$env) }
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("vegan")

doubs <- get_data("doubs")
keep  <- rowSums(doubs$fish) > 0
Y     <- doubs$fish[keep, colSums(doubs$fish[keep, ]) > 0]
X     <- scale(doubs$env[keep, ])

## ---- Canonical correlation: symmetric association between two tables -------
set.seed(1998)                     # CCorA and anova.cca below both permute
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

with_fig("14_rda", { ordiplot(rda, scaling = 2, type = "text",
                              main = "RDA: fish constrained by environment") })

cat("\n[14_CanonicalCorrelation] canonical correlation is symmetric (two-way",
    "association); RDA is directional (Y ~ X), the multivariate regression.\n")
