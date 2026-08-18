## ============================================================================
##  MS_LJMR :: 14_CorrespondenceAnalysis.R — Correspondence Analysis (CA)
##
##  Original authors: Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 14_CorrespondaceAnalysis.R. Uses frozen Doubs fish counts.
## ============================================================================
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("vegan")

doubs <- get_data("doubs")
fish  <- doubs$fish
## drop any empty rows/cols (CA cannot handle them)
fish  <- fish[rowSums(fish) > 0, colSums(fish) > 0]

## ---- CA in vegan -----------------------------------------------------------
ca <- vegan::cca(fish)
ev <- ca$CA$eig
cat("Inertia (eigenvalues) and proportion per axis:\n")
print(round(rbind(inertia = ev[1:4], prop = (ev / sum(ev))[1:4]), 4))

## ---- biplot (scaling 1: distances among sites preserved) -------------------
with_fig("14_ca_biplot", {
  plot(ca, scaling = 1, main = "Correspondence analysis (scaling 1: sites)")
})

## ---- detrending removes the arch (horseshoe) effect ------------------------
dca <- vegan::decorana(fish)
with_fig("14_ca_vs_dca", {
  op <- par(mfrow = c(1, 2)); on.exit(par(op))
  plot(ca,  main = "CA (arch effect)")
  plot(dca, main = "Detrended CA")
})

cat("\n[14_CorrespondenceAnalysis] CA = SVD of a chi-square-standardized table;",
    "ordinates sites and species together. Watch for the arch on long gradients.\n")
