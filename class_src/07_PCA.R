## ============================================================================
##  MS_LJMR :: 07_PCA.R  — Principal Component Analysis (I–III, PCoA folded in)
##
##  Original authors: Laura Jiménez & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from Soberon & Jimenez 09_PCAI / 10_PCAII / 11_PCAIII_Biplots and
##  12_PCOAI, per the Fall-2026 schedule (PCoA now folded into the PCA arc).
##
##  What changed vs. the originals
##  ------------------------------
##  * No setwd(); data via get_data(); no x11() (uses open_dev()).
##  * princomp / prcomp / SVD are shown to compute the SAME object; the
##    equivalence "PCA = spectral decomposition of the correlation matrix" is
##    demonstrated numerically.
##  * A clean biplot, a scree plot with the Kaiser rule, and a PCoA section that
##    shows metric MDS on a Euclidean distance reproduces PCA scores.
##
##  Run:  source("R/00_utils.R"); source("R/07_PCA.R")
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/07_PCA.R by scripts/make_class_src.py --
##  edit R/07_PCA.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("MASS")
#install.packages("vegan")
#install.packages("cluster")
library(MASS)
library(vegan)
library(cluster)

# Working directory -- point this at the folder holding the data files
#   BiodivCountries.csv, CitiesEurope.csv, NeotomaMorphoEnvir.csv, gen.10000.vcf  (optional), speciesCrawley3.csv
setwd("YOUR/DIRECTORY")

# Built-in data sets used below
data(iris)

## ############################################################################
##  PART A -- WHAT A PRINCIPAL COMPONENT IS            (schedule: PCA I, Oct 1)
## ############################################################################

##  The running data set is Neotoma: 615 woodrat specimens, 101 variables --
##  skull landmarks in four blocks plus 19 bioclimatic variables at the
##  collection locality.
m <- read.csv("NeotomaMorphoEnvir.csv", stringsAsFactors = TRUE)
cat("Neotoma table:", nrow(m), "specimens x", ncol(m), "variables\n")

## column blocks, selected by NAME rather than by position
dorsalN <- grep("^d[0-9]+$",   names(m), value = TRUE)   # 14
ventralN<- grep("^v[0-9]+$",   names(m), value = TRUE)   # 28
mandN   <- grep("^m[0-9]+$",   names(m), value = TRUE)   # 14
latN    <- grep("^L[0-9]+$",   names(m), value = TRUE)   # 20
bioN    <- grep("^BIO[0-9]+$", names(m), value = TRUE)   # 19
cat("blocks -- dorsal", length(dorsalN), "ventral", length(ventralN),
    "mandibular", length(mandN), "lateral", length(latN),
    "bioclim", length(bioN), "\n")

## ============================================================================
##  A1. Centering and scaling, done by hand once
## ============================================================================
pair <- c("BIO1", "BIO7")            # mean annual temperature, temp range
par(mfrow = c(1, 1))
plot(m[, pair], pch = 19, cex = .4, col = "grey40",
     main = "Two variables, as recorded")

##  Two things are wrong with that picture for our purposes: the units differ,
##  and the origin is arbitrary. Fix the origin first -- subtract the mean of
##  each column, which moves the cloud onto the centroid.
mu  <- colMeans(m[, pair])
cat("\ncolumn means:", round(mu, 3), "\n")
mus <- matrix(rep(mu, times = nrow(m)), ncol = 2, byrow = TRUE)
m2  <- as.matrix(m[, pair] - mus)

## Then fix the units -- divide by each column's standard deviation.
sds <- sapply(m[, pair], sd)
cat("column sds  :", round(sds, 3), "\n")
m3  <- cbind(m2[, 1] / sds[1], m2[, 2] / sds[2])

## Both steps at once: that is all scale() does.
m4  <- scale(m[, pair])
cat("hand-rolled z-transform equals scale():",
    isTRUE(all.equal(unname(m3), unname(m4), check.attributes = FALSE)), "\n")

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 13, height = 4.6)
op <- par(mfrow = c(1, 3))
plot(m[, pair], pch = 19, cex = .3, col = "grey40", main = "raw")
plot(m2, pch = 19, cex = .3, col = "grey40", main = "centered")
abline(v = 0, h = 0)
plot(m3, pch = 19, cex = .3, col = "grey40", main = "centered + scaled")
abline(v = 0, h = 0)
par(op)

## ---- covariance and correlation are matrix products --------------------------
##  Divide by n - 1, not n, because the mean was estimated from the same data.
n <- nrow(m)
cat("\n(t(centered) %*% centered) / (n-1)  ==  cov():\n")
print(round((t(m2) %*% m2) / (n - 1), 3))
print(round(cov(m[, pair]), 3))

cat("\nthe same product on STANDARDIZED data gives the correlation matrix:\n")
print(round((t(m4) %*% m4) / (n - 1), 3))
print(round(cor(m[, pair]), 3))

cat("\ncorrelations among four bioclim variables:\n")
print(round(cor(m[, c("BIO1", "BIO2", "BIO4", "BIO7")]), 3))

## ============================================================================
##  A2. princomp(): the parts of a PCA
## ============================================================================
bio4 <- as.matrix(m[, c("BIO1", "BIO8", "BIO15", "BIO17")])

if (requireNamespace("car", quietly = TRUE)) {
  par(mfrow = c(1, 1))
  car::scatterplotMatrix(bio4, regLine = list(col = "#8c2d3a"), smooth = FALSE,
                         col = "grey30", pch = 19,
                         main = "Four bioclimatic variables")

} else message("Skipped: scatterplot matrix of the bioclim block.  Install with: install.packages(c(\"car\"))")
## cor = FALSE means "work from the COVARIANCE matrix", i.e. do not standardize
pcbiov <- princomp(bio4, cor = FALSE, scores = TRUE)
str(pcbiov)

## ---- $sdev: square roots of the eigenvalues ---------------------------------
eigens <- pcbiov$sdev^2
cat("\neigenvalues:\n"); print(round(eigens, 3))
cat("proportion of variance:\n"); print(round(eigens / sum(eigens), 4))
print(summary(pcbiov))

## Note how brutally BIO17 (precipitation, in millimetres) dominates: on the
## covariance matrix, the variable with the largest units wins. That is the
## argument for cor = TRUE, and the reason most PCAs are run on correlations.

op <- par(mfrow = c(1, 2))
screeplot(pcbiov, type = "l", main = "scree (line)")
screeplot(pcbiov, type = "barplot", main = "scree (bars)")
par(op)

## ---- $loadings: the coefficients that define the axes ------------------------
##  Pielou considered these the most important part of a PCA: they are the
##  columns of the rotation matrix, and the eigenvectors of the spectral
##  decomposition of cov (or cor).
print(pcbiov$loadings)

## print.loadings() blanks small values for readability, which is unhelpful when
## you want the numbers. unclass() gets them all.
lds <- round(unclass(pcbiov$loadings), 3)
print(lds)
cat("\nSo  Comp.1 = ",
    paste(sprintf("%+.3f x %s", lds[, 1], rownames(lds)), collapse = " "), "\n")

## ---- $scores: the data in the rotated space ----------------------------------
cat("\nscores:", nrow(pcbiov$scores), "x", ncol(pcbiov$scores), "\n")
if (requireNamespace("car", quietly = TRUE)) {
  par(mfrow = c(1, 1))
  car::scatterplotMatrix(pcbiov$scores, regLine = FALSE, smooth = FALSE,
                         col = "grey30", pch = 19,
                         main = "The same specimens, in PC space")

}
## The scores are UNCORRELATED by construction -- that is the whole point:
cat("correlations among scores (should be ~0 off the diagonal):\n")
print(round(cor(pcbiov$scores), 10))

## ============================================================================
##  A3. princomp vs prcomp, and what center/scale actually do
## ============================================================================
##  princomp() eigen-decomposes the covariance/correlation matrix.
##  prcomp() takes the SVD of the data matrix itself -- numerically steadier,
##  and mathematically the same object. The vocabulary differs:
##      princomp$loadings  ==  prcomp$rotation
##      princomp$scores    ==  prcomp$x
pc_uu <- prcomp(bio4, retx = TRUE, center = FALSE, scale. = FALSE)
pc_cu <- prcomp(bio4, retx = TRUE, center = TRUE,  scale. = FALSE)
pc_us <- prcomp(bio4, retx = TRUE, center = FALSE, scale. = TRUE)
pc_cs <- prcomp(bio4, retx = TRUE, center = TRUE,  scale. = TRUE)
print(summary(pc_uu))
print(summary(pc_cs))

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 10, height = 9)
op <- par(mfrow = c(2, 2))
plot(pc_uu$x[, 1:2], pch = 19, cex = .3, main = "uncentered, unscaled")
plot(pc_cu$x[, 1:2], pch = 19, cex = .3, main = "centered, unscaled")
plot(pc_us$x[, 1:2], pch = 19, cex = .3, main = "uncentered, scaled")
plot(pc_cs$x[, 1:2], pch = 19, cex = .3, main = "centered, scaled")
par(op)

## Without centering, PC1 simply points at the centroid: the first component is
## spent describing where the cloud is rather than how it is shaped.
cat("\nprcomp rotation (centered + scaled):\n")
print(round(unclass(pc_cs$rotation), 3))

cat("\nprincomp(cor=TRUE) and prcomp(scale.=TRUE) agree up to sign:",
    isTRUE(all.equal(abs(unclass(princomp(bio4, cor = TRUE)$loadings)),
                     abs(unclass(pc_cs$rotation)), check.attributes = FALSE)), "\n")

## ============================================================================
##  A4. How many components? Kaiser and the broken stick
## ============================================================================
##  evplot() is Francois Gillet's function (2012, GPL-2), used throughout
##  Numerical Ecology in R. It draws the eigenvalues against two null models:
##  the average eigenvalue (Kaiser-Guttman) and the broken-stick distribution --
##  the expected share if total variance were split at random among the axes.
evplot <- function(ev, main = "") {
  n   <- length(ev)
  bsm <- data.frame(j = seq_len(n), p = 0)
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

pc_bio <- prcomp(as.matrix(m[, bioN]), center = TRUE, scale. = TRUE)
par(mfrow = c(1, 1))
evplot(pc_bio$sdev^2, "(19 bioclim variables)")

## ############################################################################
##  PART B -- PCA IN PRACTICE                         (schedule: PCA II, Oct 6)
## ############################################################################

## ============================================================================
##  B1. Where the specimens come from
## ============================================================================
spp  <- sort(unique(as.character(m$sp)))
cols <- rainbow(length(spp))
par(mfrow = c(1, 1))
if (requireNamespace("maps", quietly = TRUE)) {
  maps::map("world", xlim = c(-140, -70), ylim = c(15, 55), col = "grey70")
} else {
  plot(m$long, m$lat, type = "n", xlab = "longitude", ylab = "latitude")
}
points(m$long, m$lat, pch = 19, cex = .5,
       col = cols[factor(as.character(m$sp), levels = spp)])
legend("bottomleft", legend = spp, col = cols, pch = 19, cex = .55, bty = "n")
title("Neotoma woodrats")

## ============================================================================
##  B2. One PCA per block of variables
## ============================================================================
blocks <- list(dorsal     = as.matrix(m[, dorsalN]),
               ventral    = as.matrix(m[, ventralN]),
               lateral    = as.matrix(m[, latN]),
               mandibular = as.matrix(m[, mandN]),
               climatic   = as.matrix(m[, bioN]))
pcs <- lapply(blocks, prcomp, retx = TRUE, center = TRUE, scale. = TRUE)

cat("\nvariance captured by the first two components of each block:\n")
print(round(sapply(pcs, function(p) summary(p)$importance[3, 2]), 3))
print(summary(pcs$dorsal))
print(summary(pcs$climatic))

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 16, height = 7.5)
op <- par(mfrow = c(2, 5), mar = c(4, 4, 3, 1))
for (nm in names(pcs))
  screeplot(pcs[[nm]], type = "l", npcs = min(10, ncol(blocks[[nm]])),
            main = nm)
for (nm in names(pcs))
  plot(pcs[[nm]]$x[, 1:2], pch = 19, cex = .4, main = nm,
       col = cols[factor(as.character(m$sp), levels = spp)])
par(op)

## ============================================================================
##  B3. Comparing two ordinations: Procrustes
## ============================================================================
##  Procrustes was the innkeeper who made guests fit the bed by stretching or
##  amputating them. The analysis is gentler: it scales, rotates and translates
##  one configuration onto another, then measures what is left over.
pdv <- vegan::procrustes(pcs$dorsal, pcs$ventral)
par(mfrow = c(1, 1))
plot(pdv, main = "615 specimens: unreadable")

## 615 arrows is not a picture. Aggregate to the 15 species first -- about 40
## specimens each -- and the comparison becomes legible.
agg_block <- function(X) {
  a <- aggregate(X, by = list(sp = m$sp), FUN = mean)
  M <- as.matrix(a[, -1]); rownames(M) <- a$sp; M
}
aggs <- lapply(blocks, agg_block)
pcsA <- lapply(aggs, prcomp, retx = TRUE, center = TRUE, scale. = TRUE)

## the first argument is the TARGET, the second is rotated onto it
pro.ld <- vegan::procrustes(pcsA$lateral, pcsA$dorsal)
pro.lv <- vegan::procrustes(pcsA$lateral, pcsA$ventral, symmetric = TRUE)

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 12, height = 6)
op <- par(mfrow = c(1, 2))
plot(pro.ld, kind = 1, type = "text", main = "lateral vs dorsal")
plot(pro.lv, kind = 1, type = "text", main = "lateral vs ventral")
par(op)

## kind = 2 plots the residual per object, with the 25% (dashed), 50% (solid)
## and 75% (dashed) quantiles marked. Objects above the top line are the ones
## the two ordinations disagree about.
# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 12, height = 5)
op <- par(mfrow = c(1, 2))
plot(pro.lv, kind = 2, main = "Lateral vs ventral")
plot(pro.ld, kind = 2, main = "Lateral vs dorsal")
par(op)

## protest() permutes to test whether the agreement is better than chance.
## A small p means the two ordinations ARE associated.
set.seed(1998)
cat("\nPROTEST, lateral vs ventral:\n"); print(vegan::protest(pcsA$lateral, pcsA$ventral))
cat("\nPROTEST, lateral vs dorsal:\n");  print(vegan::protest(pcsA$lateral, pcsA$dorsal))

## ============================================================================
##  B4. PCA as fitting an ellipsoid
## ============================================================================
##  The covariance matrix is symmetric and positive definite, so its spectral
##  decomposition defines an ellipsoid: eigenvectors give the axis directions,
##  eigenvalues their squared half-lengths. The 95% ellipsoid is the level set
##  of the Mahalanobis distance at qchisq(0.95, df).
bio3   <- as.matrix(m[, c("BIO1", "BIO3", "BIO6")])   # a fairly Gaussian trio
mu3    <- colMeans(bio3)
sigma  <- cov(bio3)
psi    <- solve(sigma) / qchisq(0.95, df = 3)
es     <- eigen(psi)
stds   <- 1 / sqrt(es$values)         # semi-axis lengths of the 95% ellipsoid
cat("\n95% ellipsoid semi-axes:", round(stds, 2), "\n")
cat("fraction of points inside:",
    round(mean(mahalanobis(bio3, mu3, sigma) <= qchisq(.95, 3)), 3), "\n")

if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
  el <- rgl::ellipse3d(sigma, centre = mu3, level = .95)
  rgl::plot3d(bio3, col = "blue", alpha = .3, size = 5)
  rgl::plot3d(el, col = "green", alpha = .1, add = TRUE)
  for (k in 1:3) {
    lo <- mu3 - es$vectors[, k] * stds[k]
    hi <- mu3 + es$vectors[, k] * stds[k]
    rgl::segments3d(x = c(lo[1], hi[1]), y = c(lo[2], hi[2]),
                    z = c(lo[3], hi[3]), lwd = 3)
  }
  rgl::aspect3d("iso")
} else message("Skipped: 3-D confidence ellipsoid with its semi-axes.  Install with: install.packages(c(\"rgl\"))")
## The fit is only as good as the multinormality. Compare a well-behaved trio
## with an awkward one -- the scree plots differ, and so does the ellipsoid.
op <- par(mfrow = c(1, 2))
screeplot(princomp(m[, c("BIO1", "BIO3", "BIO6")],  cor = TRUE),
          type = "barplot", main = "BIO1/3/6 (ellipsoidal)")
screeplot(princomp(m[, c("BIO12", "BIO8", "BIO9")], cor = TRUE),
          type = "barplot", main = "BIO12/8/9 (less so)")
par(op)

## ============================================================================
##  B5. Numerical proof: PCA IS the spectral decomposition
## ============================================================================
bioS <- scale(as.matrix(m[, bioN]))
co   <- cor(bioS)
eig  <- eigen(co)
U    <- eig$vectors                 # the rotation matrix
Y    <- t(U) %*% t(bioS)            # scores, computed by hand
pc1  <- princomp(bioS, cor = TRUE, scores = TRUE)

par(mfrow = c(1, 1))
plot(t(Y)[, 1], pc1$scores[, 1], pch = 19, cex = .4, col = "#1f3b73",
     xlab = "t(U) %*% t(X), first row", ylab = "princomp scores, PC1",
     main = "Hand-rolled scores vs princomp")
abline(0, 1, col = "#8c2d3a", lwd = 2)
cat("\ncorrelation between the two routes:",
    round(abs(cor(t(Y)[, 1], pc1$scores[, 1])), 6),
    " (sign is arbitrary in any eigen-decomposition)\n")

## ############################################################################
##  PART C -- BIPLOTS                            (schedule: PCA III, Oct 8)
## ############################################################################

irisS <- scale(iris[, 1:4])
colnames(irisS) <- c("SepL", "SepW", "PetL", "PetW")
sp    <- factor(iris[, 5])

##  vegan::rda() with no constraints IS a PCA. Its vocabulary comes from
##  community ecology: "sites" are objects, "species" are variables.
pcaI <- vegan::rda(irisS)
cat("\nloadings live in $CA$v:\n"); print(round(pcaI$CA$v[, 1:2], 3))

parts <- vegan::scores(pcaI)
cat("\nsite scores (objects):\n");   print(head(parts$sites, 3))
cat("species scores (variables):\n"); print(parts$species)

## ---- a biplot, built by hand -------------------------------------------------
##  A biplot shows two things at once: where the objects fall, and how each
##  variable contributes to the axes. The variable arrows run from the origin to
##  the coordinates in scores()$species.
colvec <- c("green4", "orange", "mediumblue")
par(mfrow = c(1, 1))
plot(parts$sites, pch = 19, col = colvec[sp],
     xlim = c(-2.5, 2.75), ylim = c(-2.5, 2.5),
     main = "Scores, with loadings drawn as lines")
for (i in seq_len(nrow(parts$species)))
  lines(matrix(c(0, 0, parts$species[i, ]), nrow = 2, byrow = TRUE), lwd = 2)
text(parts$species, labels = rownames(parts$species), cex = .9, pos = 4)

## ---- and the same thing from the built-in function ---------------------------
par(mfrow = c(1, 1))
biplot(pcaI, display = "species", xlim = c(-2, 3), ylim = c(-2, 2),
              main = "Iris PCA biplot")
points(pcaI, col = colvec[sp], pch = 20)
## group ellipses and convex hulls, both from vegan
vegan::ordiellipse(pcaI, groups = sp, conf = .90, lwd = 2)
vegan::ordihull(pcaI, groups = sp, lwd = 2, lty = 3)
legend("topright", legend = levels(sp), fill = colvec, bty = "n", cex = .8)

## Is one group's ellipse really bigger than another's? Permutation test:
set.seed(1998)
print(vegan::ordiareatest(pcaI, groups = sp, area = "ellipse", permutations = 999))

## ============================================================================
##  C2. Biplot of the biodiversity table
## ============================================================================
bioC2 <- read.csv("BiodivCountries.csv", stringsAsFactors = TRUE)
richN <- c("AmphRich", "Rept_rich", "BirdRich", "MamsRich",
           "DensAmphRich", "DensRept_rich", "DensBirdRich", "DensMamsRich")
bioC  <- na.omit(bioC2[, c("RegionCode", richN)])
cat("\ncountries per region:\n"); print(table(bioC$RegionCode))

## CAR and NAM have very few countries; drop them so the biplot is not driven
## by two points.
keep      <- !bioC$RegionCode %in% c("CAR", "NAM")
bioCsmall <- droplevels(bioC[keep, ])
bioS      <- scale(bioCsmall[, -1])
colnames(bioS) <- c("AmpR","ReptR","BirdR","MamR","AmpD","RepD","BirdD","MamD")
rgs  <- factor(bioCsmall$RegionCode)

pcaB <- vegan::rda(bioS)
par(mfrow = c(1, 1))
biplot(pcaB, display = "species", xlim = c(-2, 2), ylim = c(-2, 2),
              main = "Countries: richness and richness density")
points(pcaB, col = as.integer(rgs), pch = 20)
legend("topright", legend = levels(rgs), col = seq_along(levels(rgs)),
       pch = 19, cex = .7, bty = "n")
par(mfrow = c(1, 1))
screeplot(pcaB, main = "Biodiversity PCA")

## ============================================================================
##  C3. Biplot of the Neotoma species means
## ============================================================================
morphSmall  <- aggregate(m[, bioN], by = list(sp = m$sp), FUN = mean)
morphSmallS <- scale(morphSmall[, -1])
rownames(morphSmallS) <- morphSmall$sp
pcaM <- vegan::rda(morphSmallS)
par(mfrow = c(1, 1))
biplot(pcaM, display = "species", main = "Neotoma species x bioclim")
points(pcaM, col = as.integer(factor(morphSmall$sp)), pch = 20, cex = 2)
legend("bottomleft", legend = morphSmall$sp,
       col = seq_len(nrow(morphSmall)), pch = 19, cex = .6, bty = "n")

## BiplotGUI offered an interactive biplot explorer, but it was archived from
## CRAN and needs Windows-only tooling. vegan's ordiplot(..., type = "n") plus
## identify() covers the same ground portably.
if (requireNamespace("BiplotGUI", quietly = TRUE)) {
  message("BiplotGUI is available: BiplotGUI::Biplots(Data = morphSmallS)")
} else message("Skipped: the interactive BiplotGUI explorer (archived on CRAN).  Install with: install.packages(c(\"BiplotGUI\"))")
## ############################################################################
##  PART D -- PRINCIPAL COORDINATE ANALYSIS      (schedule: folded into Oct 8)
## ############################################################################
##  PCoA (Gower 1966), a.k.a. classical or metric MDS, takes a DISTANCE matrix
##  and returns Euclidean coordinates reproducing those distances as closely as
##  possible. PCA assumes Euclidean distance on the variables; PCoA accepts any
##  distance, which is what lets it ordinate mixed-type or non-Euclidean data.
##  Contrast with NMDS (session 09), which preserves only the RANK ORDER of the
##  distances rather than their values.

## ============================================================================
##  D1. Cities of Europe: distances you already know the answer to
## ============================================================================
eur <- read.csv("CitiesEurope.csv", stringsAsFactors = TRUE)
cat("\n", nrow(eur), "European cities\n")

##  Great-circle distance from longitude/latitude (haversine). terra::distance()
##  with lonlat = TRUE does the same thing; this version needs no GIS stack.
haversine <- function(lon, lat, R = 6371) {
  p <- cbind(lon, lat) * pi / 180
  n <- nrow(p); D <- matrix(0, n, n)
  for (i in seq_len(n)) {
    dlon <- p[, 1] - p[i, 1]; dlat <- p[, 2] - p[i, 2]
    a <- sin(dlat / 2)^2 + cos(p[i, 2]) * cos(p[, 2]) * sin(dlon / 2)^2
    D[i, ] <- 2 * R * asin(pmin(1, sqrt(a)))
  }
  as.dist(D)
}
eudist <- haversine(eur$Long, eur$Lat)
attr(eudist, "Labels") <- as.character(eur$City)
cat("Madrid-Berlin great-circle distance:",
    round(as.matrix(eudist)[which(eur$City == "Madrid"),
                            which(eur$City == "Berlin")]), "km\n")

par(mfrow = c(1, 1))
if (requireNamespace("maps", quietly = TRUE)) maps::map("world", xlim = range(eur$Long) + c(-5, 5),
                               ylim = range(eur$Lat) + c(-5, 5), col = "grey70") else plot(eur$Long, eur$Lat, type = "n")
points(eur$Long, eur$Lat, pch = 16, col = "red")
text(eur$Long, eur$Lat, eur$City, cex = .7, pos = 3)
title("The cities, in geography")

## ---- cmdscale: PCoA in base R -----------------------------------------------
eupcoa <- cmdscale(eudist, eig = TRUE)
cat("\nfirst eigenvalues:\n"); print(round(eupcoa$eig[1:6], 1))
cat("first coordinates:\n");   print(round(head(eupcoa$points), 2))

par(mfrow = c(1, 1))
plot(eupcoa$eig / sum(eupcoa$eig), type = "b", pch = 19,
     xlab = "Eigenvalue", ylab = "Variance proportion", las = 1,
     main = "Two axes carry essentially everything")
abline(h = 0, lty = 2, col = "red")

## Broken-stick comparison, written out longhand (this is what evplot draws).
bstick <- function(ev) {
  n <- length(ev); b <- numeric(n); b[1] <- 0
  for (i in 2:n) b[i] <- b[i - 1] + 1 / (1 + n - i)
  100 * b / n
}
par(mfrow = c(1, 1))
brs <- t(cbind(100 * eupcoa$eig / sum(eupcoa$eig), rev(bstick(eupcoa$eig))))
barplot(brs, beside = TRUE, col = c("bisque", 2), las = 2, main = "Broken stick")
legend("topright", legend = c("eigenvalue", "broken stick"),
       fill = c("bisque", 2))

## ---- the reconstruction, and why the signs do not matter ---------------------
##  asp = 1 is essential: the axes are in the same units, so the aspect ratio
##  must be 1 or the distances are misrepresented.
# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 12, height = 6)
op <- par(mfrow = c(1, 2))
plot(eupcoa$points[, 1:2], xlab = "Axis 1", ylab = "Axis 2", pch = 16,
     col = "red", asp = 1, main = "PCoA of great-circle distances")
text(eupcoa$points[, 1:2], rownames(eupcoa$points), cex = .7, pos = 3)
plot(eupcoa$points[, 1:2] * -1, xlab = "Axis 1", ylab = "Axis 2", pch = 16,
     col = "red", asp = 1, main = "Both axes flipped: same map")
text(eupcoa$points[, 1:2] * -1, rownames(eupcoa$points), cex = .7, pos = 3)
par(op)

## How faithful is it? Correlate the reconstructed distances with the originals.
cat("\ncorrelation between original and recovered distances:",
    round(cor(as.vector(eudist), as.vector(dist(eupcoa$points[, 1:2]))), 4), "\n")

## ============================================================================
##  D2. Mixed-type data: Gower distance
## ============================================================================
##  Gower's coefficient handles numeric and categorical columns together. For
##  each variable j and pair (i, h) it defines a similarity S_jih:
##    numeric     : 1 - |x_i - x_h| / range(var_j)
##    categorical : 1 if the categories match, 0 otherwise
##    dichotomous : 1 if both present, and the comparison is EXCLUDED when both
##                  are absent (the double-zero problem)
##  with a validity weight d_jih that is 0 for excluded comparisons, so
##    G(i,h) = sum_j S_jih d_jih / sum_j d_jih
##  G is a similarity in [0, 1]; 1 - G is the dissimilarity.
biodata <- read.csv("speciesCrawley3.csv", stringsAsFactors = TRUE)
biodata$Species <- as.factor(biodata$Species)
biodata$Soil    <- as.factor(biodata$Soil)
biodata[, 2:5]  <- scale(biodata[, 2:5])
str(biodata)

biodist <- cluster::daisy(biodata[, -1], metric = "gower")
attr(biodist, "Labels") <- as.character(biodata$Species)

##  Gower distance is not Euclidean, so some eigenvalues come out NEGATIVE:
##  there is no Euclidean configuration that reproduces these distances exactly.
biopcoa <- cmdscale(biodist, k = 9, eig = TRUE)
cat("\nnegative eigenvalues:", sum(biopcoa$eig < 0), "of", length(biopcoa$eig), "\n")
cat("first eigenvalues:\n"); print(round(biopcoa$eig[1:8], 4))

pos <- sum(biopcoa$eig > 0)
op <- par(mfrow = c(1, 2))
plot(biopcoa$eig / sum(biopcoa$eig), type = "b", pch = 19, las = 1,
     xlab = "Eigenvalue", ylab = "Variance proportion", main = "all")
abline(h = 0, lty = 2, col = "red")
plot((biopcoa$eig / sum(biopcoa$eig))[1:pos], type = "b", pch = 19, las = 1,
     xlab = "Eigenvalue", ylab = "Variance proportion", main = "positive only")
abline(h = 0, lty = 2, col = "red")
par(op)

bioev4 <- biopcoa$points[, 1:4]
colnames(bioev4) <- paste("Axis", 1:4)
op <- par(mfrow = c(1, 2))
pairs(bioev4, pch = 16, cex = .5,
      col = rainbow(nlevels(biodata$Species))[biodata$Species])
par(op)

## Colour by SOIL instead of species and the structure appears: the first axes
## are picking up the categorical variable, not the taxonomy.
par(mfrow = c(1, 1))
pairs(bioev4, pch = 16, cex = .6,
      col = c("red", "blue", "green")[biodata$Soil])

## ---- ape::pcoa: eigenvalues with corrections and a broken stick --------------
if (requireNamespace("ape", quietly = TRUE)) {
  pcoaa <- ape::pcoa(biodist)
  par(mfrow = c(1, 1))
  plot(pcoaa$values$Relative_eig[1:10], type = "b", pch = 19,
       ylab = "Variance proportion", xlab = "Eigenvalue",
       main = "ape::pcoa, eigenvalues vs broken stick")
  grid(col = "lightgray")
  lines(pcoaa$values$Broken_stick[1:10], type = "b", col = "red")
  legend("topright", c("relative eigenvalue", "broken stick"),
         col = c("black", "red"), lty = 1, pch = 19, bty = "n")


  ## A PCoA biplot has to be built from scratch: the loadings are not in the
  ## distance matrix. ape::biplot.pcoa (Legendre) recovers approximate loadings
  ## by regressing the original variables on the standardized coordinates.
  ## Categorical variables have no loadings, so Soil is left out.
  par(mfrow = c(1, 1))
  ape::biplot.pcoa(pcoaa, biodata[, 2:5], dir.axis2 = -1)

} else message("Skipped: ape::pcoa and its biplot.  Install with: install.packages(c(\"ape\"))")
## ============================================================================
##  D3. Genomic data: PCA and PCoA side by side
## ============================================================================
##  Data and the PCA recipe courtesy of Ben J. Wiens. A VCF of 10000 SNPs is
##  read into a genlight object; PCA runs on the allele counts, PCoA on a
##  distance matrix built from them.
vcf_path <- "gen.10000.vcf"
if (file.exists(vcf_path) && requireNamespace("vcfR", quietly = TRUE) && requireNamespace("adegenet", quietly = TRUE) && requireNamespace("ade4", quietly = TRUE) && requireNamespace("cluster", quietly = TRUE)) {
  vcf <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  gl  <- vcfR::vcfR2genlight(vcf)
  ids <- colnames(vcfR::extract.gt(vcf))
  adegenet::pop(gl) <- factor(sort(rep(0:(length(ids) / 20 - 1), 20)))

  genopca   <- adegenet::glPca(gl, nf = 20)
  genodist  <- dist(gl)                                    # Euclidean
  genodistg <- cluster::daisy(as.matrix(gl), metric = "gower")
  genopcoa  <- ade4::dudi.pco(genodist,  scannf = FALSE, nf = 20)
  genopcoag <- ade4::dudi.pco(ade4::cailliez(genodistg), scannf = FALSE, nf = 20)

  colpal <- colorRampPalette(c("#d73027", "#fc8d59", "#fee090", "#e0f3f8",
                               "#91bfdb", "#4575b4"))
  cols_g <- colpal(nlevels(adegenet::pop(gl)))
  # Wide figure: widen the Plot pane, or open a sized device first --
  #   dev.new(width = 15, height = 5.5)
  op <- par(mfrow = c(1, 3), cex = 1.1)
  ade4::s.class(genopca$scores, adegenet::pop(gl), col = cols_g,
                axesell = FALSE, grid = FALSE, cstar = 0, pch = 16)
  title("PCA: PC1 and PC2")
  ade4::s.class(genopcoa$li,  adegenet::pop(gl), col = cols_g,
                axesell = FALSE, grid = FALSE, cstar = 0, pch = 16)
  title("PCoA (Euclidean)")
  ade4::s.class(genopcoag$li, adegenet::pop(gl), col = cols_g,
                axesell = FALSE, grid = FALSE, cstar = 0, pch = 16)
  title("PCoA (Gower)")
  par(op)
} else {
  message("Skipped: the genomic PCA/PCoA comparison (needs gen.10000.vcf).  Install with: install.packages(c(\"vcfR\", \"adegenet\"))")
}

## ============================================================================
##  D4. The identity: PCoA on Euclidean distance IS PCA
## ============================================================================
Xs <- scale(iris[, 1:4])
pco <- cmdscale(dist(Xs), k = 2, eig = TRUE)
pca <- prcomp(Xs)
cat("\ncorrelation of PCoA axis 1 with PCA PC1:",
    round(abs(cor(pco$points[, 1], pca$x[, 1])), 6), "\n")
cat("correlation of PCoA axis 2 with PCA PC2:",
    round(abs(cor(pco$points[, 2], pca$x[, 2])), 6), "\n")

op <- par(mfrow = c(1, 2))
plot(pca$x[, 1:2], pch = 19, cex = .6, col = colvec[sp], main = "PCA scores")
plot(pco$points,   pch = 19, cex = .6, col = colvec[sp],
     main = "PCoA on Euclidean distance", xlab = "Axis 1", ylab = "Axis 2")
par(op)

## the book's scree figure
par(mfrow = c(1, 1))
ev <- pca$sdev^2
barplot(ev, names.arg = paste0("PC", seq_along(ev)), col = "bisque",
        ylab = "eigenvalue", main = "Iris PCA - scree plot")
abline(h = 1, lty = 2, col = "#8c2d3a")
legend("topright", "Kaiser (eig = 1)", lty = 2, col = "#8c2d3a", bty = "n")

cat("\n[07_PCA] PCA rotates variables into uncorrelated components; the biplot",
    "shows objects and variables together; PCoA generalises the whole thing to",
    "any distance, and coincides with PCA when that distance is Euclidean.\n")
