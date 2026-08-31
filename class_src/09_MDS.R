## ============================================================================
##  MS_LJMR :: 09_MDS.R — (non-metric) Multidimensional Scaling
##
##  Original authors: Marlon E. Cobos & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 15_MultidimensionalScaling.R. Uses mtcars (built-in).
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/09_MDS.R by scripts/make_class_src.py --
##  edit R/09_MDS.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("vegan")
library(vegan)

# Working directory -- point this at the folder holding the data files
#   BiodivCountries.csv, ButterfliesQRoo2.csv
setwd("YOUR/DIRECTORY")

# Built-in data sets used below
data(mtcars)

## ############################################################################
##  PART A -- NON-METRIC MDS                              (schedule: MDS I)
## ############################################################################
##  MDS is dimension reduction WITHOUT an eigen- or singular-value
##  decomposition. The recipe is:
##    1. compute a distance matrix between objects in the full space;
##    2. choose how many dimensions you want (usually 2 or 3);
##    3. search for a configuration of points in that many dimensions whose
##       distances match the originals as closely as possible.
##  Step 3 minimises a quantity called STRESS by iterative search, not by
##  solving an eigenproblem. There is no guaranteed solution and no unique one:
##  different starts can land on different configurations.
cars <- mtcars
cat("mtcars:", nrow(cars), "cars x", ncol(cars), "features\n")
print(head(cars))

## ============================================================================
##  A1. metaMDS: the workhorse
## ============================================================================
##  vegan::metaMDS() makes many decisions for you (random starts, rotation to
##  principal axes, rescaling). What you must choose is the DISTANCE:
##
##    "manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski",
##    "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
##    "binomial", "chao", "cao", "mahalanobis"
##
##  The default is Bray-Curtis with k = 2. Note vegan expects non-negative data.
set.seed(1998)
md2 <- vegan::metaMDS(cars, trace = FALSE)          # Bray-Curtis, the default
cat("\nstress (Bray-Curtis):", round(md2$stress, 4), "\n")

##  The fitted object is large. Four parts matter:
##    $points   the configuration -- the equivalent of PCA scores
##    $species  variable scores, used to draw a biplot
##    $stress   how good the fit is
##    $dist     the distances, for Mantel or Procrustes tests
cat("\nconfiguration (points):\n");  print(head(md2$points))
cat("\nvariable scores (species):\n"); print(head(md2$species))

par(mfrow = c(1, 1))
plot(md2$points, pch = 20, asp = 1, main = "nMDS configuration (Bray-Curtis)")
text(md2$points, labels = rownames(cars), cex = .5, pos = 4)

##  In an interactive session identify() lets you click a point to label it:
##      plot(md2$points, pch = 20); identify(md2$points)
##  orditorp() does the decluttering automatically -- it writes labels where
##  there is room and plots a point where there is not.
par(mfrow = c(1, 1))
plot(md2, main = "orditorp: labels where they fit")
vegan::orditorp(md2, display = "sites")
vegan::orditorp(md2, display = "species", col = "red")

if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
  md3D <- vegan::metaMDS(cars, k = 3, trace = FALSE)
  rgl::plot3d(md3D$points, type = "s", size = 1, col = 3)
  rgl::text3d(md3D$points, texts = rownames(cars), cex = .8, pos = 4)
} else message("Skipped: three-dimensional nMDS you can rotate.  Install with: install.packages(c(\"rgl\"))")
## ============================================================================
##  A2. Stress, and the Shepard diagram
## ============================================================================
##  Stress is the square root of the normalised sum of squared differences
##  between the original distances and the distances in the reduced space.
##  Rules of thumb:  < 0.05 excellent | < 0.1 good | < 0.2 usable | > 0.2 not
##  BUT stress is an optimisation score, not a probability. A low stress on a
##  badly chosen distance is still a bad ordination -- which is exactly what the
##  Euclidean run below shows.
set.seed(1998)
md3 <- vegan::metaMDS(cars, distance = "gower",     trace = FALSE)
md4 <- vegan::metaMDS(cars, distance = "euclidean", trace = FALSE,
                      autotransform = FALSE)
cat("\nstress -- bray:", round(md2$stress, 4),
    " gower:", round(md3$stress, 4),
    " euclidean:", round(md4$stress, 4), "\n")

##  The Shepard diagram plots original distance against ordination distance for
##  every pair. A tight, monotone cloud means the map is faithful. The two R^2
##  reported are the non-metric fit (1 - stress^2) and the linear fit, the
##  squared correlation between fitted values and ordination distances.
op <- par(mfrow = c(1, 2))
plot(md3$points, pch = 19, col = "#1f3b73", asp = 1,
     main = "nMDS configuration (Gower)", xlab = "dim 1", ylab = "dim 2")
text(md3$points, labels = rownames(cars), cex = .5, pos = 3)
vegan::stressplot(md3, main = "Shepard diagram")
par(op)

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 11, height = 5)
op <- par(mfrow = c(1, 2))
plot(md4, main = "Euclidean: low stress ...")
vegan::stressplot(md4, main = "... but a curved Shepard plot")
par(op)

##  Euclidean distance is wrong for this table -- the columns are on wildly
##  different scales (displacement in hundreds, number of gears in single
##  digits) -- and the Shepard plot shows the curvature that the stress value
##  alone hides.

## ============================================================================
##  A3. Goodness of fit, object by object
## ============================================================================
##  goodness() returns the cumulative proportion of inertia accounted for, per
##  site or per species. Drawing it as bubble size shows WHICH objects the map
##  places badly, rather than giving one number for the whole configuration.
gof3 <- vegan::goodness(md3)
gof2 <- vegan::goodness(md2)

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 12, height = 6)
op <- par(mfrow = c(1, 2))
plot(md3, display = "sites", type = "n", main = "Gower: fit per car")
points(md3, display = "sites", cex = 6 * gof3 / max(gof3), col = "seagreen")
vegan::orditorp(md3, display = "sites")
plot(md2, display = "sites", type = "n", main = "Bray-Curtis: fit per car")
points(md2, display = "sites", cex = 6 * gof2 / max(gof2), col = "seagreen")
vegan::orditorp(md2, display = "sites")
par(op)
cat("\nlarger bubbles are POORLY fitted objects\n")

## ============================================================================
##  A4. Does the distance choice change the answer? Procrustes
## ============================================================================
##  A Procrustes test compares two ordinations after optimally scaling,
##  rotating and translating one onto the other -- so it compares the TOPOLOGY
##  of the configurations rather than raw coordinates.
set.seed(1998)
p23 <- vegan::protest(md2, md3)
p24 <- vegan::protest(md2, md4)
cat("\nBray vs Gower:\n");     print(p23)
cat("\nBray vs Euclidean:\n"); print(p24)

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 12, height = 6)
op <- par(mfrow = c(1, 2))
plot(p23, main = "Bray vs. Gower")
plot(p24, main = "Bray vs. Euclidean")
par(op)

## The arrows show which cars the two ordinations disagree about; in an
## interactive session identify() names them.

## ############################################################################
##  PART B -- BIPLOTS FOR MDS                             (schedule: MDS II)
## ############################################################################
##  The original session used BiplotGUI (Lagrange et al. 2009), which was
##  archived from CRAN and only ever ran on Windows. Everything it demonstrated
##  -- variable axes through the configuration, per-object and per-variable
##  quality of fit, group structure -- can be built with vegan and base
##  graphics, which is what follows.
if (requireNamespace("BiplotGUI", quietly = TRUE)) {
  message("BiplotGUI is installed: try BiplotGUI::Biplots(Data = scale(cars))")
} else message("Skipped: the interactive BiplotGUI explorer (archived on CRAN).  Install with: install.packages(c(\"BiplotGUI\"))")
## ============================================================================
##  B1. Fitting variable axes onto an nMDS configuration
## ============================================================================
##  nMDS gives no loadings -- the configuration is not a linear map of the
##  variables. envfit() recovers the direction of steepest increase of each
##  variable across the map and tests it by permutation, which is the honest
##  way to get biplot arrows out of a non-metric ordination.
set.seed(1998)
fit_vars <- vegan::envfit(md3, cars, permutations = 999)
print(fit_vars)

par(mfrow = c(1, 1))
plot(md3, display = "sites", type = "n", main = "nMDS with fitted variables")
points(md3, display = "sites", pch = 19, col = "grey50")
vegan::orditorp(md3, display = "sites", cex = .7)
plot(fit_vars, col = "#8c2d3a")

##  ordisurf() goes further: instead of one straight arrow it fits a smooth
##  surface for a variable over the ordination, so a non-monotone response is
##  visible rather than averaged into a direction.
par(mfrow = c(1, 1))
plot(md3, display = "sites", type = "n", main = "Smooth surface of mpg")
points(md3, display = "sites", pch = 19, col = "grey60")
## assigned, not printed: ordisurf() returns the fitted GAM, and letting it
## auto-print would dump a model summary into the middle of the figures
srf <- vegan::ordisurf(md3, cars$mpg, add = TRUE, col = "#1f3b73")
vegan::orditorp(md3, display = "sites", cex = .7)

## ============================================================================
##  B2. Groups on the map
## ============================================================================
##  The BiplotGUI examples grouped cases (region, furniture type). The same
##  question here: do cars split by number of cylinders?
cyl <- factor(cars$cyl)
par(mfrow = c(1, 1))
plot(md3, display = "sites", type = "n", main = "Grouped by cylinder count")
points(md3, display = "sites", pch = 19, col = c("#8c2d3a","#1f3b73","#2a7f7f")[cyl])
vegan::ordiellipse(md3, groups = cyl, conf = .90, lwd = 2)
vegan::ordihull(md3, groups = cyl, lty = 3)
legend("topright", legend = levels(cyl), fill = c("#8c2d3a","#1f3b73","#2a7f7f"),
       bty = "n", title = "cylinders")

set.seed(1998)
cat("\nDoes cylinder count explain the configuration?\n")
print(vegan::adonis2(vegan::vegdist(cars, method = "gower") ~ cyl))

## ============================================================================
##  B3. The biodiversity and butterfly tables, ordinated the same way
## ============================================================================
biodiv2 <- read.csv("BiodivCountries.csv", stringsAsFactors = TRUE)
richN   <- c("AmphRich","Rept_rich","BirdRich","MamsRich",
             "DensAmphRich","DensRept_rich","DensBirdRich","DensMamsRich")
bd      <- na.omit(biodiv2[, c("RegionCode", richN)])
set.seed(1998)
md_bd <- vegan::metaMDS(bd[, -1], distance = "gower", trace = FALSE)
cat("\nbiodiversity nMDS stress:", round(md_bd$stress, 4), "\n")

par(mfrow = c(1, 1))
plot(md_bd, display = "sites", type = "n",
     main = "Countries by richness, grouped by region")
points(md_bd, display = "sites", pch = 19, cex = .7,
       col = as.integer(bd$RegionCode))
vegan::ordiellipse(md_bd, groups = bd$RegionCode, conf = .9, lwd = 1.5)
legend("topright", legend = levels(bd$RegionCode),
       col = seq_along(levels(bd$RegionCode)), pch = 19, cex = .7, bty = "n")

## Quintana Roo: sites x wing patterns, as in sessions 04 and 08
qroo4  <- read.csv("ButterfliesQRoo2.csv", stringsAsFactors = TRUE)
sitesQ <- c("HD", "SD", "GA", "YA", "MA", "OA", "PF")
qroo2  <- aggregate(qroo4[, sitesQ], by = list(Pattern = qroo4$Pattern), FUN = mean)
QRoo   <- as.matrix(t(qroo2[, -1]))
## short labels, as in the original
vars2  <- abbreviate(as.character(qroo2$Pattern), 8)
colnames(QRoo) <- vars2
print(round(head(QRoo), 2))

set.seed(1998)
md_qr <- vegan::metaMDS(QRoo, distance = "bray", trace = FALSE)
par(mfrow = c(1, 1))
plot(md_qr, type = "n", main = "Quintana Roo sites and wing patterns")
vegan::orditorp(md_qr, display = "sites",   col = "#1f3b73", cex = 1)
vegan::orditorp(md_qr, display = "species", col = "#8c2d3a", cex = .8)
cat("\nbutterfly nMDS stress:", round(md_qr$stress, 4),
    " (7 sites in 2 dimensions -- almost no constraint, so treat with care)\n")

cat("\n[09_MDS] stress is an optimisation score, not a probability; a tight,",
    "monotone Shepard diagram means the low-D map is faithful; and nMDS has no",
    "loadings, so biplot arrows must be fitted after the fact.\n")
