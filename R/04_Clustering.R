## ============================================================================
##  MS_LJMR :: 04_Clustering.R — k-means & hierarchical clustering
##
##  Original authors: Jorge Soberón & Laura Jiménez
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 05_kmeans_clustering_QRoo_Iris.R and 06_Hierarchical_Clustering.R
##
##  Improvements vs. originals
##  --------------------------
##  * No setwd()/x11(); data via get_data() (real Butterflies/Biodiv file or
##    seeded fallback). Windows-only rgl 3-D calls replaced by base graphics.
##  * The k-means "how many clusters?" search is wrapped in a reusable helper
##    and the elbow is drawn cleanly.
##  * Hierarchical section compares linkage methods side by side.
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
##  Files to keep next to this script: BiodivCountries.csv; ButterfliesQRoo2.csv
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
#   biodiv      = read.csv("BiodivCountries.csv", stringsAsFactors = TRUE),
#   butterflies = read.csv("ButterfliesQRoo2.csv", stringsAsFactors = TRUE),
#   iris        = { utils::data("iris");   iris },
#   mtcars      = { utils::data("mtcars"); mtcars },
#   stop("unknown data set: ", name))
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("MASS", "vegan")

## ############################################################################
##  PART A -- k-MEANS
## ############################################################################

## ============================================================================
##  A1. Butterflies of Quintana Roo
## ============================================================================
##  De la Maza, R. & J. Soberon (1998) "Morphological grouping of Mexican
##  butterflies in relation to habitat association." Biodiversity &
##  Conservation 7: 927-944.
##
##  Seven localities along a post-hurricane successional gradient:
##    HD = Highly Damaged (by a hurricane)    MA = Medium-age Acahual
##    SD = Slightly Damaged                   OA = Old Acahual
##    GA = Garden                             PF = Primary Forest
##    YA = Young Acahual
##
##  The question is NOT which species live where -- it is which WING PATTERNS
##  live where. Pattern is a functional trait, so the analysis asks whether
##  habitat structure selects for appearance.
but <- get_data("butterflies")
cat("Butterflies table:", nrow(but), "species x", ncol(but), "columns\n")
print(head(but))
sites <- c("HD", "SD", "GA", "YA", "MA", "OA", "PF")

## ---- exploratory: counts per pattern, one site at a time --------------------
## First attempt: the pattern names are long, so horizontal labels collide.
## las = 2 turns them perpendicular, and mar buys room at the bottom for them.
with_fig("04_but_boxplot_raw", {
  par(mar = c(7.5, 4, 2, 0.5))
  boxplot(but$HD ~ but$Pattern, las = 2, main = "Highly Damaged", xlab = "")
}, width = 8, height = 5.5)

## Better: order the patterns by their median count instead of alphabetically.
## reorder() is the single most useful function for making a categorical axis
## informative -- alphabetical order carries no information at all.
noHD <- with(but, reorder(Pattern, HD, median, na.rm = TRUE))
noPF <- with(but, reorder(Pattern, PF, median, na.rm = TRUE))

with_fig("04_but_boxplot_ordered", {
  par(mfrow = c(1, 2), mar = c(7.5, 4, 2, 0.5))
  boxplot(but$HD ~ noHD, las = 2, main = "Highly Damaged", xlab = "")
  boxplot(but$PF ~ noPF, las = 2, main = "Primary Forest",  xlab = "")
}, width = 12, height = 5.5)

## ---- from species x sites to sites x patterns -------------------------------
agg    <- aggregate(but[, sites], by = list(Pattern = but$Pattern), FUN = sum)
pnames <- as.character(agg[, 1])
## kmeans() clusters ROWS, and we want to cluster SITES, so transpose.
tagg <- as.matrix(t(agg[, -1]))
colnames(tagg) <- pnames
cat("\nSites x patterns matrix:", nrow(tagg), "x", ncol(tagg), "\n")
print(tagg[, 1:min(6, ncol(tagg))])

## ---- how many clusters? -----------------------------------------------------
set.seed(1998)
kmsites5 <- kmeans(tagg, 5)
kmsites3 <- kmeans(tagg, 3)
cat("\nk = 5 on 7 sites (almost as many clusters as objects):\n")
print(kmsites5$cluster)
cat("\nk = 3:\n"); print(kmsites3$cluster)
cat("\nWithin-cluster sums of squares:\n")
cat("  k = 5:", round(kmsites5$withinss, 1), "\n")
cat("  k = 3:", round(kmsites3$withinss, 1), "\n")

## str() on the fitted object is worth doing once: cluster, centers, withinss,
## tot.withinss, size and iter are all there.
str(kmsites3)

if (interactive() && has_pkg("rgl")) {
  rgl::plot3d(tagg[, 4:6], type = "s", size = 1.5, col = kmsites3$cluster)
  rgl::text3d(tagg[, 4:6], text = rownames(tagg), adj = 2)
} else {
  skip_note("rotatable 3-D view of the site clusters", "rgl")
}

## ============================================================================
##  A2. Choosing k on iris: the elbow
## ============================================================================
iris <- get_data("iris")
Xi   <- as.matrix(iris[, 1:4])

##  k-means minimises the total within-cluster sum of squares. At k = 1 that is
##  just the total sum of squares, which equals sum(variances) * (n - 1) -- the
##  (n - 1) because var() is the unbiased estimator.
mss    <- numeric(10)
mss[1] <- sum(apply(Xi, 2, var)) * (nrow(Xi) - 1)
set.seed(1998)
for (i in 2:10) mss[i] <- sum(kmeans(Xi, i, nstart = 10)$withinss)
cat("\nTotal within-cluster SS for k = 1..10:\n"); print(round(mss, 1))

with_fig("04_kmeans_elbow", {
  par(mai = c(1, 1, 0.5, 0.5))
  plot(mss, type = "b", pch = 19, col = "#1f3b73",
       xlab = "number of clusters k", ylab = "total within-cluster SS",
       main = "Accumulated sums of squares")
  grid(nx = 7, ny = 7)
})

## The curve drops steeply to 3 and then flattens: adding a fourth cluster buys
## little. Note it never increases -- more clusters ALWAYS fit better, which is
## exactly why you cannot choose k by minimising this number.

## ---- k = 3, and does it recover the species? --------------------------------
set.seed(1998)
fit3 <- kmeans(Xi, 3, nstart = 25)
sp   <- factor(iris$Species)

pal_sp <- c("#8c2d3a", "#1f3b73", "#2a7f7f")     # species
pal_cl <- c("orange", "violet", "grey40")        # clusters -- deliberately different

with_fig("04_species_vs_clusters", {
  par(mfrow = c(1, 2))
  plot(Xi[, 1:2], col = pal_sp[sp], pch = 19, main = "Species (known labels)")
  legend("topright", legend = levels(sp), fill = pal_sp, bty = "n", cex = .8)
  plot(Xi[, 1:2], col = pal_cl[factor(fit3$cluster)], pch = 19,
       main = "k-means clusters (k = 3)")
}, width = 10, height = 5)

with_fig("04_kmeans_pairs", {
  pairs(Xi, col = pal_cl[factor(fit3$cluster)], pch = 19,
        main = "k-means (k = 3) assignments")
})

with_fig("04_species_pairs", {
  pairs(Xi, col = pal_sp[sp], pch = 19, main = "Species")
})

cat("\nObjects per species:\n");  print(table(iris$Species))
cat("\nObjects per cluster:\n");  print(table(fit3$cluster))
ct <- table(cluster = fit3$cluster, species = iris$Species)
cat("\nCross-tabulation:\n"); print(ct)
cat("\nAgreement (best assignment) =",
    round(sum(apply(ct, 1, max)) / nrow(iris), 3), "\n")

## setosa is recovered perfectly; versicolor and virginica bleed into each
## other, exactly as the Q-mode heatmap in session 01 predicted.

## ---- would four clusters be better? -----------------------------------------
set.seed(1998)
fit4 <- kmeans(Xi, 4, nstart = 25)
cat("\nk = 4 cross-tabulation:\n")
print(table(cluster = fit4$cluster, species = iris$Species))

with_fig("04_kmeans_pairs_k4", {
  pairs(Xi, col = c("orange", "violet", "grey40", "skyblue")[factor(fit4$cluster)],
        pch = 19, main = "k-means (k = 4)")
})

## k = 4 splits one real species in two. A better SS is not a better biology.

## ############################################################################
##  PART B -- HIERARCHICAL CLUSTERING
## ############################################################################

## ============================================================================
##  B1. Three views of the same regions
## ============================================================================
##  BiodivCountries.csv describes most countries of the world. We cluster world
##  REGIONS three times over, using three different blocks of variables, and ask
##  whether the three pictures agree.
biodiv <- get_data("biodiv")

## Columns are selected BY NAME. The original script used positions
## (divJ <- c(3, 28:31, 34:36) and so on), which silently break if a column is
## ever inserted -- and the file has 62 of them.
divN <- c("AmphRich", "Rept_rich", "BirdRich", "MamsRich",
          "CountCentPlantDiv", "CountEcoregions", "CountHotspots")
govN <- c("VA_Mean", "VA_SD", "PS_Mean", "PS_SD", "GE_Mean", "GE_SD",
          "RL_Mean", "RL_SD", "CC_Mean", "CC_SD")
capN <- c("PapersCapita", "GBIF_OwnToTotal", "CountHerbaria", "log10NmbSpc",
          "InternetPenetration", "GDRE_Mean10Years", "logGDP2010", "logGDP2011")

## d = diversity, g = governance, c = capacity -- consistent naming pays off
## once there are three parallel analyses.
md <- na.omit(biodiv[, c("RegionCode", divN)])
mg <- na.omit(biodiv[, c("RegionCode", govN)])
mc <- na.omit(biodiv[, c("RegionCode", capN)])
cat("\nComplete cases -- diversity:", nrow(md), " governance:", nrow(mg),
    " capacity:", nrow(mc), "\n")

## aggregate countries into regions
region_means <- function(d) {
  a <- aggregate(d[, -1], by = list(Region = d$RegionCode), FUN = mean)
  M <- as.matrix(a[, -1]); rownames(M) <- a$Region; M
}
mda <- region_means(md); mga <- region_means(mg); mca <- region_means(mc)
options(digits = 3)
cat("\nRegional means, diversity block:\n"); print(head(mda))

## ---- distances, and are the three views the same? ---------------------------
dd <- vegan::vegdist(mda, method = "euclidean")
dg <- vegan::vegdist(mga, method = "euclidean")
dc <- vegan::vegdist(mca, method = "euclidean")
str(dd)
cat("\nSmallest diversity distance between two regions:", round(min(dd), 3), "\n")

## A Mantel test correlates two distance matrices, permuting rows and columns
## together to get a null distribution -- the standard way to ask whether two
## descriptions of the same objects agree.
cat("\nMantel, diversity vs governance:\n"); print(vegan::mantel(dd, dg))
cat("\nMantel, diversity vs capacity:\n");   print(vegan::mantel(dd, dc))
cat("\nMantel, governance vs capacity:\n");  print(vegan::mantel(dg, dc))

## ---- trees, unstandardized then standardized --------------------------------
td <- hclust(dd, method = "ward.D")
tc <- hclust(dc, method = "ward.D")
tg <- hclust(dg, method = "ward.D")

with_fig("04_hclust_three_views_raw", {
  par(mfrow = c(1, 3))
  plot(td, xlab = "Diversity",  main = "", sub = "")
  plot(tc, xlab = "Capacity",   main = "", sub = "")
  plot(tg, xlab = "Governance", main = "", sub = "")
}, width = 13, height = 5)

## Those dendrograms are dominated by whichever variable has the largest range
## -- bird richness is in the hundreds, a governance index is around zero. So:
mds <- scale(mda); mcs <- scale(mca); mgs <- scale(mga)
dds <- vegan::vegdist(mds, method = "euclidean")
dcs <- vegan::vegdist(mcs, method = "euclidean")
dgs <- vegan::vegdist(mgs, method = "euclidean")

with_fig("04_hclust_three_views_std", {
  par(mfrow = c(1, 3))
  plot(hclust(dds, method = "ward.D"), xlab = "Diversity",  main = "", sub = "")
  plot(hclust(dcs, method = "ward.D"), xlab = "Capacity",   main = "", sub = "")
  plot(hclust(dgs, method = "ward.D"), xlab = "Governance", main = "", sub = "")
}, width = 13, height = 5)

## ============================================================================
##  B2. Linkage is a modelling choice, not a detail
## ============================================================================
##  single   = nearest neighbour -> chaining, long straggly clusters
##  complete = furthest neighbour -> compact, equal-diameter clusters
##  average  = UPGMA, a compromise
##  ward.D   = minimise the increase in within-cluster SS (the k-means criterion
##             made hierarchical)
methods <- c(single = "single", complete = "complete",
             average = "average", ward = "ward.D")

with_fig("04_hclust_linkages", {
  par(mfrow = c(2, 2), mar = c(2, 4, 2, 1))
  for (nm in names(methods))
    plot(hclust(dds, method = methods[nm]), main = nm, xlab = "", sub = "")
})

## the same four on UNSTANDARDIZED distances, for comparison
with_fig("04_hclust_linkages_raw", {
  par(mfrow = c(2, 2), mar = c(2, 4, 2, 1))
  for (nm in names(methods))
    plot(hclust(dd, method = methods[nm]), main = paste(nm, "(raw)"),
         xlab = "", sub = "")
})

## ---- Ward tree cut into groups ----------------------------------------------
tree <- hclust(dds, method = "ward.D2")
grp  <- cutree(tree, k = 3)
cat("\nWard clusters of world regions (k = 3):\n"); print(grp)
with_fig("04_hclust_ward", {
  plot(tree, main = "Ward clustering of world regions", xlab = "", sub = "")
  rect.hclust(tree, k = 3, border = c("#1f3b73", "#8c2d3a", "#2a7f7f"))
})

## ============================================================================
##  B3. How much of the tree should you believe?
## ============================================================================
##  A dendrogram always returns a tree, even from noise. pvclust bootstraps the
##  VARIABLES and reports two support values at each node:
##    BP (green) = bootstrap probability, the fraction of replicates containing
##                 that cluster
##    AU  (red)  = approximately unbiased p-value, a multiscale-bootstrap
##                 correction of BP; AU >= 95 is the usual "believe it" line
##  Note pvclust clusters the COLUMNS of what you give it, hence the transpose.
if (has_pkg("pvclust")) {
  set.seed(1998)
  nb <- getOption("msljmr.nboot", 500)     # raise to 5000 for a final figure
  dcb <- pvclust::pvclust(t(mda), method.hclust = "ward.D",
                          use.cor = "all.obs", nboot = nb, quiet = TRUE)
  with_fig("04_pvclust_biodiv", {
    plot(dcb, cex.axis = 1.2, main = "Diversity: AU / BP support")
    pvclust::pvrect(dcb, alpha = 0.95)
  })
  print(dcb)
} else skip_note("bootstrap support for the dendrogram", "pvclust")

## ============================================================================
##  B4. The butterflies again, hierarchically
## ============================================================================
mbut <- scale(tagg)
dbut <- vegan::vegdist(mbut, method = "euclidean")
sbut <- hclust(dbut, method = "ward.D")
with_fig("04_hclust_butterflies", {
  plot(sbut, xlab = "Quintana Roo sites", sub = "",
       main = "Sites clustered by wing-pattern composition", cex.axis = 1.2)
})

## Does the tree recover the successional gradient HD -> PF?
cat("\nSite order along the Ward tree:", paste(sbut$labels[sbut$order],
                                               collapse = " -> "), "\n")

if (has_pkg("pvclust")) {
  set.seed(1998)
  pvbut <- pvclust::pvclust(t(tagg), method.hclust = "ward.D",
                            use.cor = "all.obs",
                            nboot = getOption("msljmr.nboot", 500), quiet = TRUE)
  with_fig("04_pvclust_butterflies", {
    plot(pvbut, xlab = "Quintana Roo sites", main = "", cex.axis = 1.2)
  })
}

## ============================================================================
##  B5. A third data set, and other ways to draw a tree
## ============================================================================
cars <- get_data("mtcars")

if (has_pkg("pvclust")) {
  set.seed(1998)
  carclus <- pvclust::pvclust(t(cars), method.hclust = "ward.D",
                              method.dist = "euclidean",
                              nboot = getOption("msljmr.nboot", 500), quiet = TRUE)
  ## the same data under a different linkage: mcquitty (WPGMA)
  carclus2 <- pvclust::pvclust(t(cars), method.hclust = "mcquitty",
                               method.dist = "euclidean",
                               nboot = getOption("msljmr.nboot", 500), quiet = TRUE)
  with_fig("04_pvclust_cars", {
    par(mfrow = c(1, 2))
    plot(carclus,  main = "mtcars, Ward",     cex = .6, cex.axis = 1)
    plot(carclus2, main = "mtcars, McQuitty", cex = .6, cex.axis = 1)
  }, width = 13, height = 6)
}

## ape draws the same hclust object in shapes base graphics will not: a fan or
## a cladogram. Nothing about the clustering changes -- only the projection of
## the tree onto the page, which is worth saying out loud after session 03.
if (has_pkg("ape")) {
  carhclus <- hclust(dist(cars), method = "mcquitty")
  apec     <- ape::as.phylo(carhclus)
  clus4    <- cutree(carhclus, 4)
  mypal    <- c("sandybrown", "steelblue3", "forestgreen", "tomato3")
  with_fig("04_ape_trees", {
    par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
    plot(apec, type = "fan",       tip.col = mypal[clus4], cex = .6,
         main = "fan")
    plot(apec, type = "cladogram", tip.col = mypal[clus4], cex = .6,
         main = "cladogram")
  }, width = 12, height = 6.5)
} else skip_note("fan and cladogram tree layouts", "ape")

cat("\n[04_Clustering] k-means needs k and gives flat groups; hierarchical",
    "gives a whole nested family and makes you choose where to cut. Both are",
    "sensitive to scaling, and neither knows what a species is.\n")
