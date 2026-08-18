## ============================================================================
##  MS_LJMR :: 05_Clustering.R — k-means & hierarchical clustering
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
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("MASS")

## ---- k-MEANS ================================================================
iris <- get_data("iris")
Xi   <- as.matrix(iris[, 1:4])

## (a) choose k with the elbow of total within-cluster sum of squares --------
wss <- function(data, kmax = 10, nstart = 25) {
  tot <- numeric(kmax)
  tot[1] <- (nrow(data) - 1) * sum(apply(data, 2, var))
  for (k in 2:kmax) tot[k] <- sum(kmeans(data, k, nstart = nstart)$withinss)
  tot
}
w <- wss(Xi)
with_fig("05_kmeans_elbow", {
  plot(w, type = "b", pch = 19, col = "#1f3b73", lwd = 2,
       xlab = "Number of clusters k", ylab = "Total within-cluster SS",
       main = "k-means elbow (iris)")
  abline(v = 3, col = "#8c2d3a", lty = 2, lwd = 1.5)
})

## (b) fit k = 3 and compare clusters to the true species ---------------------
set.seed(1998)
km3 <- kmeans(Xi, 3, nstart = 25)
ct  <- table(cluster = km3$cluster, species = iris$Species)
cat("k-means (k=3) vs species:\n"); print(ct)
cat("\nAgreement (max assignment) =",
    round(sum(apply(ct, 1, max)) / nrow(iris), 3), "\n")

with_fig("05_kmeans_pairs", {
  pairs(Xi, col = c("#1f3b73", "#8c2d3a", "#2a7f7f")[km3$cluster], pch = 19,
        main = "k-means (k = 3) assignments")
})

## ---- HIERARCHICAL CLUSTERING ===============================================
## Aggregate the biodiversity table by region, standardize, cluster.
biodiv <- get_data("biodiv")
## Four richness counts and their four area-corrected (density) counterparts.
## These are the real BiodivCountries.csv column names; the seeded fallback in
## 00_utils.R emits the same ones, so this line works either way.
num    <- biodiv[, c("AmphRich","Rept_rich","BirdRich","MamsRich",
                     "DensAmphRich","DensRept_rich","DensBirdRich","DensMamsRich")]
## The real table is incomplete (3-15 missing values per column), and a single
## NA would propagate through mean() -> scale() -> dist() and abort hclust().
agg    <- aggregate(num, by = list(Region = biodiv$RegionCode),
                    FUN = function(x) mean(x, na.rm = TRUE))
M      <- scale(as.matrix(agg[, -1]))
rownames(M) <- agg$Region
Dreg   <- dist(M)                               # Euclidean on standardized means

## (a) compare four linkage rules --------------------------------------------
methods <- c(single = "single", complete = "complete",
             average = "average", ward = "ward.D2")
with_fig("05_hclust_linkages", {
  op <- par(mfrow = c(2, 2), mar = c(2, 4, 2, 1)); on.exit(par(op))
  for (nm in names(methods))
    plot(hclust(Dreg, method = methods[nm]), main = nm, xlab = "", sub = "")
})

## (b) Ward tree + cut into groups -------------------------------------------
tree <- hclust(Dreg, method = "ward.D2")
grp  <- cutree(tree, k = 3)
cat("\nWard clusters of world regions (k = 3):\n"); print(grp)
with_fig("05_hclust_ward", {
  plot(tree, main = "Ward clustering of world regions", xlab = "", sub = "")
  rect.hclust(tree, k = 3, border = c("#1f3b73", "#8c2d3a", "#2a7f7f"))
})

cat("\n[05_Clustering] done.\n")
