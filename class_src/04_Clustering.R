## ============================================================================
##  MS_LJMR :: 04_Clustering.R — k-means & hierarchical clustering
##
##  Original authors: Jorge Soberón & Laura Jiménez
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 05_kmeans_clustering_QRoo_Iris.R and 06_Hierarchical_Clustering.R
## ============================================================================

##  PLAIN CLASSROOM EDITION
## ============================================================================

# R packages required
# install.packages(c("vegan", "pvclust", "ggplot2"))
library(vegan)
library(ggplot2)

# Built-in dataset used for elbow & ground-truth validation
data(iris)

## ############################################################################
##  PART A -- k-MEANS CLUSTERING
## ############################################################################

## ============================================================================
##  A1. Butterflies of Quintana Roo
## ============================================================================
##  De la Maza, R. & J. Soberon (1998) "Morphological grouping of Mexican
##  butterflies in relation to habitat association." Biodiversity &
##  Conservation 7: 927-944.
##
##  Seven localities along a post-hurricane successional gradient:
##    HD = Highly Damaged (hurricane)        MA = Medium-age Acahual
##    SD = Slightly Damaged                  OA = Old Acahual
##    GA = Garden                            PF = Primary Forest
##    YA = Young Acahual
##
##  The question: Does habitat structure select for wing patterns?
##  We cluster the 7 localities based on functional wing-pattern frequencies.

# Working directory: set to folder containing ButterfliesQRoo2.csv
# setwd("YOUR/DIRECTORY")
dfile <- if (file.exists("ButterfliesQRoo2.csv")) "ButterfliesQRoo2.csv" else file.path("data", "ButterfliesQRoo2.csv")
but   <- read.csv(dfile, stringsAsFactors = TRUE)
sites <- c("HD", "SD", "GA", "YA", "MA", "OA", "PF")
cat("Butterflies table:", nrow(but), "species x", ncol(but), "columns\n")

## Aggregate species counts by Pattern -> transpose to get sites x patterns
agg    <- aggregate(but[, sites], by = list(Pattern = but$Pattern), FUN = sum)
pnames <- as.character(agg[, 1])
tagg   <- as.matrix(t(agg[, -1]))
colnames(tagg) <- pnames

cat("\nSites x patterns matrix:", nrow(tagg), "sites x", ncol(tagg), "patterns\n")
print(tagg[, 1:min(6, ncol(tagg))])

## Fit k-means with k = 3
## Note: kmeans() clusters ROWS (sites) using Euclidean distance in pattern space.
## Always use nstart > 1 to avoid bad local minima!
set.seed(1998)
kmsites3 <- kmeans(tagg, centers = 3, nstart = 25)

cat("\nk = 3 site clusters (Disturbed / Intermediate Acahual / Forest):\n")
print(kmsites3$cluster)
cat("Within-cluster sums of squares:", round(kmsites3$withinss, 1), "\n")

## ============================================================================
##  A2. Choosing k & Ground-Truth Validation (iris)
## ============================================================================
##  With 7 sites, k = 3 was guided by ecology. But how do we choose k objectively
##  in larger data? And does k-means actually recover biological reality?
Xi <- as.matrix(iris[, 1:4])

## The Elbow Method:
## k-means minimizes total within-cluster sum of squares (WSS).
## More clusters ALWAYS decrease WSS, so we look for where the curve "elbows".
mss    <- numeric(10)
mss[1] <- sum(apply(Xi, 2, var)) * (nrow(Xi) - 1)  # Total SS at k = 1
set.seed(1998)
for (i in 2:10) mss[i] <- sum(kmeans(Xi, i, nstart = 25)$withinss)

op <- par(mai = c(1, 1, 0.5, 0.5))
plot(1:10, mss, type = "b", pch = 19, col = "#1f3b73",
     xlab = "number of clusters k", ylab = "total within-cluster SS",
     main = "k-means Elbow Criterion")
grid(nx = 7, ny = 7)
abline(v = 3, lty = 2, col = "#8c2d3a")
par(op)

## Fit k = 3 and evaluate against true species labels
set.seed(1998)
fit3   <- kmeans(Xi, 3, nstart = 25)
sp     <- factor(iris$Species)
pal_sp <- c("#8c2d3a", "#1f3b73", "#2a7f7f")
pal_cl <- c("orange", "violet", "grey40")

op <- par(mfrow = c(1, 2))
plot(Xi[, 1:2], col = pal_sp[sp], pch = 19, main = "True Species (known)")
legend("topright", legend = levels(sp), fill = pal_sp, bty = "n", cex = 0.8)
plot(Xi[, 1:2], col = pal_cl[factor(fit3$cluster)], pch = 19,
     main = "k-means Clusters (k = 3)")
legend("topright", legend = levels(factor(fit3$cluster)), 
       fill = pal_cl , bty = "n", cex = 0.8)
par(op)

## Confusion matrix and agreement metric
ct <- table(cluster = fit3$cluster, species = iris$Species)
cat("\nCross-tabulation (k = 3 vs. Species):\n")
print(ct)
cat("Agreement (best assignment) =",
    round(sum(apply(ct, 1, max)) / nrow(iris), 3), "\n")
## setosa is recovered perfectly; versicolor and virginica overlap at the boundary.

## Pitfall: Would k = 4 be better?
fit4 <- kmeans(Xi, 4, nstart = 25)
cat("\nCross-tabulation (k = 4 vs. Species):\n")
print(table(cluster = fit4$cluster, species = iris$Species))
## k = 4 splits one real species in two. Lower WSS does not mean better biology!

## ############################################################################
##  PART B -- HIERARCHICAL CLUSTERING
## ############################################################################

## ============================================================================
##  B1. The Butterflies Hierarchically: Distances & Scaling
## ============================================================================
##  Agglomerative clustering builds a nested tree from a distance matrix.
##  Because abundant wing patterns would dominate Euclidean distance, we SCALE:
mbut_raw <- tagg
mbut_std <- scale(tagg)

dbut_raw <- vegan::vegdist(mbut_raw, method = "euclidean")
dbut_std <- vegan::vegdist(mbut_std, method = "euclidean")

## Compare trees: Raw vs Standardized (Ward's method)
op <- par(mfrow = c(1, 2))
plot(hclust(dbut_raw, method = "ward.D2"), xlab = "", sub = "",
     main = "Raw (dominated by common patterns)")
plot(hclust(dbut_std, method = "ward.D2"), xlab = "", sub = "",
     main = "Standardized (functional composition)")
par(op)

sbut <- hclust(dbut_std, method = "ward.D2")

## ============================================================================
##  B2. Linkage Rules & Tree Cutting
## ============================================================================
##  Linkage defines the distance between two clusters:
##    single   = nearest neighbour -> chaining
##    complete = furthest neighbour -> compact, equal diameter
##    average  = UPGMA -> intermediate compromise
##    ward.D2  = minimum increase in within-cluster SS (k-means analogue)
methods <- c(single = "single", complete = "complete",
             average = "average", "Ward's SS" = "ward.D2")

op <- par(mfrow = c(2, 2), mar = c(2, 4, 2, 1))
for (nm in names(methods))
  plot(hclust(dbut_std, method = methods[nm]), main = nm, xlab = "", sub = "")
par(op)

## Cutting the tree into discrete groups (k = 3)
grp <- cutree(sbut, k = 3)
cat("\nWard clusters of sites (k = 3):\n")
print(grp)

par(mfrow = c(1, 1))
plot(sbut, main = "Ward clustering of Quintana Roo sites", xlab = "", sub = "")
rect.hclust(sbut, k = 3, border = c("#1f3b73", "#8c2d3a", "#2a7f7f"))

## ============================================================================
##  B3. Seeing the Biology: Functional Wing Patterns Across Clusters
## ============================================================================
##  Why do the sites group the way they do? We plot the functional wing-pattern
##  profiles across the successional gradient.
##
##  NB the panels below are the A PRIORI successional stages, not cutree()'s
##  output -- the k = 3 cut puts YA alone and leaves MA with the disturbed sites:
site_order <- sites
cat("\nk = 3 cut:", paste(cutree(sbut, 3)[sites], collapse = " "),
    " vs a priori stage: 1 1 1 2 2 3 3\n")
cluster_levels <- c("Disturbed", "Acahual", "Mature Forest")
site_clusters <- c(
  HD = "Disturbed", SD = "Disturbed", GA = "Disturbed",
  YA = "Acahual", MA = "Acahual",
  OA = "Mature Forest", PF = "Mature Forest"
)

site_info <- data.frame(
  Site    = factor(site_order, levels = site_order),
  Cluster = factor(site_clusters[site_order], levels = cluster_levels)
)

df_long <- expand.grid(Site = factor(site_order, levels = site_order),
                       Pattern = colnames(tagg))
df_long$Z_score <- mapply(function(s, p) mbut_std[as.character(s), as.character(p)],
                          df_long$Site, df_long$Pattern)
df_long$Count   <- mapply(function(s, p) tagg[as.character(s), as.character(p)],
                          df_long$Site, df_long$Pattern)
df_long <- merge(df_long, site_info, by = "Site")

# Order patterns along the successional gradient (early -> mature)
site_ranks <- setNames(1:7, site_order)
pattern_score <- sapply(colnames(tagg), function(p) sum(tagg[, p] * site_ranks) / sum(tagg[, p]))
df_long$Pattern <- factor(df_long$Pattern, levels = names(sort(pattern_score)))

p_but <- ggplot(df_long, aes(x = Site, y = Pattern)) +
  geom_point(aes(size = Count, color = Z_score)) +
  scale_size_area(max_size = 11, breaks = c(10, 100, 300, 600), name = "Specimens\n(Count)") +
  scale_color_gradient2(low = "#2a7f7f", mid = "#e0e0e0", high = "#8c2d3a", midpoint = 0,
                        name = "Standardized\nAbundance (Z)") +
  facet_grid(~ Cluster, scales = "free_x", space = "free_x")   # a priori stage +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 11, color = "#1f3b73"),
    axis.text.x = element_text(face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 13, color = "#1f3b73"),
    plot.subtitle = element_text(size = 10, color = "#444444")
  ) +
  labs(
    title = "Functional Wing Patterns Across Forest Succession",
    subtitle = "De la Maza & Soberon (1998): Canopy closure selects for dark/reflective patterns; open ground selects for sand/shrub",
    x = "Successional Gradient: Disturbed (HD, SD, GA) -> Acahual (YA, MA) -> Mature Forest (OA, PF)",
    y = "Functional Wing Pattern"
  )
print(p_but)

## ============================================================================
##  B3b. Does k-means recover the gradient? Raw vs standardized
## ============================================================================
##  A partition tracks a gradient only if its groups are CONTIGUOUS blocks along
##  it. Check both scalings, for k = 2 and 3.
contig <- function(g) all(tapply(succ[names(g)], g,
                                 function(r) max(r) - min(r) + 1 == length(r)))
km <- do.call(rbind, lapply(2:3, function(k) {
  set.seed(1998); gr <- kmeans(mbut_raw, k, nstart = 50)$cluster[sites]
  set.seed(1998); gs <- kmeans(mbut_std, k, nstart = 50)$cluster[sites]
  rbind(data.frame(k = paste("k =", k), Scaling = "Raw counts",   Site = sites,
                   Cluster = factor(gr), ok = contig(gr)),
        data.frame(k = paste("k =", k), Scaling = "Standardized", Site = sites,
                   Cluster = factor(gs), ok = contig(gs)))
}))
km$Site <- factor(km$Site, levels = sites)
print(unique(km[, c("k", "Scaling", "ok")]), row.names = FALSE)

print(ggplot(km, aes(Site, Scaling, fill = Cluster)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = Cluster), colour = "white", fontface = "bold") +
  facet_wrap(~ k, ncol = 1) +
  scale_fill_manual(values = c("#8c2d3a", "#1f3b73", "#2a7f7f")) +
  theme_bw(base_size = 12) + theme(panel.grid = element_blank()) +
  labs(title = "k-means along the successional gradient",
       subtitle = "Contiguous blocks = the partition tracks succession",
       x = "HD (disturbed) -> PF (primary forest)", y = NULL))

##  The two scalings succeed at DIFFERENT k. Standardized gives the clean k = 2
##  split (mature vs disturbed), the one the bootstrap supports at 100%. Raw
##  counts look better at k = 3, but they sort sites by TOTAL CATCH -- 239
##  individuals at HD to 2521 at PF, Spearman 0.93 -- abundance, not wing
##  pattern. Scaling decides which gradient you find.
cat("\ntotal specimens per site:\n"); print(rowSums(mbut_raw)[sites])

## ============================================================================
##  B4. Tree Reliability: pvclust Bootstrap Support
## ============================================================================
##  A dendrogram ALWAYS draws branches, even from random noise.
##  pvclust bootstraps variables (columns) and reports two values:
##    BP (green) = bootstrap probability (raw resampling frequency)
##    AU (red)   = approximately unbiased p-value (multiscale correction; AU >= 95% = supported)
##  pvclust clusters COLUMNS, so we pass t(mbut_std).
if (requireNamespace("pvclust", quietly = TRUE)) {
  set.seed(1998)
  cat("\nRunning pvclust bootstrap on butterfly sites (nboot = 300)...\n")
  pvbut <- pvclust::pvclust(t(mbut_std), method.hclust = "ward.D2",
                            method.dist = "euclidean", nboot = 300, quiet = TRUE)
  par(mfrow = c(1, 1))
  plot(pvbut, main = "Quintana Roo sites: AU / BP support", cex.axis = 1.1)
  pvclust::pvrect(pvbut, alpha = 0.95)
} else {
  message("Skipped pvclust: install with install.packages(\"pvclust\")")
}

cat("\n[04_Clustering] Summary:",
    "\n  * k-means partitions into flat groups; needs k specified in advance.",
    "\n  * Hierarchical builds a nested family; cutting height determines groups.",
    "\n  * Both depend critically on scaling and distance metrics.",
    "\n  * Neither knows what a true biological entity is -- always validate!\n")
