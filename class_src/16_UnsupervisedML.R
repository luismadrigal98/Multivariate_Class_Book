## ============================================================================
##  MS_LJMR :: 16_UnsupervisedML.R — Unsupervised Machine Learning
##
##  Author: Luis J. Madrigal-Roca  (new session for the MS_LJMR edition)
##
##  NEW session (revised Fall-2026 schedule).
##  Readings: ISLR Ch. 12; McInnes et al. 2018 (UMAP); van der Maaten & Hinton
##  2008 (t-SNE).
##
##  Idea: PCA is a *linear* map that preserves global variance. Modern manifold
##  learners (t-SNE, UMAP) are *non-linear* and preserve local neighbourhoods,
##  so they separate clusters that PCA leaves entangled. We show all three on a
##  gene-expression matrix (leukemia subtypes).
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/16_UnsupervisedML.R by scripts/make_class_src.py --
##  edit R/16_UnsupervisedML.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("Rtsne")
#install.packages("uwot")
library(Rtsne)
library(uwot)
library(RColorBrewer)

# Working directory
#   Point this at the folder that holds this script and its data/ subfolder.
#   Every file name below is resolved relative to it, so the script and its
#   data have to travel together -- or at least stay in step.
#
#   This session reads, from data/: demo/pam.csv, leukemiaExpressionSubset.rds
setwd("YOUR/DIRECTORY")

set.seed(1998)

## ############################################################################
##  PART A -- DIMENSION REDUCTION: PCA vs t-SNE vs UMAP
## ############################################################################

## ---- Data: genes x samples -> samples x genes ------------------------------
expr  <- readRDS("data/leukemiaExpressionSubset.rds")                    # matrix, rows = genes
X     <- t(as.matrix(expr))                      # rows = samples
## Sample names look like "ALL_GSM330151.CEL": the subtype is the field before
## the UNDERSCORE. Splitting on the period instead leaves the accession in place
## and gives every sample a class of its own -- see the purity check below.
type  <- factor(sub("_.*$", "", rownames(X)))    # ALL / AML / CLL / CML / NoL
pal   <- brewer.pal(n = 5, name = "Dark2")
cols  <- pal[type]
cat("samples per subtype:\n"); print(table(type))

## ---- 1. PCA (linear baseline) ----------------------------------------------
pc <- prcomp(X, center = TRUE, scale. = TRUE)
ve <- pc$sdev^2 / sum(pc$sdev^2)
cat("\nPCA variance explained (PC1-3):", round(ve[1:3], 3), "\n")

## ---- 2. t-SNE (local, probabilistic) ---------------------------------------
ts <- Rtsne::Rtsne(X, dims = 2, perplexity = 15, pca = TRUE, verbose = FALSE)

## ---- 3. UMAP (local, graph-based; usually faster & keeps more structure) ---
um <- uwot::umap(X, n_neighbors = 15, min_dist = 0.1, n_components = 2)

## ---- Compare the three embeddings ------------------------------------------
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
plot(pc$x[, 1:2], col = cols, pch = 19, main = "PCA",
     xlab = sprintf("PC1 (%.0f%%)", 100 * ve[1]),
     ylab = sprintf("PC2 (%.0f%%)", 100 * ve[2]))
plot(ts$Y, col = cols, pch = 19, main = "t-SNE", xlab = "dim 1", ylab = "dim 2")
plot(um,   col = cols, pch = 19, main = "UMAP",  xlab = "dim 1", ylab = "dim 2")
legend("topright", levels(type), col = pal, pch = 19, bty = "n", cex = .9)
par(op)

## ---- Quantitative check: neighbourhood purity ------------------------------
##  For each sample, what fraction of its k nearest neighbours in the map share
##  its subtype? A number, rather than an impression of the picture.
purity <- function(Y, lab, k = 10) {
  D <- as.matrix(dist(Y)); n <- nrow(Y); s <- 0
  for (i in 1:n) {
    nn <- order(D[i, ])[2:(k + 1)]
    s  <- s + mean(lab[nn] == lab[i])
  }
  s / n
}
cat("k-NN label purity  PCA:", round(purity(pc$x[, 1:2], type), 3),
    " t-SNE:", round(purity(ts$Y, type), 3),
    " UMAP:",  round(purity(um, type), 3), "\n")
cat("  (chance level with", nlevels(type), "equal subtypes =",
    round(1 / nlevels(type), 3), ")\n")

## ---- Pitfall: the layout depends on perplexity -----------------------------
##  Perplexity is roughly "how many neighbours each point should care about".
##  It is not a nuisance parameter to leave at its default: the same data can
##  look like one cloud or several islands depending on it.
##
##  It is also bounded by the sample size -- Rtsne requires
##  perplexity < (n - 1) / 3, so with 60 samples anything from 20 up is refused.
##  Small data sets simply cannot be asked to care about many neighbours.
cat("\nlargest perplexity allowed for", nrow(X), "samples:",
    floor((nrow(X) - 1) / 3 - 1e-9), "\n")
op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (pp in c(5, 12, 18)) {
  set.seed(1998)
  y <- Rtsne::Rtsne(X, dims = 2, perplexity = pp, pca = TRUE, verbose = FALSE)$Y
  plot(y, col = cols, pch = 19, main = paste("perplexity =", pp),
       xlab = "dim 1", ylab = "dim 2")
}
par(op)

## ############################################################################
##  PART B -- ASSOCIATION RULES: unsupervised learning without a projection
## ############################################################################
##  Everything above answers "which objects are alike?" by putting them in a
##  low-dimensional space. Association-rule mining asks a different unsupervised
##  question entirely -- "which things occur TOGETHER?" -- and answers it with
##  no distance, no projection and no embedding.
##
##  The data is a set of TRANSACTIONS, each a basket of ITEMS. A rule
##  {A, B} => {C} is scored by three numbers:
##
##    support(X)        = fraction of transactions containing X
##    confidence(A=>C)  = support(A and C) / support(A)
##                        i.e. P(C | A) -- how often the rule holds when it applies
##    lift(A=>C)        = confidence(A=>C) / support(C)
##                        i.e. how much MORE likely C is given A than at random.
##                        lift = 1 means independence; lift > 1 is association.
##
##  Our transactions are the presence-absence matrix of session 21: each SITE is
##  a basket and each SPECIES an item, so a rule reads "where these species are
##  present, that one usually is too".
if (requireNamespace("arules", quietly = TRUE)) {
  library(arules)
  pamdf <- read.csv("data/demo/pam.csv")
  pam   <- as.matrix(pamdf[, -1]) == 1            # logical: sites x species
  tr    <- as(pam, "transactions")
  cat("\n"); print(tr)
  cat("mean species per site:", round(mean(rowSums(pam)), 1), "\n")

  ## apriori searches for itemsets above a support floor, then turns them into
  ## rules above a confidence floor. Both floors are choices: set them too low
  ## and you drown, too high and you see only the obvious.
  rules <- apriori(tr,
                   parameter = list(supp = 0.25, conf = 0.9, minlen = 2, maxlen = 3),
                   control   = list(verbose = FALSE))
  cat("rules at support >= 0.25, confidence >= 0.9:", length(rules), "\n")

  ## Most of those are redundant: a longer rule whose shorter version is at
  ## least as good tells you nothing extra.
  rules <- rules[!is.redundant(rules)]
  cat("after dropping redundant rules:", length(rules), "\n")

  cat("\nstrongest associations, ranked by LIFT:\n")
  inspect(head(sort(rules, by = "lift"), 6))

  ## ---- why lift and not confidence ----------------------------------------
  ##  Confidence alone rewards common consequents: if a species sits in 90% of
  ##  sites, almost anything predicts it with 90% confidence while telling you
  ##  nothing. Lift divides that out. Compare the two rankings:
  cat("\nsame rules ranked by CONFIDENCE instead -- note the lift column:\n")
  inspect(head(sort(rules, by = "confidence"), 6))

  top <- sort(rules, by = "lift")[1]
  cat("\nmost-lifted rule, read out:  lift =", round(quality(top)$lift, 2),
      "means the consequent is that many times more likely at sites holding the",
      "antecedent than at a site picked at random.\n")

  ## ---- the standard diagnostic scatter -------------------------------------
  q <- quality(rules)
  par(mfrow = c(1, 1))
  plot(q$support, q$confidence, pch = 19, cex = .5,
       col = c("blue", "red"),
       xlab = "support", ylab = "confidence")
  legend("bottomright", legend = c("low lift", "high lift"), bty = "n",
         pch = 19, col = c("blue", "red"))

  ## ---- is a rule more than co-occurrence by chance? ------------------------
  ##  Nothing above is a significance test. With 120 species there are many
  ##  thousands of candidate rules, so some will look strong by chance. Fisher's
  ##  exact test on the 2x2 table, Bonferroni-corrected, is the honest check.
  cat("\nFisher test on the top rule's 2x2 table (Bonferroni over all rules):\n")
  lhs_i <- which(as(lhs(top), "matrix")[1, ])   # lhs() already is an itemMatrix
  rhs_i <- which(as(rhs(top), "matrix")[1, ])
  a <- rowSums(pam[, lhs_i, drop = FALSE]) == length(lhs_i)
  b <- pam[, rhs_i]
  ft <- fisher.test(table(antecedent = a, consequent = b))
  cat("  raw p =", signif(ft$p.value, 3),
      "  Bonferroni p =", signif(min(1, ft$p.value * length(rules)), 3), "\n")
} else {
  message("Skipped: the association-rule section.  Install with: install.packages(c(\"arules\"))")
}

cat("\n[16_UnsupervisedML] non-linear maps raise local purity but distort global",
    "geometry; association rules find co-occurrence directly, and lift -- not",
    "confidence -- is the number that says whether a pattern is surprising.\n")
