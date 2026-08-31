## ============================================================================
##  MS_LJMR :: 12_DiscriminantAnalysis.R — Linear & Quadratic Discriminant
##
##  Original authors: Marlon E. Cobos & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 18_DiscriminantAnalysis.R. Uses frozen taxon data.
## ============================================================================

##  PLAIN CLASSROOM EDITION
##  Generated from R/12_DiscriminantAnalysis.R by scripts/make_class_src.py --
##  edit R/12_DiscriminantAnalysis.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("MASS")
#install.packages("vegan")
library(MASS)
library(vegan)

# Working directory -- point this at the folder holding the data files
#   taxon.csv
setwd("YOUR/DIRECTORY")

# Built-in data sets used below
data(iris)

## ============================================================================
##  1. THE PROBLEM
## ============================================================================
##  Discriminant analysis predicts DISCRETE classes from CONTINUOUS predictors.
##  The linear case assumes every group has the same covariance matrix, and then
##  finds the axes maximising the ratio of between-group to within-group
##  variance -- Fisher's criterion. Unlike PCA, it is told the labels, and that
##  is exactly why it separates better.
taxon <- read.csv("taxon.csv", stringsAsFactors = TRUE)           # Crawley, "The R Book"
taxon$Taxon <- factor(taxon$Taxon)
cat("taxon table:", nrow(taxon), "x", ncol(taxon), "\n"); print(head(taxon))
cat("\nspecimens per taxon:\n"); print(table(taxon$Taxon))

## ============================================================================
##  2. FIRST, TEST THE ASSUMPTION
## ============================================================================
##  LDA assumes homogeneous within-group dispersion. Anderson's (2006)
##  permutational test -- vegan::betadisper() plus permutest() -- is the
##  distance-based multivariate analogue of Levene's test. H0 is that the groups
##  have the SAME dispersion, so a large p is what lets you proceed.
dtax   <- vegan::vegdist(taxon[, -1])          # Bray-Curtis distances
groups <- taxon$Taxon
mod    <- vegan::betadisper(dtax, groups)
set.seed(1998)
print(anova(mod))
print(vegan::permutest(mod))

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 11, height = 5.5)
op <- par(mfrow = c(1, 2))
plot(mod, main = "Dispersion of each taxon")
boxplot(mod, main = "Distance to group centroid")
par(op)

##  The same check on a classic ecological table: vegan's varespec, 24 lichen
##  pasture plots split into grazed and ungrazed.
if (requireNamespace("vegan", quietly = TRUE)) {
  e <- new.env(); utils::data("varespec", package = "vegan", envir = e)
  dis <- vegan::vegdist(e$varespec)
  gr  <- factor(c(rep(1, 16), rep(2, 8)), labels = c("grazed", "ungrazed"))
  modv <- vegan::betadisper(dis, gr)
  set.seed(1998)
  print(vegan::permutest(modv))
  par(mfrow = c(1, 1))
  plot(modv, main = "varespec: grazed vs ungrazed dispersion")

}

## ============================================================================
##  3. LINEAR DISCRIMINANT ANALYSIS
## ============================================================================
options(digits = 3)
model <- MASS::lda(Taxon ~ ., data = taxon)     # CV = FALSE by default
print(model)

##  Read three things off that output:
##    Prior probabilities -- with no argument, the observed group frequencies
##    Coefficients        -- the weights defining each discriminant axis; large
##                           absolute values mark the discriminating variables
##    Proportion of trace -- how much of the between-group separation each axis
##                           carries (here LD1 dominates)
lev  <- as.integer(taxon$Taxon)
par(mfrow = c(1, 1))
plot(model, col = lev)

##  dimen = 1 collapses onto the first axis and draws histogram + density per
##  group -- the clearest picture of what "separation" means.
par(mfrow = c(1, 1))
plot(model, dimen = 1, type = "both")

##  $svd holds the singular values, so squaring gives a scree plot for LDA.
par(mfrow = c(1, 1))
plot(model$svd^2 / sum(model$svd^2), type = "b", pch = 19, col = "#1f3b73",
     xlab = "discriminant axis", ylab = "proportion of trace",
     main = "How much separation per axis")

## ---- posterior probabilities -------------------------------------------------
##  CV = TRUE runs leave-one-out cross-validation and returns, for every
##  specimen, the predicted class and the posterior probability of each class.
model2 <- MASS::lda(Taxon ~ ., data = taxon, CV = TRUE)
cat("\nposterior probabilities (first rows):\n")
print(round(head(model2$posterior), 3))
cat("\nleave-one-out confusion matrix:\n")
print(table(observed = taxon$Taxon, predicted = model2$class))
cat("LOO accuracy:", round(mean(model2$class == taxon$Taxon), 3), "\n")

preds <- predict(model)$class
cat("\nresubstitution -- predicted vs actual counts:\n")
print(rbind(predicted = table(preds), actual = table(taxon$Taxon)))

## ============================================================================
##  4. TRAIN AND TEST
## ============================================================================
##  Classifying the data you fitted on is not evidence. Hold part of it back.
set.seed(123)
train  <- sample(seq_len(nrow(taxon)), 80)
dTrain <- taxon[train, ]
dTest  <- taxon[-train, ]
cat("\ntraining set:\n"); print(table(dTrain$Taxon))
cat("test set:\n");       print(table(dTest$Taxon))

modelA  <- MASS::lda(Taxon ~ ., data = dTrain, CV = TRUE)   # LOO within train
modelA1 <- MASS::lda(Taxon ~ ., data = dTrain, CV = FALSE)  # a usable model

cat("\nconfusion matrix inside the training set (leave-one-out):\n")
print(table(observed = dTrain$Taxon, predicted = modelA$class))

predTest <- predict(modelA1, dTest)
cat("\nconfusion matrix on the HELD-OUT test set:\n")
print(table(observed = dTest$Taxon, predicted = predTest$class))
cat("hold-out accuracy:", round(mean(predTest$class == dTest$Taxon), 3), "\n")

predTrain <- predict(modelA1)$x
# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 11, height = 5.5)
op <- par(mfrow = c(1, 2))
plot(predTrain[, 1:2], col = as.integer(dTrain$Taxon), pch = 19,
     main = "training set in LD space", xlab = "LD1", ylab = "LD2")
plot(predTest$x[, 1:2], col = as.integer(dTest$Taxon), pch = 19,
     main = "test set, projected", xlab = "LD1", ylab = "LD2")
text(predTest$x[, 1:2], labels = predTest$class, cex = .6, pos = 3)
par(op)

if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
  rgl::plot3d(predTrain, col = as.integer(dTrain$Taxon), type = "s", size = 1.5)
  rgl::text3d(predTrain, text = modelA$class, cex = 1.25, adj = c(0, 0))
} else message("Skipped: rotatable 3-D view of the discriminant space.  Install with: install.packages(c(\"rgl\"))")
## ============================================================================
##  5. LDA VERSUS PCA ON THE SAME DATA
## ============================================================================
##  PCA is not told the groups, so it maximises total variance and the taxa
##  overlap. LDA is told, and separates them. Neither is "better" -- they answer
##  different questions.
pca <- prcomp(taxon[, -1], scale. = TRUE)
lda_all <- predict(MASS::lda(Taxon ~ ., data = taxon))$x

## (verify/make_figures.py writes its own lda_vs_pca.png for the book; this one
## is the R-side twin, kept under the session prefix so the two never collide)
# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 11, height = 5.5)
op <- par(mfrow = c(1, 2))
plot(pca$x[, 1:2], col = lev, pch = 19, main = "PCA: unsupervised",
     xlab = "PC1", ylab = "PC2")
plot(lda_all[, 1:2], col = lev, pch = 19, main = "LDA: supervised",
     xlab = "LD1", ylab = "LD2")
legend("topright", legend = levels(taxon$Taxon), col = 1:4, pch = 19,
       bty = "n", cex = .8)
par(op)

## ============================================================================
##  6. THE SAME MACHINERY ON IRIS
## ============================================================================

diris  <- vegan::vegdist(iris[, -5])
modi   <- vegan::betadisper(diris, iris$Species)
set.seed(1998)
print(vegan::permutest(modi))

set.seed(1998)
itrain <- sample(seq_len(nrow(iris)), 80)
modI   <- MASS::lda(Species ~ ., data = iris[itrain, ])
iTest  <- iris[-itrain, ]
pI     <- predict(modI, iTest)
cat("\niris hold-out confusion matrix:\n")
print(table(observed = iTest$Species, predicted = pI$class))
cat("accuracy:", round(mean(pI$class == iTest$Species), 3), "\n")

# Wide figure: widen the Plot pane, or open a sized device first --
#   dev.new(width = 11, height = 5.5)
op <- par(mfrow = c(1, 2))
plot(pI$x, col = as.integer(factor(iTest$Species)), pch = 19,
     main = "iris test set in LD space")
plot(modI, dimen = 1, type = "both")
par(op)

## the book's figure
par(mfrow = c(1, 1))
plot(lda_all[, 1:2], col = lev, pch = 19,
     xlab = "LD1", ylab = "LD2", main = "Linear discriminant space")

## ============================================================================
##  7. WHEN THE COVARIANCES ARE NOT EQUAL: QDA
## ============================================================================
##  Quadratic discriminant analysis drops the equal-covariance assumption and
##  fits one covariance matrix per group. The decision boundaries become
##  quadratic surfaces, at the cost of many more parameters -- it needs more
##  data per group than LDA does.
options(digits = 4)
modelQ <- MASS::qda(Taxon ~ ., data = taxon, CV = TRUE)
cat("\nQDA leave-one-out confusion matrix:\n")
print(table(observed = taxon$Taxon, predicted = modelQ$class))
cat("QDA LOO accuracy:", round(mean(modelQ$class == taxon$Taxon), 3),
    " vs LDA:", round(mean(model2$class == taxon$Taxon), 3), "\n")

cat("\n[12_DiscriminantAnalysis] test the equal-dispersion assumption before",
    "trusting LDA, always report a held-out confusion matrix, and remember the",
    "separation LDA shows you was handed to it in the labels.\n")
