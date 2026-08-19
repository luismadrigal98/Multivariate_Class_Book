## ============================================================================
##  MS_LJMR :: 18_DiscriminantAnalysis.R — Linear & Quadratic Discriminant
##
##  Original authors: Marlon E. Cobos & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 18_DiscriminantAnalysis.R. Uses frozen taxon data.
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
##  Files to keep next to this script: taxon.csv
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) read.csv("taxon.csv", stringsAsFactors = TRUE)
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("MASS", "vegan")

taxon <- get_data("taxon")
taxon$Taxon <- factor(taxon$Taxon)

## ---- check the LDA assumption of equal covariances (Anderson's test) -------
d   <- vegan::vegdist(taxon[, -1], method = "euclidean")
bd  <- vegan::betadisper(d, taxon$Taxon)
set.seed(1998)                     # permutest() is a permutation test: seed it
cat("Homogeneity of dispersions (permutation test):\n")
print(vegan::permutest(bd)$tab[1, ])

## ---- train/test split ------------------------------------------------------
set.seed(123)
tr  <- sample(nrow(taxon), 2/3 * nrow(taxon))
train <- taxon[tr, ]; test <- taxon[-tr, ]

## ---- LDA -------------------------------------------------------------------
model <- MASS::lda(Taxon ~ ., data = train)
cat("\nProportion of trace (discriminant axes):\n")
print(round(model$svd^2 / sum(model$svd^2), 3))
pred  <- predict(model, test)$class
cat("\nConfusion matrix (held-out test):\n")
print(table(true = test$Taxon, predicted = pred))
cat("Hold-out accuracy:", round(mean(pred == test$Taxon), 3), "\n")

with_fig("18_lda", { plot(model, col = as.integer(train$Taxon)) })

## ---- QDA (unequal covariances) ---------------------------------------------
modelQ <- MASS::qda(Taxon ~ ., data = train)
predQ  <- predict(modelQ, test)$class
cat("\nQDA hold-out accuracy:", round(mean(predQ == test$Taxon), 3), "\n")

cat("\n[18_DiscriminantAnalysis] LDA maximises between/within variance under",
    "equal-covariance; test that assumption, else use QDA.\n")
