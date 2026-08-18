## ============================================================================
##  MS_LJMR :: 18_DiscriminantAnalysis.R — Linear & Quadratic Discriminant
##
##  Original authors: Marlon E. Cobos & Jorge Soberón
##                    (BIOL 943 Multivariate Methods, University of Kansas)
##  Revised for the MS_LJMR edition by: Luis J. Madrigal-Roca
##
##  Revised from 18_DiscriminantAnalysis.R. Uses frozen taxon data.
## ============================================================================
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("MASS", "vegan")

taxon <- get_data("taxon")
taxon$Taxon <- factor(taxon$Taxon)

## ---- check the LDA assumption of equal covariances (Anderson's test) -------
d   <- vegan::vegdist(taxon[, -1], method = "euclidean")
bd  <- vegan::betadisper(d, taxon$Taxon)
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
