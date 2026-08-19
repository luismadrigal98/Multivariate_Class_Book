## ============================================================================
##  MS_LJMR :: 25_SupervisedML.R — Supervised Machine Learning
##
##  Author: Luis J. Madrigal-Roca  (adapted from his ML course, UnitCART & UnitCARET)
##
##  NEW session (revised Fall-2026 schedule). Integrates Luis Madrigal-Roca's
##  ML-course units:  UnitCART (rpart / bagging / random forest / boosting) and
##  UnitCARET (unified training & tuning with caret).
##  Readings: James et al. 2021 (ISLR) Ch. 4 & 8; Cutler et al. 2007.
##
##  Storyline
##  ---------
##  1. A single decision tree (CART): interpretable but high-variance.
##  2. Ensembles reduce variance/bias: bagging -> random forest -> boosting.
##  3. caret gives one interface for resampling, tuning and honest comparison.
##
##  Everything is self-contained: iris (built-in) is the running example so it
##  reproduces without external data; swap in Butterflies / Limenitis via
##  get_data() for the ecological version used in class.
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
##  Files to keep next to this script: none — the data sets used here ship with R itself
##
# need <- function(...) invisible(lapply(c(...), function(p) {
#   if (!requireNamespace(p, quietly = TRUE))
#     stop("This session needs: install.packages(\"", p, "\")", call. = FALSE)
#   library(p, character.only = TRUE)
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) { utils::data("iris");   iris }
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("rpart", "rpart.plot", "randomForest", "gbm", "caret")
set.seed(1998)                                   # reproducibility (as in ML units)

## ---- Data & train/test split -----------------------------------------------
iris  <- get_data("iris")
idx   <- caret::createDataPartition(iris$Species, p = 0.7, list = FALSE)
train <- iris[idx, ]; test <- iris[-idx, ]

## ---- 1. A single CART tree -------------------------------------------------
fit_cart <- rpart::rpart(Species ~ ., data = train, method = "class",
                         control = rpart::rpart.control(cp = 0.01))
with_fig("25_cart_tree", {
  rpart.plot::rpart.plot(fit_cart, box.palette = "BuGn", main = "CART decision tree")
})
acc <- function(model, newdata, truth = newdata$Species, type = "class")
  mean(predict(model, newdata, type = type) == truth)
cat("Single CART test accuracy:", round(acc(fit_cart, test), 3), "\n")

## ---- 2. Random forest (bagging of de-correlated trees) ---------------------
fit_rf <- randomForest::randomForest(Species ~ ., data = train,
                                     ntree = 500, importance = TRUE)
cat("Random-forest test accuracy:",
    round(mean(predict(fit_rf, test) == test$Species), 3), "\n")
with_fig("25_rf_importance", {
  randomForest::varImpPlot(fit_rf, main = "Random-forest variable importance")
})

## ---- 3. caret: one interface, honest resampling, tuning --------------------
ctrl <- caret::trainControl(method = "cv", number = 5)
grids <- list(
  rpart = caret::train(Species ~ ., data = train, method = "rpart",
                       trControl = ctrl, tuneLength = 6),
  rf    = caret::train(Species ~ ., data = train, method = "rf",
                       trControl = ctrl, tuneLength = 3),
  gbm   = caret::train(Species ~ ., data = train, method = "gbm",
                       trControl = ctrl, verbose = FALSE))

## compare the tuned models on identical CV folds
res <- caret::resamples(grids)
cat("\nCross-validated accuracy (caret resamples):\n")
print(summary(res)$statistics$Accuracy[, c("Mean", "Median", "Max.")])

with_fig("25_caret_compare", {
  dotplot(res, metric = "Accuracy", main = "Tuned model comparison (5-fold CV)")
})

## ---- 4. Confusion matrix of the best model on the held-out test set --------
best <- grids[[ which.max(vapply(grids, function(g) max(g$results$Accuracy), 0)) ]]
cm   <- caret::confusionMatrix(predict(best, test), test$Species)
cat("\nBest model:", best$method, "\n"); print(cm$table)
cat("Test accuracy:", round(cm$overall["Accuracy"], 3), "\n")

cat("\n[25_SupervisedML] done.\n")
