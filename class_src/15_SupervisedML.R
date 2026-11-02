## ============================================================================
##  MS_LJMR :: 15_SupervisedML.R — Supervised Machine Learning
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

##  PLAIN CLASSROOM EDITION
##  Generated from R/15_SupervisedML.R by scripts/make_class_src.py --
##  edit R/15_SupervisedML.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("rpart")
#install.packages("rpart.plot")
#install.packages("randomForest")
#install.packages("gbm")
#install.packages("caret")
#install.packages("mlbench")
library(rpart)
library(rpart.plot)
library(randomForest)
library(gbm)
library(caret)
library(mlbench)

# Built-in data sets used below
data(iris)
data(Vehicle)

set.seed(1998)                                   # reproducibility (as in ML units)

## ---- Data & train/test split -----------------------------------------------

idx   <- caret::createDataPartition(iris$Species, p = 0.7, list = FALSE)
train <- iris[idx, ]; test <- iris[-idx, ]

## ---- 1. A single CART tree -------------------------------------------------
fit_cart <- rpart::rpart(Species ~ ., data = train, method = "class",
                         control = rpart::rpart.control(cp = 0.01))
par(mfrow = c(1, 1))
rpart.plot::rpart.plot(fit_cart, box.palette = "BuGn", main = "CART decision tree")
acc <- function(model, newdata, truth = newdata$Species, type = "class")
  mean(predict(model, newdata, type = type) == truth)
cat("Single CART test accuracy:", round(acc(fit_cart, test), 3), "\n")

## ---- 2. Random forest (bagging of de-correlated trees) ---------------------
fit_rf <- randomForest::randomForest(Species ~ ., data = train,
                                     ntree = 500, importance = TRUE)
cat("Random-forest test accuracy:",
    round(mean(predict(fit_rf, test) == test$Species), 3), "\n")
par(mfrow = c(1, 1))
randomForest::varImpPlot(fit_rf, main = "Random-forest variable importance")

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

par(mfrow = c(1, 1))
## dotplot() returns a lattice object; inside a function it must be printed
## explicitly or nothing is drawn and no file is written.
print(dotplot(res, metric = "Accuracy",
              main = "Tuned model comparison (5-fold CV)"))

## ---- 4. Confusion matrix of the best model on the held-out test set --------
best <- grids[[ which.max(vapply(grids, function(g) max(g$results$Accuracy), 0)) ]]
cm   <- caret::confusionMatrix(predict(best, test), test$Species)
cat("\nBest model:", best$method, "\n"); print(cm$table)
cat("Test accuracy:", round(cm$overall["Accuracy"], 3), "\n")

## ---- 5. MLBENCH dataset

## ---- Data & train/test split -----------------------------------------------
idx   <- caret::createDataPartition(Vehicle$Class, p = 0.7, list = FALSE)
train <- Vehicle[idx, ]; test <- Vehicle[-idx, ]

## ---- 1. A single CART tree -------------------------------------------------
fit_cart <- rpart::rpart(Class ~ ., data = train, method = "class",
                         control = rpart::rpart.control(cp = 0.01))
rpart.plot::rpart.plot(fit_cart, box.palette = "BuGn", main = "CART decision tree")

acc <- function(model, newdata, truth = newdata$Class, type = "class")
  mean(predict(model, newdata, type = type) == truth)
cat("Single CART test accuracy:", round(acc(fit_cart, test), 3), "\n")

## ---- 2. Random forest (bagging of de-correlated trees) ---------------------
fit_rf <- randomForest::randomForest(Class ~ ., data = train,
                                     ntree = 500, importance = TRUE)
cat("Random-forest test accuracy:",
    round(mean(predict(fit_rf, test) == test$Class), 3), "\n")
randomForest::varImpPlot(fit_rf, n.var = 10, main = "Random-forest variable importance")

## ---- 3. caret: one interface, honest resampling, tuning --------------------
## Seed the CV folds or set seed right before each train() call for exact fold parity
ctrl <- caret::trainControl(method = "cv", number = 5)

set.seed(1998)
fit_cv_rpart <- caret::train(Class ~ ., data = train, method = "rpart",
                             trControl = ctrl, tuneLength = 6)
set.seed(1998)
fit_cv_bag   <- caret::train(Class ~ ., data = train, method = "treebag",
                             trControl = ctrl)
set.seed(1998)
fit_cv_rf    <- caret::train(Class ~ ., data = train, method = "rf",
                             trControl = ctrl, tuneLength = 4)
set.seed(1998)
fit_cv_gbm   <- caret::train(Class ~ ., data = train, method = "gbm",
                             trControl = ctrl, verbose = FALSE)

grids <- list(rpart = fit_cv_rpart, treebag = fit_cv_bag, rf = fit_cv_rf, gbm = fit_cv_gbm)

print(grids)

res2 <- caret::resamples(grids)

print(dotplot(res2, metric = "Accuracy",
              main = "Tuned model comparison (5-fold CV)"))

cat("\n[15_SupervisedML] done.\n")