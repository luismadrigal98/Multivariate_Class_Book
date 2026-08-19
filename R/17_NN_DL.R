## ============================================================================
##  MS_LJMR :: 17_NN_DL.R — Neural Networks & Deep Learning
##
##  Author: Luis J. Madrigal-Roca  (adapted from his ML course, UnitNN)
##
##  NEW session (revised Fall-2026 schedule). Integrates Luis Madrigal-Roca's
##  UnitNN (keras/tensorflow: dense nets, CNNs, transfer learning).
##  Reading: Borowiec et al. 2022 (Deep learning for biologists).
##
##  Two tracks
##  ----------
##  A. A SELF-CONTAINED dense network with the base-R package `nnet`, so the
##     session runs on any machine with no Python/keras/CUDA. We use it to teach
##     the core ideas: layers, weights, non-linear decision surfaces, training.
##  B. A keras/tensorflow TEMPLATE (dense + a small CNN) mirroring UnitNN, for
##     students with a GPU. It is guarded by requireNamespace() so sourcing this
##     file never fails when keras is absent — it only prints the recipe.
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
#   suppressPackageStartupMessages(library(p, character.only = TRUE))
# }))
# with_fig <- function(name, expr, ...) invisible(force(expr))   # draw on screen
# get_data <- function(name) { utils::data("iris");   iris }
## ---------------------------------------------------------------------------
if (!exists("get_data")) source(file.path("R", "00_utils.R"))
need("nnet")
set.seed(1998)

## ===========================================================================
##  TRACK A — a runnable multilayer perceptron with nnet
## ===========================================================================
iris  <- get_data("iris")
Xs    <- scale(as.matrix(iris[, 1:4]))
Y     <- class.ind(iris$Species)                 # one-hot targets
idx   <- sample(nrow(iris), 0.7 * nrow(iris))

## one hidden layer of 8 units, softmax output, weight decay for regularization
net <- nnet::nnet(Xs[idx, ], Y[idx, ], size = 8, softmax = TRUE,
                  decay = 5e-3, maxit = 400, trace = FALSE)

pred <- predict(net, Xs[-idx, ], type = "class")
cat("nnet test accuracy:",
    round(mean(pred == iris$Species[-idx]), 3), "\n")
cat("Network architecture: 4 inputs -> 8 hidden (logistic) -> 3 softmax;",
    length(net$wts), "weights total.\n")

## visualise the non-linear decision surface on the two petal variables
petal <- scale(as.matrix(iris[, 3:4]))
net2  <- nnet::nnet(petal, Y, size = 8, softmax = TRUE, decay = 5e-3,
                    maxit = 400, trace = FALSE)
gx <- seq(min(petal[,1]), max(petal[,1]), length = 200)
gy <- seq(min(petal[,2]), max(petal[,2]), length = 200)
grid <- expand.grid(PetalLength = gx, PetalWidth = gy)
zz   <- max.col(predict(net2, as.matrix(grid)))
with_fig("17_nn_boundary", {
  image(gx, gy, matrix(zz, 200, 200),
        col = c("#dbe2f0", "#f0dbe0", "#dbeeee"),
        xlab = "Petal length (z)", ylab = "Petal width (z)",
        main = "Neural-network decision surface (nnet)")
  points(petal, col = c("#1f3b73","#8c2d3a","#2a7f7f")[iris$Species], pch = 19)
})

## ===========================================================================
##  TRACK B — keras/tensorflow template (UnitNN style). Runs only if installed.
## ===========================================================================
if (requireNamespace("keras", quietly = TRUE)) {
  library(keras)
  ## --- dense classifier ---------------------------------------------------
  model <- keras_model_sequential() |>
    layer_dense(units = 16, activation = "relu", input_shape = 4) |>
    layer_dropout(0.2) |>
    layer_dense(units = 8, activation = "relu") |>
    layer_dense(units = 3, activation = "softmax")
  model |> compile(optimizer = "adam",
                   loss = "categorical_crossentropy", metrics = "accuracy")
  history <- model |> fit(Xs[idx, ], Y[idx, ], epochs = 60, batch_size = 16,
                          validation_split = 0.2, verbose = 0)
  plot(history)                                   # training/validation curves

  ## --- small CNN template (image track from UnitNN) -----------------------
  ## cnn <- keras_model_sequential() |>
  ##   layer_conv_2d(32, c(3,3), activation="relu", input_shape=c(64,64,3)) |>
  ##   layer_max_pooling_2d(c(2,2)) |>
  ##   layer_conv_2d(64, c(3,3), activation="relu") |>
  ##   layer_max_pooling_2d(c(2,2)) |>
  ##   layer_flatten() |> layer_dense(128, activation="relu") |>
  ##   layer_dense(n_classes, activation="softmax")
} else {
  message("keras not installed — Track B is a template only. ",
          "See UnitNN in the ML course for the full GPU workflow.")
}

cat("\n[17_NN_DL] done.\n")
