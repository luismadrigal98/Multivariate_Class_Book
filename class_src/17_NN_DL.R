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

##  PLAIN CLASSROOM EDITION
##  Generated from R/17_NN_DL.R by scripts/make_class_src.py --
##  edit R/17_NN_DL.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("nnet")
#BiocManager::install("EBImage")
library(nnet)
library(EBImage)

# Built-in data sets used below
data(iris)

set.seed(1998)

## ===========================================================================
##  TRACK A — a runnable multilayer perceptron with nnet
## ===========================================================================

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
par(mfrow = c(1, 1))
image(gx, gy, matrix(zz, 200, 200),
      col = c("#dbe2f0", "#f0dbe0", "#dbeeee"),
      xlab = "Petal length (z)", ylab = "Petal width (z)",
      main = "Neural-network decision surface (nnet)")
points(petal, col = c("#1f3b73","#8c2d3a","#2a7f7f")[iris$Species], pch = 19)

## ===========================================================================
##  TRACK B — keras/tensorflow template (UnitNN style). Runs only if installed.
## ===========================================================================
if (requireNamespace("keras3", quietly = TRUE)) {
  library(keras3)
  ## --- dense classifier ---------------------------------------------------
  model <- keras_model_sequential() |>
    layer_dense(units = 16, activation = "relu") |>
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

## ===========================================================================
##  TRACK C — keras/tensorflow template (UnitNN style). Serious example
## ===========================================================================

image_main_dir <- "data/Flower_data_set"

preprocess_images <- function(base_path, target_width = 320, 
                              target_height = 240) {
  # Get all directories (flower categories)
  categories <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)
  
  # Initialize lists to store images and labels
  images <- list()
  labels <- vector()
  file_paths <- vector()
  
  # Loop through each category
  for (category in categories) {
    cat("Processing category:", category, "\n")
    
    # Get all image files in the category
    img_files <- list.files(file.path(base_path, category), 
                            pattern = "\\.jpg|\\.jpeg|\\.png|\\.gif",
                            full.names = TRUE,
                            recursive = FALSE)
    
    # Process each image
    for (img_path in img_files) {
      tryCatch({
        # Read image
        img <- readImage(img_path)
        
        # Ensure image has 3 channels (RGB)
        if (length(dim(img)) == 2) {
          # If grayscale, convert to RGB by replicating the channel
          img <- array(img, dim = c(dim(img), 3))
        } else if (dim(img)[3] == 4) {
          # If RGBA, remove alpha channel
          img <- img[,,1:3]
        }
        
        # Resize image
        img <- resize(img, w = target_width, h = target_height)
        
        # Normalize pixel values to [0,1]
        img <- normalize(img)
        
        # Store image and label
        images[[length(images) + 1]] <- img
        labels <- c(labels, category)
        file_paths <- c(file_paths, img_path)
        
      }, error = function(e) {
        warning(sprintf("Error processing %s: %s", img_path, e$message))
      })
    }
  }
  
  # Convert images list to array
  image_array <- array(unlist(images), 
                       dim = c(length(images), target_height, target_width, 3))
  
  # Convert labels to factor and then to numeric
  labels <- as.numeric(factor(labels)) - 1
  
  return(list(
    images = image_array,
    labels = labels,
    categories = categories,
    file_paths = file_paths
  ))
}



# Process all images
cat("Loading and preprocessing images...\n")
Data_flowers <- preprocess_images(image_main_dir, 160, 120)

split_data <- function(data, test_size = 0.2, seed = 42) {
  set.seed(seed)
  
  # Create index for splitting
  total_samples <- length(data$labels)
  indices <- sample(total_samples)
  test_samples <- round(total_samples * test_size)
  
  # Split data
  train_indices <- indices[1:(total_samples - test_samples)]
  test_indices <- indices[(total_samples - test_samples + 1):total_samples]
  
  # Create training and testing sets
  train_data <- list(
    images = data$images[train_indices,,,],
    labels = data$labels[train_indices],
    path = data$file_paths[train_indices]
  )
  
  test_data <- list(
    images = data$images[test_indices,,,],
    labels = data$labels[test_indices],
    path = data$file_paths[test_indices]
  )
  
  return(list(
    train = train_data,
    test = test_data
  ))
}

# Split data
cat("Splitting data into training and testing sets...\n")
split <- split_data(Data_flowers)

# Prepare hyperparameters

callbacks_list <- list(
  callback_early_stopping(
    monitor = "val_loss",
    patience = 10,
    restore_best_weights = TRUE
  ),
  callback_reduce_lr_on_plateau(
    monitor = "val_loss",
    factor = 0.5,
    patience = 5,
    min_lr = 1e-6
  ),
  callback_model_checkpoint(
    filepath = "best_model.keras",
    monitor = "val_accuracy",
    save_best_only = TRUE
  )
)

lr_schedule <- learning_rate_schedule_exponential_decay(
  initial_learning_rate = 0.001,
  decay_steps = 1000,
  decay_rate = 0.8
)

create_model_flower_CNN2 <- function(num_classes) {
  model <- keras_model_sequential() |>
    
    # Full model ("The Beast")
    
    # First Convolutional Block
    layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu') |>
    layer_batch_normalization() |>
    layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu') |>
    layer_max_pooling_2d(pool_size = c(2,2)) |>
    layer_dropout(0.25) |>
    
    # # Second Convolutional Block
    # layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') |>
    # layer_batch_normalization() |>
    # layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') |>
    # layer_max_pooling_2d(pool_size = c(2,2)) |>
    # layer_dropout(0.25) |>
    # 
    # # Third Convolutional Block
    # layer_conv_2d(filters = 128, kernel_size = c(3,3), activation = 'relu') |>
    # layer_batch_normalization() |>
    # layer_conv_2d(filters = 128, kernel_size = c(3,3), activation = 'relu') |>
    # layer_max_pooling_2d(pool_size = c(2,2)) |>
    # layer_dropout(0.25) |>
    
    # Dense layers
    layer_flatten() |>
    layer_dense(units = 64, activation = 'relu') |>
    layer_batch_normalization() |>
    layer_dropout(0.5) |>
    layer_dense(units = num_classes, activation = 'softmax')
  
  model |> compile(
    optimizer = optimizer_adam(learning_rate = 0.001),
    loss = 'sparse_categorical_crossentropy',
    metrics = c('accuracy')
  )
  
  return(model)
}

model_CNN <- create_model_flower_CNN2(length(unique(Data_flowers$labels)))

model_fitter <- function(model, X, Y_train, 
                         validation_prop = 0.2, epochs = 100, batch_size = 32,
                         callbacks_list, ...) 
{
  #' Fit a keras model
  #' 
  #' @param model A keras model
  #' @param X A matrix of predictors
  #' @param Y_train A matrix or list of responses
  #' @param callbacks_list A list of callbacks
  #' @param validation_prop Proportion of data to use for validation
  #' @param epochs Number of training epochs
  #' @param batch_size Batch size for training
  #' @param ... Additional arguments passed to keras::fit()
  #' @return A fitted keras model
  #' ___________________________________________________________________________
  
  history <- model %>% fit(
    x = X,
    y = Y_train,
    validation_split = validation_prop,
    callbacks = callbacks_list,
    epochs = epochs,
    batch_size = batch_size,
    verbose = 1,
    ...
  )
  
  return(history)
}

fitted_model <- model_fitter(model_CNN, X = split$train$images, 
                             Y_train = split$train$labels, 
                             callbacks_list = callbacks_list, epochs = 100, 
                             batch_size = 32)

# Evaluate the model

X_test <- split$test$images
Y_test <- split$test$labels
labels <- Data_flowers$categories

model_evaluator <- function(model, x, y, labels = NULL, predict = F, ...)
{
  #' Wrapper of the evaluate function in keras to assess the performance of the
  #' models using unseen data.
  #' 
  #' @param model A keras model.
  #' @param x An array of predictors (testing slot)
  #' @param y An array of target values (testing slot)
  #' @param labels A vector with the labels of the target variable
  #' @param predict Whether to do a probabilistic prediction in classification
  #' problems
  #' 
  #' @return A list with evaluation info.
  #' ___________________________________________________________________________
  
  evaluation_results <- list()
  
  evaluation_results[['overall_evaluation']] <- model |> 
    keras3::evaluate(x, y, ...)
  
  if (predict) 
  {
    # Check that the labels are provided
    if (is.null(labels))
    {
      stop("Labels must be provied for making a probabilistic prediction")
    }
    
    # Get probabilistic predictions
    prob_predictions <- model |> predict(x)
    
    # Create a data frame with probabilities and corresponding labels
    prob_df <- as.data.frame(prob_predictions)
    names(prob_df) <- labels
    
    evaluation_results[['probabilistic_prediction']] <- prob_df
  } else 
  {
    evaluation_results[['probabilistic_prediction']] <- NULL
  }
  
  return(evaluation_results)
} 

evaluation <- model_evaluator(model_CNN, x = X_test, y = Y_test, 
                              labels = labels, predict = T)

heatmap(x = as.matrix(evaluation$probabilistic_prediction), 
        keep.dendro = F, )

cat("\n[17_NN_DL] done.\n")
