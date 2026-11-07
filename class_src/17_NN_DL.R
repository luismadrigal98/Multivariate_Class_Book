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

##  HAND-OWNED -- this file is maintained by hand and scripts/make_class_src.py
##  will NOT regenerate it. Fixes made in R/ do not reach it automatically; apply
##  them here too. Remove this line to hand the file back to the generator.
##
##  Hand-owned because: Track C, the full flower-CNN training run from the ML
##  course, exists only in this class edition.
##  PLAIN CLASSROOM EDITION
##  Generated from R/17_NN_DL.R by scripts/make_class_src.py --
##  edit R/17_NN_DL.R and regenerate; changes made here will be overwritten.
## ============================================================================

# R packages required
#install.packages("nnet")
#BiocManager::install("EBImage")
library(nnet)
library(EBImage)
library(ggplot2)
library(dplyr)
library(tidyr)

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


##  WHY THIS IS OFF BY DEFAULT
##  Track C is a full training run, not a demonstration: it decodes 2746 JPEGs
##  (~1.2 GB in memory) and then trains a CNN for up to 30 epochs. Measured on a
##  62 GB workstation that is ~20 s to load and ~20 min to train, with the R
##  console blocked the whole time; on a teaching laptop it is much longer and
##  can exhaust memory. Until now it ran every time this file was sourced, even
##  when you only wanted Tracks A and B -- which is what makes the session look
##  frozen mid-class.
run_track_c <- FALSE          # set TRUE when you actually want to train

image_main_dir <- "data/Flower_data_set"

## ---- function definitions (nothing runs yet) --------------------------------

preprocess_images <- function(base_path, target_width = 320, 
                              target_height = 240) {
  # Get all directories (flower categories)
  categories <- list.dirs(base_path, full.names = FALSE, recursive = FALSE)
  
  # Pre-size the containers. Growing them with c() inside the loop copies the
  # whole vector on every iteration, and the loop runs once per image.
  all_files <- unlist(lapply(categories, function(cg)
    list.files(file.path(base_path, cg),
               pattern = "\\.jpg|\\.jpeg|\\.png|\\.gif", full.names = TRUE)))
  n_tot      <- length(all_files)
  images     <- vector("list", n_tot)
  labels     <- character(n_tot)
  file_paths <- character(n_tot)
  k          <- 0L
  
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
        k             <- k + 1L
        images[[k]]   <- img
        labels[k]     <- category
        file_paths[k] <- img_path
        if (k %% 250 == 0) cat("  ", k, "/", n_tot, "images\n")
        
      }, error = function(e) {
        warning(sprintf("Error processing %s: %s", img_path, e$message))
      })
    }
  }
  
  # drop any slots left empty by unreadable files
  if (k < n_tot) {
    images <- images[seq_len(k)]; labels <- labels[seq_len(k)]
    file_paths <- file_paths[seq_len(k)]
  }

  # Convert images list to array.
  #
  # NOTE: the obvious one-liner is WRONG --
  #     array(unlist(images), dim = c(n, target_height, target_width, 3))
  # unlist() lays each image out in its own column-major order, but the array
  # is filled with the FIRST index varying fastest, i.e. across IMAGES. The
  # result interleaves every picture into every slice, so each "image" the CNN
  # sees is a mixture of all of them. (Check it on three 1-colour images: every
  # slice comes out striped.) EBImage also returns (width, height, channel),
  # so each image needs transposing to (height, width, channel) as well.
  n <- length(images)
  image_array <- array(0, dim = c(n, target_height, target_width, 3))
  for (i in seq_len(n))
    image_array[i, , , ] <- aperm(images[[i]], c(2, 1, 3))
  rm(images)                      # ~1.2 GB; drop it before the caller copies
  
  # Convert labels to factor and then to numeric
  labels <- as.numeric(factor(labels)) - 1
  
  return(list(
    images = image_array,
    labels = labels,
    categories = categories,
    file_paths = file_paths
  ))
}

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

create_model_flower_CNN2 <- function(num_classes) {
  model <- keras_model_sequential() |>
    
    # Full model ("The Beast")
    
    # First Convolutional Block
    layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu') |>
    layer_batch_normalization() |>
    layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu') |>
    layer_max_pooling_2d(pool_size = c(2,2)) |>
    layer_dropout(0.25) |>

    # Second Convolutional Block
    layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') |>
    layer_batch_normalization() |>
    layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') |>
    layer_max_pooling_2d(pool_size = c(2,2)) |>
    layer_dropout(0.25) |>

    # Third Convolutional Block
    layer_conv_2d(filters = 128, kernel_size = c(3,3), activation = 'relu') |>
    layer_batch_normalization() |>
    layer_conv_2d(filters = 128, kernel_size = c(3,3), activation = 'relu') |>
    layer_max_pooling_2d(pool_size = c(2,2)) |>
    layer_dropout(0.25) |>

    # Leaving the two blocks above commented out does NOT make the model
    # lighter -- it makes it 5.4x HEAVIER (9.28M parameters against 1.73M),
    # because layer_flatten() then sees a 58x78x32 feature map instead of a
    # pooled-down one, and the dense layer after it carries that cost.
    
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
    ## verbose = 1 draws a per-BATCH progress bar with carriage returns -- 55
    ## redraws per epoch. RStudio's console handles that badly and the IDE can
    ## stop responding. verbose = 2 prints one tidy line per epoch instead.
    verbose = 2,
    ...
  )
  
  return(history)
}

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

    # ---- confusion matrix -------------------------------------------------
    # The highest-probability class is the prediction. Both vectors are made
    # factors over the SAME levels so the table stays square even when a class
    # is never predicted -- otherwise diag() silently misaligns.
    pred_class <- factor(labels[max.col(prob_predictions, ties.method = "first")],
                         levels = labels)
    obs_class  <- factor(labels[y + 1], levels = labels)
    cm <- table(Observed = obs_class, Predicted = pred_class)

    evaluation_results[['confusion_matrix']] <- cm
    evaluation_results[['accuracy']]         <- sum(diag(cm)) / sum(cm)
    evaluation_results[['per_class']] <- data.frame(
      class     = labels,
      n         = as.integer(rowSums(cm)),
      recall    = round(diag(cm) / pmax(1, rowSums(cm)), 3),  # of the true X, how many found
      precision = round(diag(cm) / pmax(1, colSums(cm)), 3),  # of those called X, how many right
      row.names = NULL)
  } else 
  {
    evaluation_results[['probabilistic_prediction']] <- NULL
  }
  
  return(evaluation_results)
} 

## ---- the expensive part, only on request ------------------------------------
if (!run_track_c) {
  message("Track C skipped. Set run_track_c <- TRUE to load the images and ",
          "train the CNN (~20 min, ~1.2 GB).")
} else if (!requireNamespace("keras3", quietly = TRUE)) {
  message("Track C needs keras3 with a working Python TensorFlow behind it.")
} else if (!dir.exists(image_main_dir)) {
  message("Track C needs the flower images in ", image_main_dir, ".")
} else {
  library(keras3)        # Track C attaches keras3 itself rather than relying
                         # on Track B having happened to run first
  cat("Loading and preprocessing images...\n")
  Data_flowers <- preprocess_images(image_main_dir, 160, 120)
  cat("Splitting data into training and testing sets...\n")
  split <- split_data(Data_flowers)
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
  model_CNN <- create_model_flower_CNN2(length(unique(Data_flowers$labels)))
  fitted_model <- model_fitter(model_CNN, X = split$train$images, 
                               Y_train = split$train$labels, 
                               callbacks_list = callbacks_list, epochs = 30, 
                               batch_size = 32)
  X_test <- split$test$images
  Y_test <- split$test$labels
  labels <- Data_flowers$categories
  evaluation <- model_evaluator(model_CNN, x = X_test, y = Y_test, 
                                labels = labels, predict = T)
  ## ---- how well is it actually doing? ------------------------------------
  cm <- evaluation$confusion_matrix
  cat("\nConfusion matrix (rows = observed, columns = predicted):\n")
  print(cm)
  cat(sprintf("\noverall accuracy: %.3f   (chance with %d equal classes = %.3f)\n",
              evaluation$accuracy, nrow(cm), 1 / nrow(cm)))
  cat("\nper class:\n"); print(evaluation$per_class, row.names = FALSE)

  ## Accuracy alone hides WHICH classes are confused for which. The off-diagonal
  ## cells are the interesting part: a mutually confused pair is a pair these
  ## features cannot separate.
  cm_df      <- as.data.frame(cm)
  cm_df$Prop <- cm_df$Freq / rep(rowSums(cm), times = ncol(cm))
  print(ggplot(cm_df, aes(x = Predicted, y = Observed, fill = Prop)) +
          geom_tile(colour = "white", linewidth = 1) +
          geom_text(aes(label = Freq), colour = "grey15", size = 4) +
          scale_fill_gradient(low = "#f0f2f5", high = "#1f3b73", limits = c(0, 1),
                              name = "proportion of\nobserved row") +
          theme_bw(base_size = 12) + theme(panel.grid = element_blank()) +
          labs(x = "Predicted class", y = "Observed class"))

  prob_df <- evaluation$probabilistic_prediction
  prob_df$Predicted_Class <- colnames(prob_df)[max.col(prob_df, ties.method = "first")]
  prob_df$Observation <- 1:nrow(prob_df)
  df_long <- prob_df |> 
    dplyr::arrange(Predicted_Class, Observation) |> 
    dplyr::mutate(Observation_Rank = row_number()) |> # Creates a clean x-axis order
    dplyr::select(-Predicted_Class, -Observation) |> 
    pivot_longer(cols = -Observation_Rank, names_to = "Class", values_to = "Probability")
  ggplot(df_long, aes(x = Observation_Rank, y = Probability, fill = Class)) +
    geom_col(width = 1) +
    scale_fill_viridis_d(option = "turbo") +
    labs(title = "Class Probabilities per Image", x = "Images (Sorted by Prediction)", 
         y = "Probability") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggplot(df_long, aes(x = Class, y = Observation_Rank, fill = Probability)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
    labs(title = "Probability Heatmap", x = "Flower Class", y = "Test Set Images") +
    theme_bw() +
    theme(axis.text.y = element_blank(), panel.grid = element_blank())
  prob_df$True_Label <- split$test$labels
  density_df <- prob_df |> 
    dplyr::select(-Predicted_Class, -Observation) |>
    pivot_longer(cols = -True_Label, names_to = "Class", 
                 values_to = "Probability")
  ggplot(density_df, aes(x = Probability, fill = Class)) +
    geom_density(alpha = 0.6, color = NA) +
    facet_wrap(~ Class, scales = "free_y") +
    labs(title = "Probability Distributions by Class", x = "Assigned Probability", 
         y = "Density") +
    theme_bw()
  cat("\n[17_NN_DL] done.\n")
}

cat("
[17_NN_DL] done.
")
