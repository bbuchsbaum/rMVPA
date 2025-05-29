# Add caret-based models to MVPAModels registry

# Random Forest model adapted from caret
MVPAModels$rf <- list(
  type = "Classification",
  library = "randomForest",
  label = "rf",
  loop = NULL,
  parameters = data.frame(
    parameters = "mtry",
    class = "numeric",
    label = "#Randomly Selected Predictors"
  ),
  
  grid = function(x, y, len = NULL, search = "grid") {
    # Adapted from caret's grid function
    if (search == "grid") {
      # Use caret's var_seq function if available, otherwise approximate
      if (requireNamespace("caret", quietly = TRUE)) {
        data.frame(mtry = caret::var_seq(p = ncol(x), 
                                         classification = is.factor(y), 
                                         len = len))
      } else {
        # Simple approximation if caret not available
        ncol_x <- ncol(x)
        if (is.factor(y)) {
          # For classification, use sqrt(p) as default
          default_mtry <- floor(sqrt(ncol_x))
          if (is.null(len) || len == 1) {
            data.frame(mtry = default_mtry)
          } else {
            mtry_seq <- unique(floor(seq(1, ncol_x, length.out = len)))
            data.frame(mtry = mtry_seq)
          }
        } else {
          # For regression, use p/3 as default
          default_mtry <- max(floor(ncol_x/3), 1)
          if (is.null(len) || len == 1) {
            data.frame(mtry = default_mtry)
          } else {
            mtry_seq <- unique(floor(seq(1, ncol_x, length.out = len)))
            data.frame(mtry = mtry_seq)
          }
        }
      }
    } else {
      # Random search
      data.frame(mtry = unique(sample(1:ncol(x), size = len, replace = TRUE)))
    }
  },
  
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    # Direct adaptation from caret
    model <- randomForest::randomForest(x, y, mtry = param$mtry, ...)
    # Store levels for proper prediction
    model$obsLevels <- lev
    model
  },
  
  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    # Direct adaptation from caret
    if (!is.null(newdata)) {
      predict(modelFit, newdata)
    } else {
      predict(modelFit)
    }
  },
  
  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    # Direct adaptation from caret
    if (!is.null(newdata)) {
      predict(modelFit, newdata, type = "prob")
    } else {
      predict(modelFit, type = "prob")
    }
  }
)

# Sparse PLS model adapted from caret
MVPAModels$spls <- list(
  type = c("Regression", "Classification"),
  library = "spls",
  label = "spls",
  loop = NULL,
  parameters = data.frame(
    parameters = c("K", "eta", "kappa"),
    class = c("numeric", "numeric", "numeric"),
    label = c("#Components", "Threshold", "Kappa")
  ),
  
  grid = function(x, y, len = NULL, search = "grid") {
    # Direct adaptation from caret
    if (search == "grid") {
      expand.grid(K = 1:min(nrow(x), ncol(x)), 
                  eta = seq(0.1, 0.9, length = len), 
                  kappa = 0.5)
    } else {
      data.frame(kappa = runif(len, min = 0, max = 0.5),
                 eta = runif(len, min = 0, max = 1),
                 K = sample(1:min(nrow(x), ncol(x)), size = len, replace = TRUE))
    }
  },
  
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    # Ensure K doesn't exceed number of observations
    param$K <- min(param$K, length(y))
    
    if (is.factor(y)) {
      # For classification, we need to implement splsda
      # Since caret uses internal function, we'll need to either:
      # 1. Use mixOmics package which has splsda, or
      # 2. Convert to numeric and use spls with post-processing
      
      # Option 2: Convert factor to numeric for spls
      y_numeric <- as.numeric(y)
      model <- spls::spls(x, y_numeric, K = param$K, eta = param$eta, kappa = param$kappa, ...)
      # Store original levels for prediction
      model$obsLevels <- lev
      model$is_classification <- TRUE
      model
    } else {
      # For regression, use spls directly
      model <- spls::spls(x, y, K = param$K, eta = param$eta, kappa = param$kappa, ...)
      model$is_classification <- FALSE
      model
    }
  },
  
  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    if (isTRUE(modelFit$is_classification)) {
      # For classification, predict numeric values and convert to classes
      pred_numeric <- predict(modelFit, newdata)
      # Round to nearest integer and convert to factor levels
      pred_class <- round(pred_numeric)
      pred_class[pred_class < 1] <- 1
      pred_class[pred_class > length(modelFit$obsLevels)] <- length(modelFit$obsLevels)
      modelFit$obsLevels[pred_class]
    } else {
      # For regression, use standard prediction
      predict(modelFit, newdata)
    }
  },
  
  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    if (isTRUE(modelFit$is_classification)) {
      # Simple probability estimation based on distance to class centers
      pred_numeric <- predict(modelFit, newdata)
      n_classes <- length(modelFit$obsLevels)
      
      # Create probability matrix
      probs <- matrix(0, nrow = nrow(newdata), ncol = n_classes)
      colnames(probs) <- modelFit$obsLevels
      
      # Simple softmax-like probability assignment
      for (i in seq_len(nrow(probs))) {
        # Distance to each class (1, 2, ..., n_classes)
        distances <- abs(pred_numeric[i] - seq_len(n_classes))
        # Convert distances to probabilities
        probs[i,] <- exp(-distances) / sum(exp(-distances))
      }
      
      probs
    } else {
      # For regression, probabilities don't make sense
      stop("Probabilities are not available for regression models")
    }
  }
)

# Alternative SPLS implementation using mixOmics if available
MVPAModels$spls_mixomics <- list(
  type = "Classification",
  library = "mixOmics",
  label = "spls_mixomics",
  loop = NULL,
  parameters = data.frame(
    parameters = c("ncomp", "keepX"),
    class = c("numeric", "numeric"),
    label = c("#Components", "#Variables to keep")
  ),
  
  grid = function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
      expand.grid(ncomp = 1:min(10, ncol(x)), 
                  keepX = floor(seq(ncol(x)/10, ncol(x), length = len)))
    } else {
      data.frame(ncomp = sample(1:min(10, ncol(x)), size = len, replace = TRUE),
                 keepX = sample(10:ncol(x), size = len, replace = TRUE))
    }
  },
  
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    # Use mixOmics splsda for classification
    keepX_vec <- rep(param$keepX, param$ncomp)
    model <- mixOmics::splsda(X = x, Y = y, ncomp = param$ncomp, keepX = keepX_vec, ...)
    model$obsLevels <- lev
    model
  },
  
  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    pred <- predict(modelFit, newdata)
    # mixOmics predict returns a list with class predictions
    as.character(pred$class$max.dist[, ncol(pred$class$max.dist)])
  },
  
  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    pred <- predict(modelFit, newdata)
    # Return the posterior probabilities
    pred$predict[,,modelFit$ncomp]
  }
)

# Add model aliases for backward compatibility
MVPAModels$sda <- MVPAModels$sda_notune
MVPAModels$glmnet <- MVPAModels$glmnet_opt
