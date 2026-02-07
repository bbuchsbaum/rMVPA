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

# Linear SVM model (e1071) adapted for MVPAModels
MVPAModels$svmLinear <- list(
  type = "Classification",
  library = "e1071",
  label = "svmLinear",
  loop = NULL,
  parameters = data.frame(
    parameters = "cost",
    class = "numeric",
    label = "Cost"
  ),
  grid = function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
      if (is.null(len) || len <= 1) {
        data.frame(cost = 1)
      } else {
        # Exponentially spaced costs around 1
        exps <- seq(-2, 2, length.out = len)
        data.frame(cost = 2^exps)
      }
    } else {
      data.frame(cost = 2^runif(len, min = -3, max = 3))
    }
  },
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    if (!requireNamespace("e1071", quietly = TRUE)) {
      stop("Package 'e1071' is required for svmLinear model.")
    }
    mod <- e1071::svm(x = x, y = y,
                      kernel = "linear",
                      cost = param$cost,
                      probability = TRUE, ...)
    mod$obsLevels <- lev
    mod
  },
  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    predict(modelFit, newdata, type = "class")
  },
  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    # e1071 returns probabilities in the attr of predict output
    pred <- predict(modelFit, newdata, probability = TRUE)
    probs <- attr(pred, "probabilities")
    # Ensure columns align with training levels
    lev <- modelFit$obsLevels
    if (!is.null(lev)) {
      # Add missing columns with 0 and reorder
      missing <- setdiff(lev, colnames(probs))
      if (length(missing) > 0) {
        for (m in missing) probs <- cbind(probs, setNames(rep(0, nrow(probs)), m))
      }
      probs <- probs[, lev, drop = FALSE]
    }
    probs
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

# XGBoost model adapted from caret
MVPAModels$xgbTree <- list(
  type = c("Regression", "Classification"),
  library = c("xgboost", "plyr"),
  label = "xgbTree",
  loop = NULL,
  parameters = data.frame(
    parameters = c("nrounds", "max_depth", "eta", "gamma", "colsample_bytree", "min_child_weight", "subsample"),
    class = c("numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"),
    label = c("# Boosting Iterations", "Max Tree Depth", "Shrinkage", "Minimum Loss Reduction", 
              "Subsample Ratio of Columns", "Minimum Sum of Instance Weight", "Subsample Percentage")
  ),
  
  grid = function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
      expand.grid(
        max_depth = seq(1, len), 
        nrounds = floor((1:len) * 50), 
        eta = c(0.3, 0.4), 
        gamma = 0, 
        colsample_bytree = c(0.6, 0.8), 
        min_child_weight = c(1), 
        subsample = seq(0.5, 1, length = len)
      )
    } else {
      data.frame(
        nrounds = floor(runif(len, min = 1, max = 1000)),
        max_depth = sample(1:10, replace = TRUE, size = len),
        eta = runif(len, min = 0.001, max = 0.6),
        gamma = runif(len, min = 0, max = 10),
        colsample_bytree = runif(len, min = 0.3, max = 0.7),
        min_child_weight = sample(0:20, size = len, replace = TRUE),
        subsample = runif(len, min = 0.25, max = 1)
      )
    }
  },
  
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    # Convert to matrix if needed
    if (!inherits(x, "xgb.DMatrix")) {
      x <- as.matrix(x)
    }
    
    if (is.factor(y)) {
      # Classification
      if (length(lev) == 2) {
        # Binary classification
        y <- ifelse(y == lev[1], 1, 0)
        if (!inherits(x, "xgb.DMatrix")) {
          x <- xgboost::xgb.DMatrix(x, label = y, missing = NA)
        } else {
          xgboost::setinfo(x, "label", y)
        }
        if (!is.null(wts)) {
          xgboost::setinfo(x, "weight", wts)
        }
        out <- xgboost::xgb.train(
          list(eta = param$eta, max_depth = param$max_depth,
               gamma = param$gamma, colsample_bytree = param$colsample_bytree,
               min_child_weight = param$min_child_weight, subsample = param$subsample),
          data = x, nrounds = param$nrounds, objective = "binary:logistic", ...
        )
      } else {
        # Multiclass classification
        y <- as.numeric(y) - 1
        if (!inherits(x, "xgb.DMatrix")) {
          x <- xgboost::xgb.DMatrix(x, label = y, missing = NA)
        } else {
          xgboost::setinfo(x, "label", y)
        }
        if (!is.null(wts)) {
          xgboost::setinfo(x, "weight", wts)
        }
        out <- xgboost::xgb.train(
          list(eta = param$eta, max_depth = param$max_depth,
               gamma = param$gamma, colsample_bytree = param$colsample_bytree,
               min_child_weight = param$min_child_weight, subsample = param$subsample),
          data = x, num_class = length(lev), nrounds = param$nrounds,
          objective = "multi:softprob", ...
        )
      }
    } else {
      # Regression
      if (!inherits(x, "xgb.DMatrix")) {
        x <- xgboost::xgb.DMatrix(x, label = y, missing = NA)
      } else {
        xgboost::setinfo(x, "label", y)
      }
      if (!is.null(wts)) {
        xgboost::setinfo(x, "weight", wts)
      }
      out <- xgboost::xgb.train(
        list(eta = param$eta, max_depth = param$max_depth,
             gamma = param$gamma, colsample_bytree = param$colsample_bytree,
             min_child_weight = param$min_child_weight, subsample = param$subsample),
        data = x, nrounds = param$nrounds, objective = "reg:squarederror", ...
      )
    }
    
    # Store problem type and levels for prediction
    out$obsLevels <- lev
    out$problemType <- ifelse(is.factor(y), "Classification", "Regression")
    
    out
  },
  
  predict = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    if (!inherits(newdata, "xgb.DMatrix")) {
      newdata <- as.matrix(newdata)
      newdata <- xgboost::xgb.DMatrix(data = newdata, missing = NA)
    }
    
    out <- predict(modelFit, newdata)
    
    if (modelFit$problemType == "Classification") {
      if (length(modelFit$obsLevels) == 2) {
        # Binary classification
        out <- ifelse(out >= 0.5, modelFit$obsLevels[1], modelFit$obsLevels[2])
      } else {
        # Multiclass classification
        out <- matrix(out, ncol = length(modelFit$obsLevels), byrow = TRUE)
        out <- modelFit$obsLevels[apply(out, 1, which.max)]
      }
    }
    
    out
  },
  
  prob = function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    if (!inherits(newdata, "xgb.DMatrix")) {
      newdata <- as.matrix(newdata)
      newdata <- xgboost::xgb.DMatrix(data = newdata, missing = NA)
    }
    
    out <- predict(modelFit, newdata)
    
    if (length(modelFit$obsLevels) == 2) {
      # Binary classification
      out <- cbind(out, 1 - out)
      colnames(out) <- modelFit$obsLevels
    } else {
      # Multiclass classification
      out <- matrix(out, ncol = length(modelFit$obsLevels), byrow = TRUE)
      colnames(out) <- modelFit$obsLevels
    }
    
    as.data.frame(out, stringsAsFactors = TRUE)
  }
)
