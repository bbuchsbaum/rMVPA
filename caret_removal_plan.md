Okay, this is an excellent and well-thought-out "surgery plan." It correctly identifies the key `caret` touchpoints and proposes very sensible, phased replacements. My additions will focus on formalizing this into a step-by-step project plan, adding a bit more detail to each replacement, considering potential edge cases, and ensuring the new helper functions are well-defined.

## Formalized Plan: Excising `caret` from `rMVPA`

**Overall Strategy:**
We will follow the phased approach outlined in the "surgery plan." Each phase will aim to replace a distinct piece of `caret` functionality, be tested thoroughly, and keep the package in a runnable state. The primary replacements will leverage `rsample` for resampling and `yardstick` for metrics, both from the `tidymodels` ecosystem, which are lightweight and modern. For tuning, we'll start with a custom loop (Option A) and consider `tune` (Option B) as a future enhancement.

**Phase 0: Preparation**
1.  **Branching:** Create a new feature branch in Git (e.g., `feature/remove-caret`).
2.  **Dependency Update (Anticipatory):** Add `rsample` and `yardstick` to `Suggests` in the `DESCRIPTION` file for now. We'll move them to `Imports` as they become strictly necessary. This allows for gradual integration.

---

**Phase 1: Decouple Model Discovery and Basic Infrastructure**

*   **Target:** `load_model()` and `glmnet_opt`'s fold creation.
*   **Goal:** Remove direct `caret` calls for model loading and the simplest resampling case.

1.  **Modify `load_model()` (R/common.R):**
    *   **`caret` call:** `if (length(caret::getModelInfo(name)) > 0) …`
    *   **Replacement:** Remove this `else if` block entirely.
    *   **New Logic:**
        ```R
        load_model <- function(name) {
          registry <- MVPAModels # Existing environment
          
          ret <- if (!is.null(registry[[name]])) {
            registry[[name]]   
          } else {
            # Removed caret fallback
            stop(paste0("Model '", name, "' not found in MVPAModels. ",
                        "Register custom models using `register_mvpa_model()` if needed."))
          }
          
          ret$label <- name # Existing logic
          ret
        }
        ```
    *   **Rationale:** `MVPAModels` already covers all tested and exemplified models. This change makes the model source explicit.

2.  **Implement `create_mvpa_folds()` Helper (New file, e.g., R/resampling_helpers.R):**
    *   **Purpose:** Replace `caret::createFolds`.
    *   **Implementation:**
        ```R
        #' Create Cross-Validation Folds
        #'
        #' Generates a list of row indices for k-fold cross-validation.
        #' Optionally stratified by a factor response.
        #'
        #' @param y A vector, typically the response variable. If a factor,
        #'        folds can be stratified.
        #' @param k Integer, the number of folds.
        #' @param stratified Logical, whether to perform stratified sampling if y is a factor.
        #' @param seed Optional integer for reproducible fold creation.
        #' @return A list of k integer vectors, each containing row indices for one fold.
        #' @importFrom rsample vfold_cv analysis
        #' @keywords internal
        create_mvpa_folds <- function(y, k = 5, stratified = TRUE, seed = NULL) {
          if (!is.null(seed)) set.seed(seed)
          n <- length(y)
          
          if (is.factor(y) && stratified) {
            # Use rsample for robust stratified k-fold CV
            # Create a dummy data frame for rsample
            df_for_rsample <- data.frame(.response_var_for_stratification = y, .indices = seq_len(n))
            folds_obj <- rsample::vfold_cv(df_for_rsample, v = k, 
                                           strata = .response_var_for_stratification, repeats = 1)
            # Extract the assessment (hold-out) indices for each fold
            # rsample returns analysis (training) and assessment (testing) sets.
            # createFolds returns indices for *each fold directly* (i.e., the test set for that fold).
            # So we extract the assessment indices.
            out_indices <- lapply(folds_obj$splits, function(split) rsample::assessment(split)$.indices)
          } else {
            # Simple random k-fold CV
            shuffled_indices <- sample(seq_len(n))
            out_indices <- split(shuffled_indices, cut(seq_along(shuffled_indices), breaks = k, labels = FALSE))
          }
          return(out_indices)
        }
        ```
    *   **Note:** `rsample::vfold_cv` returns `rset` objects. We extract the *assessment* indices to mimic `caret::createFolds` which returns hold-out indices for each fold.

3.  **Patch `glmnet_opt` Model (R/classifiers.R):**
    *   **`caret` call:** `foldid <- caret::createFolds(y, k=5, list=FALSE)`
    *   **Replacement:**
        ```R
        # Inside MVPAModels$glmnet_opt$fit function:
        # ...
        # foldid_list <- create_mvpa_folds(y, k=5, seed=1234) # seed for reproducibility
        # foldid <- integer(length(y))
        # for (i in seq_along(foldid_list)) {
        #   foldid[foldid_list[[i]]] <- i
        # }
        # ...
        ```
    *   **Note:** `glmnet::cv.glmnet` expects `foldid` to be a vector where `foldid[i]` is the fold assignment for observation `i`. `create_mvpa_folds` returns a list of indices *per fold*. The conversion loop maps this back. Or, `cv.glmnet` can create its own folds if `foldid` is NULL. For simplicity and direct replacement, we can just pass `foldid = NULL` to `cv.glmnet` and let it handle its internal fold creation. This avoids needing `create_mvpa_folds` here.
        *   **Decision:** Simpler to let `cv.glmnet` create its own folds by default if `foldid` is not critical to replicate exactly. The `glmnet_opt` model uses `epsgo`, which in turn calls `tune.glmnet.interval`, which *then* calls `cv.glmnet`. The `foldid` is passed down.
            *   *Correction:* The `tune.glmnet.interval` function within `c060` (used by `glmnet_opt`) expects `foldid` to be pre-generated. So, the `create_mvpa_folds` and mapping is necessary if we want to keep that structure.
            *   **Final approach for `glmnet_opt`:**
                ```R
                # Inside MVPAModels$glmnet_opt$fit:
                # ...
                # fam <- if (length(lev) > 2) "multinomial" else "binomial"

                # Generate fold assignments for cv.glmnet
                fold_indices_list <- create_mvpa_folds(y, k = 5, seed = 1234) # Using our helper
                foldid_vector <- integer(length(y))
                for(f_idx in seq_along(fold_indices_list)) {
                    foldid_vector[fold_indices_list[[f_idx]]] <- f_idx
                }
                
                fit <- epsgo(Q.func="tune.glmnet.interval",
                             bounds=bounds,
                             parms.coding="none",
                             seed=1234,
                             show="none",
                             fminlower=-100,
                             x=x, y=y, family=fam,
                             type.min="lambda.1se",
                             foldid=foldid_vector, # Pass the generated fold assignments
                             type.measure="mse")
                # ...
                ```

4.  **Update Prediction Methods (R/model_fit.R):**
    *   **Context:** `predict.class_model_fit` and `predict.regression_model_fit` currently use `object$model$prob` and `object$model$predict`. `object$model` is the caret model *specification*. We need them to use `object$fit`, which is the actual *trained model object*.
    *   **Change in `predict.class_model_fit`:**
        ```R
        # Original:
        # probs <- object$model$prob(object$fit, mat)
        # New (conceptual - depends on how underlying model's predict method works for probs):
        # Assume object$fit is the direct output of e.g. sda::sda() or other model training.
        # We need a consistent way to get probabilities. If model is from MVPAModels,
        # the model spec itself provides the 'prob' function.
        # So, this logic might already be correct if MVPAModels are primary.
        # The key is that object$model refers to the MVPAModels spec, not caret's.
        
        # Reviewing MVPAModels entries:
        # MVPAModels$sda_notune$prob calls predict(modelFit, ..., verbose=FALSE)$posterior
        # MVPAModels$corclass$prob calls prob_corsimFit(modelFit, ...)
        # This seems fine if object$model refers to an MVPAModels entry.
        # The user's "surgery plan" step 1 already makes MVPAModels primary.
        # THEREFORE: This step might require no code change IF load_model is already updated
        # and MVPAModels structure is consistent.
        
        # Ensure `object$fit` is what the model's `prob` function expects.
        # Example: if modelFit is sda::sda output, then object$model$prob is a wrapper.
        # This means `object$model$prob(object$fit, ...)` is the correct pattern.
        # No change needed here if `load_model` is fixed first.
        ```
    *   **Change in `predict.regression_model_fit`:** Similar to above, ensure `object$model$predict(object$fit, ...)` is the correct pattern. This seems to be the case.

5.  **Remove `load_caret_libs()` (R/model_fit.R):**
    *   This function is now obsolete. Delete it.

6.  **Testing - Phase 1:**
    *   Run `R CMD check`.
    *   Run existing tests. Focus on models from `MVPAModels`, especially `glmnet_opt`.
    *   Ensure examples using `load_model` still work for `MVPAModels`.

---

**Phase 2: Resampling Control and Performance Metrics**

*   **Target:** `get_control()`, `mclass_summary()`, `compute_performance.mvpa_model()` and its helpers.
*   **Goal:** Replace `caret::trainControl` and `caret` summary functions with `yardstick`.

1.  **Rewrite `get_control()` (R/model_fit.R):**
    *   **`caret` call:** `caret::trainControl(...)`, `caret::twoClassSummary`, `mclass_summary`.
    *   **Replacement:**
        ```R
        # R/model_fit.R
        get_control <- function(y, nreps) { # nreps is for tune_model bootstrap, not used here
          is_class <- is.factor(y)
          
          metric_name <- if (is_class) {
            if (nlevels(y) == 2) "roc_auc" else "mn_log_loss" # Or "accuracy" / "kap"
          } else {
            "rmse" # Or "rsq"
          }
          
          # The 'control' object is now simpler, primarily carrying the metric name.
          # Cross-validation for tuning will be handled by the new tune_model.
          list(
            metric = metric_name,
            # No longer need classProbs, returnData, etc. from caret::trainControl
            # as these are handled differently now.
            # We can add other control parameters if our new tune_model needs them.
            number = nreps # Keep nreps if tune_model will use it for bootstrap.
          )
        }
        ```

2.  **Remove/Replace `mclass_summary()` (R/model_fit.R):**
    *   This was a custom summary function for `caret::trainControl`.
    *   It calculated a mean AUC across classes.
    *   **Replacement:** This specific calculation can be done with `yardstick` if needed, e.g., by computing AUC per class and averaging. `yardstick::roc_auc_vec` handles multiclass by default using handTill or one-vs-all. We can use `mn_log_loss` for multiclass tuning or overall accuracy.
    *   **Action:** Delete `mclass_summary`. Performance evaluation will be more direct.

3.  **Update Performance Calculation (R/performance.R and R/mvpa_model.R):**
    *   Modify `compute_performance.mvpa_model` in `R/mvpa_model.R`. It calls helpers like `get_multiclass_perf`.
    *   Modify `binary_perf()`, `multiclass_perf()`, `performance.regression_result()` in `R/performance.R`.
    *   **Key `yardstick` functions to use:**
        *   `binary_perf`:
            *   Accuracy: `yardstick::accuracy_vec(truth = observed, estimate = predicted)`
            *   AUC: `yardstick::roc_auc_vec(truth = observed, estimate = probs[,2], event_level = "second")` (assuming `probs[,2]` is prob of positive class, and positive class is second level).
        *   `multiclass_perf`:
            *   Accuracy: `yardstick::accuracy_vec(truth = observed, estimate = predicted)`
            *   AUC (average one-vs-all): `yardstick::roc_auc_vec(truth = observed, estimate = probs, estimator="macro")`
            *   Per-class AUC (if `class_metrics=TRUE`): Can loop `roc_auc_vec` with `event_level` set for each class or use a grouped summary.
        *   `performance.regression_result`:
            *   R2: `yardstick::rsq_vec(truth = observed, estimate = predicted)`
            *   RMSE: `yardstick::rmse_vec(truth = observed, estimate = predicted)`
            *   Spearman correlation needs `stats::cor`.
    *   **Example `binary_perf` rewrite (R/performance.R):**
        ```R
        # R/performance.R
        binary_perf <- function(observed, predicted, probs) {
          # Ensure observed is a factor with levels in a consistent order
          # and probs columns match this order.
          lvls <- levels(observed)
          if (length(lvls) != 2) stop("binary_perf expects 2 levels in observed.")
          
          # Assuming probs has columns named after levels(observed)
          # and we want AUC for the *second* level.
          # If not, this needs adjustment.
          prob_positive_class <- probs[, lvls[2]] 
          
          res_acc <- yardstick::accuracy_vec(truth = observed, estimate = factor(predicted, levels=lvls))
          res_auc <- tryCatch(
             yardstick::roc_auc_vec(truth = observed, estimate = prob_positive_class, event_level = "second"),
             error = function(e) NA_real_ # Handle cases where AUC cannot be computed
          )

          c(Accuracy = res_acc, AUC = res_auc) # yardstick returns numeric, ensure it's named
        }
        ```
    *   Update `get_multiclass_perf`, `get_binary_perf`, `get_regression_perf` in `R/mvpa_model.R` to correctly call the new performance functions.

4.  **Testing - Phase 2:**
    *   `R CMD check`.
    *   Verify performance metrics match (or are reasonably close to) previous `caret`-based calculations for a few key models. Differences might arise from default settings in `yardstick` vs. `caret`.

---

**Phase 3: Hyperparameter Tuning**

*   **Target:** `tune_model()` (R/model_fit.R).
*   **Goal:** Replace `caret::train()` with a custom tuning loop (Option A from user plan).

1.  **Rewrite `tune_model()` (R/model_fit.R):**
    *   **`caret` call:** `caret::train(...)`
    *   **Replacement (Custom Loop - Option A):**
        ```R
        # R/model_fit.R
        # @param mspec A model specification (e.g., from mvpa_model).
        # @param x Training data matrix (samples x features).
        # @param y Training response vector.
        # @param wts Optional class weights.
        # @param param_grid A data.frame where each row is a parameter combination.
        # @param nreps Number of bootstrap or CV folds for tuning.
        # @importFrom rsample bootstraps # or vfold_cv if using CV
        # @importFrom yardstick metric_set rmse roc_auc accuracy mn_log_loss rsq
        # @importFrom purrr map_dfr
        tune_model <- function(mspec, x, y, wts, param_grid, nreps = 10) {
          
          # Determine the metric from mspec$model$problemType or a new control object
          # This re-uses the logic from the updated get_control()
          temp_control <- get_control(y, nreps) # nreps passed but might not be used by get_control
          metric_to_optimize <- temp_control$metric 
          
          # Define the yardstick metric function based on metric_to_optimize
          # And define if higher is better for this metric
          higher_is_better <- TRUE
          metric_fn <- switch(metric_to_optimize,
            "roc_auc" = yardstick::roc_auc_vec,
            "mn_log_loss" = { higher_is_better <- FALSE; yardstick::mn_log_loss_vec },
            "accuracy" = yardstick::accuracy_vec,
            "rmse" = { higher_is_better <- FALSE; yardstick::rmse_vec },
            "rsq" = yardstick::rsq_vec,
            stop("Unsupported metric for tuning: ", metric_to_optimize)
          )
          
          # Create resamples (e.g., bootstrap or k-fold CV for tuning)
          # Using bootstraps as in the original caret::trainControl example
          # Ensure 'y' is a vector, not a matrix, for rsample stratification
          y_vector <- if(is.matrix(y) && ncol(y) == 1) y[,1] else y
          
          # Create a data frame for rsample
          # Ensure x is a data.frame for rsample functions
          df_for_rsample <- as.data.frame(x)
          df_for_rsample$.response_var_for_stratification <- y_vector # Add y for stratification
          
          resamples <- rsample::bootstraps(df_for_rsample, times = nreps, 
                                          strata = if(is.factor(y_vector)) .response_var_for_stratification else NULL)

          # Loop over each parameter combination in param_grid
          tuning_results <- purrr::map_dfr(seq_len(nrow(param_grid)), function(param_idx) {
            current_params <- param_grid[param_idx, , drop = FALSE]
            
            # Loop over each resample
            fold_metrics <- purrr::map_dbl(resamples$splits, function(split) {
              train_data_fold <- rsample::analysis(split)
              test_data_fold  <- rsample::assessment(split)
              
              # Separate predictors and response
              y_train_fold <- train_data_fold$.response_var_for_stratification
              x_train_fold <- train_data_fold[, !(names(train_data_fold) %in% ".response_var_for_stratification"), drop = FALSE]
              
              y_test_fold  <- test_data_fold$.response_var_for_stratification
              x_test_fold  <- test_data_fold[, !(names(test_data_fold) %in% ".response_var_for_stratification"), drop = FALSE]

              # Fit model on this training fold with current_params
              # fit_model expects x to be a matrix and y a vector/factor
              # mspec$model is the entry from MVPAModels
              # mspec already contains the 'model' element which is the list spec (e.g. MVPAModels$sda_notune)
              
              # The 'fit_model.mvpa_model' calls mspec$model$fit(...)
              # Ensure 'mspec' passed to fit_model.mvpa_model is the outer 'mspec'
              # or that the inner one has the right 'model' element.
              
              # fit_model.mvpa_model expects 'obj' to be the mvpa_model spec.
              # It then calls obj$model$fit.
              # 'param' argument for fit_model.mvpa_model needs to be current_params.
              
              # Simplification: fit_model.mvpa_model itself might not be the right one to call here,
              # as it's for the *final* fit. We need to call mspec$model$fit directly for tuning.
              
              # Assuming mspec$model has $fit, $predict, $prob
              # lev is levels(y), classProbs is typically TRUE for classification tuning
              
              # Ensure x_train_fold is matrix
              x_train_fold_mat <- as.matrix(x_train_fold)
              
              # Correct way to call the underlying model's fit function:
              # fit_object is the actual raw model (e.g. sda object, glmnet object)
              fit_object <- mspec$model$fit(x_train_fold_mat, y_train_fold, 
                                            wts = wts, # Pass weights if available
                                            param = current_params, # current_params for this iteration
                                            lev = levels(y_vector), 
                                            classProbs = is.factor(y_vector)) 

              # Predict on this test fold
              x_test_fold_mat <- as.matrix(x_test_fold)

              if (metric_to_optimize == "roc_auc" || metric_to_optimize == "mn_log_loss") {
                 # Need probabilities
                 preds_probs <- mspec$model$prob(fit_object, x_test_fold_mat)
                 # Ensure column names match levels of y_test_fold if factor
                 if (is.factor(y_test_fold) && is.matrix(preds_probs)) {
                    colnames(preds_probs) <- levels(y_test_fold)
                 }
                 # Calculate metric
                 val <- metric_fn(truth = y_test_fold, estimate = preds_probs, 
                                  event_level = if(metric_to_optimize=="roc_auc") "second" else NULL)
              } else {
                 # Need class predictions or numeric predictions
                 preds <- mspec$model$predict(fit_object, x_test_fold_mat)
                 if (is.factor(y_test_fold)) {
                    preds <- factor(preds, levels = levels(y_test_fold))
                 }
                 # Calculate metric
                 val <- metric_fn(truth = y_test_fold, estimate = preds)
              }
              return(val)
            }) # End loop over resamples
            
            # Aggregate metric across folds for this parameter set
            tibble::tibble(.param_id = param_idx, aggregated_metric = mean(fold_metrics, na.rm = TRUE))
          }) # End loop over param_grid
          
          # Select best parameters
          best_idx <- if (higher_is_better) {
            which.max(tuning_results$aggregated_metric)
          } else {
            which.min(tuning_results$aggregated_metric)
          }
          
          if (length(best_idx) == 0) { # All metrics were NA
            warning("tune_model: All parameter combinations resulted in NA performance. Returning first parameter set.")
            best_idx <- 1
          }
          
          best_tune <- param_grid[tuning_results$.param_id[best_idx], , drop = FALSE]
          return(best_tune)
        }
        ```
    *   **Dependencies:** Add `rsample` and `yardstick` to `Imports` in `DESCRIPTION`. `purrr` is already an import.
    *   **Considerations:**
        *   This loop is sequential. For parallelization, `furrr::future_map_dfr` can replace `purrr::map_dfr` for the outer loop (parameters) or inner loop (resamples).
        *   Error handling within the loops is basic. More robust error catching might be needed.
        *   This replaces `caret::train`'s sophisticated tuning logic with a simpler grid search. For many `MVPAModels`, this is sufficient as their grids are small.

2.  **Integrate `tune_model()` into `train_model.mvpa_model` (R/model_fit.R):**
    *   The existing `train_model.mvpa_model` calls `tune_grid()` to get `param`.
    *   Then, if `param` has >1 row, it should call the *new* `tune_model()` to get `best_param`.
    *   Then, it calls `fit_model()` with `best_param`.
    *   **Modification in `train_model.mvpa_model`:**
        ```R
        # R/model_fit.R (inside train_model.mvpa_model)
        # ... after feature selection ...
        # train_dat is now the feature-selected data matrix
        
        # Get initial parameter grid (could be single row or multiple)
        current_param_grid <- tune_grid(obj, train_dat, y, len = obj$tune_reps) # obj$tune_reps might be used by model$grid
        
        best_param <- if (!is.vector(current_param_grid) && !is.null(nrow(current_param_grid)) && nrow(current_param_grid) > 1) {
          # If multiple parameter sets, call the new tune_model
          bp <- tune_model(obj, train_dat, y, wts, current_param_grid, nreps = obj$tune_reps)
          futile.logger::flog.debug("Best tuning parameters: %s", 
                                   paste(capture.output(print(bp)), collapse="\n"))
          bp
        } else {
          # If only one parameter set, use it directly
          current_param_grid
        }
        
        # Fit final model with best_param
        fit <- fit_model(obj, train_dat, y, wts=wts, param=best_param, classProbs=TRUE)
        model_fit(obj$model, y, fit, mtype, best_param, indices, feature_mask)
        # ...
        ```

3.  **Testing - Phase 3:**
    *   `R CMD check`.
    *   Thoroughly test models that use `tune_grid` (e.g., `sda` with its default grid, `rf`).
    *   Compare results to ensure tuning selects reasonable parameters. This is the most complex phase to validate.

---

**Phase 4: Cleanup and Finalization**

1.  **Update `DESCRIPTION`:**
    *   Remove `caret` from `Imports`.
    *   Ensure `rsample`, `yardstick` are in `Imports`. (`purrr`, `dplyr`, `tibble` are already there).
2.  **Documentation:**
    *   Update any documentation, vignettes, or examples that refer to `caret` models not in `MVPAModels` or `caret` specific tuning options.
    *   Document the `MVPAModels` structure and the `register_mvpa_model()` helper.
3.  **Final `R CMD check --as-cran`:** Ensure everything passes.

---

**New Helper Functions to Implement:**

1.  **`register_mvpa_model(name, model_spec)` (e.g., in R/common.R or R/classifiers.R):**
    *   Allows users to add their own models to the `MVPAModels` environment.
    *   `model_spec` must be a list with `type`, `library`, `label`, `parameters`, `grid`, `fit`, `predict`, `prob` elements matching the `MVPAModels` convention.
    ```R
    #' Register a Custom MVPA Model
    #'
    #' Adds a user-defined model specification to the rMVPA model registry (MVPAModels).
    #'
    #' @param name A character string, the unique name for the model.
    #' @param model_spec A list containing the model specification. It must include
    #'   elements: `type` ("Classification" or "Regression"), `library` (character vector
    #'   of required packages), `label` (character, usually same as name), `parameters`
    #'   (data.frame of tunable parameters: name, class, label), `grid` (function to
    #'   generate tuning grid), `fit` (function), `predict` (function), and `prob`
    #'   (function, for classification).
    #' @export
    #' @examples
    #' \dontrun{
    #' my_model_spec <- list(
    #'   type = "Classification", library = "e1071", label = "my_svm",
    #'   parameters = data.frame(parameter = "cost", class = "numeric", label = "Cost"),
    #'   grid = function(x, y, len = NULL) data.frame(cost = 10^(-3:3)),
    #'   fit = function(x, y, wts, param, lev, ...) e1071::svm(x, y, cost = param$cost, ...),
    #'   predict = function(modelFit, newdata, ...) predict(modelFit, newdata),
    #'   prob = function(modelFit, newdata, ...) {
    #'     # For SVM with probability=TRUE, need to handle attr(predict(...), "probabilities")
    #'     pred_obj <- predict(modelFit, newdata, probability = TRUE)
    #'     attr(pred_obj, "probabilities")
    #'   }
    #' )
    #' register_mvpa_model("my_svm", my_model_spec)
    #' loaded_model <- load_model("my_svm")
    #' }
    register_mvpa_model <- function(name, model_spec) {
      # Basic validation of model_spec structure
      required_elements <- c("type", "library", "label", "parameters", "grid", "fit", "predict", "prob")
      if (!all(required_elements %in% names(model_spec))) {
        stop("model_spec is missing one or more required elements: ", 
             paste(setdiff(required_elements, names(model_spec)), collapse=", "))
      }
      if (!is.data.frame(model_spec$parameters) || 
          !all(c("parameter", "class", "label") %in% names(model_spec$parameters))) {
        stop("'model_spec$parameters' must be a data.frame with columns: parameter, class, label.")
      }
      # Assign to MVPAModels environment
      MVPAModels[[name]] <- model_spec
      invisible(NULL)
    }
    ```

2.  **`create_mvpa_folds()`:** Already detailed in Phase 1.

---

**Benefits Gained (Reiterated):**
*   **Reduced Dependencies:** Significant drop in installation footprint and complexity.
*   **Modern Tooling:** Leverages `rsample` and `yardstick` for robust and maintainable resampling and metrics.
*   **Explicit Model Handling:** Clearer distinction between internal `MVPAModels` and user-provided ones.
*   **Control:** More direct control over the tuning loop.
*   **Future Flexibility:** Easier to integrate other tuning backends (like `tune::tune_grid`) or advanced resampling strategies later without being tied to `caret`'s large ecosystem.

This formalized plan builds directly on your excellent analysis and provides a concrete path forward. The main effort will be in thoroughly implementing and testing the new `tune_model` function and ensuring metric calculations align.This is an excellent, detailed "surgery plan." It's well-structured, pragmatic, and correctly identifies the key `caret` dependencies and proposes sound replacements. The phased approach is spot on for managing complexity and ensuring the package remains runnable.

Here's a formalized path forward, building upon your plan, with some additional details, considerations, and helper function definitions.

## Formalized Plan: Excising `caret` from `rMVPA`

**Overall Strategy:**
We will adopt the phased approach outlined. The core idea is to replace `caret`'s model management, resampling, and tuning infrastructure with a combination of in-house solutions and lightweight, focused `tidymodels` packages (`rsample` for resampling, `yardstick` for metrics). This will significantly reduce dependencies and provide more direct control.

**Phase 0: Preparation**
1.  **Branching:** Create a new feature branch in Git (e.g., `feature/remove-caret`).
2.  **Dependency Management:**
    *   In `DESCRIPTION`:
        *   Add `rsample` and `yardstick` to `Imports:`.
        *   `caret` will remain in `Imports:` until Phase 4.
    *   Add `#' @importFrom rsample ...` and `#' @importFrom yardstick ...` in relevant R files as functions are used.

---

**Phase 1: Decouple Model Discovery and Basic Infrastructure**

*   **Goal:** Remove `caret` from model loading and the simplest resampling case used in `glmnet_opt`. Ensure core model fitting and prediction pathways are independent of `caret`'s *model specification* objects.

1.  **Modify `load_model()` (R/common.R):**
    *   **`caret` call:** `if (length(caret::getModelInfo(name)) > 0) …`
    *   **Action:** Remove this `else if` block.
    *   **New Logic (as per your plan):**
        ```R
        load_model <- function(name) {
          registry <- MVPAModels # Existing environment
          
          ret <- if (!is.null(registry[[name]])) {
            registry[[name]]   
          } else {
            stop(paste0("Model '", name, "' not found in MVPAModels. ",
                        "Register custom models using `register_mvpa_model()` if needed."))
          }
          
          ret$label <- name # Existing logic
          ret
        }
        ```

2.  **Implement `create_mvpa_folds()` Helper (New file: `R/resampling_utils.R` or similar):**
    *   **Purpose:** Replace `caret::createFolds`.
    *   **Implementation (using `rsample` for robustness and stratification):**
        ```R
        #' Create Cross-Validation Folds
        #'
        #' Generates a list of row indices for k-fold cross-validation.
        #' Can perform stratified sampling if y is a factor.
        #'
        #' @param y A vector, typically the response variable.
        #' @param k Integer, the number of folds.
        #' @param list Logical, if TRUE, return a list of indices for each fold.
        #'        If FALSE, return a vector of fold assignments for each observation.
        #'        (Mimicking caret's `list` argument).
        #' @param seed Optional integer for reproducible fold creation.
        #' @return If `list=TRUE`, a list of k integer vectors. If `list=FALSE`, an integer
        #'         vector of fold assignments.
        #' @importFrom rsample vfold_cv
        #' @keywords internal
        create_mvpa_folds <- function(y, k = 5, list = TRUE, seed = NULL) {
          if (!is.null(seed)) set.seed(seed)
          n <- length(y)
          
          # Create a dummy data frame for rsample
          df_for_rsample <- data.frame(.indices = seq_len(n))
          strata_arg <- NULL
          if (is.factor(y) && k < n) { # Stratification possible and meaningful
             df_for_rsample$.response_var_for_stratification <- y
             strata_arg <- ".response_var_for_stratification"
          }
          
          folds_obj <- rsample::vfold_cv(df_for_rsample, v = k, strata = strata_arg, repeats = 1)
          
          if (list) {
            # Extract assessment (hold-out) indices for each fold
            out_indices <- lapply(folds_obj$splits, function(split) rsample::assessment(split)$.indices)
          } else {
            # Create a vector of fold assignments
            out_indices <- integer(n)
            for (i in seq_along(folds_obj$splits)) {
              out_indices[rsample::assessment(folds_obj$splits[[i]])$.indices] <- i
            }
          }
          return(out_indices)
        }
        ```

3.  **Patch `glmnet_opt` Model (R/classifiers.R):**
    *   **`caret` call:** `foldid <- caret::createFolds(y,k=5,list=FALSE)`
    *   **Replacement (inside `MVPAModels$glmnet_opt$fit`):**
        ```R
        # ...
        # fam <- if (length(lev) > 2) "multinomial" else "binomial" # Already exists

        # Generate fold assignments for cv.glmnet
        foldid_vector <- create_mvpa_folds(y, k = 5, list = FALSE, seed = 1234) 
        
        fit <- epsgo(Q.func="tune.glmnet.interval",
                     # ... other epsgo args ...
                     foldid = foldid_vector, # Pass the generated fold assignments
                     # ...
                     )
        # ...
        ```

4.  **Verify Prediction Methods (R/model_fit.R):**
    *   `predict.class_model_fit` and `predict.regression_model_fit` use `object$model$prob` and `object$model$predict` respectively.
    *   **Action:** Confirm that after Step 1 (modifying `load_model`), `object$model` correctly refers to an entry from `MVPAModels`. The structure of `MVPAModels` entries (which include `fit`, `predict`, `prob` functions) means these prediction methods should already be calling the correct internal functions defined within `MVPAModels` rather than relying on `caret` model objects. No code change should be needed here *if* `MVPAModels` specs are self-contained.

5.  **Remove `load_caret_libs()` (R/model_fit.R):**
    *   Delete this function.

6.  **Testing - Phase 1:**
    *   Run `R CMD check`.
    *   Execute existing tests. Pay close attention to tests involving `load_model` (ensure it only loads from `MVPAModels`) and any tests directly or indirectly using `glmnet_opt`.

---

**Phase 2: Resampling Control and Performance Metrics**

*   **Goal:** Replace `caret::trainControl` and `caret`'s summary functions with simpler internal logic and `yardstick` for metrics.

1.  **Rewrite `get_control()` (R/model_fit.R):**
    *   **`caret` call:** `caret::trainControl(...)`, `caret::twoClassSummary`, `mclass_summary`.
    *   **Replacement (as per your plan):**
        ```R
        # R/model_fit.R
        get_control <- function(y, nreps) { # nreps for bootstrap tuning
          is_class <- is.factor(y)
          
          # Determine primary metric for tuning and optimization
          metric_name <- if (is_class) {
            if (nlevels(y) == 2) "roc_auc" else "accuracy" # Defaulting to accuracy for multiclass
          } else {
            "rmse" 
          }
          
          list(
            metric = metric_name, 
            number = nreps # For bootstrap reps in tuning, if applicable
          )
        }
        ```

2.  **Remove `mclass_summary()` (R/model_fit.R):**
    *   Delete this function. Its functionality (mean AUC for multiclass) will be replaced by direct `yardstick` calls or simpler aggregation if needed.

3.  **Update Performance Calculation (R/performance.R and R/mvpa_model.R):**
    *   **File: R/performance.R**
        *   **Modify `binary_perf()`:**
            ```R
            # @importFrom yardstick accuracy roc_auc_vec
            binary_perf <- function(observed, predicted, probs) {
              lvls <- levels(observed)
              if (length(lvls) != 2) stop("binary_perf expects 2 levels in observed.")
              
              # Ensure predicted is a factor with the same levels as observed
              predicted_factor <- factor(predicted, levels = lvls)
              
              # Assuming probs has columns named after levels(observed) or in the same order.
              # And positive class is the second level.
              prob_positive_class <- if (ncol(probs) == 2) probs[, lvls[2]] else probs[,1] # Adapt if probs is single col

              res_acc <- yardstick::accuracy_vec(truth = observed, estimate = predicted_factor)
              res_auc <- tryCatch(
                 yardstick::roc_auc_vec(truth = observed, estimate = prob_positive_class, event_level = "second"),
                 error = function(e) NA_real_
              )
              c(Accuracy = res_acc, AUC = res_auc)
            }
            ```
        *   **Modify `multiclass_perf()`:**
            ```R
            # @importFrom yardstick accuracy roc_auc_vec
            multiclass_perf <- function(observed, predicted, probs, class_metrics=FALSE) {
              lvls <- levels(observed)
              predicted_factor <- factor(predicted, levels = lvls)

              acc <- yardstick::accuracy_vec(truth = observed, estimate = predicted_factor)
              
              # Default multiclass AUC (Hand-Till generalization, macro-averaged)
              # Ensure probs colnames match levels(observed)
              colnames(probs) <- lvls 
              avg_auc <- tryCatch(
                yardstick::roc_auc_vec(truth = observed, estimate = probs, estimator = "macro"), 
                error = function(e) NA_real_
              )
              
              metrics <- c(Accuracy = acc, AUC = avg_auc)
              
              if (class_metrics) {
                # Per-class AUC (one-vs-all)
                auc_per_class <- sapply(lvls, function(lvl) {
                  tryCatch(
                    yardstick::roc_auc_vec(truth = observed, estimate = probs, event_level = "second", case_weights = NULL,
                                          # One-vs-all by focusing on prob of current level vs others
                                          .roc_auc_ovr_impl_custom(truth = observed, estimate = probs[,lvl], event_level=lvl)),
                    error = function(e) NA_real_
                  )
                })
                names(auc_per_class) <- paste0("AUC_", lvls)
                metrics <- c(metrics, auc_per_class)
              }
              metrics
            }
            # Helper for per-class AUC if yardstick doesn't directly support it easily for multiclass one-vs-all from matrix
            .roc_auc_ovr_impl_custom <- function(truth, estimate, event_level) {
                 positive_class_prob <- estimate # This estimate column IS the prob of event_level
                 binary_truth <- factor(ifelse(truth == event_level, event_level, "other_"), levels=c("other_", event_level))
                 yardstick::roc_auc_vec(truth = binary_truth, estimate = positive_class_prob, event_level = "second")
            }
            ```
        *   **Modify `performance.regression_result()`:**
            ```R
            # @importFrom yardstick rsq_vec rmse_vec
            # @importFrom stats cor
            performance.regression_result <- function(x, split_list,...) {
              # ... (split_list logic remains) ...
              obs <- x$observed
              pred <- x$predicted
              
              res_rsq <- yardstick::rsq_vec(truth = obs, estimate = pred)
              res_rmse <- yardstick::rmse_vec(truth = obs, estimate = pred)
              res_spearcor <- tryCatch(stats::cor(obs, pred, method="spearman", use="pairwise.complete.obs"), error = function(e) NA_real_)
              
              c(R2=res_rsq, RMSE=res_rmse, spearcor=res_spearcor)
            }
            ```
    *   **File: R/mvpa_model.R:**
        *   The `get_multiclass_perf`, `get_binary_perf`, `get_regression_perf` wrappers around the `performance.*` S3 methods seem fine as they mainly handle `split_list`. Ensure they correctly call the updated S3 methods.

4.  **Testing - Phase 2:**
    *   `R CMD check`.
    *   Verify performance metrics. This is critical. Create small, deterministic test cases where expected Accuracy/AUC/RMSE/R2 can be manually calculated or are known, and compare `rMVPA`'s output.

---

**Phase 3: Hyperparameter Tuning**

*   **Goal:** Replace `caret::train()` with a custom tuning loop.

1.  **Rewrite `tune_model()` (R/model_fit.R):**
    *   **`caret` call:** `caret::train(...)`
    *   **Replacement (Custom Loop - Option A):**
        ```R
        # R/model_fit.R
        # @param mspec The full mvpa_model specification.
        # @param x Training data matrix (samples x features). This is *after* any
        #        feature selection done in train_model.mvpa_model.
        # @param y Training response vector.
        # @param wts Optional class weights.
        # @param param_grid A data.frame where each row is a parameter combination.
        # @param nreps Number of resamples (e.g., bootstrap iterations or CV folds) for tuning.
        # @importFrom rsample bootstraps # or vfold_cv
        # @importFrom yardstick metric_set rmse roc_auc accuracy mn_log_loss rsq
        # @importFrom purrr map_dfr
        # @importFrom stats predict
        tune_model <- function(mspec, x, y, wts, param_grid, nreps = 10) {
          
          # Get metric to optimize from the mspec's control object (derived via get_control)
          # The 'model' element within mspec is the list from MVPAModels
          control_obj <- get_control(y, nreps) 
          metric_to_optimize <- control_obj$metric 
          
          higher_is_better <- TRUE # Default for AUC, Accuracy, Rsq
          metric_fn_yardstick <- switch(metric_to_optimize,
            "roc_auc"     = yardstick::roc_auc_vec,
            "mn_log_loss" = { higher_is_better <- FALSE; yardstick::mn_log_loss_vec },
            "accuracy"    = yardstick::accuracy_vec,
            "rmse"        = { higher_is_better <- FALSE; yardstick::rmse_vec },
            "rsq"         = yardstick::rsq_vec,
            stop("Unsupported metric for tuning: ", metric_to_optimize)
          )
          
          y_vector <- if(is.matrix(y) && ncol(y) == 1) y[,1] else y
          
          df_for_rsample <- as.data.frame(x)
          df_for_rsample$.response_var_for_stratification <- y_vector
          
          # Using bootstraps for tuning, consistent with original caret::trainControl default for this package
          resamples_obj <- rsample::bootstraps(df_for_rsample, times = nreps, 
                                            strata = if(is.factor(y_vector)) .response_var_for_stratification else NULL)

          tuning_metrics <- purrr::map_dfr(seq_len(nrow(param_grid)), .f = function(param_idx) {
            current_params_df <- param_grid[param_idx, , drop = FALSE]
            
            # Performance over resamples for this parameter set
            resample_perf <- purrr::map_dbl(resamples_obj$splits, .f = function(split) {
              train_df_fold <- rsample::analysis(split)
              test_df_fold  <- rsample::assessment(split)
              
              y_train_fold <- train_df_fold$.response_var_for_stratification
              x_train_fold <- as.matrix(train_df_fold[, !(names(train_df_fold) %in% ".response_var_for_stratification"), drop = FALSE])
              
              y_test_fold  <- test_df_fold$.response_var_for_stratification
              x_test_fold  <- as.matrix(test_df_fold[, !(names(test_df_fold) %in% ".response_var_for_stratification"), drop = FALSE])

              # mspec$model is the list from MVPAModels (e.g., MVPAModels$sda_notune)
              # Call its $fit element
              # The `param` argument to the model's fit function should be the current set
              fit_obj <- mspec$model$fit(x_train_fold, y_train_fold, 
                                         wts = wts, # Pass weights if available
                                         param = current_params_df, 
                                         lev = levels(y_vector), 
                                         classProbs = is.factor(y_vector)) # Pass classProbs

              # Predict and evaluate
              metric_value <- NA_real_
              if (metric_to_optimize == "roc_auc" || metric_to_optimize == "mn_log_loss") {
                 preds_probs <- mspec$model$prob(fit_obj, x_test_fold)
                 if (is.factor(y_test_fold) && is.matrix(preds_probs) && !is.null(levels(y_test_fold))) {
                    colnames(preds_probs) <- levels(y_test_fold)
                 }
                 metric_value <- metric_fn_yardstick(truth = y_test_fold, estimate = preds_probs, 
                                                     event_level = if(metric_to_optimize=="roc_auc" && nlevels(y_test_fold) == 2) "second" else NULL)
              } else {
                 preds <- mspec$model$predict(fit_obj, x_test_fold)
                 if (is.factor(y_test_fold) && !is.null(levels(y_test_fold))) {
                    preds <- factor(preds, levels = levels(y_test_fold))
                 }
                 metric_value <- metric_fn_yardstick(truth = y_test_fold, estimate = preds)
              }
              metric_value
            }) # End loop over resamples_obj$splits
            
            # Return mean performance for this parameter set
            data.frame(.param_id = param_idx, mean_metric = mean(resample_perf, na.rm = TRUE))
          }) # End loop over param_grid
          
          best_row_idx <- if (higher_is_better) {
            which.max(tuning_metrics$mean_metric)
          } else {
            which.min(tuning_metrics$mean_metric)
          }
          
          if (length(best_row_idx) == 0 || is.na(tuning_metrics$mean_metric[best_row_idx])) { 
            warning("tune_model: All parameter combinations resulted in NA performance. Returning first parameter set.")
            best_row_idx <- 1
          }
          
          best_param_id <- tuning_metrics$.param_id[best_row_idx]
          best_tune_params <- param_grid[best_param_id, , drop = FALSE]
          
          return(best_tune_params)
        }
        ```

2.  **Integrate `tune_model()` into `train_model.mvpa_model` (R/model_fit.R):**
    *   The logic outlined in your "surgery plan" (Phase 3, Step 2) is correct.
        ```R
        # R/model_fit.R (inside train_model.mvpa_model)
        # ... after feature selection ...
        # train_dat is feature-selected data matrix

        # Get initial parameter grid. obj is mspec. tune_grid.mvpa_model is called.
        current_param_grid <- tune_grid(obj, train_dat, y, len = obj$tune_reps) 
        
        best_param <- if (!is.vector(current_param_grid) && !is.null(nrow(current_param_grid)) && nrow(current_param_grid) > 1) {
          bp <- tune_model(obj, train_dat, y, wts, current_param_grid, nreps = obj$tune_reps) # obj is mspec
          futile.logger::flog.debug("Best tuning parameters: %s", 
                                   paste(capture.output(print(bp)), collapse="\n"))
          bp
        } else {
          current_param_grid
        }
        
        # Fit final model with best_param
        # obj refers to the mspec. fit_model.mvpa_model calls obj$model$fit
        fit <- fit_model(obj, train_dat, y, wts=wts, param=best_param, classProbs = is.factor(y))
        model_fit(obj$model, y, fit, mtype, best_param, indices, feature_mask)
        # ...
        ```

3.  **Testing - Phase 3:**
    *   `R CMD check`.
    *   Crucial: Test models with tuning grids (e.g., `sda` with its default grid, `rf` if used, `corclass` with `robust=TRUE/FALSE`). Compare parameter selection if possible or at least ensure it runs and selects *a* parameter set.

---

**Phase 4: Cleanup and Finalization**

1.  **Update `DESCRIPTION`:**
    *   Remove `caret` and `c060` (if `glmnet_opt`'s `epsgo` is replaced or `c060` is re-evaluated for direct use) from `Imports`.
    *   Ensure `rsample`, `yardstick`, `rlang` (already there), `purrr` (already there) are in `Imports`. If `glmnet_opt` still uses `epsgo` from `c060`, then `c060` needs to stay. (The plan above *keeps* `epsgo` for `glmnet_opt` but replaces its internal `createFolds`).
2.  **Documentation:**
    *   Update vignettes, examples, and function documentation that might refer to `caret` or imply `caret` models can be used directly.
    *   Clearly document the `MVPAModels` structure and `register_mvpa_model()`.
3.  **Final `R CMD check --as-cran`.**

---

**Helper Function: `register_mvpa_model` (R/common.R or R/classifiers.R)**
(As detailed in your Phase 3, copied here for completeness)
```R
    #' Register a Custom MVPA Model
    #'
    #' Adds a user-defined model specification to the rMVPA model registry (MVPAModels).
    #'
    #' @param name A character string, the unique name for the model.
    #' @param model_spec A list containing the model specification. It must include
    #'   elements: `type` ("Classification" or "Regression"), `library` (character vector
    #'   of required packages for the *model itself*, not for rMVPA's wrappers), 
    #'   `label` (character, usually same as name), `parameters`
    #'   (data.frame of tunable parameters: parameter, class, label), `grid` (function to
    #'   generate tuning grid, takes x, y, len args), `fit` (function), `predict` (function), 
    #'   and `prob` (function for classification, takes modelFit, newdata; should return matrix/df with colnames as levels).
    #' @export
    #' @examples
    #' \dontrun{
    #' # Example of how a user might define an e1071 SVM spec
    #' my_svm_spec <- list(
    #'   type = "Classification", library = "e1071", label = "my_svm",
    #'   parameters = data.frame(parameter = "cost", class = "numeric", label = "Cost (C)"),
    #'   # grid should return a data.frame with columns matching 'parameter' names in 'parameters'
    #'   grid = function(x, y, len = NULL) { 
    #'      data.frame(cost = if (is.null(len) || len == 1) 1 else 10^seq(-2, 2, length.out = len))
    #'   },
    #'   # fit function receives: x, y, wts (weights), param (current params from grid), 
    #'   # lev (levels of y), last (unused), weights (unused), classProbs (unused by e1071::svm)
    #'   fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    #'      e1071::svm(x, y, cost = param$cost, probability = TRUE, ...) # Ensure probability=TRUE for prob
    #'   },
    #'   # predict function receives: modelFit (output of $fit), newdata
    #'   predict = function(modelFit, newdata, ...) {
    #'      predict(modelFit, newdata, ...)
    #'   },
    #'   # prob function receives: modelFit, newdata
    #'   # Should return a matrix/df with columns named as in levels(y)
    #'   prob = function(modelFit, newdata, ...) {
    #'     pred_obj <- predict(modelFit, newdata, probability = TRUE)
    #'     attr(pred_obj, "probabilities") 
    #'   }
    #' )
    #' register_mvpa_model("my_svm", my_svm_spec)
    #' # Now load_model("my_svm") would work.
    #' }
    register_mvpa_model <- function(name, model_spec) {
      required_elements <- c("type", "library", "label", "parameters", "grid", "fit", "predict", "prob")
      if (!all(required_elements %in% names(model_spec))) {
        stop("model_spec is missing one or more required elements: ", 
             paste(setdiff(required_elements, names(model_spec)), collapse=", "))
      }
      if (!is.data.frame(model_spec$parameters) || 
          !all(c("parameter", "class", "label") %in% names(model_spec$parameters))) {
        stop("'model_spec$parameters' must be a data.frame with columns: parameter, class, label.")
      }
      MVPAModels[[name]] <- model_spec
      invisible(NULL)
    }
```
Okay, here's a ticketed list of items for the `caret` refactoring, based on the formalized plan. Each ticket includes the targeted file(s), a description of the task, and acceptance criteria (Definition of Done - DoD).

---
**Project: Refactor rMVPA - Remove `caret` Dependency**

**Preamble:**
*   All work should be done on a dedicated feature branch (e.g., `feature/remove-caret`).
*   After each phase (or significant ticket), run `R CMD check` and relevant existing tests to ensure the package remains in a runnable state.
*   Incrementally update `NAMESPACE` with `@importFrom` directives as new functions from `rsample` and `yardstick` are used.

---

**Phase 0: Preparation**

1.  **[ ] Ticket #001: Setup - Branching and Initial Dependency Declaration**
    *   **File(s):** `DESCRIPTION`
    *   **Task:**
        1.  Create a new Git feature branch (e.g., `feature/remove-caret`).
        2.  Add `rsample` and `yardstick` to the `Suggests:` field in the `DESCRIPTION` file. (They will be moved to `Imports:` later as their functions are directly used).
    *   **DoD:**
        *   New Git branch created and checked out.
        *   `DESCRIPTION` file updated with `rsample` and `yardstick` in `Suggests:`.
        *   Package installs and loads without error.

---

**Phase 1: Decouple Model Discovery and Basic Infrastructure**

2.  **[ ] Ticket #002: Modify `load_model()` - Remove `caret` Fallback**
    *   **File(s):** `R/common.R`
    *   **Task:** Remove the `else if (length(caret::getModelInfo(name)) > 0) …` block from `load_model()`. Update the error message for unknown models.
    *   **DoD:**
        *   `load_model()` no longer calls `caret::getModelInfo()`.
        *   `load_model()` successfully loads models defined in `MVPAModels`.
        *   `load_model()` throws an error for models not in `MVPAModels`.
        *   Relevant unit tests for `load_model` pass.

3.  **[ ] Ticket #003: Implement `create_mvpa_folds()` Helper**
    *   **File(s):** New file (e.g., `R/resampling_utils.R`)
    *   **Task:** Implement the `create_mvpa_folds(y, k, list, seed)` function using `rsample::vfold_cv` for stratified k-fold CV (if `y` is factor) or simple random k-fold CV.
    *   **DoD:**
        *   `create_mvpa_folds()` function exists and is exported (or internally available if preferred).
        *   Unit tests for `create_mvpa_folds()` cover:
            *   Correct number of folds returned.
            *   `list=TRUE` returns a list of index vectors.
            *   `list=FALSE` returns a vector of fold assignments.
            *   Stratification works correctly for factor `y`.
            *   `seed` argument ensures reproducibility.
            *   Handles edge cases (e.g., `k` > `n`, `k` = `n`).

4.  **[ ] Ticket #004: Patch `glmnet_opt` Model for Fold Creation**
    *   **File(s):** `R/classifiers.R` (within `MVPAModels$glmnet_opt$fit`)
    *   **Task:** Replace the `caret::createFolds()` call with a call to the new `create_mvpa_folds(y, k=5, list=FALSE, seed=1234)` to generate the `foldid_vector`.
    *   **DoD:**
        *   `MVPAModels$glmnet_opt$fit` no longer calls `caret::createFolds`.
        *   `glmnet_opt` model trains successfully using folds generated by `create_mvpa_folds`.
        *   Existing tests for `glmnet_opt` (if any specifically test fold generation) pass.

5.  **[ ] Ticket #005: Verify Prediction Method Calls**
    *   **File(s):** `R/model_fit.R` (`predict.class_model_fit`, `predict.regression_model_fit`)
    *   **Task:** Review these prediction methods. Confirm that `object$model$prob` and `object$model$predict` correctly refer to functions within the `MVPAModels` specifications after `load_model` changes, not `caret` model objects.
    *   **DoD:**
        *   Code review confirms correct dispatch to `MVPAModels` functions.
        *   No code changes are required *if* the assumption holds.
        *   Existing prediction tests continue to pass.

6.  **[ ] Ticket #006: Remove `load_caret_libs()` Function**
    *   **File(s):** `R/model_fit.R`
    *   **Task:** Delete the `load_caret_libs()` function.
    *   **DoD:**
        *   Function `load_caret_libs()` is removed.
        *   Package compiles and loads without error.

7.  **[ ] Ticket #007: Phase 1 Integration Testing**
    *   **File(s):** N/A (Testing activity)
    *   **Task:** Run `R CMD check`. Execute all existing unit tests. Pay special attention to model loading, `glmnet_opt` training, and general prediction pathways.
    *   **DoD:**
        *   `R CMD check` passes with no new errors/warnings related to these changes.
        *   All existing relevant unit tests pass.

---

**Phase 2: Resampling Control and Performance Metrics**

8.  **[ ] Ticket #008: Rewrite `get_control()`**
    *   **File(s):** `R/model_fit.R`
    *   **Task:** Rewrite `get_control(y, nreps)` to return a simple list containing the `metric` name (e.g., "roc_auc", "accuracy", "rmse") and `number` (for `nreps`), removing `caret::trainControl` logic.
    *   **DoD:**
        *   `get_control()` no longer calls `caret::trainControl`.
        *   Function returns the expected list structure based on `y`.

9.  **[ ] Ticket #009: Remove `mclass_summary()` Function**
    *   **File(s):** `R/model_fit.R`
    *   **Task:** Delete the `mclass_summary()` function.
    *   **DoD:**
        *   Function `mclass_summary()` is removed.
        *   Package compiles and loads.

10. **[ ] Ticket #010: Update `binary_perf()` using `yardstick`**
    *   **File(s):** `R/performance.R`
    *   **Task:** Rewrite `binary_perf(observed, predicted, probs)` to use `yardstick::accuracy_vec` and `yardstick::roc_auc_vec`.
    *   **DoD:**
        *   `binary_perf()` uses `yardstick` functions.
        *   Returns named vector `c(Accuracy = ..., AUC = ...)`.
        *   Metrics are validated against known examples or previous `caret` outputs.

11. **[ ] Ticket #011: Update `multiclass_perf()` using `yardstick`**
    *   **File(s):** `R/performance.R`
    *   **Task:** Rewrite `multiclass_perf(observed, predicted, probs, class_metrics)` to use `yardstick::accuracy_vec` and `yardstick::roc_auc_vec` (with `estimator="macro"` for average AUC). Implement per-class AUC if `class_metrics=TRUE`.
    *   **DoD:**
        *   `multiclass_perf()` uses `yardstick` functions.
        *   Returns named vector including `Accuracy`, `AUC`, and per-class AUCs if requested.
        *   Metrics validated.

12. **[ ] Ticket #012: Update `performance.regression_result()` using `yardstick`**
    *   **File(s):** `R/performance.R`
    *   **Task:** Rewrite `performance.regression_result(x, split_list, ...)` to use `yardstick::rsq_vec`, `yardstick::rmse_vec`, and `stats::cor` for Spearman correlation.
    *   **DoD:**
        *   `performance.regression_result()` uses `yardstick` and `stats::cor`.
        *   Returns named vector `c(R2=..., RMSE=..., spearcor=...)`.
        *   Metrics validated.

13. **[ ] Ticket #013: Update Performance Wrappers in `mvpa_model.R`**
    *   **File(s):** `R/mvpa_model.R`
    *   **Task:** Ensure `get_multiclass_perf()`, `get_binary_perf()`, `get_regression_perf()` correctly call the updated S3 performance methods.
    *   **DoD:**
        *   Wrapper functions correctly dispatch to the new `yardstick`-based performance calculators.
        *   `compute_performance.mvpa_model` functions as expected.

14. **[ ] Ticket #014: Phase 2 Integration Testing**
    *   **File(s):** N/A (Testing activity)
    *   **Task:** Run `R CMD check`. Execute all existing unit tests, focusing on those that calculate and report performance metrics. Compare outputs with previous versions if possible.
    *   **DoD:**
        *   `R CMD check` passes.
        *   Performance metrics are correctly calculated and reported.

---

**Phase 3: Hyperparameter Tuning**

15. **[ ] Ticket #015: Rewrite `tune_model()` - Custom Tuning Loop**
    *   **File(s):** `R/model_fit.R`
    *   **Task:** Implement the new `tune_model(mspec, x, y, wts, param_grid, nreps)` function. This will involve:
        1.  Determining the optimization metric and direction (higher better/lower better) from `get_control()`.
        2.  Setting up resamples (e.g., `rsample::bootstraps`).
        3.  Looping through `param_grid`.
        4.  For each parameter set, looping through resamples:
            *   Fitting the model (`mspec$model$fit`) on the analysis set.
            *   Predicting (`mspec$model$predict` or `mspec$model$prob`) on the assessment set.
            *   Calculating the chosen metric using `yardstick`.
        5.  Averaging the metric across resamples for each parameter set.
        6.  Selecting the `param_grid` row that yields the best average metric.
    *   **DoD:**
        *   `tune_model()` function implemented, no longer calls `caret::train`.
        *   Uses `rsample` for resampling and `yardstick` for metric calculation.
        *   Correctly identifies and returns the best parameter set from `param_grid`.
        *   Unit tests for `tune_model` with simple mock models and grids verify its logic.

16. **[ ] Ticket #016: Integrate new `tune_model()` into `train_model.mvpa_model`**
    *   **File(s):** `R/model_fit.R` (`train_model.mvpa_model`)
    *   **Task:** Modify `train_model.mvpa_model` so that if the `current_param_grid` (obtained from `tune_grid(obj, ...)`) has more than one row, it calls the new `tune_model()` to determine `best_param`. Otherwise, `current_param_grid` is used as `best_param`.
    *   **DoD:**
        *   `train_model.mvpa_model` correctly calls the new `tune_model` when appropriate.
        *   The final model is fitted using the `best_param` determined by `tune_model` or the single provided parameter set.
        *   Existing tests for models involving tuning (e.g., `sda` with its default grid, `corclass` with multiple options) pass and select reasonable parameters.

17. **[ ] Ticket #017: Phase 3 Integration Testing**
    *   **File(s):** N/A (Testing activity)
    *   **Task:** Run `R CMD check`. Thoroughly test models with tuning grids. Verify that the tuning process completes and reasonable parameters are selected. This is a critical validation step.
    *   **DoD:**
        *   `R CMD check` passes.
        *   Models with hyperparameter tuning (`tune_grid` in `mvpa_model` call results in multiple rows) train successfully.

---

**Phase 4: Cleanup and Finalization**

18. **[ ] Ticket #018: Update `DESCRIPTION` File - Final Dependencies**
    *   **File(s):** `DESCRIPTION`
    *   **Task:**
        1.  Remove `caret` from `Imports:`.
        2.  If `glmnet_opt` still uses `epsgo` from `c060`, ensure `c060` is in `Imports:`. Otherwise, remove it if `epsgo` usage is also refactored out (outside current scope).
        3.  Ensure `rsample` and `yardstick` are in `Imports:`.
    *   **DoD:**
        *   `DESCRIPTION` file accurately reflects the new dependencies.
        *   Package installs correctly with the updated dependencies.

19. **[ ] Ticket #019: Implement `register_mvpa_model()` Helper**
    *   **File(s):** `R/common.R` or `R/classifiers.R`
    *   **Task:** Implement the `register_mvpa_model(name, model_spec)` function to allow users to add models to `MVPAModels`. Include basic validation for the `model_spec` structure.
    *   **DoD:**
        *   `register_mvpa_model()` function implemented and exported.
        *   A simple test case demonstrates adding a custom model specification and then loading it with `load_model()`.

20. **[ ] Ticket #020: Documentation Update**
    *   **File(s):** All Rd files, vignettes, README.
    *   **Task:**
        1.  Remove any mentions of `caret` models being loadable via `load_model` if they are not part of `MVPAModels`.
        2.  Update documentation for functions whose behavior or dependencies changed (e.g., `mvpa_model` regarding tuning if internals changed significantly).
        3.  Document the `MVPAModels` structure and the new `register_mvpa_model()` function.
        4.  Ensure examples are updated and do not rely on implicit `caret` behavior.
    *   **DoD:**
        *   All user-facing documentation is accurate and reflects the removal of `caret`.
        *   Examples are functional with the refactored code.

21. **[ ] Ticket #021: Final Comprehensive Check and Merge Preparation**
    *   **File(s):** Entire codebase.
    *   **Task:**
        1.  Perform a final `R CMD check --as-cran`.
        2.  Review all changes for any lingering `caret` references or unintended consequences.
        3.  Ensure all unit tests pass.
        4.  Consider building vignettes to check for issues.
    *   **DoD:**
        *   `R CMD check --as-cran` is clean.
        *   All tests pass.
        *   Code review complete.
        *   Branch is ready to be merged into the main development line.

---

This ticketed list provides a systematic way to track progress and ensure all aspects of the refactoring are addressed.