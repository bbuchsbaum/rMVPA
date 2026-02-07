pkgname <- "rMVPA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rMVPA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MVPAModels")
### * MVPAModels

flush(stderr()); flush(stdout())

### Name: MVPAModels
### Title: Pre-defined MVPA Classification Models
### Aliases: MVPAModels
### Keywords: datasets

### ** Examples

# Load simple SDA classifier
model <- load_model("sda_notune")

# Load correlation classifier
model <- load_model("corclass")




cleanEx()
nameEx("add_interaction_contrasts")
### * add_interaction_contrasts

flush(stderr()); flush(stdout())

### Name: add_interaction_contrasts
### Title: Add Interaction Contrasts to an msreve_design
### Aliases: add_interaction_contrasts

### ** Examples

## Not run: 
##D # Example with non-overlapping contrasts (zero interaction)
##D C1 <- matrix(c(1,-1,0,0, 0,0,1,-1), nrow=4, 
##D              dimnames=list(NULL, c("A","B")))
##D # A compares conditions 1 vs 2, B compares 3 vs 4
##D # Their interaction will be zero and skipped
##D 
##D # Example with overlapping contrasts (non-zero interaction)  
##D C2 <- matrix(c(1,1,-1,-1, 1,-1,1,-1), nrow=4,
##D              dimnames=list(NULL, c("Main1","Main2")))
##D # These contrasts overlap and will produce a meaningful interaction
## End(Not run)



cleanEx()
nameEx("average_labels")
### * average_labels

flush(stderr()); flush(stdout())

### Name: average_labels
### Title: Average NeuroVec Data by Labels
### Aliases: average_labels

### ** Examples

## Not run: 
##D # Basic averaging
##D averaged <- average_labels(scandat, condition_labels, mask)
##D 
##D # With z-score normalization of each volume
##D averaged_norm <- average_labels(scandat, condition_labels, mask, 
##D                                 normalize = "z", normalize_by = "volume")
##D 
##D # Scale to unit norm for RSA
##D averaged_unit <- average_labels(scandat, condition_labels, mask,
##D                                 normalize = "unit")
##D                                 
##D # Get just the data matrix
##D data_mat <- average_labels(scandat, condition_labels, mask,
##D                            return_matrix = TRUE)
## End(Not run)




cleanEx()
nameEx("balance_partitions")
### * balance_partitions

flush(stderr()); flush(stdout())

### Name: balance_partitions
### Title: Balance Cross-Validation Partitions
### Aliases: balance_partitions balance_partitions.default
###   balance_partitions.blocked_cross_validation
###   balance_partitions.kfold_cross_validation
###   balance_partitions.twofold_blocked_cross_validation
###   balance_partitions.bootstrap_blocked_cross_validation
###   balance_partitions.sequential_blocked_cross_validation
###   balance_partitions.custom_cross_validation

### ** Examples

# Create an imbalanced dataset design (more class 'b')
design_df <- data.frame(condition = factor(rep(c("a", "b", "b"), 20)),
                       block = rep(1:6, each = 10))
des <- mvpa_design(design_df, y_train = ~ condition, block_var = ~ block)

# Create standard blocked partitions (likely unbalanced)
cval_unbalanced <- blocked_cross_validation(des$block_var)
print("Unbalanced Counts (Example Fold 1 Train):")
print(table(des$y_train[unlist(crossval_samples(cval_unbalanced,
          design_df, des$y_train)$train[[1]]$idx)]))

# Balance using sub-sampling (default)
cval_sub <- balance_partitions(cval_unbalanced, des, seed = 1)
print(cval_sub)
print("Subsample Balanced Counts (Example Fold 1 Train):")
print(table(crossval_samples(cval_sub, design_df, des$y_train)$ytrain[[1]]))

# Balance using over-sampling
cval_over <- balance_partitions(cval_unbalanced, des, method = "oversample", seed = 2)
print(cval_over)
print("Oversample Balanced Counts (Example Fold 1 Train):")
print(table(crossval_samples(cval_over, design_df, des$y_train)$ytrain[[1]]))




cleanEx()
nameEx("banded_ridge_da")
### * banded_ridge_da

flush(stderr()); flush(stdout())

### Name: banded_ridge_da
### Title: Convenience wrapper: build a grouped-ridge domain-adaptation
###   model from matrices
### Aliases: banded_ridge_da grouped_ridge_da

### ** Examples

## Not run: 
##D ms <- grouped_ridge_da(
##D   dataset = dset,
##D   X_train = X_enc,
##D   spec = blocks(low = 100, mid = 100, high = 100, sem = 100),
##D   gamma = gamma,
##D   block_var_test = recall_runs,
##D   mode = "stacked",
##D   lambdas = c(low = 10, mid = 10, high = 10, sem = 10),
##D   alpha_recall = 0.2
##D )
## End(Not run)



cleanEx()
nameEx("banded_ridge_da_model")
### * banded_ridge_da_model

flush(stderr()); flush(stdout())

### Name: banded_ridge_da_model
### Title: Grouped (banded) ridge domain-adaptation model (continuous
###   predictors → brain)
### Aliases: banded_ridge_da_model grouped_ridge_da_model

### ** Examples

## Not run: 
##D # Build encoding predictors and declare feature sets
##D fs_enc <- feature_sets(X_enc, blocks(low = 100, mid = 100, high = 100, sem = 100))
##D 
##D # Build recall predictors from a soft alignment posterior gamma
##D fs_rec <- expected_features(fs_enc, gamma, drop_null = TRUE, renormalize = FALSE)
##D 
##D # Create design and model spec
##D des <- feature_sets_design(fs_enc, fs_rec, block_var_test = recall_runs)
##D ms <- grouped_ridge_da_model(
##D   dataset = dset,
##D   design = des,
##D   mode = "coupled",
##D   lambdas = c(low = 10, mid = 10, high = 10, sem = 10),
##D   alpha_recall = 0.2,
##D   rho = 5,
##D   compute_delta_r2 = TRUE,
##D   return_diagnostics = TRUE
##D )
##D 
##D res <- run_regional(ms, region_mask)
## End(Not run)



cleanEx()
nameEx("blocks")
### * blocks

flush(stderr()); flush(stdout())

### Name: blocks
### Title: Define consecutive column blocks as feature sets
### Aliases: blocks

### ** Examples

spec <- blocks(low = 100, mid = 100, high = 100, sem = 100)



cleanEx()
nameEx("by_set")
### * by_set

flush(stderr()); flush(stdout())

### Name: by_set
### Title: Define feature sets by a per-column set label
### Aliases: by_set

### ** Examples

X <- matrix(rnorm(10 * 6), 10, 6)
set <- rep(c("audio", "vision"), each = 3)
fs <- feature_sets(X, by_set(set, order = c("vision", "audio")))
fs



cleanEx()
nameEx("category_rdm")
### * category_rdm

flush(stderr()); flush(stdout())

### Name: category_rdm
### Title: Create Hypothesis RDM from Category Structure
### Aliases: category_rdm

### ** Examples

# Create category structure
categories <- c(cat = "animal", dog = "animal", bird = "animal",
               car = "vehicle", plane = "vehicle", boat = "vehicle")

# Create category-based RDM
rdm <- category_rdm(categories)

# Custom similarity values
rdm <- category_rdm(categories, 
                   within_category_sim = 0.9,
                   between_category_sim = 0.1)




cleanEx()
nameEx("classification_result")
### * classification_result

flush(stderr()); flush(stdout())

### Name: classification_result
### Title: Create a 'classification_result' instance
### Aliases: classification_result

### ** Examples

# A vector of observed values
yobs <- factor(rep(letters[1:4], 5))

# Predicted probabilities
probs <- data.frame(a = runif(1:20), b = runif(1:20), c = runif(1:20), d = runif(1:20))
probs <- sweep(probs, 1, rowSums(probs), "/")

# Get the max probability per row and use this to determine the predicted class
maxcol <- max.col(probs)
predicted <- levels(yobs)[maxcol]

# Construct a classification result
cres <- classification_result(yobs, predicted, probs)

# Compute default performance measures (Accuracy, AUC)
performance(cres)



cleanEx()
nameEx("combine_prediction_tables")
### * combine_prediction_tables

flush(stderr()); flush(stdout())

### Name: combine_prediction_tables
### Title: Combine prediction tables
### Aliases: combine_prediction_tables

### ** Examples

# Create example prediction tables
observed = factor(sample(letters[1:2], 10, replace = TRUE))
predtab1 <- data.frame(.rownum = 1:10,
                       roinum = rep(1, 10),
                       observed = observed,
                       prob_A = runif(10),
                       prob_B = runif(10))
predtab2 <- data.frame(.rownum = 1:10,
                       roinum = rep(2, 10),
                       observed = observed,
                       prob_A = runif(10),
                       prob_B = runif(10))

# Combine the tables
combined_table <- combine_prediction_tables(list(predtab1, predtab2))



cleanEx()
nameEx("compute_performance-methods")
### * compute_performance-methods

flush(stderr()); flush(stdout())

### Name: compute_performance
### Title: Compute Performance for an Object
### Aliases: compute_performance compute_performance.mvpa_model
### Keywords: internal

### ** Examples

cres <- binary_classification_result(
  observed  = factor(c("a", "b")),
  predicted = factor(c("a", "b")),
  probs     = matrix(c(0.8, 0.2, 0.3, 0.7), ncol = 2,
                     dimnames = list(NULL, c("a", "b")))
)
dummy <- list(performance = performance)
class(dummy) <- "mvpa_model"
compute_performance(dummy, cres)



cleanEx()
nameEx("contrast_rsa_model")
### * contrast_rsa_model

flush(stderr()); flush(stdout())

### Name: contrast_rsa_model
### Title: Constructor for contrast_rsa_model
### Aliases: contrast_rsa_model

### ** Examples

# --- Minimal Setup ---
# 1. Create dummy data and an mvpa_dataset

  # Dummy data: 16 samples, 10 voxels, 4 conditions, 2 runs
  set.seed(123)
  n_samples <- 16
  n_voxels <- 10
  n_conditions <- 4 # condA, condB, condC, condD
  n_runs <- 2

  dummy_array <- array(rnorm(n_voxels * n_samples), c(n_voxels, 1, 1, n_samples))
  dummy_space <- neuroim2::NeuroSpace(c(n_voxels, 1, 1, n_samples))
  dummy_sl_vec <- neuroim2::NeuroVec(dummy_array, dummy_space)

  dummy_mask <- neuroim2::NeuroVol(
    array(1, c(n_voxels, 1, 1)),
    neuroim2::NeuroSpace(c(n_voxels, 1, 1))
  )

  condition_labels <- factor(
    rep(paste0("cond", LETTERS[1:n_conditions]),
        each = n_samples / n_conditions)
  )
  run_labels <- factor(rep(1:n_runs, each = n_samples / n_runs))

  # Create mvpa_dataset (without Y and block_var)
  mvpa_dat <- rMVPA::mvpa_dataset(
    train_data = dummy_sl_vec,
    mask = dummy_mask
  )

  # Create mvpa_design
  mvpa_des <- rMVPA::mvpa_design(
    train_design = data.frame(condition = condition_labels, run = run_labels),
    y_train = ~condition,
    block_var = ~run
  )

  condition_levels <- levels(rMVPA::y_train(mvpa_des))
  K <- length(condition_levels)

  C_mat <- matrix(0, nrow = K, ncol = 2)
  rownames(C_mat) <- condition_levels
  C_mat["condA", 1] <- 1; C_mat["condB", 1] <- 1
  C_mat["condC", 1] <- -1; C_mat["condD", 1] <- -1
  C_mat["condA", 2] <- 1; C_mat["condB", 2] <- -1
  C_mat <- base::scale(C_mat, center = TRUE, scale = FALSE)
  colnames(C_mat) <- c("AB_vs_CD", "A_vs_B")

  msreve_des <- rMVPA::msreve_design(
    mvpa_design = mvpa_des, # Use mvpa_des here
    contrast_matrix = C_mat
  )

  # --- Example 1: Basic contrast_rsa_model ---
  model_basic <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des
  )
  print(model_basic)

  # --- Example 1b: Requesting multiple metrics ---
  model_multi_metric <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des,
    output_metric = c("beta_delta", "recon_score", "beta_only")
  )
  print(model_multi_metric)

  # --- Example 2: Using L2_norm for U_hat and normalize_delta ---
  model_l2_norm_delta <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des,
    estimation_method = "L2_norm",
    normalize_delta = TRUE,
    output_metric = "beta_delta_norm"
  )
  print(model_l2_norm_delta)

  # --- Example 3: Ridge Regression (HKB) ---
  model_ridge <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des,
    regression_type = "ridge_hkb"
  )
  print(model_ridge)

  # --- Example 4: Reconstruction Score Output ---
  model_recon <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des,
    output_metric = "recon_score"
  )
  print(model_recon)

  # --- Example 5: Composite Score Output ---
  C_mat_ortho <- rMVPA::orthogonalize_contrasts(C_mat)
  msreve_des_ortho <- rMVPA::msreve_design(
      mvpa_design = mvpa_des, # Use mvpa_des here
      contrast_matrix = C_mat_ortho
  )
  print(paste("Is contrast matrix orthonormal:", attr(msreve_des_ortho, "is_orthonormal")))

  model_composite <- contrast_rsa_model(
    dataset = mvpa_dat,
    design = msreve_des_ortho,
    output_metric = "composite",
    normalize_delta = TRUE 
  )
  print(model_composite)

  # --- Example 6: Crossnobis estimation_method ---
  # This only shows setting the method. Actual training would require passing
  # a pre-computed whitening_matrix_W to compute_crossvalidated_means_sl,
  # which is called by train_model.contrast_rsa_model.
  model_crossnobis <- contrast_rsa_model(
      dataset = mvpa_dat,
      design = msreve_des,
      estimation_method = "crossnobis"
  )
  print(model_crossnobis)



cleanEx()
nameEx("contrasts")
### * contrasts

flush(stderr()); flush(stdout())

### Name: contrasts
### Title: Generate Contrast Matrices
### Aliases: contrasts

### ** Examples

labs <- c("faces","animals","plants","tools",
          "vehicles","furniture","buildings","food")

# 1) Mini-DSL: 2x2 Factorial (Animacy x Size) + Interaction, Orthonormal
C1 <- contrasts(
        labels = labs,
        spec   = ~ anim( faces + animals + plants + food ~ . )
                 + size( faces + animals + tools + furniture ~ . )
                 + anim:size,
        orth   = TRUE)
print(colnames(C1))
print(attr(C1, "source"))
print(round(crossprod(C1), 5))

# 2) Mini-DSL: One-vs-rest, Centered, Unit Length (L2)
C2 <- contrasts(labels = labs,
                spec   = ~ faces( faces ~ . ) + tools( tools ~ . ),
                scale = "l2")
print(round(colSums(C2^2), 5)) # Should be 1

# 3) Metadata + Formula: Centered, Scaled (SD)
meta <- tibble::tribble(
  ~label,      ~anim, ~size,
  "faces",        1,    0,
  "animals",      1,    0,
  "plants",       1,    1,
  "tools",        0,    0,
  "vehicles",     0,    1,
  "furniture",    0,    0,
  "buildings",    0,    1,
  "food",         1,    1)
# Note: labels argument is ignored here, order comes from meta$label
# Also note: This function masks stats::contrasts
C3 <- contrasts(metadata = meta,
                spec     = ~ anim + size + anim:size,
                scale    = "sd")
print(round(colMeans(C3), 5)) # Should be 0
print(round(apply(C3, 2, sd), 5)) # Should be 1



cleanEx()
nameEx("create_dist")
### * create_dist

flush(stderr()); flush(stdout())

### Name: create_dist
### Title: Create a Distance Function Object
### Aliases: create_dist

### ** Examples

# Create a Euclidean distance function object
dist_obj_euc <- create_dist("euclidean")

# Create a correlation distance function object with a specified correlation method
dist_obj_cor <- create_dist("cordist", method="spearman")




cleanEx()
nameEx("cross_validation")
### * cross_validation

flush(stderr()); flush(stdout())

### Name: bootstrap_blocked_cross_validation
### Title: bootstrap_blocked_cross_validation
### Aliases: bootstrap_blocked_cross_validation blocked_cross_validation
###   sequential_blocked_cross_validation custom_cross_validation

### ** Examples

block_var <- rep(1:5, each=50)
weights <- runif(length(block_var))
weights[1] = 0
cval <- bootstrap_blocked_cross_validation(block_var, weights=weights)
X <- matrix(rnorm(length(block_var) * 10), length(block_var), 10)
y <- rep(letters[1:5], length.out=length(block_var))

sam <- crossval_samples(cval, as.data.frame(X), y)
block_var <- rep(1:5, each=50)
cval <- blocked_cross_validation(block_var)
X <- matrix(rnorm(length(block_var) * 10), length(block_var), 10)
y <- rep(letters[1:5], length.out=length(block_var))

sam <- crossval_samples(cval, as.data.frame(X), y)
block_var <- rep(1:5, each=50)
nfolds <- 2
nreps <- 4
cval <- sequential_blocked_cross_validation(block_var, nfolds, nreps)
X <- matrix(rnorm(length(block_var) * 10), length(block_var), 10)
y <- rep(letters[1:5], length.out=length(block_var))

sam <- crossval_samples(cval, as.data.frame(X), y)
sample_set <- list(
  list(train = 1:80, test = 81:100),
  list(train = 1:60, test = 61:100),
  list(train = 1:40, test = 41:100)
)
cval <- custom_cross_validation(sample_set)
X <- matrix(rnorm(100 * 10), 100, 10)
y <- rep(letters[1:4], length.out=100)

sam <- crossval_samples(cval, as.data.frame(X), y)



cleanEx()
nameEx("crossv_block")
### * crossv_block

flush(stderr()); flush(stdout())

### Name: crossv_block
### Title: Block Cross-Validation Data Preparation
### Aliases: crossv_block

### ** Examples

X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
y <- rep(letters[1:4], 25)
block_var <- rep(1:4, each = 25)
cv <- crossv_block(X, y, block_var)



cleanEx()
nameEx("crossv_k")
### * crossv_k

flush(stderr()); flush(stdout())

### Name: crossv_k
### Title: K-fold Cross-Validation Data Preparation
### Aliases: crossv_k

### ** Examples

data <- iris[,-5]
y <- iris$Species
result <- crossv_k(data, y, k = 5)



cleanEx()
nameEx("crossval_samples")
### * crossval_samples

flush(stderr()); flush(stdout())

### Name: crossval_samples
### Title: Cross-validation samples
### Aliases: crossval_samples
###   crossval_samples.sequential_blocked_cross_validation
###   crossval_samples.kfold_cross_validation
###   crossval_samples.blocked_cross_validation
###   crossval_samples.bootstrap_blocked_cross_validation
###   crossval_samples.custom_cross_validation
###   crossval_samples.twofold_blocked_cross_validation
###   crossval_samples.mvpa_model

### ** Examples

cval <- kfold_cross_validation(len = 20, nfolds = 4)
dat  <- as.data.frame(matrix(rnorm(20 * 2), 20, 2))
y    <- factor(rep(letters[1:4], 5))
crossval_samples(cval, dat, y)



cleanEx()
nameEx("custom_performance")
### * custom_performance

flush(stderr()); flush(stdout())

### Name: custom_performance
### Title: Apply Custom Performance Metric to Prediction Result
### Aliases: custom_performance

### ** Examples

cres <- binary_classification_result(
  observed  = factor(c("A", "B")),
  predicted = factor(c("A", "A")),
  probs = matrix(c(0.9, 0.1,
                   0.6, 0.4),
                 ncol = 2, byrow = TRUE,
                 dimnames = list(NULL, c("A", "B")))
)
acc_fun <- function(x) c(acc = mean(x$observed == x$predicted))
custom_performance(cres, acc_fun)



cleanEx()
nameEx("distance-constructors")
### * distance-constructors

flush(stderr()); flush(stdout())

### Name: cordist
### Title: Distance Function Constructors
### Aliases: cordist mahadist eucdist euclidean robustmahadist pcadist

### ** Examples

dist_obj_1 <- cordist(method="pearson")
dist_obj_2 <- mahadist()
dist_obj_3 <- eucdist()
dist_obj_4 <- robustmahadist()
dist_obj_5 <- pcadist(ncomp=2, dist_method="cosine")




cleanEx()
nameEx("expected_features")
### * expected_features

flush(stderr()); flush(stdout())

### Name: expected_features
### Title: Build expected-domain features from a soft alignment matrix
### Aliases: expected_features

### ** Examples

X <- matrix(rnorm(10 * 5), 10, 5)
fs_enc <- feature_sets(X, blocks(a = 2, b = 3))
gamma <- matrix(runif(6 * 10), 6, 10)
gamma <- gamma / rowSums(gamma)
fs_rec <- expected_features(fs_enc, gamma, drop_null = FALSE, renormalize = TRUE)



cleanEx()
nameEx("feature_selector")
### * feature_selector

flush(stderr()); flush(stdout())

### Name: feature_selector
### Title: Create a feature selection specification
### Aliases: feature_selector

### ** Examples

fsel <- feature_selector("FTest", "top_k", 1000)
fsel <- feature_selector("FTest", "top_p", .1)
class(fsel) == "FTest"



cleanEx()
nameEx("feature_sets")
### * feature_sets

flush(stderr()); flush(stdout())

### Name: feature_sets
### Title: Feature sets: grouped predictor matrices
### Aliases: feature_sets

### ** Examples

X <- matrix(rnorm(20 * 8), 20, 8)
fs <- feature_sets(X, blocks(low = 3, sem = 5))
fs
# 1) Matrix input + blocks()
X <- matrix(rnorm(20 * 8), 20, 8)
fs <- feature_sets(X, blocks(low = 3, sem = 5))

# 2) List input (already split per set)
Xlist <- list(low = X[, 1:3], sem = X[, 4:8])
fs2 <- feature_sets(Xlist)



cleanEx()
nameEx("feature_sets_design")
### * feature_sets_design

flush(stderr()); flush(stdout())

### Name: feature_sets_design
### Title: Feature-sets design (mvpa_design extension for continuous
###   regression)
### Aliases: feature_sets_design

### ** Examples

# Train predictors (TR × features), split into named sets:
X_enc <- matrix(rnorm(20 * 8), 20, 8)
fs_enc <- feature_sets(X_enc, blocks(low = 3, sem = 5))

# Test predictors (TR × features), for example from a soft alignment:
gamma <- matrix(runif(10 * 20), 10, 20)
gamma <- gamma / rowSums(gamma)
fs_rec <- expected_features(fs_enc, gamma, drop_null = FALSE, renormalize = TRUE)

des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1:2, each = 5))



cleanEx()
nameEx("fit_model-methods")
### * fit_model-methods

flush(stderr()); flush(stdout())

### Name: fit_model
### Title: Fit Model
### Aliases: fit_model fit_model.mvpa_model

### ** Examples





cleanEx()
nameEx("format_result")
### * format_result

flush(stderr()); flush(stdout())

### Name: format_result
### Title: Format Result Object
### Aliases: format_result format_result.mvpa_model
###   format_result.feature_rsa_model

### ** Examples




cleanEx()
nameEx("gen_sample_dataset")
### * gen_sample_dataset

flush(stderr()); flush(stdout())

### Name: gen_sample_dataset
### Title: Generate Sample Dataset for MVPA Analysis
### Aliases: gen_sample_dataset

### ** Examples

# Generate categorical image dataset
dataset <- gen_sample_dataset(
  D = c(10,10,10),
  nobs = 100,
  response_type = "categorical",
  data_mode = "image",
  blocks = 3,
  nlevels = 2
)

# Generate continuous surface dataset
surf_data <- gen_sample_dataset(
  D = 1000,  # number of vertices
  nobs = 50,
  response_type = "continuous",
  data_mode = "surface"
)

# Generate dataset with external test set
test_dataset <- gen_sample_dataset(
  D = c(8,8,8),
  nobs = 80,
  response_type = "categorical",
  nlevels = 3,
  external_test = TRUE
)




cleanEx()
nameEx("get_nfolds")
### * get_nfolds

flush(stderr()); flush(stdout())

### Name: get_nfolds
### Title: Get the Number of Folds
### Aliases: get_nfolds

### ** Examples

cval <- kfold_cross_validation(len = 20, nfolds = 4)
get_nfolds(cval)



cleanEx()
nameEx("group_means")
### * group_means

flush(stderr()); flush(stdout())

### Name: group_means
### Title: Compute Group Means of a Matrix
### Aliases: group_means

### ** Examples

# Create a random matrix
data <- matrix(rnorm(100 * 100), 100, 100)

# Define a grouping variable
groups <- factor(rep(1:5, each = 20))

# Calculate group means for each row
row_means <- group_means(data, margin = 1, group = groups)

# Calculate group means for each column
col_means <- group_means(data, margin = 2, group = groups)



cleanEx()
nameEx("has_crossval-methods")
### * has_crossval-methods

flush(stderr()); flush(stdout())

### Name: has_crossval
### Title: Cross-Validation Availability
### Aliases: has_crossval

### ** Examples

ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10)
cval <- blocked_cross_validation(ds$design$block_var)
mdl <- load_model("sda_notune")
mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                    "classification", crossval = cval)
has_crossval(mspec)



cleanEx()
nameEx("has_test_set-methods")
### * has_test_set-methods

flush(stderr()); flush(stdout())

### Name: has_test_set
### Title: Test Set Availability
### Aliases: has_test_set has_test_set.mvpa_design

### ** Examples

ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10, external_test = TRUE)
has_test_set(ds$design)



cleanEx()
nameEx("hrfdecoder_design")
### * hrfdecoder_design

flush(stderr()); flush(stdout())

### Name: hrfdecoder_design
### Title: Construct an hrfdecoder design (explicit mvpa_design extension)
### Aliases: hrfdecoder_design

### ** Examples

## Not run: 
##D library(fmridesign)
##D library(fmrihrf)
##D 
##D # Create event table
##D events_df <- data.frame(
##D   onset = c(5, 15, 35, 45),
##D   condition = factor(c("A", "B", "A", "B")),
##D   run = c(1, 1, 2, 2)
##D )
##D 
##D # Define sampling frame (TR structure)
##D sframe <- sampling_frame(blocklens = c(60, 60), TR = 2)
##D 
##D # Build event model
##D evmod <- event_model(
##D   onset ~ hrf(condition, basis = "spmg1"),
##D   data = events_df,
##D   block = ~run,
##D   sampling_frame = sframe
##D )
##D 
##D # Create hrfdecoder design
##D block_var <- rep(1:2, each = 60)  # 60 TRs per run
##D design <- hrfdecoder_design(
##D   event_model = evmod,
##D   events = events_df,
##D   block_var = block_var
##D )
## End(Not run)



cleanEx()
nameEx("hrfdecoder_model")
### * hrfdecoder_model

flush(stderr()); flush(stdout())

### Name: hrfdecoder_model
### Title: hrfdecoder Model (continuous-time MVPA)
### Aliases: hrfdecoder_model

### ** Examples

## Not run: 
##D library(fmridesign)
##D library(fmrihrf)
##D library(hrfdecoder)
##D 
##D # 1. Create event table
##D events_df <- data.frame(
##D   onset = seq(10, 290, by = 20),  # Events every 20 seconds
##D   condition = rep(c("A", "B", "C"), length.out = 15),
##D   run = rep(1:3, each = 5)
##D )
##D 
##D # 2. Define temporal structure
##D sframe <- sampling_frame(blocklens = c(100, 100, 100), TR = 2)
##D 
##D # 3. Build event model
##D evmod <- event_model(
##D   onset ~ hrf(condition, basis = "spmg3"),
##D   data = events_df,
##D   block = ~run,
##D   sampling_frame = sframe
##D )
##D 
##D # 4. Create rMVPA dataset
##D mask <- neuroim2::NeuroVol(...)
##D fmri_data <- neuroim2::NeuroVec(...)  # 300 TRs
##D dset <- mvpa_dataset(train_data = fmri_data, mask = mask)
##D 
##D # 5. Create hrfdecoder design
##D block_var <- rep(1:3, each = 100)
##D design <- hrfdecoder_design(
##D   event_model = evmod,
##D   events = events_df,
##D   block_var = block_var
##D )
##D 
##D # 6. Specify model with custom hyperparameters
##D mspec <- hrfdecoder_model(
##D   dataset = dset,
##D   design = design,
##D   lambda_W = 10,
##D   lambda_HRF = 1,
##D   lambda_smooth = 5,
##D   basis = fmrihrf::spmg1(),
##D   window = c(4, 8),
##D   max_iter = 15
##D )
##D 
##D # 7. Run searchlight analysis
##D results <- run_searchlight(mspec, radius = 8, method = "randomized", niter = 4)
## End(Not run)



cleanEx()
nameEx("kfold_cross_validation")
### * kfold_cross_validation

flush(stderr()); flush(stdout())

### Name: kfold_cross_validation
### Title: kfold_cross_validation
### Aliases: kfold_cross_validation

### ** Examples

cval <- kfold_cross_validation(len=100, nfolds=10)
sample_data <- as.data.frame(matrix(rnorm(100*10), 100, 10))
sample_y <- rep(letters[1:5], 20)
samples <- crossval_samples(cval, data = sample_data, y = sample_y)
stopifnot(nrow(samples) == 10)



cleanEx()
nameEx("load_model")
### * load_model

flush(stderr()); flush(stdout())

### Name: load_model
### Title: Load a Pre-defined MVPA Model
### Aliases: load_model

### ** Examples

# Load custom MVPA model
model <- load_model("sda_notune")

# Load correlation classifier with parameter tuning options
corr_model <- load_model("corclass")
print(corr_model$parameters)  # View tunable parameters




cleanEx()
nameEx("make_feature_contrasts")
### * make_feature_contrasts

flush(stderr()); flush(stdout())

### Name: make_feature_contrasts
### Title: Generate Contrasts from a Feature Matrix (Optional PCA)
### Aliases: make_feature_contrasts

### ** Examples

# Example feature matrix (4 conditions, 5 features)
feat_mat <- matrix(rnorm(20), nrow = 4,
                   dimnames = list(paste0("Cond", 1:4), paste0("F", 1:5)))

# Use raw features (first 3)
C_raw <- make_feature_contrasts(feat_mat[, 1:3], use_pca = FALSE, prefix="RawFeat_")
print(C_raw)

# Use PCA, selecting top 2 PCs
C_pca <- make_feature_contrasts(feat_mat, use_pca = TRUE, n_pcs = 2, prefix="PCA_")
print(C_pca)

# Use PCA, selecting >= 80% variance explained
C_pca_pve <- make_feature_contrasts(feat_mat, use_pca = TRUE, pve = 0.8, prefix="PCA_")
print(C_pca_pve)

# Reorder based on labels
C_pca_reorder <- make_feature_contrasts(feat_mat, labels=c("Cond3", "Cond1", "Cond4", "Cond2"),
                                      use_pca = TRUE, n_pcs = 2, prefix="PCA_")
print(C_pca_reorder)




cleanEx()
nameEx("manova_design")
### * manova_design

flush(stderr()); flush(stdout())

### Name: manova_design
### Title: Create a MANOVA Design
### Aliases: manova_design

### ** Examples

# Create simple dissimilarity matrices
dissimilarity_matrix_y  <- matrix(rnorm(9), nrow = 3)
dissimilarity_matrix_x1 <- matrix(rnorm(9), nrow = 3)
dissimilarity_matrix_x2 <- matrix(rnorm(9), nrow = 3)

# Create a MANOVA design
formula   <- y ~ x1 + x2
data_list <- list(
  y  = dissimilarity_matrix_y,
  x1 = dissimilarity_matrix_x1,
  x2 = dissimilarity_matrix_x2
)
manova_design_obj <- manova_design(formula, data_list)



cleanEx()
nameEx("manova_model")
### * manova_model

flush(stderr()); flush(stdout())

### Name: manova_model
### Title: Create a MANOVA Model
### Aliases: manova_model

### ** Examples

## Not run: 
##D # Create a MANOVA model using gen_sample_dataset
##D dset <- gen_sample_dataset(D = c(5, 5, 5), nobs = 50, nlevels = 3, blocks = 3)
##D 
##D # Create dissimilarity matrices for MANOVA design
##D formula <- y ~ x1 + x2
##D data_list <- list(
##D   y = matrix(rnorm(9), nrow = 3),
##D   x1 = matrix(rnorm(9), nrow = 3),
##D   x2 = matrix(rnorm(9), nrow = 3)
##D )
##D design <- manova_design(formula, data_list)
##D manova_model_obj <- manova_model(dset$dataset, design)
## End(Not run)



cleanEx()
nameEx("merge_results")
### * merge_results

flush(stderr()); flush(stdout())

### Name: merge_results
### Title: Merge Multiple Results
### Aliases: merge_results

### ** Examples





cleanEx()
nameEx("msreve_design")
### * msreve_design

flush(stderr()); flush(stdout())

### Name: msreve_design
### Title: Constructor for msreve_design
### Aliases: msreve_design

### ** Examples

# Assume 'mvpa_des_obj' is a pre-existing mvpa_design object
# e.g. from mvpa_design(data=my_data_frame, formula = ~ condition_labels + run_labels,
#                       block_var = "run_labels")
# Let\'s say mvpa_des_obj implies 6 conditions based on unique(my_data_frame$condition_labels)
K <- 6 # Number of conditions
Q <- 2 # Number of contrasts

# Example contrast matrix (K x Q)
C_mat <- matrix(c(
 # C1: Cond 1,2,3 vs 4,5,6
  1,  1,  1, -1, -1, -1,
 # C2: Cond 1,2 vs 3 (and 0 for 4,5,6 for simplicity here)
  1,  1, -2,  0,  0,  0
), nrow = K, ncol = Q, byrow = FALSE)
colnames(C_mat) <- c("GroupComparison", "SubComparison")

# if (inherits(mvpa_des_obj, "mvpa_design")) {
#  design_obj <- msreve_design(mvpa_des_obj, C_mat, name="example_msreve")
#  print(design_obj)
# }
#
# # Automatically add pairwise interactions
# design_obj_int <- msreve_design(mvpa_des_obj, C_mat,
#                                include_interactions = TRUE)
# colnames(design_obj_int$contrast_matrix)



cleanEx()
nameEx("mvpa_dataset")
### * mvpa_dataset

flush(stderr()); flush(stdout())

### Name: mvpa_dataset
### Title: Create an MVPA Dataset Object
### Aliases: mvpa_dataset

### ** Examples

# Use gen_sample_dataset helper to create a simple dataset
sample_data <- gen_sample_dataset(c(5, 5, 5), nobs = 100, blocks = 4)
dataset <- sample_data$dataset

# Access components
print(dim(dataset$train_data))
print(sum(dataset$mask > 0))




cleanEx()
nameEx("mvpa_design")
### * mvpa_design

flush(stderr()); flush(stdout())

### Name: mvpa_design
### Title: Create an MVPA Design Object
### Aliases: mvpa_design

### ** Examples

# Basic design with only training data
train_df <- data.frame(condition = rep(c("A", "B"), each = 50),
                       block = rep(1:5, each = 20),
                       group = rep(c("Group1", "Group2"), 50))
design <- mvpa_design(train_df, y_train = ~ condition)

# Design with test data and blocking variable
test_df <- data.frame(condition = rep(c("A", "B"), each = 25))
design_with_test <- mvpa_design(
  train_df, 
  test_df, 
  y_train = ~ condition, 
  y_test = ~ condition,
  block_var = ~ block
)

# Design with split_by factor
design_split <- mvpa_design(
  train_df, 
  y_train = ~ condition,
  split_by = ~ group
)




cleanEx()
nameEx("mvpa_model")
### * mvpa_model

flush(stderr()); flush(stdout())

### Name: mvpa_model
### Title: Create an MVPA Model
### Aliases: mvpa_model

### ** Examples


mod <- load_model("sda")
arr_data <- array(rnorm(6*6*6*100), c(6,6,6,100))
sp <- neuroim2::NeuroSpace(c(6,6,6,100))
traindat <- neuroim2::NeuroVec(arr_data, sp)
mask <- neuroim2::LogicalNeuroVol(array(rnorm(6*6*6)>-.2, c(6,6,6)), neuroim2::NeuroSpace(c(6,6,6)))

mvdset <- mvpa_dataset(traindat,mask=mask)
design <- data.frame(fac=rep(letters[1:4], 25), block=rep(1:10, each=10))
cval <- blocked_cross_validation(design$block)
mvdes <- mvpa_design(design, y_train = ~ fac, block_var = ~ block)

custom_perf <- function(result) {
  c(accuracy=sum(result$observed == result$predicted)/length(result$observed))
}
mvpmod <- mvpa_model(mod, dataset=mvdset, design=mvdes, crossval=cval, performance=custom_perf)



cleanEx()
nameEx("mvpa_surface_dataset")
### * mvpa_surface_dataset

flush(stderr()); flush(stdout())

### Name: mvpa_surface_dataset
### Title: Create a Surface-Based MVPA Dataset Object
### Aliases: mvpa_surface_dataset

### ** Examples

## Not run: 
##D # Create surface dataset with automatic mask
##D train_surf <- NeuroSurfaceVector(geometry, data)
##D dataset <- mvpa_surface_dataset(train_surf, name="lh")
##D 
##D # Create dataset with test data and custom mask
##D test_surf <- NeuroSurfaceVector(geometry, test_data)
##D mask <- numeric(length(nodes(geometry)))
##D mask[roi_indices] <- 1
##D dataset <- mvpa_surface_dataset(train_surf, test_surf, mask, name="rh")
## End(Not run)




cleanEx()
nameEx("mvpa_sysinfo")
### * mvpa_sysinfo

flush(stderr()); flush(stdout())

### Name: mvpa_sysinfo
### Title: Report System and Package Information for rMVPA
### Aliases: mvpa_sysinfo

### ** Examples

## Not run: 
##D # Display system information in the console
##D mvpa_sysinfo()
##D 
##D # Capture the information in a variable
##D sys_info <- mvpa_sysinfo()
##D print(sys_info$r_version)
##D print(sys_info$dependencies$rsample)
## End(Not run)



cleanEx()
nameEx("nmf_preprocess_maps")
### * nmf_preprocess_maps

flush(stderr()); flush(stdout())

### Name: nmf_preprocess_maps
### Title: Preprocess Maps for Spatial NMF
### Aliases: nmf_preprocess_maps

### ** Examples

## Not run: 
##D # For AUC-0.5 maps (chance-centered)
##D prepped <- nmf_preprocess_maps(auc_maps, method = "auc")
##D result <- spatial_nmf_maps(prepped$maps, mask = mask, k = 5)
##D 
##D # For raw AUC maps
##D prepped <- nmf_preprocess_maps(auc_maps, method = "auc_raw")
##D 
##D # For z-score maps with small positive floor
##D prepped <- nmf_preprocess_maps(zmaps, method = "shift", min_val = 0.01)
## End(Not run)




cleanEx()
nameEx("orthogonalize_contrasts")
### * orthogonalize_contrasts

flush(stderr()); flush(stdout())

### Name: orthogonalize_contrasts
### Title: Orthogonalize a Contrast Matrix
### Aliases: orthogonalize_contrasts

### ** Examples

K <- 6 # Number of conditions
Q <- 2 # Number of contrasts
C_orig <- matrix(c( 1,  1,  1, -1, -1, -1,  # Contrast 1
                    1, -1,  0,  1, -1,  0), # Contrast 2 (not orthogonal to C1)
                 nrow=K, ncol=Q)
colnames(C_orig) <- c("MainEffect", "InteractionLike")
C_ortho <- orthogonalize_contrasts(C_orig)
# print(round(crossprod(C_ortho), 10)) # Should be close to identity matrix

# Example with a rank-deficient matrix (3rd contrast is sum of first two)
C_rank_def <- cbind(C_orig, C_orig[,1] + C_orig[,2])
colnames(C_rank_def) <- c("C1", "C2", "C3_dependent")
C_ortho_def <- orthogonalize_contrasts(C_rank_def)
# print(round(crossprod(C_ortho_def), 10))
# The 3rd column of C_ortho_def will be zeros.



cleanEx()
nameEx("performance-methods")
### * performance-methods

flush(stderr()); flush(stdout())

### Name: performance
### Title: Compute Performance Metrics
### Aliases: performance

### ** Examples

cres <- binary_classification_result(
  observed  = factor(c("a", "b")),
  predicted = factor(c("a", "b")),
  probs     = matrix(c(0.8, 0.2, 0.3, 0.7), ncol = 2,
                     dimnames = list(NULL, c("a", "b")))
)
performance(cres)



cleanEx()
nameEx("predicted_class")
### * predicted_class

flush(stderr()); flush(stdout())

### Name: predicted_class
### Title: Calculate the Predicted Class from Probability Matrix
### Aliases: predicted_class

### ** Examples

prob <- matrix(c(0.2, 0.8,
                 0.6, 0.4),
               nrow = 2, byrow = TRUE,
               dimnames = list(NULL, c("A", "B")))
predicted_class(prob)



cleanEx()
nameEx("prep_regional")
### * prep_regional

flush(stderr()); flush(stdout())

### Name: prep_regional
### Title: Prepare regional data for MVPA analysis
### Aliases: prep_regional

### ** Examples

# Create example data
sample_data <- gen_sample_dataset(c(5, 5, 5), nobs = 100, blocks = 4)

# Create a simple region mask with 3 ROIs
mask_vol <- sample_data$dataset$mask
region_mask <- neuroim2::NeuroVol(
  sample(1:3, size = sum(mask_vol > 0), replace = TRUE),
  space = neuroim2::space(mask_vol),
  indices = which(mask_vol > 0)
)

# Create a basic model spec
model_spec <- list(dataset = sample_data$dataset)

# Prepare regional data
regional_data <- prep_regional(model_spec, region_mask)



cleanEx()
nameEx("process_roi-methods")
### * process_roi-methods

flush(stderr()); flush(stdout())

### Name: process_roi
### Title: Process ROI
### Aliases: process_roi process_roi.default
###   process_roi.custom_internal_model_spec

### ** Examples





cleanEx()
nameEx("regional_mvpa_result")
### * regional_mvpa_result

flush(stderr()); flush(stdout())

### Name: regional_mvpa_result
### Title: Create a 'regional_mvpa_result' instance
### Aliases: regional_mvpa_result

### ** Examples

# Create example inputs
model_spec <- list(dataset = "Example dataset")
performance_table <- data.frame(accuracy = c(0.8, 0.85))
prediction_table <- data.frame(observed = factor(rep(letters[1:2], 5)),
                                predicted = factor(rep(letters[1:2], 5)))
vol_results <- list(vol1 = "Example vol_result 1", vol2 = "Example vol_result 2")
fits <- list(fit1 = "Example fit 1", fit2 = "Example fit 2")

# Construct a regional_mvpa_result
regional_result <- regional_mvpa_result(model_spec, performance_table,
                                        prediction_table, vol_results, fits = fits)



cleanEx()
nameEx("register_mvpa_model")
### * register_mvpa_model

flush(stderr()); flush(stdout())

### Name: register_mvpa_model
### Title: Register a Custom MVPA Model
### Aliases: register_mvpa_model

### ** Examples

## Not run: 
##D # Example of how a user might define an e1071 SVM spec
##D my_svm_spec <- list(
##D   type = "Classification", library = "e1071", label = "my_svm",
##D   parameters = data.frame(parameter = "cost", class = "numeric", label = "Cost (C)"),
##D   # grid should return a data.frame with columns matching 'parameter' names in 'parameters'
##D   grid = function(x, y, len = NULL) { 
##D      data.frame(cost = if (is.null(len) || len == 1) 1 else 10^seq(-2, 2, length.out = len))
##D   },
##D   # fit function receives: x, y, wts (weights), param (current params from grid), 
##D   # lev (levels of y), last (unused), weights (unused), classProbs (unused by e1071::svm)
##D   fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
##D      e1071::svm(x, y, cost = param$cost, probability = TRUE, ...) # Ensure probability=TRUE for prob
##D   },
##D   # predict function receives: modelFit (output of $fit), newdata
##D   predict = function(modelFit, newdata, ...) {
##D      predict(modelFit, newdata, ...)
##D   },
##D   # prob function receives: modelFit, newdata
##D   # Should return a matrix/df with columns named as in levels(y)
##D   prob = function(modelFit, newdata, ...) {
##D     pred_obj <- predict(modelFit, newdata, probability = TRUE)
##D     attr(pred_obj, "probabilities") 
##D   }
##D )
##D register_mvpa_model("my_svm", my_svm_spec)
##D # Now load_model("my_svm") would work.
## End(Not run)



cleanEx()
nameEx("rsa_design")
### * rsa_design

flush(stderr()); flush(stdout())

### Name: rsa_design
### Title: Construct a design for an RSA (Representational Similarity
###   Analysis) model
### Aliases: rsa_design

### ** Examples

dismat <- dist(matrix(rnorm(100*100), 100, 100))
rdes <- rsa_design(~ dismat, list(dismat=dismat))



cleanEx()
nameEx("rsa_model")
### * rsa_model

flush(stderr()); flush(stdout())

### Name: rsa_model
### Title: Construct an RSA (Representational Similarity Analysis) model
### Aliases: rsa_model

### ** Examples

# Create a random MVPA dataset (image data)
arr  <- array(rnorm(100 * 5), c(5, 5, 4, 5))   # 5 voxels x 5 voxels x 4 slices x 5 observations
sp   <- neuroim2::NeuroSpace(c(5, 5, 4, 5))
vec  <- neuroim2::NeuroVec(arr, sp)
mask <- neuroim2::LogicalNeuroVol(array(1, c(5, 5, 4)), neuroim2::NeuroSpace(c(5, 5, 4)))
mvpa_data <- mvpa_dataset(train_data = vec, mask = mask)

# Create two random RDMs (distance matrices) over the 5 observations
data_mat  <- matrix(rnorm(5 * 10), 5, 10)
dismat1   <- dist(data_mat)
dismat2   <- dist(matrix(rnorm(5 * 10), 5, 10))
rdes <- rsa_design(~ dismat1 + dismat2,
                   list(dismat1 = dismat1, dismat2 = dismat2))

# Create an RSA model with standard 'lm' (returns t-values):
rsa_mod <- rsa_model(mvpa_data, rdes, regtype = "lm")

# Create an RSA model enforcing non-negativity for dismat2 only:
# Requires the 'glmnet' package to be installed
# rsa_mod_nneg <- rsa_model(mvpa_data, rdes, regtype="lm",
#                          nneg = list(dismat2 = TRUE))

# Create an RSA model using 'lm' but returning semi-partial correlations:
rsa_mod_sp <- rsa_model(mvpa_data, rdes, regtype = "lm",
                        semipartial = TRUE)

# Train the model using a trial-by-feature matrix
fit_params <- train_model(rsa_mod_sp, data_mat, y = NULL, indices = NULL)
# 'fit_params' = named vector of semi-partial correlations for each predictor




cleanEx()
nameEx("run_custom_regional")
### * run_custom_regional

flush(stderr()); flush(stdout())

### Name: run_custom_regional
### Title: Run a Custom Analysis Function Regionally
### Aliases: run_custom_regional

### ** Examples

# Generate sample dataset
dset_info <- gen_sample_dataset(D = c(8,8,8), nobs = 50, nlevels = 2)
dataset_obj <- dset_info$dataset
design_obj <- dset_info$design # Not used by custom_func here, but needed for setup

# Create a region mask with 3 ROIs
mask_arr <- array(0, dim(dataset_obj$mask))
mask_arr[1:4, 1:4, 1:4] <- 1
mask_arr[5:8, 1:4, 1:4] <- 2
mask_arr[1:4, 5:8, 5:8] <- 3
region_mask_vol <- neuroim2::NeuroVol(mask_arr, neuroim2::space(dataset_obj$mask))

# Define a custom function: calculate mean and sd for each ROI
my_roi_stats <- function(roi_data, roi_info) {
  # roi_data is samples x features matrix
  # roi_info$id is the region number
  # roi_info$indices are the feature indices
  mean_signal <- mean(roi_data, na.rm = TRUE)
  sd_signal <- sd(roi_data, na.rm = TRUE)
  num_features <- ncol(roi_data)
  list(
    roi_id = roi_info$id, # Can include id if desired, or rely on output table
    mean_signal = mean_signal,
    sd_signal = sd_signal,
    n_features = num_features
  )
}

# Run the custom regional analysis



cleanEx()
nameEx("run_custom_searchlight")
### * run_custom_searchlight

flush(stderr()); flush(stdout())

### Name: run_custom_searchlight
### Title: Run a Custom Analysis Function in a Searchlight
### Aliases: run_custom_searchlight

### ** Examples

# Generate sample dataset
dset_info <- gen_sample_dataset(D = c(10, 10, 10), nobs = 30, nlevels = 2)
dataset_obj <- dset_info$dataset

# Define a custom function: calculate mean and sd within the sphere
my_sl_stats <- function(sl_data, sl_info) {
  # sl_data is samples x features_in_sphere matrix
  # sl_info contains center_index, indices, etc.
  mean_signal <- mean(sl_data, na.rm = TRUE)
  sd_signal <- sd(sl_data, na.rm = TRUE)
  n_features <- ncol(sl_data)
  list(
    mean_signal = mean_signal,
    sd_signal = sd_signal,
    n_vox_in_sphere = n_features
  )
}

# Run the custom searchlight (standard method)




cleanEx()
nameEx("run_future-methods")
### * run_future-methods

flush(stderr()); flush(stdout())

### Name: run_future
### Title: Run Future
### Aliases: run_future run_future.default
### Keywords: internal

### ** Examples

frame <- tibble::tibble(
  .id = 1:2,
  rnum = c("roi1", "roi2"),
  roi = list(1:3, 4:5),
  size = c(3, 2)
)
mod_spec <- list(process_roi = function(mod_spec, roi, rnum, ...) {
  tibble::tibble(
    result = list(mean(roi)),
    indices = list(seq_along(roi)),
    performance = list(NULL),
    id = rnum
  )
})
run_future(mod_spec, frame, NULL)




cleanEx()
nameEx("run_regional-methods")
### * run_regional-methods

flush(stderr()); flush(stdout())

### Name: run_regional
### Title: Region of Interest Based MVPA Analysis
### Aliases: run_regional run_regional_base run_regional.default
###   run_regional.mvpa_model run_regional.rsa_model
###   run_regional.vector_rsa_model run_regional.banded_ridge_da_model

### ** Examples





cleanEx()
nameEx("run_searchlight")
### * run_searchlight

flush(stderr()); flush(stdout())

### Name: run_searchlight
### Title: Run Searchlight Analysis
### Aliases: run_searchlight

### ** Examples





cleanEx()
nameEx("run_searchlight.contrast_rsa_model")
### * run_searchlight.contrast_rsa_model

flush(stderr()); flush(stdout())

### Name: run_searchlight.contrast_rsa_model
### Title: Run Searchlight Analysis for Contrast RSA Model
### Aliases: run_searchlight.contrast_rsa_model

### ** Examples

# Assuming 'spec' is a valid contrast_rsa_model object
# Standard (recommended) method:
# results <- run_searchlight(spec, radius = 8, method = "standard")
# plot(results$results[[1]]) # Plot the map for the first contrast

# Assuming contrast_rsa_model examples have run and objects like
# 'mvpa_dat', 'msreve_des', 'model_basic', 'model_recon' are available.
# This requires the setup from contrast_rsa_model examples.

if (requireNamespace("neuroim2", quietly = TRUE) && 
    requireNamespace("rMVPA", quietly = TRUE) &&
    exists("model_basic") && inherits(model_basic, "contrast_rsa_model") &&
    exists("model_recon") && inherits(model_recon, "contrast_rsa_model")) {

  # --- Example 1: Run searchlight with basic model ---
  # Use a very small radius for quick example run.
  # Actual searchlight analyses would use a more appropriate radius (e.g., 3-4 voxels).
  # With dummy data, results won't be meaningful; focus is on execution.

  message("Running searchlight example 1 (basic model, radius=1)... May take a moment.")
  sl_results_basic <- tryCatch({
    run_searchlight(model_basic, radius = 1, method = "standard")
  }, error = function(e) {
    message("Searchlight (basic model) example failed: ", e$message)
    NULL
  })
  if (!is.null(sl_results_basic)) {
    print(sl_results_basic)
  }

  # --- Example 2: Run searchlight with recon_score output ---
  message("Running searchlight example 2 (recon_score model, radius=1)... May take a moment.")
  sl_results_recon <- tryCatch({
    run_searchlight(model_recon, radius = 1, method = "standard")
  }, error = function(e) {
    message("Searchlight (recon_score model) example failed: ", e$message)
    NULL
  })
  if (!is.null(sl_results_recon)) {
    print(sl_results_recon)
  }

  # Note on Crossnobis with searchlight:
  # To run a searchlight with 'estimation_method = "crossnobis"' from 'model_crossnobis',
  # the 'whitening_matrix_W' needs to be passed through the searchlight machinery
  # to 'compute_crossvalidated_means_sl'. This typically involves passing it via
  # the `...` argument of `run_searchlight` and ensuring `mvpa_iterate` and
  # `train_model` propagate it. This advanced usage is not shown here as it
  # requires modification to the general `mvpa_iterate` or a custom processor.

} else {
 message("Skipping run_searchlight.contrast_rsa_model example execution here.")
 message("It can be time-consuming and depends on prior setup.")
}




cleanEx()
nameEx("save_results")
### * save_results

flush(stderr()); flush(stdout())

### Name: save_results
### Title: Save MVPA Results to Disk
### Aliases: save_results

### ** Examples





cleanEx()
nameEx("second_order_similarity")
### * second_order_similarity

flush(stderr()); flush(stdout())

### Name: second_order_similarity
### Title: Compute Second-Order Similarity Scores
### Aliases: second_order_similarity

### ** Examples

# Suppose we have X (10x5), a reference D (10x10), block var, and a correlation distfun:
X <- matrix(rnorm(50), 10, 5)
D <- matrix(runif(100), 10, 10)
block <- rep(1:2, each=5)
dist_obj <- cordist(method="pearson")
scores <- second_order_similarity(dist_obj, X, D, block, method="spearman")




cleanEx()
nameEx("select_features-methods")
### * select_features-methods

flush(stderr()); flush(stdout())

### Name: select_features
### Title: Select Features
### Aliases: select_features select_features.catscore select_features.FTest

### ** Examples

fsel <- feature_selector("FTest", "top_k", 2)
coords <- rbind(c(1,1,1), c(2,2,2), c(3,3,3))
space <- neuroim2::NeuroSpace(c(10,10,10))
roi_data <- matrix(rnorm(100*3), 100, 3)
ROI <- neuroim2::ROIVec(space, coords=coords, roi_data)
Y <- factor(rep(c("a", "b"), each=50))
featureMask <- select_features(fsel, neuroim2::values(ROI), Y)
sum(featureMask) == 2

fsel2 <- feature_selector("FTest", "top_p", .1)
featureMask <- select_features(fsel2, neuroim2::values(ROI), Y)




cleanEx()
nameEx("set_log_level")
### * set_log_level

flush(stderr()); flush(stdout())

### Name: set_log_level
### Title: Set rMVPA Logging Level
### Aliases: set_log_level

### ** Examples

## Not run: 
##D   rMVPA::set_log_level("DEBUG")
##D   rMVPA::set_log_level("WARN")
## End(Not run)




cleanEx()
nameEx("strip_dataset-methods")
### * strip_dataset-methods

flush(stderr()); flush(stdout())

### Name: strip_dataset
### Title: Strip Dataset from Model Specification
### Aliases: strip_dataset strip_dataset.default strip_dataset.mvpa_model

### ** Examples




cleanEx()
nameEx("summarize_remap_items")
### * summarize_remap_items

flush(stderr()); flush(stdout())

### Name: summarize_remap_items
### Title: Summarize REMAP item-level residuals for an ROI
### Aliases: summarize_remap_items

### ** Examples

## Not run: 
##D res <- run_regional(model, region_mask, return_fits = TRUE)
##D items_tbl <- summarize_remap_items(res, roi = 1)
## End(Not run)



cleanEx()
nameEx("summarize_remap_roi")
### * summarize_remap_roi

flush(stderr()); flush(stdout())

### Name: summarize_remap_roi
### Title: Summarize REMAP diagnostics at the ROI level
### Aliases: summarize_remap_roi

### ** Examples

## Not run: 
##D res <- run_regional(model, region_mask, return_fits = TRUE)
##D summarize_remap_roi(res)
## End(Not run)



cleanEx()
nameEx("table_to_rdm")
### * table_to_rdm

flush(stderr()); flush(stdout())

### Name: table_to_rdm
### Title: Convert Similarity Table to RDM
### Aliases: table_to_rdm

### ** Examples

# Create a similarity table based on theoretical predictions
sim_table <- data.frame(
  label1 = c("cat", "cat", "dog", "car"),
  label2 = c("dog", "bird", "bird", "plane"),
  similarity = c(0.7, 0.5, 0.6, 0.8)
)

# Create RDM for specific conditions
conditions <- c("cat", "dog", "bird", "car", "plane", "boat")
rdm <- table_to_rdm(sim_table, conditions)

# Get as matrix instead of dist
rdm_matrix <- table_to_rdm(sim_table, conditions, as_dist = FALSE)

# Use in RSA design
## Not run: 
##D rsa_des <- rsa_design(~ theoretical_rdm,
##D                      data = list(theoretical_rdm = rdm))
## End(Not run)




cleanEx()
nameEx("temporal")
### * temporal

flush(stderr()); flush(stdout())

### Name: temporal
### Title: Temporal RDM wrapper for formula usage
### Aliases: temporal

### ** Examples

## Not run: 
##D # Use directly in RSA formula
##D rdes <- rsa_design(
##D   ~ task_rdm + temporal(trial_index, block=run, kernel="adjacent", width=2),
##D   data = list(task_rdm = task_rdm, trial_index = 1:100, run = run_ids),
##D   block_var = ~ run
##D )
## End(Not run)




cleanEx()
nameEx("temporal_nuisance_for_msreve")
### * temporal_nuisance_for_msreve

flush(stderr()); flush(stdout())

### Name: temporal_nuisance_for_msreve
### Title: Temporal nuisance RDM at condition level (for MS-ReVE)
### Aliases: temporal_nuisance_for_msreve

### ** Examples

## Not run: 
##D # Create temporal nuisance for MS-ReVE
##D temp_K <- temporal_nuisance_for_msreve(
##D   mvpa_design = mvpa_des,
##D   time_idx = seq_len(nrow(mvpa_des$train_design)),
##D   reduce = "min",
##D   kernel = "exp", 
##D   lambda = 3,
##D   within_blocks_only = TRUE
##D )
##D 
##D # Use in msreve_design
##D msreve_des <- msreve_design(
##D   mvpa_design = mvpa_des,
##D   contrast_matrix = C_mat,
##D   nuisance_rdms = list(temp_decay = temp_K)
##D )
## End(Not run)




cleanEx()
nameEx("temporal_rdm")
### * temporal_rdm

flush(stderr()); flush(stdout())

### Name: temporal_rdm
### Title: Temporal/ordinal nuisance RDM (trial-level)
### Aliases: temporal_rdm

### ** Examples

# Create temporal RDM for 20 trials across 4 runs
trial_index <- 1:20
run_labels <- rep(1:4, each = 5)

# Exponential decay kernel within runs only
temp_rdm <- temporal_rdm(trial_index, block = run_labels, 
                         kernel = "exp", lambda = 2,
                         within_blocks_only = TRUE)

# Use in RSA design
## Not run: 
##D rdes <- rsa_design(~ task_rdm + temporal_rdm(trial_idx, block=run, kernel="adjacent"),
##D                    data = list(task_rdm = my_task_rdm,
##D                               trial_idx = seq_len(n_trials),
##D                               run = run_ids),
##D                    block_var = ~ run,
##D                    keep_intra_run = TRUE)
## End(Not run)




cleanEx()
nameEx("test_design-methods")
### * test_design-methods

flush(stderr()); flush(stdout())

### Name: test_design
### Title: Test Design Extraction
### Aliases: test_design test_design.mvpa_design

### ** Examples

ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10, external_test = TRUE)
test_design(ds$design)



cleanEx()
nameEx("train_indices")
### * train_indices

flush(stderr()); flush(stdout())

### Name: train_indices
### Title: Get Training Indices for a Fold
### Aliases: train_indices

### ** Examples

cval <- kfold_cross_validation(len = 20, nfolds = 4)
train_indices(cval, 1)



cleanEx()
nameEx("train_model")
### * train_model

flush(stderr()); flush(stdout())

### Name: train_model
### Title: Train a classification, regression, or representational model.
### Aliases: train_model train_model.manova_model train_model.mvpa_model
###   train_model.rsa_model train_model.vector_rsa_model
###   train_model.feature_rsa_model train_model.contrast_rsa_model
### Keywords: internal

### ** Examples

# This example shows the structure of the returned list but doesn't actually run the function
# For a multi-metric model: output_metric = c("beta_delta", "recon_score", "beta_only")



cleanEx()
nameEx("transform_contrasts")
### * transform_contrasts

flush(stderr()); flush(stdout())

### Name: transform_contrasts
### Title: Apply Transformations to an Existing Contrast Matrix
### Aliases: transform_contrasts

### ** Examples

C_manual <- matrix(c( 1, -1,
                      1,  1,
                      0,  0,
                      0,  0), nrow = 4, byrow = TRUE,
                   dimnames = list(paste0("Cond", 1:4), c("MainA", "MainB")))

# Center and make orthonormal
C_transformed <- transform_contrasts(C_manual, orth = TRUE)
print(C_transformed)
print(attr(C_transformed, "source"))

# Center and scale to unit L2 norm
C_l2 <- transform_contrasts(C_manual, scale = "l2")
print(round(colSums(C_l2^2), 5))




cleanEx()
nameEx("tune_grid-methods")
### * tune_grid-methods

flush(stderr()); flush(stdout())

### Name: tune_grid
### Title: Extract Tuning Grid
### Aliases: tune_grid

### ** Examples




cleanEx()
nameEx("twofold_blocked_cross_validation")
### * twofold_blocked_cross_validation

flush(stderr()); flush(stdout())

### Name: twofold_blocked_cross_validation
### Title: twofold_blocked_cross_validation
### Aliases: twofold_blocked_cross_validation

### ** Examples

blockvar <- rep(1:5, each=10)
nreps <- 5
cval <- twofold_blocked_cross_validation(blockvar, nreps=nreps)
samples <- crossval_samples(cval, as.data.frame(matrix(rnorm(50*50),50,50)), y=rep(letters[1:5],10))
stopifnot(nrow(samples) == nreps)



cleanEx()
nameEx("y_test-methods")
### * y_test-methods

flush(stderr()); flush(stdout())

### Name: y_test
### Title: Test Labels/Response Extraction
### Aliases: y_test y_test.mvpa_design y_test.mvpa_model

### ** Examples

ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10, external_test = TRUE)
y_test(ds$design)



cleanEx()
nameEx("y_train-methods")
### * y_train-methods

flush(stderr()); flush(stdout())

### Name: y_train
### Title: Training Labels/Response Extraction
### Aliases: y_train y_train.mvpa_design y_train.mvpa_model
###   y_train.feature_rsa_model y_train.feature_rsa_design

### ** Examples

ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10)
y_train(ds$design)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
