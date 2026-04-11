# Create an MVPA Model

Create an MVPA model based on a classification or regression model from
the MVPAModels registry.

## Usage

``` r
mvpa_model(
  model,
  dataset,
  design,
  model_type = c("classification", "regression"),
  crossval = NULL,
  feature_selector = NULL,
  tune_grid = NULL,
  tune_reps = 15,
  performance = NULL,
  class_metrics = FALSE,
  compute_performance = TRUE,
  return_predictions = TRUE,
  return_fits = FALSE
)
```

## Arguments

- model:

  A character string naming a model from the MVPAModels registry, or a
  custom model specification list.

- dataset:

  An \`mvpa_dataset\` instance.

- design:

  An \`mvpa_design\` instance.

- model_type:

  A character string indicating the problem type: "classification" or
  "regression".

- crossval:

  An optional \`cross_validation\` instance.

- feature_selector:

  An optional \`feature_selector\` instance.

- tune_grid:

  An optional parameter tuning grid as a \`data.frame\`.

- tune_reps:

  The number of replications used during parameter tuning. Only relevant
  if \`tune_grid\` is supplied.

- performance:

  An optional custom function for computing performance metrics.

- class_metrics:

  A logical flag indicating whether to compute performance metrics for
  each class.

- compute_performance:

  A `logical` indicating whether to compute and store performance
  measures for each voxel set (defaults to TRUE).

- return_predictions:

  A `logical` indicating whether to return row-wise predictions for each
  voxel set (defaults to TRUE).

- return_fits:

  A `logical` indicating whether to return the model fit for each voxel
  set (defaults to FALSE).

## Value

An \`mvpa_model\` object (a list with class 'mvpa_model').

## Details

If \`performance\` is supplied, it must be a function that takes one
argument and returns a named list of scalar values. The argument the
function takes is a class deriving from \`classification_result\`
appropriate for the problem at hand. See example below.

## Examples

``` r
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
print(mvpmod)
#> mvpa_model object. 
#> model:  sda 
#> model type:  classification 
#> 
#>  Blocked Cross-Validation 
#> 
#> - Dataset Information 
#>   - Observations:  100 
#>   - Number of Folds:  10 
#> - Block Information 
#>   - Total Blocks:  10 
#>   - Mean Block Size:  10  (SD:  0 ) 
#>   - Block Sizes:  1: 10, 2: 10, 3: 10, 4: 10, 5: 10, 6: 10, 7: 10, 8: 10, 9: 10, 10: 10 
#> 
#> 
#>  MVPA Dataset 
#> 
#> - Training Data 
#>   - Dimensions:  6 x 6 x 6 x 100 observations 
#>   - Type:  DenseNeuroVec 
#> - Test Data 
#>   -  None 
#> - Mask Information 
#>   - Areas:  TRUE : 118 
#>   - Active voxels/vertices:  118 
#> 
#> 
#>  MVPA Design 
#> 
#> - Training Data 
#>   - Observations:  100 
#>   - Response Type:  Factor
#>   - Levels:  a, b, c, d 
#>   - Class Distribution:  a: 25, b: 25, c: 25, d: 25 
#> - Test Data 
#>   -  None 
#> - Structure 
#>   - Blocking:  Present
#>   - Number of Blocks:  10 
#>   - Mean Block Size:  10  (SD:  0 ) 
#>   - Split Groups:  None 
#> 
```
