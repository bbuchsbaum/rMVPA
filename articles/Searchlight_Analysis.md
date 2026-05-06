# Searchlight Analysis

## Searchlight Analysis

### Simulated data

We first create a small volumetric dataset for demonstration. The call
below generates a 4D image with 6×6×6 spatial voxels and 80 observations
along the time dimension. Observations are grouped into 4 blocks (20
trials each), and the response factor `y` has two levels. The return
value contains both the data (`mvpa_dataset`) and the design
(`mvpa_design`).

``` r


dataset <- gen_sample_dataset(D=c(6,6,6), nobs = 80, blocks=4, nlevels=2)
print(dataset)
#> $dataset
#> 
#>  MVPA Dataset 
#> 
#> - Training Data 
#>   - Dimensions:  6 x 6 x 6 x 80 observations 
#>   - Type:  DenseNeuroVec 
#> - Test Data 
#>   -  None 
#> - Mask Information 
#>   - Areas:  TRUE : 216 
#>   - Active voxels/vertices:  216 
#> 
#> 
#> $design
#> 
#>  MVPA Design 
#> 
#> - Training Data 
#>   - Observations:  80 
#>   - Response Type:  Factor
#>   - Levels:  a, b 
#>   - Class Distribution:  a: 40, b: 40 
#> - Test Data 
#>   -  None 
#> - Structure 
#>   - Blocking:  Present
#>   - Number of Blocks:  4 
#>   - Mean Block Size:  20  (SD:  0 ) 
#>   - Split Groups:  None
```

### Cross‑validation

Because trials are organized into runs (blocks), we use a blocked
cross‑validation scheme: train on `k−1` blocks and test on the held‑out
block. This leave‑one‑group‑out strategy respects temporal correlations
within runs and yields a more realistic estimate of generalization.

``` r

block <- dataset$design$block_var
crossval <- blocked_cross_validation(block)
crossval
#> 
#>  Blocked Cross-Validation 
#> 
#> - Dataset Information 
#>   - Observations:  80 
#>   - Number of Folds:  4 
#> - Block Information 
#>   - Total Blocks:  4 
#>   - Mean Block Size:  20  (SD:  0 ) 
#>   - Block Sizes:  1: 20, 2: 20, 3: 20, 4: 20
```

### Model specification

We use the shrinkage discriminant analysis model (`sda_notune`), which
estimates its shrinkage parameter from the training folds. See the [sda
package](https://cran.r-project.org/web/packages/sda/index.html) for
details.

``` r

sda_model <- load_model("sda_notune") 
model <- mvpa_model(model=sda_model, dataset=dataset$dataset, design=dataset$design, crossval=crossval)
model
#> mvpa_model object. 
#> model:  sda_notune 
#> model type:  classification 
#> 
#>  Blocked Cross-Validation 
#> 
#> - Dataset Information 
#>   - Observations:  80 
#>   - Number of Folds:  4 
#> - Block Information 
#>   - Total Blocks:  4 
#>   - Mean Block Size:  20  (SD:  0 ) 
#>   - Block Sizes:  1: 20, 2: 20, 3: 20, 4: 20 
#> 
#> 
#>  MVPA Dataset 
#> 
#> - Training Data 
#>   - Dimensions:  6 x 6 x 6 x 80 observations 
#>   - Type:  DenseNeuroVec 
#> - Test Data 
#>   -  None 
#> - Mask Information 
#>   - Areas:  TRUE : 216 
#>   - Active voxels/vertices:  216 
#> 
#> 
#>  MVPA Design 
#> 
#> - Training Data 
#>   - Observations:  80 
#>   - Response Type:  Factor
#>   - Levels:  a, b 
#>   - Class Distribution:  a: 40, b: 40 
#> - Test Data 
#>   -  None 
#> - Structure 
#>   - Blocking:  Present
#>   - Number of Blocks:  4 
#>   - Mean Block Size:  20  (SD:  0 ) 
#>   - Split Groups:  None
```

### Standard searchlight

[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
returns image volumes with performance metrics at each sphere center.
For two‑class problems we report cross‑validated accuracy and AUC
(centered at 0 by subtracting 0.5). The `radius` is in millimeters;
`method = "standard"` evaluates one model per sphere centered on each
voxel.

``` r

result <- run_searchlight(model, radius=4, method="standard")
result
```

## Randomized searchlight

The randomized variant samples non‑overlapping spheres and propagates
each sphere’s score to all voxels it covers. Repeating this for `niter`
exhaustive passes yields, at each voxel, the average performance across
the spheres that included it. This can improve spatial localization and
often reduces computation, since the number of models scales with
`nvoxels / radius × niter` rather than the total number of voxels.

``` r

result <- run_searchlight(model, radius=4, method="randomized", niter=8)
result
```

## Using different classifiers

Any classifier in the rMVPA registry can be used in a searchlight. You
can also register your own with
[`register_mvpa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/register_mvpa_model.md).
Two robust options that ship with rMVPA are `hdrda` and `pca_lda`:

``` r

# High‑Dimensional Regularized Discriminant Analysis
hdrda <- load_model("hdrda")
model_hdrda <- mvpa_model(model = hdrda, dataset = dataset$dataset, design = dataset$design,
                          crossval = crossval)
res_hdrda <- run_searchlight(model_hdrda, radius = 4, method = "randomized", niter = 2)
res_hdrda
```

``` r

# PCA + LDA pipeline
pca_lda <- load_model("pca_lda")
model_pca_lda <- mvpa_model(model = pca_lda, dataset = dataset$dataset, design = dataset$design,
                            crossval = crossval)
res_pca_lda <- run_searchlight(model_pca_lda, radius = 4, method = "randomized", niter = 2)
res_pca_lda
```

## Domain‑adaptive cross‑decoding (REMAP‑RRR)

REMAP‑RRR models memory as perception plus a low‑rank correction in a
jointly whitened space. It learns a residual map Δ on the prototype
residuals (Y_w − X_w) ≈ X_w Δ, builds predicted memory templates Ŷ_w =
X_w + λ X_w Δ (with λ tuned on training items), and classifies held‑out
memory trials by correlation to Ŷ_w. Provide an external test set in
your dataset/design (same pattern as localizer→WM), and optionally a
`link_by` column shared by train/test designs to pair items across
domains.

``` r

# Synthetic source→target setup
toy <- gen_sample_dataset(D = c(6,6,6), nobs = 120, nlevels = 4, blocks = 3, external_test = TRUE)
regionMask <- neuroim2::NeuroVol(sample(1:5, size = length(toy$dataset$mask), replace = TRUE),
                                 neuroim2::space(toy$dataset$mask))

mspec <- remap_rrr_model(
  dataset = toy$dataset,
  design  = toy$design,
  link_by = NULL,            # defaults to class-wise pairing
  rank    = 0,               # identity fallback (no rrpack required)
  leave_one_key_out = TRUE
)
remap_res <- run_regional(mspec, regionMask)
names(remap_res$vol_results)  # includes remap_improv, delta_frob_mean, lambda_mean
```

Set `rank = "auto"` or a positive integer to enable reduced‑rank mapping
via `rrpack`.

## See also

- [`vignette("Regional_Analysis")`](http://bbuchsbaum.github.io/rMVPA/articles/Regional_Analysis.md)
  – the ROI-based counterpart to searchlight analysis
- [`vignette("CrossValidation")`](http://bbuchsbaum.github.io/rMVPA/articles/CrossValidation.md)
  – cross-validation strategies for fMRI
- [`vignette("CustomAnalyses")`](http://bbuchsbaum.github.io/rMVPA/articles/CustomAnalyses.md)
  – plug in your own analysis functions
