# Temporal/ordinal nuisance RDM (trial-level)

Creates a temporal or ordinal nuisance representational dissimilarity
matrix (RDM) for use in RSA analyses. This function generates various
kernels to model temporal proximity effects in neuroimaging data and can
return either similarity- or distance-like outputs.

## Usage

``` r
temporal_rdm(
  index,
  block = NULL,
  kernel = c("adjacent", "boxcar", "linear", "poly", "exp", "gauss"),
  width = 1L,
  power = 1,
  lambda = NULL,
  sigma = NULL,
  within_blocks_only = TRUE,
  wrap = FALSE,
  units = c("auto", "sec", "TR", "index"),
  TR = NULL,
  metric = c("similarity", "distance"),
  normalize = c("rank", "z", "none"),
  as_dist = TRUE
)
```

## Arguments

- index:

  numeric or integer vector representing trial order or time (length N
  observations)

- block:

  optional vector (length N) of run/block identifiers

- kernel:

  character string specifying the kernel type, one of:

  adjacent

  :   Binary kernel for immediate neighbors within specified width

  boxcar

  :   Binary kernel including diagonal up to specified width

  linear

  :   Linear distance based on lag

  poly

  :   Polynomial distance with specified power

  exp

  :   Exponential similarity with specified lambda

  gauss

  :   Gaussian similarity with specified sigma

- width:

  integer window for "adjacent"/"boxcar" kernels (default 1)

- power:

  exponent for "poly" kernel (default 1)

- lambda:

  decay constant for "exp" kernel (in units; default chosen
  automatically if NULL)

- sigma:

  standard deviation for "gauss" kernel (in units; default chosen
  automatically if NULL)

- within_blocks_only:

  logical; if TRUE, zero nuisance across blocks (default TRUE)

- wrap:

  logical; if TRUE, treat index as circular (default FALSE)

- units:

  one of `"auto"`, `"sec"`, `"TR"`, or `"index"`; used to choose
  sensible defaults for decay parameters

- TR:

  repetition time in seconds (optional, but recommended when
  `units="sec"` or `units="TR"`)

- metric:

  return type: `"distance"` (increasing with lag) or `"similarity"`
  (decreasing with lag). Default `"distance"`.

- normalize:

  character string specifying normalization method:

  rank

  :   Rank transform (ties averaged)

  z

  :   Z-score normalization

  none

  :   No normalization

- as_dist:

  logical; if TRUE return a `dist` object, otherwise return matrix
  (default TRUE)

## Value

A `dist` object or symmetric matrix (N x N) with 0 on the diagonal,
representing temporal/ordinal relationships between observations

## Details

This function creates temporal nuisance RDMs for modeling carry-over
effects, scanner drift, or other temporal confounds in fMRI data. The
resulting RDM can be included as a nuisance regressor in RSA models to
account for temporal proximity effects while preserving statistical
power. By default, it returns a distance-like quantity that increases
with temporal separation (`metric="distance"`).

## Examples

``` r
# Create temporal RDM for 20 trials across 4 runs
trial_index <- 1:20
run_labels <- rep(1:4, each = 5)

# Exponential decay kernel within runs only
temp_rdm <- temporal_rdm(trial_index, block = run_labels, 
                         kernel = "exp", lambda = 2,
                         within_blocks_only = TRUE)

# Use in RSA design
if (FALSE) { # \dontrun{
rdes <- rsa_design(~ task_rdm + temporal_rdm(trial_idx, block=run, kernel="adjacent"),
                   data = list(task_rdm = my_task_rdm,
                              trial_idx = seq_len(n_trials),
                              run = run_ids),
                   block_var = ~ run,
                   keep_intra_run = TRUE)
} # }
```
