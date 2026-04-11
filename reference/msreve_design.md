# Constructor for msreve_design

Creates an msreve_design object, which encapsulates the necessary design
information for a Multi-Dimensional Signed Representational Voxel
Encoding (MS-ReVE) analysis.

## Usage

``` r
msreve_design(
  mvpa_design,
  contrast_matrix,
  name = "msreve_design_01",
  include_interactions = FALSE,
  nuisance_rdms = NULL
)
```

## Arguments

- mvpa_design:

  An object of class `mvpa_design`, containing information about
  conditions, blocks, and cross-validation.

- contrast_matrix:

  A numeric matrix (`K x Q`) where `K` is the number of conditions and
  `Q` is the number of contrasts. Each column represents a contrast
  vector. It is highly recommended that columns are named to identify
  the contrasts.

- name:

  An optional character string to name the design.

- include_interactions:

  Logical. If TRUE, automatically add pairwise interaction contrasts
  using
  [`add_interaction_contrasts`](http://bbuchsbaum.github.io/rMVPA/reference/add_interaction_contrasts.md).

- nuisance_rdms:

  Optional named list of K x K matrices or `dist` objects representing
  nuisance RDMs to be included as additional predictors in the MS-ReVE
  regression. These are typically temporal or spatial nuisance patterns
  that should be accounted for but are not of primary interest.

## Value

An object of class `msreve_design`, which is a list containing:

- mvpa_design:

  The input `mvpa_design` object.

- contrast_matrix:

  The input `contrast_matrix`.

- name:

  The name of the design.

## Examples

``` r
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
```
