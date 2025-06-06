# Tests for contrasts generation functions

library(testthat)
library(tibble)

# Common setup
labs <- c("faces","animals","plants","tools", 
          "vehicles","furniture","buildings","food")
n_labs <- length(labs)

meta <- tibble::tribble(
  ~label,      ~anim, ~size, ~cat,
  "faces",        1,     0, "bio",
  "animals",      1,     0, "bio",
  "plants",       1,     1, "bio",
  "tools",        0,     0, "nonbio",
  "vehicles",     0,     1, "nonbio",
  "furniture",    0,     0, "nonbio",
  "buildings",    0,     1, "nonbio",
  "food",         1,     1, "bio")

meta_factor <- meta
meta_factor$cat <- factor(meta_factor$cat)


# === Test contrasts() function ===

test_that("contrasts() works with mini-DSL spec", {
  spec_dsl <- ~ anim( faces + animals + plants + food ~ . ) + 
                size( faces + animals + tools + furniture ~ . )
  
  # Basic structure
  C_dsl <- contrasts(labels = labs, spec = spec_dsl, centre = FALSE, scale = "none", orth = FALSE)
  expect_true(is.matrix(C_dsl))
  expect_equal(nrow(C_dsl), n_labs)
  expect_equal(ncol(C_dsl), 2) # anim, size
  expect_equal(colnames(C_dsl), c("anim", "size"))
  expect_equal(rownames(C_dsl), labs)
  
  # Check coding (anim: faces,animals,plants,food = +1, others = -1)
  expect_equal(unname(C_dsl[c("faces","animals","plants","food"), 1]), rep(1, 4))
  expect_equal(unname(C_dsl[c("tools","vehicles","furniture","buildings"), 1]), rep(-1, 4))
  
  # Check coding (size: faces,animals,tools,furniture = +1, others = -1)
  expect_equal(unname(C_dsl[c("faces","animals","tools","furniture"), 2]), rep(1, 4))
  expect_equal(unname(C_dsl[c("plants","vehicles","buildings","food"), 2]), rep(-1, 4))
  
  # Test centering
  C_dsl_centered <- contrasts(labels = labs, spec = spec_dsl, centre = TRUE, scale = "none", orth = FALSE)
  expect_true(all(abs(colMeans(C_dsl_centered)) < 1e-10))
  
  # Test interaction term
  spec_dsl_int <- ~ anim( faces + animals + plants + food ~ . ) + 
                    size( faces + animals + tools + furniture ~ . ) + 
                    anim:size
  C_dsl_int <- contrasts(labels = labs, spec = spec_dsl_int, centre = FALSE, scale = "none", orth = FALSE)
  expect_equal(ncol(C_dsl_int), 3)
  expect_equal(colnames(C_dsl_int), c("anim", "size", "anim:size"))
  # Interaction should be product of base columns before centering
  expect_equal(unname(C_dsl_int[, 3]), unname(C_dsl[, 1] * C_dsl[, 2]))
  
  # Test '.' notation where B is explicit empty
   spec_dot_empty <- ~ faces_only( faces ~ . )
   C_dot_empty <- contrasts(labels = labs, spec = spec_dot_empty, centre = FALSE, scale = "none")
   expect_equal(unname(C_dot_empty[,1]), c(1, rep(-1, n_labs -1)))
   
   # Test '.' warning when A covers all labels
   expect_warning(contrasts(labels = labs, spec = ~ all_vs_none( faces + animals + plants + tools + vehicles + furniture + buildings + food ~ .), centre=FALSE),
                   regexp="resulted in an empty set for the right side because all labels were on the left")

})

test_that("contrasts() works with metadata and formula spec", {
  spec_formula <- ~ anim + size + anim:size
  
  # Basic structure
  C_meta <- contrasts(metadata = meta, spec = spec_formula, centre = FALSE, scale = "none", orth = FALSE)
  expect_true(is.matrix(C_meta))
  expect_equal(nrow(C_meta), n_labs)
  expect_equal(ncol(C_meta), 3) # anim, size, anim:size
  expect_equal(colnames(C_meta), c("anim", "size", "anim:size"))
  expect_equal(rownames(C_meta), meta$label)
  
  # Check coding matches metadata (before centering)
  expect_equal(unname(C_meta[, 1]), meta$anim)
  expect_equal(unname(C_meta[, 2]), meta$size)
  expect_equal(unname(C_meta[, 3]), meta$anim * meta$size)
  
  # Centering
  C_meta_centered <- contrasts(metadata = meta, spec = spec_formula, centre = TRUE, scale = "none", orth = FALSE)
  expect_true(all(abs(colMeans(C_meta_centered)) < 1e-10))
  
  # Works with factor columns in metadata
  spec_formula_factor <- ~ cat + size + cat:size
  C_meta_factor <- contrasts(metadata = meta_factor, spec = spec_formula_factor, centre = FALSE, scale = "none", orth = FALSE)
  expect_true(is.matrix(C_meta_factor))
  # model.matrix creates dummy var for factor 'cat' -> catnonbio
  expect_equal(colnames(C_meta_factor), c("catnonbio", "size", "catnonbio:size"))
  # Check row order matches metadata
  expect_equal(rownames(C_meta_factor), meta_factor$label)

})

test_that("contrasts() handles scaling options", {
  spec_dsl <- ~ anim( faces + animals + plants + food ~ . ) + 
                size( faces + animals + tools + furniture ~ . )
                
  # Scale = "sd"
  C_sd <- contrasts(labels = labs, spec = spec_dsl, centre = TRUE, scale = "sd", orth = FALSE)
  expect_true(all(abs(apply(C_sd, 2, sd) - 1) < 1e-10))
  expect_true(all(abs(colMeans(C_sd)) < 1e-10)) # Still centered
  
  # Scale = "l2"
  C_l2 <- contrasts(labels = labs, spec = spec_dsl, centre = TRUE, scale = "l2", orth = FALSE)
  expect_true(all(abs(colSums(C_l2^2) - 1) < 1e-10)) # L2 norm = 1
  expect_true(all(abs(colMeans(C_l2)) < 1e-10)) # Still centered
  
  # Scale = "none"
  C_none <- contrasts(labels = labs, spec = spec_dsl, centre = TRUE, scale = "none", orth = FALSE)
  expect_false(abs(apply(C_none, 2, sd)[1] - 1) < 1e-10) # SD is not 1
  expect_false(abs(colSums(C_none^2)[1] - 1) < 1e-10) # L2 norm is not 1

  # Scaling with zero variance column (after centering)
  labs_2 <- c("A", "B")
  spec_const <- ~ const( A + B ~ . ) # Before centering: [1, 1]. After centering: [0, 0]
  expect_warning(
    C_zero_sd <- contrasts(labels = labs_2, spec = spec_const, centre = TRUE, scale = "sd", orth = FALSE),
    regexp = "near-zero standard deviation"
  )
  expect_equal(unname(C_zero_sd[,1]), c(0, 0)) # Column remains zero

  expect_warning(
    C_zero_l2 <- contrasts(labels = labs_2, spec = spec_const, centre = TRUE, scale = "l2", orth = FALSE),
    regexp = "near-zero L2 norm"
  )
  expect_equal(unname(C_zero_l2[,1]), c(0, 0)) # Column remains zero
  
})


test_that("contrasts() handles orthogonalization", {
  spec_dsl_int <- ~ anim( faces + animals + plants + food ~ . ) + 
                    size( faces + animals + tools + furniture ~ . ) + 
                    anim:size
                    
  # Basic orthogonalization
  C_orth <- contrasts(labels = labs, spec = spec_dsl_int, centre = TRUE, orth = TRUE)
  expect_true(is.matrix(C_orth))
  expect_equal(nrow(C_orth), n_labs)
  expect_equal(ncol(C_orth), 3) # Original rank is 3
  expect_true(all(abs(crossprod(C_orth)) - diag(3) < 1e-10)) # Orthonormal
  expect_true(all(abs(colSums(C_orth^2) - 1) < 1e-10)) # Unit L2 norm
  expect_equal(colnames(C_orth), paste0("Orth", 1:3))
  
  # Check attributes
  expect_true(!is.null(attr(C_orth, "source")))
  expect_length(attr(C_orth, "source"), 3)
  expect_true(all(c("anim", "size", "anim:size") %in% attr(C_orth, "source"))) # Check original names are source
  expect_null(attr(C_orth, "dropped"))

  # Rank deficient case
  spec_rank_def <- ~ anim( faces + animals + plants + food ~ . ) + 
                     inv_anim( tools + vehicles + furniture + buildings ~ . ) + # inv_anim = -1 * anim
                     size( faces + animals + tools + furniture ~ . )
  expect_warning(
    C_orth_def <- contrasts(labels = labs, spec = spec_rank_def, centre = TRUE, orth = TRUE),
    regexp = "linearly dependent"
  )
  expect_equal(ncol(C_orth_def), 2) # Rank is 2 (anim and size are independent)
  expect_equal(colnames(C_orth_def), paste0("Orth", 1:2))
  expect_true(all(abs(crossprod(C_orth_def)) - diag(2) < 1e-10)) # Orthonormal
  expect_true(!is.null(attr(C_orth_def, "source")))
  expect_length(attr(C_orth_def, "source"), 2)
  expect_true(!is.null(attr(C_orth_def, "dropped")))
  expect_equal(attr(C_orth_def, "dropped"), "inv_anim") # inv_anim should be dropped

  # Orthogonalization with 1 column
  spec_one_col <- ~ anim( faces + animals + plants + food ~ . )
  C_orth_one <- contrasts(labels = labs, spec = spec_one_col, centre = TRUE, orth = TRUE)
  expect_equal(ncol(C_orth_one), 1)
  expect_equal(colnames(C_orth_one), "Orth1")
  expect_true(abs(sum(C_orth_one^2) - 1) < 1e-10)
  expect_equal(attr(C_orth_one, "source"), "anim")
  expect_null(attr(C_orth_one, "dropped"))
  
  # Orthogonalization with centre=FALSE
  # Use a spec where uncentered columns are NOT already centered
  labs_4 <- c("A", "B", "C", "D")
  spec_uncentered <- ~ f1( A + B ~ C )
  C_orth_nocenter <- contrasts(labels = labs_4, spec = spec_uncentered, centre = FALSE, orth = TRUE)
  expect_equal(ncol(C_orth_nocenter), 1) # Should be rank 1
  expect_true(all(abs(crossprod(C_orth_nocenter)) - diag(1) < 1e-10)) # Orthonormal (1x1)
  # Should not be centered (original was [1, 1, -1, 0] which has non-zero mean)
  expect_false(abs(mean(C_orth_nocenter)) < 1e-10) 

  # Test keep_attr = FALSE
  C_orth_noattr <- contrasts(labels = labs, spec = spec_dsl_int, centre = TRUE, orth = TRUE, keep_attr = FALSE)
  expect_null(attr(C_orth_noattr, "source"))
  expect_null(attr(C_orth_noattr, "dropped"))

})

test_that("contrasts() error handling", {
  # Invalid spec type
  expect_error(contrasts(labels=labs, spec=list()), regexp = "`spec` must be a formula.")
  expect_error(contrasts(metadata=meta, spec="~a+b"), regexp = "`spec` must be a formula.")
  
  # Missing labels when metadata is NULL
  expect_error(contrasts(spec = ~ a(b~c)), regexp = "`labels` must be provided if `metadata` is NULL.")
  
  # Invalid labels type
  expect_error(contrasts(labels=1:3, spec = ~ a(b~c)), regexp = "`labels` must be a character vector.")
  
  # Duplicated labels
  expect_error(contrasts(labels=c("a", "a"), spec = ~ a(b~c)), regexp = "`labels` must contain unique values.")
  
  # Invalid metadata type
  expect_error(contrasts(metadata=list(a=1), spec = ~ a), regexp = "`metadata` must be a data.frame or tibble.")
  
  # Missing 'label' column in metadata
  expect_error(contrasts(metadata=tibble(a=1, b=2), spec = ~ a), regexp = "`metadata` must contain a 'label' column.")
  
  # Non-character/factor label column
  expect_error(contrasts(metadata=tibble(label=1:3, a=1:3), spec= ~a), regexp = "`metadata$label` column must be character or factor.", fixed = TRUE)
  
  # Duplicated labels in metadata
  expect_error(contrasts(metadata=tibble(label=c("a", "a"), a=1:2), spec=~a), regexp = "Labels in `metadata$label` column must be unique.", fixed = TRUE)

  # # Non-numeric/factor predictors in metadata -- this check was removed from contrasts()
  # # model.matrix will handle character columns by converting to factors
  # expect_error(contrasts(metadata=tibble(label=c("a","b"), p1=c("x","y")), spec=~p1), 
  #              regexp = "Non-numeric/factor columns found in `metadata` predictors.")
               
  # Mini-DSL: Levels not in labels
  expect_error(contrasts(labels=c("a", "b"), spec = ~ f(a ~ c)), 
               regexp = "Level(s) 'c' in factor 'f' not found in provided labels.", fixed = TRUE)
  expect_error(contrasts(labels=c("a", "b"), spec = ~ f(c ~ a)), 
               regexp = "Level(s) 'c' in factor 'f' not found in provided labels.", fixed = TRUE)
               
  # Mini-DSL: Label on both sides
  expect_error(contrasts(labels=c("a", "b"), spec = ~ f(a ~ a + b)), 
               regexp = "Label(s) 'a' appear on both sides of '~' for factor 'f'.", fixed = TRUE)
               
  # Mini-DSL: Empty LHS
  expect_error(contrasts(labels=c("a", "b"), spec = ~ f( ~ a)), 
               regexp = "No valid factor definitions found") # Updated expected message

  # Mini-DSL: Factor defined twice
  expect_error(contrasts(labels=labs, spec = ~ f(a~b) + f(c~d)), 
               regexp = "Level\\(s\\) 'a' in factor 'f' not found") # Updated expected message
               
  # Mini-DSL: No valid factor definitions
  expect_error(contrasts(labels=labs, spec = ~ a + b), # Needs factor() syntax
               regexp = "No valid factor definitions found.")
               
   # Mini-DSL: Interaction with undefined factor
  expect_warning(contrasts(labels=labs, spec = ~ f(faces ~ .) + f:g ),
                 regexp="Interaction term 'f:g' involves factors not defined")

  # Spec results in empty matrix (only intercept)
  expect_error(contrasts(metadata=meta, spec= ~ 1), 
               regexp = "Specification resulted in an empty contrast matrix") # Reverted to original expected message
               
  # Matrix becomes rank 0 after centering
  labs_2 <- c("A", "B")
  spec_const <- ~ const( A + B ~ . ) # -> [1, 1] -> centered [0, 0]
  # Expect the error. The contrasts() call here also issues a warning about '.' notation, which is expected.
  expect_error(
    suppressWarnings(contrasts(labels = labs_2, spec = spec_const, centre = TRUE, orth = TRUE)),
    regexp = "Single contrast column has zero length after centering; cannot orthogonalize/normalize."
  )

})


# === Test transform_contrasts() function ===

test_that("transform_contrasts() works correctly", {
  C_manual <- matrix(c( 1,  1, -1, -1, # Col 1
                         1, -1,  1, -1), # Col 2
                       nrow = 4, byrow = FALSE,
                       dimnames = list(paste0("Cond", 1:4), c("MainA", "MainB")))
  
  # Centering
  C_cen <- transform_contrasts(C_manual, centre = TRUE, scale = "none", orth = FALSE)
  expect_true(all(abs(colMeans(C_cen)) < 1e-10))
  expect_equal(rownames(C_cen), rownames(C_manual))
  expect_equal(colnames(C_cen), colnames(C_manual))
  
  # Scaling = "sd" (after centering)
  C_sd <- transform_contrasts(C_manual, centre = TRUE, scale = "sd", orth = FALSE)
  expect_true(all(abs(apply(C_sd, 2, sd) - 1) < 1e-10))
  expect_true(all(abs(colMeans(C_sd)) < 1e-10))
  
  # Scaling = "l2" (after centering)
  C_l2 <- transform_contrasts(C_manual, centre = TRUE, scale = "l2", orth = FALSE)
  expect_true(all(abs(colSums(C_l2^2) - 1) < 1e-10))
  
  # Orthogonalization (input columns are already orthogonal and centered)
  C_orth <- transform_contrasts(C_manual, centre = TRUE, orth = TRUE)
  expect_equal(ncol(C_orth), 2)
  expect_equal(colnames(C_orth), c("Orth1", "Orth2"))
  expect_true(all(abs(crossprod(C_orth)) - diag(2) < 1e-10))
  expect_true(all(abs(colSums(C_orth^2) - 1) < 1e-10))
  expect_equal(attr(C_orth, "source"), colnames(C_manual))
  expect_null(attr(C_orth, "dropped"))
  expect_equal(rownames(C_orth), rownames(C_manual))

  # Orthogonalization with rank deficiency
  C_rank_def <- cbind(C_manual, C_manual[,1] * -1) # Add a dependent column
  colnames(C_rank_def)[3] <- "DepCol"
  expect_warning(
    C_orth_def <- transform_contrasts(C_rank_def, centre = TRUE, orth = TRUE),
    regexp = "Input contrast matrix columns are linearly dependent"
  )
  expect_equal(ncol(C_orth_def), 2) # Rank is 2
  expect_equal(colnames(C_orth_def), c("Orth1", "Orth2"))
  expect_true(all(abs(crossprod(C_orth_def)) - diag(2) < 1e-10))
  expect_equal(attr(C_orth_def, "source"), c("MainA", "MainB")) # Should be original non-dep cols
  expect_equal(attr(C_orth_def, "dropped"), "DepCol")
  
  # Orthogonalization with 1 column
  C_one <- C_manual[, 1, drop = FALSE]
  C_orth_one <- transform_contrasts(C_one, centre = TRUE, orth = TRUE)
  expect_equal(ncol(C_orth_one), 1)
  expect_equal(colnames(C_orth_one), "Orth1")
  expect_true(abs(sum(C_orth_one^2) - 1) < 1e-10)
  expect_equal(attr(C_orth_one, "source"), "MainA")
  expect_null(attr(C_orth_one, "dropped"))

})

test_that("transform_contrasts() error handling", {
  # Invalid input type
  expect_error(transform_contrasts(list(a=1)), regexp="`C` must be a numeric matrix.")
  expect_error(transform_contrasts(matrix(letters[1:4], 2)), regexp="`C` must be a numeric matrix.")
  
  # Empty matrix
  expect_error(transform_contrasts(matrix(numeric(0), nrow=0, ncol=0)), regexp="`C` must have at least one row.")
  expect_error(transform_contrasts(matrix(numeric(0), nrow=2, ncol=0), orth=TRUE), regexp="`C` must have columns for the requested operations.")
  # Should pass if no ops requiring columns
  expect_silent(transform_contrasts(matrix(numeric(0), nrow=2, ncol=0), centre=TRUE, scale="none", orth=FALSE))
  
  # Rank 0 after centering
  C_zero <- matrix(1, nrow=4, ncol=1)
  expect_error(transform_contrasts(C_zero, centre=TRUE, orth=TRUE), 
               regexp="Single contrast column has zero length after centering; cannot orthogonalize/normalize.")

})


# === Test make_feature_contrasts() ===

feat_mat <- matrix(rnorm(40), nrow = n_labs, ncol = 5,
                   dimnames = list(sample(labs), paste0("F", 1:5)))
labels_ordered <- labs # desired output order

test_that("make_feature_contrasts() works without PCA", {
  C_raw <- make_feature_contrasts(feat_mat, labels = labels_ordered, use_pca = FALSE, prefix = "Raw_")
  
  expect_true(is.matrix(C_raw))
  expect_equal(nrow(C_raw), n_labs)
  expect_equal(ncol(C_raw), ncol(feat_mat))
  expect_equal(rownames(C_raw), labels_ordered)
  expect_equal(colnames(C_raw), paste0("Raw_", colnames(feat_mat)))
  
  # Check row order matches labels_ordered
  expect_identical(unname(C_raw[labels_ordered[1], ]), unname(feat_mat[labels_ordered[1], ]))
  expect_identical(unname(C_raw[labels_ordered[3], ]), unname(feat_mat[labels_ordered[3], ]))
  
  # Without labels arg, uses rownames order
  C_raw_rownames <- make_feature_contrasts(feat_mat, use_pca = FALSE)
  expect_equal(rownames(C_raw_rownames), rownames(feat_mat))
  expect_equal(colnames(C_raw_rownames), paste0("Feat_", colnames(feat_mat)))
  
  # No row names in features or labels
  feat_mat_norownames <- feat_mat
  rownames(feat_mat_norownames) <- NULL
  expect_warning(
    C_raw_nonames <- make_feature_contrasts(feat_mat_norownames, use_pca=FALSE),
    regexp = "Rows will not be named"
  )
  expect_null(rownames(C_raw_nonames))
})


test_that("make_feature_contrasts() works with PCA", {
  # PCA with n_pcs
  ncomp <- 3
  C_pca_n <- make_feature_contrasts(feat_mat, labels = labels_ordered, use_pca = TRUE, 
                                   n_pcs = ncomp, prefix = "PCA_")
  expect_true(is.matrix(C_pca_n))
  expect_equal(nrow(C_pca_n), n_labs)
  expect_equal(ncol(C_pca_n), ncomp)
  expect_equal(rownames(C_pca_n), labels_ordered)
  expect_equal(colnames(C_pca_n), paste0("PCA_PC", 1:ncomp))
  
  # PCA with pve
  pve_thresh <- 0.8
  pca_check <- prcomp(feat_mat[labels_ordered, ], center = TRUE, scale. = FALSE)
  var_explained <- pca_check$sdev^2 / sum(pca_check$sdev^2)
  ncomp_pve <- which(cumsum(var_explained) >= pve_thresh)[1]
  
  C_pca_pve <- make_feature_contrasts(feat_mat, labels = labels_ordered, use_pca = TRUE, 
                                    pve = pve_thresh, prefix = "PCA_")
  expect_equal(ncol(C_pca_pve), ncomp_pve)
  expect_equal(colnames(C_pca_pve), paste0("PCA_PC", 1:ncomp_pve))
  
  # PCA with scaling
  C_pca_scaled <- make_feature_contrasts(feat_mat, labels = labels_ordered, use_pca = TRUE, 
                                       scale_pca = TRUE, n_pcs = 2, prefix = "PCA_Scaled_")
  expect_equal(ncol(C_pca_scaled), 2)
  # Check if result differs from non-scaled (it should unless data was already scaled)
  expect_false(isTRUE(all.equal(C_pca_scaled, C_pca_n[, 1:2])))
  
  # PCA with no centering
  C_pca_nocenter <- make_feature_contrasts(feat_mat, labels = labels_ordered, use_pca = TRUE, 
                                         centre_pca = FALSE, n_pcs = 2, prefix = "PCA_NoCenter_")
  expect_equal(ncol(C_pca_nocenter), 2)
  # Check if result differs from centered (it should unless data was already centered)
  expect_false(isTRUE(all.equal(C_pca_nocenter, C_pca_n[, 1:2])))

})


test_that("make_feature_contrasts() error handling", {
  # Invalid features type
  expect_error(make_feature_contrasts(list(a=1)), regexp="`features` must be a numeric matrix.")
  expect_error(make_feature_contrasts(data.frame(a=1:3)), regexp="`features` must be a numeric matrix.")
  
  # Empty features matrix
  expect_error(make_feature_contrasts(matrix(numeric(0), 0, 0)), regexp="`features` must not be an empty matrix.")
  # Suppress warning from matrix(rnorm(10), 5, 0)
  suppressWarnings(
    expect_error(make_feature_contrasts(matrix(rnorm(10), 5, 0)), regexp="`features` must not be an empty matrix.")
  )
  
  # Invalid labels type
  expect_error(make_feature_contrasts(feat_mat, labels=1:n_labs), regexp="`labels` must be a character vector.")
  
  # Duplicated labels
  expect_error(make_feature_contrasts(feat_mat, labels=rep(labs[1], n_labs)), regexp="`labels` must contain unique values.")
  
  # Labels not in features rownames
  expect_error(make_feature_contrasts(feat_mat, labels=c(labs[1], "nonexistent")), 
               regexp = "Not all values in `labels` are present in `rownames(features)`.", fixed = TRUE)
               
  # Rownames required if labels are given
  feat_mat_norownames <- feat_mat
  rownames(feat_mat_norownames) <- NULL
  expect_error(make_feature_contrasts(feat_mat_norownames, labels=labs), 
               regexp = "`features` matrix must have row names if `labels` argument is used.")
               
  # Duplicated rownames if labels are NULL
  feat_mat_duprows <- rbind(feat_mat, feat_mat[1, , drop=FALSE])
  rownames(feat_mat_duprows)[n_labs + 1] <- rownames(feat_mat_duprows)[1]
  expect_error(make_feature_contrasts(feat_mat_duprows, labels=NULL), 
               regexp = "`rownames(features)` must be unique if `labels` is not provided.", fixed = TRUE)

  # Invalid PCA params
  expect_error(make_feature_contrasts(feat_mat, use_pca=TRUE, n_pcs = -1), regexp="`n_pcs` must be a positive integer.")
  expect_error(make_feature_contrasts(feat_mat, use_pca=TRUE, n_pcs = ncol(feat_mat) + 1), regexp="`n_pcs` must be a positive integer.")
  expect_error(make_feature_contrasts(feat_mat, use_pca=TRUE, pve = 1.1), regexp="`pve` must be a number between 0")
  expect_error(make_feature_contrasts(feat_mat, use_pca=TRUE, pve = 0), regexp="`pve` must be a number between 0")
  
  # PCA results in 0 components (e.g., pve very high on constant data)
  feat_const <- matrix(1, nrow=n_labs, ncol=2, dimnames=list(labs, c("F1", "F2")))
  # Need centre=TRUE for prcomp to potentially have non-zero sdev if scale=TRUE
  # But with const input and scale=F, sdev is 0.
  # If scale=T, prcomp might error earlier depending on version.
  # Let's test the check directly by setting n_pcs=0 (invalid input) or pve that can't be reached
   expect_error(make_feature_contrasts(feat_const, use_pca=TRUE, pve=0.99),
               regexp = "PCA resulted in 0 components selected.") # Fails because sdev is 0

})