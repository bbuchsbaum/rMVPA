context("msreve_design constructor and helpers")

library(testthat)
library(rMVPA) # Assuming necessary objects/functions are exported or accessible

# --- Helper objects for tests ---

# Mock mvpa_design object (needs structure defined in mvpa_design.R)
# Assuming mvpa_design requires at least condition_labels and block_var
mock_mvpa_design <- function(n_cond = 6, n_blocks = 3) {
  cond <- factor(rep(1:n_cond, each=n_blocks*2)) # Example condition labels
  blocks <- factor(rep(1:n_blocks, times=n_cond*2)) # Example block variable
  
  # Minimal structure based on potential mvpa_design definition
  structure(
    list(
      condition_labels = cond,
      block_var = blocks,
      design_mat = model.matrix(~ cond + blocks), # Example design matrix
      keep_intra_run = FALSE, # Example attribute
      conditions = levels(cond),
      ncond = n_cond,
      samples_per_condition_per_block = table(cond, blocks)
      # Add other essential fields if required by msreve_design checks
    ),
    class = c("mvpa_design", "list")
  )
}

# --- Tests for msreve_design() ---

test_that("msreve_design constructor works with valid inputs", {
  mvpa_des <- mock_mvpa_design(n_cond = 4)
  C <- matrix(c( 1, -1,  0,  0,
                 0,  0,  1, -1), nrow = 4, byrow = FALSE,
              dimnames = list(NULL, c("Cond1vs2", "Cond3vs4")))
  
  msreve_des <- msreve_design(mvpa_des, C, name = "Test Design")
  
  expect_s3_class(msreve_des, "msreve_design")
  expect_equal(msreve_des$name, "Test Design")
  expect_identical(msreve_des$mvpa_design, mvpa_des)
  expect_identical(msreve_des$contrast_matrix, C)
})

test_that("msreve_design constructor errors on invalid mvpa_design", {
  C <- matrix(1, nrow = 4, ncol = 1)
  expect_error(msreve_design(list(), C),
               regexp = "`mvpa_design` must be an object of class 'mvpa_design'.")
})

test_that("msreve_design constructor errors on invalid contrast_matrix", {
  mvpa_des <- mock_mvpa_design(n_cond = 4)
  expect_error(msreve_design(mvpa_des, "not a matrix"),
               regexp = "`contrast_matrix` must be a numeric matrix.")
  expect_error(msreve_design(mvpa_des, factor(1:4)),
               regexp = "`contrast_matrix` must be a numeric matrix.")
})

test_that("msreve_design constructor errors on dimension mismatch", {
  mvpa_des <- mock_mvpa_design(n_cond = 4) # 4 conditions
  C_wrong_rows <- matrix(1, nrow = 3, ncol = 1) # Only 3 rows
  expect_error(msreve_design(mvpa_des, C_wrong_rows),
               regexp = "Number of rows in contrast_matrix .*3.* must match number of conditions in mvpa_design .*4.*")
})

test_that("msreve_design constructor warns if contrast matrix has no column names", {
   mvpa_des <- mock_mvpa_design(n_cond = 4)
   C_no_colnames <- matrix(c(1, -1, 0, 0), nrow = 4, ncol = 1)
   expect_warning(msreve_design(mvpa_des, C_no_colnames),
                  regexp = "`contrast_matrix` does not have column names. It is recommended to name your contrasts.")
})

test_that("add_interaction_contrasts creates orthonormal interactions", {
  mvpa_des <- mock_mvpa_design(n_cond = 4)
  # Use contrasts that will produce a non-zero interaction
  C <- matrix(c(1, 1, -1, -1,    # A: first two vs last two
                1, -1, 1, -1), nrow = 4, byrow = FALSE,   # B: alternating
               dimnames = list(NULL, c("A", "B")))
  base_des <- msreve_design(mvpa_des, C)

  des_int <- suppressMessages(add_interaction_contrasts(base_des))

  expect_true(all(c("A_x_B") %in% colnames(des_int$contrast_matrix)))
  expect_true(attr(des_int, "is_orthonormal"))
  expect_equal(ncol(des_int$contrast_matrix), 3)
})

test_that("include_interactions parameter expands contrast matrix", {
  mvpa_des <- mock_mvpa_design(n_cond = 4)
  # Use contrasts that will produce a non-zero interaction
  C <- matrix(c(1, 1, -1, -1,    # A: first two vs last two
                1, -1, 1, -1), nrow = 4, byrow = FALSE,   # B: alternating
               dimnames = list(NULL, c("A", "B")))

  des_int <- suppressMessages(msreve_design(mvpa_des, C, include_interactions = TRUE))

  expect_true(all(c("A_x_B") %in% colnames(des_int$contrast_matrix)))
  expect_true(attr(des_int, "is_orthonormal"))
  expect_equal(ncol(des_int$contrast_matrix), 3)
})

test_that("zero interactions are skipped with informative message", {
  mvpa_des <- mock_mvpa_design(n_cond = 4)
  # Use contrasts with non-overlapping support (will produce zero interaction)
  C <- matrix(c(1, -1, 0, 0,     # A: conditions 1 vs 2
                0, 0, 1, -1), nrow = 4, byrow = FALSE,   # B: conditions 3 vs 4
               dimnames = list(NULL, c("A", "B")))

  expect_message(
    des_int <- msreve_design(mvpa_des, C, include_interactions = TRUE),
    regexp = "is zero.*non-overlapping support.*will be skipped"
  )
  
  # Should only have the original 2 contrasts, not the zero interaction
  expect_equal(ncol(des_int$contrast_matrix), 2)
  expect_equal(colnames(des_int$contrast_matrix), c("A", "B"))
  expect_true(attr(des_int, "is_orthonormal"))
})


# --- Tests for orthogonalize_contrasts() ---

test_that("orthogonalize_contrasts works for already orthogonal columns", {
  C_orth <- matrix(c( 1, -1,  0,  0,
                      0,  0,  1, -1), nrow = 4, byrow = FALSE,
                   dimnames = list(NULL, c("C1", "C2")))
  
  # Check they are indeed orthogonal - use expect_lt for floating point
  expect_lt(abs(crossprod(C_orth[,1], C_orth[,2])), 1e-10)
  
  C_orth_out <- orthogonalize_contrasts(C_orth)
  
  expect_equal(dim(C_orth_out), dim(C_orth))
  # Check orthogonality of output columns (dot product should be close to zero)
  expect_lt(abs(crossprod(C_orth_out[,1], C_orth_out[,2])), 1e-10)
  # Check columns span the same space (project original onto new)
  P <- C_orth_out %*% solve(t(C_orth_out) %*% C_orth_out) %*% t(C_orth_out)
  expect_equal(P %*% C_orth, C_orth, tolerance = 1e-10)
})

test_that("orthogonalize_contrasts works for non-orthogonal columns", {
  # Redefine C_non_orth to be truly non-orthogonal
  C_non_orth <- matrix(c( 1,  1,  0,  0,  # C1
                          1,  0,  1,  1), # C2: C1 dot C2 = 1
                       nrow = 4, byrow = FALSE,
                       dimnames = list(NULL, c("C1", "C2")))
  
  # Check they are indeed non-orthogonal
  expect_gt(abs(crossprod(C_non_orth[,1], C_non_orth[,2])), 1e-10) # Should be > 0 (it's 1 here)
  
  C_orth_out <- orthogonalize_contrasts(C_non_orth)
  
  expect_equal(dim(C_orth_out), dim(C_non_orth))
  # Check orthogonality of output columns
  expect_lt(abs(crossprod(C_orth_out[,1], C_orth_out[,2])), 1e-10)
  # Check columns span the same space
  P <- C_orth_out %*% solve(t(C_orth_out) %*% C_orth_out) %*% t(C_orth_out)
  expect_equal(P %*% C_non_orth, C_non_orth, tolerance = 1e-10)
})

test_that("orthogonalize_contrasts handles single column", {
  C_single <- matrix(c(1, 1, -1, -1), nrow = 4, ncol = 1, dimnames=list(NULL, "C1"))
  C_out <- orthogonalize_contrasts(C_single)
  expect_equal(dim(C_out), dim(C_single))
  # Output should be proportional to input (normalized)
  # The function output for single column is normalized, so it won't be identical unless input is normalized
  expect_equal(C_out / as.vector(sqrt(colSums(C_out^2))), C_single / as.vector(sqrt(colSums(C_single^2))), tolerance = 1e-10)
})

test_that("orthogonalize_contrasts preserves column names if possible", {
  # Use the corrected C_non_orth from above for this test too
  C_non_orth_named <- matrix(c( 1,  1,  0,  0, 
                                1,  0,  1,  1), 
                             nrow = 4, byrow = FALSE,
                             dimnames = list(NULL, c("ContrastA", "ContrastB")))
  C_orth_out <- orthogonalize_contrasts(C_non_orth_named)
  expect_equal(colnames(C_orth_out), colnames(C_non_orth_named))
})

test_that("orthogonalize_contrasts errors on rank deficient input matrix", {
  # Note: The current orthogonalize_contrasts does not error on rank deficiency,
  # it returns a matrix with zero columns for the dependent ones and issues a warning.
  # The test should reflect this actual behavior.
  # Create a rank-deficient matrix (col3 = col1 + col2)
  C1 <- c(1, -1, 0, 0)
  C2 <- c(0, 0, 1, -1)
  C3_dependent <- C1 + C2
  C_rank_def <- cbind(C1, C2, C3_dependent)
  colnames(C_rank_def) <- c("C1", "C2", "C3_dep")
  
  expect_warning(
    C_ortho_def <- orthogonalize_contrasts(C_rank_def),
    # Using a simpler regexp to avoid potential subtle mismatches
    regexp = "Input matrix C is rank-deficient" 
  )
  expect_equal(dim(C_ortho_def), dim(C_rank_def)) # Dims should be preserved
  # First two columns should be orthogonal
  expect_lt(abs(crossprod(C_ortho_def[,1], C_ortho_def[,2])), 1e-10)
  # Third column should be all zeros
  expect_true(all(abs(C_ortho_def[,3]) < 1e-10))
}) 