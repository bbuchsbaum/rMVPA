# Mock mvpa_design for cross-validation tests (shared by the contrast RSA
# and cv-means tests). Textually identical copies previously lived in
# test-contrast_rsa_additional.R, test-contrast_rsa_metrics.R and
# test-contrast_rsa_model.R.

mock_mvpa_design_cv <- function(n_samples = 24, n_cond = 4, n_blocks = 3, Y_labels=NULL) {
  if (is.null(Y_labels)) {
    cond_labels <- factor(rep(paste0("Cond", 1:n_cond), each = n_samples / n_cond))
  } else {
    cond_labels <- factor(Y_labels)
    n_cond <- length(levels(cond_labels))
  }
  block_var <- factor(rep(1:n_blocks, length.out = n_samples))
  conditions <- levels(cond_labels)
  design_matrix <- matrix(0, nrow = n_samples, ncol = 1)

  structure(
    list(
      Y = cond_labels,
      block_var = block_var,
      conditions = conditions,
      ncond = n_cond,
      design_matrix = design_matrix,
      samples = 1:n_samples
    ),
    class = c("mvpa_design", "list")
  )
}
