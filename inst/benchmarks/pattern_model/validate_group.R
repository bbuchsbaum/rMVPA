# Repeated sufficient-statistic experiments with known sampling and population
# covariances. These deliberately do not reuse the confirmation regression.
# Run from the repository root: Rscript --vanilla inst/benchmarks/pattern_model/validate_group.R
devtools::load_all(quiet = TRUE)
set.seed(860)
s <- 24; p <- 2000; r <- 2
within <- matrix(c(.05, .02, .02, .08), r)
between <- diag(c(.10, .15))
basis <- structure(list(matrix = diag(r), target_ids = c("y1", "y2"),
                        type = "continuous", basis_id = "simulation"), class = "pattern_basis")
subjects <- lapply(seq_len(s), function(j) {
  estimates <- matrix(rnorm(p*r), p, r) %*% chol(within + between)
  structure(list(estimate = estimates, covariance = NULL,
    covariance_factor = within, residual_variance = rep(1, p),
    basis = basis, feature_ids = paste0("v", seq_len(p)),
    nuisance = matrix(1, 10, 1, dimnames = list(NULL, "intercept")),
    provenance = list(subject_id = paste0("s", j),
      observation_ids = paste0("s", j, ":confirm:", 1:10),
      discovery_ids = paste0("s", j, ":discovery:", 1:10),
      preprocessing_id = "simulation-original-units")), class = "pattern_confirmation")
})
elapsed <- system.time(group <- pattern_group(subjects, basis))
print(data.frame(subjects = s, rank = r, null_experiments = p,
                 mean_p = mean(group$omnibus$p),
                 rejection_05 = mean(group$omnibus$p < .05),
                 mean_heterogeneity_trace = mean(group$heterogeneity$trace),
                 true_heterogeneity_trace = sum(diag(between))))
cat("Elapsed seconds:", elapsed[["elapsed"]], "; result bytes:", as.numeric(object.size(group)), "\n")
stopifnot(abs(mean(group$omnibus$p) - .5) < .035,
          abs(mean(group$omnibus$p < .05) - .05) < .025,
          abs(mean(group$heterogeneity$trace) - sum(diag(between))) < .025)
# Heterogeneous within-subject precision is an approximate-inference stress
# case, rather than another claim of exact Hotelling calibration.
set.seed(861)
for (j in seq_len(s)) {
  V <- within * (j/s * 4)
  subjects[[j]]$estimate <- matrix(rnorm(p*r), p, r) %*% chol(V + between)
  subjects[[j]]$covariance_factor <- V
}
group_heterogeneous <- pattern_group(subjects, basis)
cat("Heterogeneous precision: mean p =", mean(group_heterogeneous$omnibus$p),
    "; rejection at .05 =", mean(group_heterogeneous$omnibus$p < .05), "\n")
print(sessionInfo())
