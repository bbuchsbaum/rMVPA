# Run from the repository root with:
# Rscript --vanilla inst/benchmarks/pattern_model/validate_confirmation.R
# Independent error columns are repeated null experiments, not brain voxels
# used to claim replication. Shared block signs retain each family's dependence.
devtools::load_all(quiet = TRUE)
setup_lm <- getFromNamespace(".pattern_lm_setup", "rMVPA")
mass_lm <- getFromNamespace(".pattern_mass_lm", "rMVPA")
wild <- getFromNamespace(".pattern_wild_test", "rMVPA")
plan <- confirmation_plan("sign_flip", n_resamples = 199, seed = 810)
results <- list()
elapsed <- system.time(for (seed in 801:840) {
  set.seed(seed)
  n <- 240; g <- 30; p <- 25
  block <- rep(seq_len(g), each = n/g)
  Tm <- matrix(rnorm(n*2), n, 2)
  # A block-level disturbance shared across its rows, but independent across
  # experiments/columns. Scores vary both within and between blocks.
  Tm <- Tm + matrix(rnorm(g*2), g, 2)[block, ]
  E <- matrix(rnorm(n*p), n, p) + matrix(rnorm(g*p), g, p)[block, ]
  s <- setup_lm(cbind(Tm, 1), 2, as.character(block))
  observed <- mass_lm(E, s)
  boot <- wild(E, s, observed, plan)
  results[[length(results) + 1L]] <- data.frame(seed = seed,
    sandwich_p = observed$p_omnibus, wild_p = boot$p_omnibus,
    wild_max_p = boot$p_omnibus_max)
})
null <- do.call(rbind, results)
print(data.frame(method = c("CR1 Wald F", "restricted block wild"),
                 experiments = nrow(null),
                 mean_p = c(mean(null$sandwich_p), mean(null$wild_p)),
                 rejection_05 = c(mean(null$sandwich_p < .05), mean(null$wild_p < .05))))
cat("Global-null family rejection fraction:", mean(vapply(results, function(x) any(x$wild_max_p < .05), logical(1))), "\n")
cat("Null simulation elapsed seconds:", elapsed[["elapsed"]], "\n")
stopifnot(abs(mean(null$sandwich_p < .05) - .05) < .04,
          abs(mean(null$wild_p < .05) - .05) < .04,
          abs(mean(null$wild_p) - .5) < .06)

set.seed(850)
n <- 400; p <- 2000; block <- rep(1:40, each = 10)
Tm <- matrix(rnorm(n*2), n, 2)
Y <- matrix(rnorm(n*p), n, p)
s <- setup_lm(cbind(Tm, 1), 2, as.character(block))
timing <- system.time(out <- mass_lm(Y, s))
cat("CR1 workload:", n, "rows,", p, "features, rank 2, 40 blocks\n")
cat("Elapsed seconds:", timing[["elapsed"]], "; retained fit summary bytes:", as.numeric(object.size(out)), "\n")
cat("No observation-by-observation or feature-by-feature covariance is formed.\n")
print(sessionInfo())
