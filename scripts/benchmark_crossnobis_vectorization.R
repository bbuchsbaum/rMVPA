# Benchmark for compute_crossnobis_distances_sl
# compares the previous loop implementation with the new vectorised version

library(microbenchmark)

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(rMVPA)
}

compute_crossnobis <- get("compute_crossnobis_distances_sl",
                          envir = asNamespace("rMVPA"))

old_impl <- function(U_folds, P_voxels) {
  K <- dim(U_folds)[1]
  pair_idx <- which(lower.tri(matrix(TRUE, nrow = K, ncol = K)), arr.ind = TRUE)
  out <- numeric(nrow(pair_idx))
  for (p in seq_len(nrow(pair_idx))) {
    i <- pair_idx[p, 1]; j <- pair_idx[p, 2]
    delta <- U_folds[i, , ] - U_folds[j, , ]
    ip <- tcrossprod(t(delta))
    diag(ip) <- 0
    out[p] <- sum(ip) / (P_voxels * dim(U_folds)[3] * (dim(U_folds)[3] - 1))
  }
  names(out) <- paste0(dimnames(U_folds)[[1]][pair_idx[, 1]],
                       "_vs_",
                       dimnames(U_folds)[[1]][pair_idx[, 2]])
  out
}

set.seed(1)
K <- 40; V <- 100; M <- 8
U_folds <- array(rnorm(K * V * M), dim = c(K, V, M),
                 dimnames = list(paste0("C", 1:K), NULL, NULL))

stopifnot(isTRUE(all.equal(old_impl(U_folds, V),
                           compute_crossnobis(U_folds, P_voxels = V),
                           tolerance = 1e-12)))

microbenchmark(
  loop = old_impl(U_folds, V),
  gram = compute_crossnobis(U_folds, P_voxels = V),
  times = 10
)
