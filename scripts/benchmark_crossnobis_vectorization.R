# Benchmark for compute_crossnobis_distances_sl
# compares the previous loop implementation with the new vectorised version

library(microbenchmark)
library(rMVPA)

old_impl <- function(U_folds, P_voxels) {
  K <- dim(U_folds)[1]
  pair_idx <- utils::combn(K, 2)
  out <- numeric(ncol(pair_idx))
  for (p in seq_len(ncol(pair_idx))) {
    i <- pair_idx[1, p]; j <- pair_idx[2, p]
    delta <- U_folds[i, , ] - U_folds[j, , ]
    ip <- tcrossprod(t(delta))
    diag(ip) <- 0
    out[p] <- sum(ip) / (P_voxels * dim(U_folds)[3] * (dim(U_folds)[3] - 1))
  }
  out
}

set.seed(1)
K <- 40; V <- 100; M <- 8
U_folds <- array(rnorm(K * V * M), dim = c(K, V, M),
                 dimnames = list(paste0("C", 1:K), NULL, NULL))

microbenchmark(
  loop = old_impl(U_folds, V),
  vectorised = compute_crossnobis_distances_sl(U_folds, P_voxels = V),
  times = 10
)

