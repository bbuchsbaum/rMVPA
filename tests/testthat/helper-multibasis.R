testthat::skip_if_not_installed("neuroim2")
library(neuroim2)

make_multibasis_fixture <- function(D = c(3, 3, 3), n_events = 4, k = 3) {
  sp3 <- neuroim2::NeuroSpace(D)
  sp4 <- neuroim2::NeuroSpace(c(D, n_events))
  mask <- neuroim2::NeuroVol(array(1L, D), sp3)

  basis_arrays <- lapply(seq_len(k), function(b) {
    arr <- array(0, c(D, n_events))
    for (e in seq_len(n_events)) {
      arr[, , , e] <- b * 100 + e
    }
    arr
  })

  basis_vecs <- lapply(basis_arrays, function(arr) neuroim2::NeuroVec(arr, sp4))

  list(
    D = D,
    n_events = n_events,
    k = k,
    mask = mask,
    basis_arrays = basis_arrays,
    basis_vecs = basis_vecs
  )
}

make_concat_vec <- function(fx, ordering = c("event_major", "basis_major")) {
  ordering <- match.arg(ordering)
  n_events <- fx$n_events
  k <- fx$k
  concat_arr <- array(0, c(fx$D, n_events * k))

  if (ordering == "event_major") {
    idx <- 1L
    for (e in seq_len(n_events)) {
      for (b in seq_len(k)) {
        concat_arr[, , , idx] <- fx$basis_arrays[[b]][, , , e]
        idx <- idx + 1L
      }
    }
  } else {
    idx <- 1L
    for (b in seq_len(k)) {
      for (e in seq_len(n_events)) {
        concat_arr[, , , idx] <- fx$basis_arrays[[b]][, , , e]
        idx <- idx + 1L
      }
    }
  }

  neuroim2::NeuroVec(concat_arr, neuroim2::NeuroSpace(c(fx$D, n_events * k)))
}

write_basis_files <- function(vec_list) {
  files <- vapply(seq_along(vec_list), function(i) {
    f <- tempfile(pattern = sprintf("basis_%02d_", i), fileext = ".nii.gz")
    neuroim2::write_vec(vec_list[[i]], f)
    f
  }, character(1))
  files
}

gen_multibasis_sample_dataset <- function(D = c(4, 4, 4), n_events = 20, k = 2,
                                           n_classes = 2, n_blocks = 4) {
  fx <- make_multibasis_fixture(D = D, n_events = n_events, k = k)
  y <- factor(rep(letters[seq_len(n_classes)], length.out = n_events))
  blocks <- rep(seq_len(n_blocks), length.out = n_events)
  design <- rMVPA::mvpa_design(
    data.frame(y = y, block = blocks),
    y_train = ~y, block_var = ~block
  )
  dset <- rMVPA::mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  cv <- rMVPA::blocked_cross_validation(design$block_var)
  list(dataset = dset, design = design, crossval = cv, fixture = fx)
}

# Like make_multibasis_fixture but with random (spatially varying) data.
# Needed for models that require spatial variance (contrast_rsa, vector_rsa, etc.)
make_random_multibasis_fixture <- function(D = c(4, 4, 4), n_events = 20, k = 2) {
  sp3 <- neuroim2::NeuroSpace(D)
  sp4 <- neuroim2::NeuroSpace(c(D, n_events))
  mask <- neuroim2::NeuroVol(array(1L, D), sp3)

  basis_arrays <- lapply(seq_len(k), function(b) {
    array(rnorm(prod(D) * n_events), c(D, n_events))
  })

  basis_vecs <- lapply(basis_arrays, function(arr) neuroim2::NeuroVec(arr, sp4))

  list(D = D, n_events = n_events, k = k, mask = mask,
       basis_arrays = basis_arrays, basis_vecs = basis_vecs)
}

gen_random_multibasis_sample_dataset <- function(D = c(4, 4, 4), n_events = 20, k = 2,
                                                  n_classes = 2, n_blocks = 4) {
  fx <- make_random_multibasis_fixture(D = D, n_events = n_events, k = k)
  # Use a properly crossed design: all conditions in every block
  reps_per_block <- ceiling(n_events / n_blocks)
  blocks <- rep(seq_len(n_blocks), each = reps_per_block)[seq_len(n_events)]
  y <- factor(rep(letters[seq_len(n_classes)], length.out = n_events))
  design <- rMVPA::mvpa_design(
    data.frame(y = y, block = blocks),
    y_train = ~y, block_var = ~block
  )
  dset <- rMVPA::mvpa_multibasis_dataset(train_data = fx$basis_vecs, mask = fx$mask)
  cv <- rMVPA::blocked_cross_validation(design$block_var)
  list(dataset = dset, design = design, crossval = cv, fixture = fx)
}
