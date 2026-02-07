## Shared helper for synthetic mediation simulations used across tests
simulate_repmed_latent <- function(K = 24,
                                   a = 1,
                                   b = 1,
                                   c = 0,
                                   sd_m = 0.04,
                                   sd_y = 0.04) {

  items <- paste0("it", seq_len(K))

  ## Latent coordinates for predictor, mediator, outcome
  z_x <- rnorm(K)
  z_m <- a * z_x + rnorm(K, sd = sd_m)
  z_y <- b * z_m + c * z_x + rnorm(K, sd = sd_y)

  X_rdm <- as.matrix(dist(z_x, method = "euclidean"))
  Y_rdm <- as.matrix(dist(z_y, method = "euclidean"))
  rownames(X_rdm) <- colnames(X_rdm) <- items
  rownames(Y_rdm) <- colnames(Y_rdm) <- items

  list(
    items = items,
    z_x   = z_x,
    z_m   = z_m,
    z_y   = z_y,
    X_rdm = X_rdm,
    Y_rdm = Y_rdm
  )
}
