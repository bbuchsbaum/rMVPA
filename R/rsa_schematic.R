#' Pipeline schematic for standard RSA
#'
#' Generates a multi-panel diagram of the standard RSA workflow on a small
#' synthetic dataset: pattern matrix \eqn{\to} neural RDM \eqn{\to} model
#' RDMs \eqn{\to} lower-triangle fit. Used by the RSA vignette to give
#' newcomers a one-glance picture of the analysis. Pure base graphics; no
#' extra dependencies. Returns the underlying matrices invisibly so callers
#' can re-use them in follow-up chunks without recomputing.
#'
#' @param file Optional path; when supplied, opens a PNG device and writes
#'   the figure there. When \code{NULL} (default), draws to the current
#'   device. The function manages \code{layout()} and \code{par()} via
#'   \code{on.exit}.
#' @param seed Integer seed for the synthetic patterns. Same seed always
#'   produces the same schematic.
#' @return Invisibly, a list with components \code{patterns},
#'   \code{neural_rdm}, \code{model_rdms}, \code{betas}, \code{rho_cat}.
#' @keywords internal
#' @noRd
#' @export
rsa_schematic_basic <- function(file = NULL, seed = 7) {
  set.seed(seed)
  K <- 8
  V <- 30
  category <- factor(rep(c("face", "object"), each = K / 2))

  cat_axis <- ifelse(category == "face", 1, -1)
  id_axis  <- scale(seq_len(K))[, 1]
  W_cat <- matrix(stats::rnorm(V), 1, V)
  W_id  <- matrix(stats::rnorm(V), 1, V) * 0.6
  patterns <- cat_axis %*% W_cat + id_axis %*% W_id +
    matrix(stats::rnorm(K * V, sd = 0.6), K, V)
  rownames(patterns) <- paste0(substr(category, 1, 1), seq_len(K))

  neural_rdm <- as.matrix(stats::dist(patterns))
  model_rdm_cat <- as.matrix(stats::dist(as.numeric(category)))
  model_rdm_id  <- as.matrix(stats::dist(seq_len(K)))

  if (!is.null(file)) {
    grDevices::png(file, width = 1200, height = 720, res = 140)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::layout(matrix(c(1, 2, 3, 4,
                            5, 6, 7, 8), 2, 4, byrow = TRUE),
                   widths = c(1.05, 0.30, 1.0, 1.0))

  .heatmap_panel(patterns, "1. Patterns (conditions x voxels)",
                 xlab = "voxel", ylab = "condition",
                 palette = grDevices::hcl.colors(48, "BluYl"))
  .arrow_panel("dist()")
  .heatmap_panel(neural_rdm, "2. Neural RDM", mark_diag = TRUE)
  .arrow_panel("vec(lower tri)")

  .heatmap_panel(model_rdm_cat, "3a. Model RDM: category", mark_diag = TRUE,
                 palette = grDevices::hcl.colors(48, "Purp", rev = TRUE))
  .heatmap_panel(model_rdm_id, "3b. Model RDM: identity", mark_diag = TRUE,
                 palette = grDevices::hcl.colors(48, "Purp", rev = TRUE))

  ny <- .lower_tri(neural_rdm)
  mx <- .lower_tri(model_rdm_cat)
  graphics::par(mar = c(3.4, 3.4, 1.8, 0.8))
  plot(mx, ny, pch = 16, col = "#3a6ea5", cex = 0.6,
       xlab = "model RDM (lower tri)", ylab = "neural RDM (lower tri)",
       main = "4. Compare lower triangles")
  graphics::abline(stats::lm(ny ~ mx), col = "firebrick", lwd = 2)
  rho <- suppressWarnings(stats::cor(mx, ny, method = "spearman"))
  graphics::legend("topleft", bty = "n", cex = 0.85,
                   legend = sprintf("Spearman rho = %.2f", rho))

  X <- cbind(category = .lower_tri(model_rdm_cat),
             identity = .lower_tri(model_rdm_id))
  y <- .lower_tri(neural_rdm)
  fit <- stats::lm(y ~ X)
  betas <- stats::coef(fit)[-1]
  names(betas) <- c("cat", "id")
  graphics::par(mar = c(3.4, 3.4, 1.8, 0.8))
  graphics::barplot(betas, col = c("#5a3786", "#c89f3a"), border = NA,
                    main = "5. Multi-RDM regression", ylab = "beta")
  graphics::abline(h = 0, col = "grey40")

  invisible(list(
    patterns = patterns,
    neural_rdm = neural_rdm,
    model_rdms = list(category = model_rdm_cat, identity = model_rdm_id),
    betas = betas,
    rho_cat = rho
  ))
}

#' Pipeline schematic for contrast / MS-ReVE RSA
#'
#' Multi-panel diagram of the \code{contrast_rsa_model()} workflow:
#' contrast matrix \eqn{C} \eqn{\to} predicted dissimilarity
#' \eqn{\Delta_q = -c_q c_q^\top} per contrast \eqn{\to} regression of the
#' neural RDM on the stack \eqn{\to} signed voxel map
#' \eqn{\beta_q \Delta_{q,v}}. Used by the Contrast RSA vignette.
#'
#' @inheritParams rsa_schematic_basic
#' @return Invisibly, a list with components \code{C}, \code{Delta},
#'   \code{neural_rdm}, \code{beta}, \code{beta_delta}.
#' @keywords internal
#' @noRd
#' @export
rsa_schematic_msreve <- function(file = NULL, seed = 11) {
  set.seed(seed)
  K <- 4
  V <- 24
  cond <- factor(c("A", "B", "C", "D"))

  C_mat <- cbind(AB_vs_CD = c( 1,  1, -1, -1) / 2,
                 A_vs_B   = c( 1, -1,  0,  0) / 2)
  rownames(C_mat) <- as.character(cond)

  Wq <- C_mat[, 1, drop = FALSE] %*%
          matrix(c(stats::rnorm(12, 1, 0.3), stats::rnorm(12, 0, 0.3)), 1, V) +
        C_mat[, 2, drop = FALSE] %*%
          matrix(c(stats::rnorm(12, 0, 0.3), stats::rnorm(12, 1, 0.3)), 1, V)
  patterns <- Wq + matrix(stats::rnorm(K * V, sd = 0.4), K, V)
  rownames(patterns) <- as.character(cond)

  neural_rdm <- as.matrix(stats::dist(patterns)) ^ 2

  delta <- function(c_q) {
    M <- outer(c_q, c_q)
    diag(M) <- 0
    -M
  }
  Delta_list <- lapply(seq_len(ncol(C_mat)), function(q) delta(C_mat[, q]))
  names(Delta_list) <- colnames(C_mat)

  yv <- .lower_tri(neural_rdm)
  Xv <- vapply(Delta_list, .lower_tri, numeric(length(yv)))
  beta <- stats::coef(stats::lm(yv ~ Xv))[-1]
  names(beta) <- colnames(C_mat)

  pc <- scale(patterns, center = TRUE, scale = FALSE)
  proj <- t(pc) %*% C_mat
  beta_delta <- sweep(proj, 2, beta, `*`)

  if (!is.null(file)) {
    grDevices::png(file, width = 1200, height = 760, res = 140)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::layout(matrix(c(1, 2, 3, 4,
                            5, 6, 7, 8), 2, 4, byrow = TRUE),
                   widths = c(0.85, 0.30, 1.0, 1.0))

  .heatmap_panel(C_mat, "1. Contrast matrix C (K x Q)",
                 palette = grDevices::hcl.colors(48, "Blue-Red 2"))
  .arrow_panel("Delta_q = c_q c_q^T")
  .heatmap_panel(Delta_list[[1]], paste0("2a. Delta for ", names(Delta_list)[1]),
                 palette = grDevices::hcl.colors(48, "Purp", rev = TRUE),
                 mark_diag = TRUE)
  .heatmap_panel(Delta_list[[2]], paste0("2b. Delta for ", names(Delta_list)[2]),
                 palette = grDevices::hcl.colors(48, "Purp", rev = TRUE),
                 mark_diag = TRUE)

  .heatmap_panel(neural_rdm, "3. Neural RDM (searchlight)", mark_diag = TRUE)

  graphics::par(mar = c(3.4, 3.4, 1.8, 0.8))
  graphics::barplot(beta, col = c("#5a3786", "#c89f3a"), border = NA,
                    main = "4. beta_q  (RDM regression)", ylab = "beta")
  graphics::abline(h = 0, col = "grey40")

  .heatmap_panel(beta_delta, "5. beta_q * Delta_{q,v}",
                 xlab = "contrast", ylab = "voxel",
                 palette = grDevices::hcl.colors(48, "Blue-Red 2"))

  graphics::par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot.new()
  graphics::plot.window(c(0, 1), c(0, 1))
  graphics::text(0.02, 0.85, "MS-ReVE summary", adj = 0, cex = 1.05, font = 2)
  graphics::text(0.02, 0.66, "neural RDM ~ sum_q beta_q * Delta_q + e",
                 adj = 0, cex = 0.9)
  graphics::text(0.02, 0.50, "voxel score: beta_q * Delta_{q,v}",
                 adj = 0, cex = 0.9)
  graphics::text(0.02, 0.34, "+  voxel pulls the contrast",
                 adj = 0, cex = 0.85, col = "#a83232")
  graphics::text(0.02, 0.22, "-  voxel opposes it",
                 adj = 0, cex = 0.85, col = "#3a6ea5")

  invisible(list(
    C = C_mat,
    Delta = Delta_list,
    neural_rdm = neural_rdm,
    beta = beta,
    beta_delta = beta_delta
  ))
}

# --- internal panel helpers ------------------------------------------------

#' @keywords internal
#' @noRd
.heatmap_panel <- function(M, title, xlab = "", ylab = "",
                           palette = grDevices::hcl.colors(48, "YlOrBr", rev = TRUE),
                           zlim = NULL, mark_diag = FALSE) {
  if (is.null(zlim)) zlim <- range(M, finite = TRUE)
  graphics::par(mar = c(2.2, 2.5, 1.8, 0.6))
  graphics::image(seq_len(ncol(M)), seq_len(nrow(M)),
                  t(M[nrow(M):1, , drop = FALSE]),
                  col = palette, zlim = zlim, axes = FALSE,
                  xlab = xlab, ylab = ylab, main = title)
  graphics::box(col = "grey30")
  if (mark_diag) {
    n <- nrow(M)
    graphics::segments(0.5, n + 0.5, n + 0.5, 0.5, col = "grey60", lty = 3)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.arrow_panel <- function(label) {
  graphics::par(mar = c(0.4, 0.4, 0.4, 0.4))
  plot.new()
  graphics::plot.window(c(0, 1), c(0, 1))
  graphics::arrows(0.05, 0.5, 0.95, 0.5, length = 0.12, lwd = 2, col = "grey25")
  graphics::text(0.5, 0.78, label, cex = 0.95, col = "grey20")
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.lower_tri <- function(M) M[lower.tri(M)]
