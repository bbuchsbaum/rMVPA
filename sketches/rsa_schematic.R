# RSA schematic prototype (exploratory)
# ----------------------------------------------------------------------------
# Generates a multi-panel figure that walks the reader through the RSA
# pipeline using a small synthetic dataset. Two variants:
#
#   rsa_schematic_basic()    -- standard RSA (model RDM ↔ neural RDM)
#   rsa_schematic_msreve()   -- contrast / MS-ReVE (Cm Cm^T ↔ neural RDM, beta_q)
#
# Pure base graphics, no extra dependencies. Run interactively:
#
#   source("sketches/rsa_schematic.R")
#   rsa_schematic_basic("sketches/rsa_basic.png")
#   rsa_schematic_msreve("sketches/rsa_msreve.png")
#
# Status: prototype. See the discussion in the related branch / PR for
# whether to fold this into a vignette as a generated figure.
# ----------------------------------------------------------------------------

# --- helpers ---------------------------------------------------------------

.heatmap_panel <- function(M, title, xlab = "", ylab = "",
                           palette = hcl.colors(48, "YlOrBr", rev = TRUE),
                           zlim = NULL, mark_diag = FALSE) {
  if (is.null(zlim)) zlim <- range(M, finite = TRUE)
  op <- par(mar = c(2.2, 2.5, 1.8, 0.6))
  on.exit(par(op), add = TRUE)
  image(seq_len(ncol(M)), seq_len(nrow(M)), t(M[nrow(M):1, , drop = FALSE]),
        col = palette, zlim = zlim, axes = FALSE, xlab = xlab, ylab = ylab,
        main = title)
  box(col = "grey30")
  if (mark_diag) {
    n <- nrow(M)
    segments(0.5, n + 0.5, n + 0.5, 0.5, col = "grey60", lty = 3)
  }
}

.arrow_panel <- function(label) {
  op <- par(mar = c(0.4, 0.4, 0.4, 0.4))
  on.exit(par(op), add = TRUE)
  plot.new(); plot.window(c(0, 1), c(0, 1))
  arrows(0.05, 0.5, 0.95, 0.5, length = 0.12, lwd = 2, col = "grey25")
  text(0.5, 0.78, label, cex = 0.95, col = "grey20")
}

.lower_tri <- function(M) M[lower.tri(M)]

# --- standard RSA schematic -------------------------------------------------

rsa_schematic_basic <- function(file = NULL, seed = 7) {
  set.seed(seed)
  K <- 8                                  # conditions
  V <- 30                                 # voxels in toy ROI
  category <- factor(rep(c("face", "object"), each = K / 2))

  # Latent neural patterns: category structure + identity + noise.
  cat_axis <- ifelse(category == "face", 1, -1)
  id_axis  <- scale(seq_len(K))[, 1]
  W_cat <- matrix(rnorm(V), 1, V)
  W_id  <- matrix(rnorm(V), 1, V) * 0.6
  patterns <- cat_axis %*% W_cat + id_axis %*% W_id + matrix(rnorm(K * V, sd = 0.6), K, V)
  rownames(patterns) <- paste0(substr(category, 1, 1), seq_len(K))

  neural_rdm <- as.matrix(dist(patterns))
  model_rdm_cat <- as.matrix(dist(as.numeric(category)))
  model_rdm_id  <- as.matrix(dist(seq_len(K)))

  # Layout: 2 rows. Row 1 = pattern matrix → neural RDM. Row 2 = model RDMs → fit.
  if (!is.null(file)) png(file, width = 1200, height = 720, res = 140)
  on.exit(if (!is.null(file)) dev.off(), add = TRUE)

  layout(matrix(c(1, 2, 3, 4,
                  5, 6, 7, 8), 2, 4, byrow = TRUE),
         widths = c(1.05, 0.30, 1.0, 1.0))

  # Row 1
  .heatmap_panel(patterns, "1. Patterns (conditions × voxels)",
                 xlab = "voxel", ylab = "condition",
                 palette = hcl.colors(48, "BluYl"))
  .arrow_panel("dist()")
  .heatmap_panel(neural_rdm, "2. Neural RDM",
                 mark_diag = TRUE)
  .arrow_panel("vec(lower △)")

  # Row 2
  .heatmap_panel(model_rdm_cat, "3a. Model RDM: category", mark_diag = TRUE,
                 palette = hcl.colors(48, "Purp", rev = TRUE))
  .heatmap_panel(model_rdm_id,  "3b. Model RDM: identity", mark_diag = TRUE,
                 palette = hcl.colors(48, "Purp", rev = TRUE))

  # Scatter
  ny <- .lower_tri(neural_rdm)
  mx <- .lower_tri(model_rdm_cat)
  op <- par(mar = c(3.4, 3.4, 1.8, 0.8))
  plot(mx, ny, pch = 16, col = "#3a6ea5", cex = 0.6,
       xlab = "model RDM (lower △)", ylab = "neural RDM (lower △)",
       main = "4. Compare lower triangles")
  abline(lm(ny ~ mx), col = "firebrick", lwd = 2)
  rho <- suppressWarnings(cor(mx, ny, method = "spearman"))
  legend("topleft", bty = "n", cex = 0.85,
         legend = sprintf("Spearman ρ = %.2f", rho))
  par(op)

  # Bar of multi-RDM regression
  X <- cbind(category = .lower_tri(model_rdm_cat),
             identity = .lower_tri(model_rdm_id))
  y <- .lower_tri(neural_rdm)
  fit <- lm(y ~ X)
  betas <- coef(fit)[-1]
  op <- par(mar = c(3.4, 3.4, 1.8, 0.8))
  bp <- barplot(betas, col = c("#5a3786", "#c89f3a"), border = NA,
                main = "5. Multi-RDM regression",
                ylab = "beta", names.arg = c("cat", "id"))
  abline(h = 0, col = "grey40")
  par(op)

  invisible(list(patterns = patterns, neural_rdm = neural_rdm,
                 model_rdms = list(category = model_rdm_cat, identity = model_rdm_id),
                 betas = betas, rho_cat = rho))
}

# --- MS-ReVE / contrast_rsa schematic ---------------------------------------

rsa_schematic_msreve <- function(file = NULL, seed = 11) {
  set.seed(seed)
  K <- 4                                       # conditions
  V <- 24                                      # voxels
  cond <- factor(c("A", "B", "C", "D"))

  # Two centred contrasts: AB_vs_CD and A_vs_B.
  C_mat <- cbind(AB_vs_CD = c( 1,  1, -1, -1) / 2,
                 A_vs_B   = c( 1, -1,  0,  0) / 2)
  rownames(C_mat) <- as.character(cond)

  # Each contrast has its own per-voxel "encoding" — we want the AB_vs_CD
  # axis to be carried by voxels 1-12 and the A_vs_B axis by voxels 13-24.
  Wq <- matrix(0, K, V)
  Wq <- C_mat[, 1, drop = FALSE] %*% matrix(c(rnorm(12, 1, 0.3), rnorm(12, 0, 0.3)), 1, V) +
        C_mat[, 2, drop = FALSE] %*% matrix(c(rnorm(12, 0, 0.3), rnorm(12, 1, 0.3)), 1, V)
  patterns <- Wq + matrix(rnorm(K * V, sd = 0.4), K, V)
  rownames(patterns) <- as.character(cond)

  neural_rdm <- as.matrix(dist(patterns))^2     # squared distances are conventional here

  # Predicted RDMs per contrast: Δ_q = c_q c_q^T (sign-flipped to dissimilarity)
  delta <- function(c_q) {
    M <- outer(c_q, c_q)
    diag(M) <- 0
    -M                                          # high cc^T => low dissim
  }
  Delta_list <- lapply(seq_len(ncol(C_mat)), function(q) delta(C_mat[, q]))
  names(Delta_list) <- colnames(C_mat)

  # Fit Δ_neural ≈ Σ β_q Δ_q on the lower triangle.
  yv <- .lower_tri(neural_rdm)
  Xv <- vapply(Delta_list, .lower_tri, numeric(length(yv)))
  beta <- coef(lm(yv ~ Xv))[-1]
  names(beta) <- colnames(C_mat)

  # Voxel projections Δ_{q,v} = pattern^T %*% c_q  (centered patterns).
  pc <- scale(patterns, center = TRUE, scale = FALSE)
  proj <- t(pc) %*% C_mat                         # V × Q
  beta_delta <- sweep(proj, 2, beta, `*`)         # V × Q signed contributions

  if (!is.null(file)) png(file, width = 1200, height = 760, res = 140)
  on.exit(if (!is.null(file)) dev.off(), add = TRUE)

  layout(matrix(c(1, 2, 3, 4,
                  5, 6, 7, 8), 2, 4, byrow = TRUE),
         widths = c(0.85, 0.30, 1.0, 1.0))

  # Row 1: contrast matrix → predicted RDMs
  .heatmap_panel(C_mat, "1. Contrast matrix C (K × Q)",
                 palette = hcl.colors(48, "Blue-Red 2"))
  .arrow_panel("Δ_q = c_q c_qᵀ")
  .heatmap_panel(Delta_list[[1]], paste0("2a. Δ for ", names(Delta_list)[1]),
                 palette = hcl.colors(48, "Purp", rev = TRUE), mark_diag = TRUE)
  .heatmap_panel(Delta_list[[2]], paste0("2b. Δ for ", names(Delta_list)[2]),
                 palette = hcl.colors(48, "Purp", rev = TRUE), mark_diag = TRUE)

  # Row 2: neural RDM, beta bar, beta_delta voxel map
  .heatmap_panel(neural_rdm, "3. Neural RDM (searchlight)",
                 mark_diag = TRUE)

  op <- par(mar = c(3.4, 3.4, 1.8, 0.8))
  bp <- barplot(beta, col = c("#5a3786", "#c89f3a"), border = NA,
                main = "4. β_q  (RDM regression)", ylab = "beta")
  abline(h = 0, col = "grey40")
  par(op)

  # Voxel map: rows = voxels, cols = contrasts
  .heatmap_panel(beta_delta, "5. β_q Δ_{q,v}  (signed voxel map)",
                 xlab = "contrast", ylab = "voxel",
                 palette = hcl.colors(48, "Blue-Red 2"))

  # Legend / equation panel
  op <- par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot.new(); plot.window(c(0, 1), c(0, 1))
  text(0.02, 0.85, "MS-ReVE summary", adj = 0, cex = 1.05, font = 2)
  text(0.02, 0.66, "Δ̂(neural)  =  Σ_q  β_q · Δ_q + ε", adj = 0, cex = 0.95)
  text(0.02, 0.50, "voxel score: β_q · Δ_{q,v}", adj = 0, cex = 0.95)
  text(0.02, 0.34, "+  voxel pulls the contrast", adj = 0, cex = 0.85, col = "#a83232")
  text(0.02, 0.22, "−  voxel opposes it",          adj = 0, cex = 0.85, col = "#3a6ea5")
  par(op)

  invisible(list(C = C_mat, Delta = Delta_list, neural_rdm = neural_rdm,
                 beta = beta, beta_delta = beta_delta))
}

# --- if run as script -------------------------------------------------------
if (sys.nframe() == 0) {
  out_dir <- "sketches"
  dir.create(out_dir, showWarnings = FALSE)
  rsa_schematic_basic(file.path(out_dir, "rsa_basic.png"))
  rsa_schematic_msreve(file.path(out_dir, "rsa_msreve.png"))
  message("wrote ", file.path(out_dir, "rsa_basic.png"),
          " and ", file.path(out_dir, "rsa_msreve.png"))
}
