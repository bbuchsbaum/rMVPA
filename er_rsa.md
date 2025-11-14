# ERA–RSA Notes (Discussion 1)

This note captures, as faithfully as possible, the proposed “Encoding–Retrieval Analysis (ERA) + Representational Geometry” approach and how it would slot into rMVPA. It summarizes the core math, intended API, per‑ROI computation, and extensible ideas discussed.

## Conceptual Overview
- First‑order (ERA): classic encoding–retrieval similarity, item‑wise E–R similarity in a naive cross‑decoding style.
- Second‑order (Geometry): similarity of representational geometries — correlation between encoding RDM and retrieval RDM, aligned by item key.
- Joint view: analyze the first‑order × second‑order space to understand whether strong ERA arises with preserved geometry, despite geometry warping, etc.

## Setup and Notation
Assume a single subject and one ROI (searchlight or region).
- Encoding trials: matrix `X_E ∈ R^{N_E × V}`
- Retrieval trials: matrix `X_R ∈ R^{N_R × V}`
- Each trial has an item key `k ∈ {1,…,K}` (e.g., `image_id`).
- We use `link_by` to align encoding and retrieval by item key (e.g., `link_by = "image_id"`).

### Encoding prototypes and retrieval patterns
Per item k:
- Encoding prototype: `e_k = mean{x_t^(E) : trial t has key k}`
- Retrieval pattern: either a single trial `r_k` (usual case) or a prototype if multiple retrieval trials exist.

Let `E ∈ R^{K × V}` be the matrix of encoding prototypes (rows are `e_k`), and `R ∈ R^{K × V}` the matrix of retrieval patterns (rows are `r_k`), aligned by key.

## First‑Order: ERA Similarity (E × R)
Define the encoding–retrieval similarity matrix:

`S^{ER}_{ij} = sim(e_i, r_j)`

- `sim` can be correlation, crossnobis, cosine, etc.
- Diagonal `S^{ER}_{ii}`: same‑item E–R similarity (standard ERA)
- Off‑diagonals: cross‑item generalization/confusions

From `S^{ER}` we derive:
- ERA accuracy: proportion of items where `S^{ER}_{ii}` is the largest in row i (1‑vs‑all decoding).
- ERA effect size: e.g., `mean(diag) − mean(off‑diag)` or `d'` between diagonal and off‑diagonal distributions.

This is “naive cross‑decoding” in RSA clothing.

## Second‑Order: Encoding vs Retrieval Geometries
Build separate RDMs for each phase:
- `D^{EE}_{ij} = d(e_i, e_j)`
- `D^{RR}_{ij} = d(r_i, r_j)`

with e.g. `d = 1 − r` or Euclidean. Then compute a second‑order ER geometry similarity:

`geom_corr = corr( vec(upper(D^{EE})), vec(upper(D^{RR})) )`

- Typically Spearman correlation of the upper triangles.
- Variants: Kendall; partial correlations controlling for model RDMs (semantics, category, DNNs, behavior, etc.).

## Outputs per ROI
- ERA metrics: accuracy, diagonal vs off‑diagonal means, delta (`diag − off`), AUC if desired.
- Geometry metrics: ER RDM correlation (`geom_corr`), optionally partials.
- Joint view: plot or model `ERA_acc` vs `geom_corr` to ask whether reinstatement is tied to preserved geometry or remapped geometry.

## rMVPA Integration: `era_rsa_model`
A new analysis type parallel to `rsa_model` / `vector_rsa_model`, specialized for dual‑phase (encoding/retrieval) data.

### API Sketch
```
era_rsa_model <- function(dataset,
                          design,
                          link_by,
                          simfun = cordist(),   # for ERA similarities
                          distfun = cordist(),  # for RDMs (often 1 - r)
                          center = TRUE,
                          shrink = TRUE,
                          ...) {
  create_model_spec(
    "era_rsa_model",
    dataset = dataset,
    design  = design,
    link_by = link_by,
    simfun  = simfun,
    distfun = distfun,
    center  = center,
    shrink  = shrink,
    compute_performance = TRUE,
    return_predictions  = TRUE,
    ...
  )
}
```
- Use `mvpa_dataset(train_data = encoding, test_data = retrieval, mask)`.
- Use `mvpa_design(train_design, test_design, y_train = ~ item, y_test = ~ item, block_var = ~ run)` and pass `link_by` (e.g., "item").
- `simfun`/`distfun` are rMVPA distance objects; for RDMs we typically want `1 − correlation`.

### Per‑ROI Computation: `process_roi.era_rsa_model`
1) Extract ROI data and designs
- `Xenc <- values(roi$train_roi)`; `Xret <- values(roi$test_roi)` (as matrices)
- `dtr <- mod_spec$design$train_design`; `dte <- mod_spec$design$test_design`
- `keys_tr <- dtr[[link_by]]`; `keys_te <- dte[[link_by]]`

2) Build prototypes
- Encoding: average repeats by key
- Retrieval: single trial or average repeats

Implementation idea:
```
E <- rowsum(Xenc, group = factor(keys_tr), reorder = TRUE) /
     as.vector(table(factor(keys_tr)))
R <- rowsum(Xret, group = factor(keys_te), reorder = TRUE) /
     pmax(1, as.vector(table(factor(keys_te))))

common <- intersect(rownames(E), rownames(R))
E <- E[common, , drop = FALSE]
R <- R[common, , drop = FALSE]
K <- nrow(E)
```

3) First‑order ERA similarity `S_ER`
- Compute similarity between every encoding prototype and every retrieval prototype, e.g. correlation:
```
S_ER <- cor(t(E), t(R))  # K × K: rows = E, cols = R
```
- Metrics:
  - `diag_vals <- diag(S_ER)`
  - `off_vals`: all off‑diagonals (upper or upper+lower)
  - `era_acc <- mean(apply(S_ER, 1, which.max) == seq_len(K))`
  - `era_delta <- mean(diag_vals) - mean(off_vals)`

4) Second‑order RDMs and their similarity
```
D_EE <- distfun(E)  # K × K
D_RR <- distfun(R)

v_EE <- D_EE[upper.tri(D_EE, diag = FALSE)]
v_RR <- D_RR[upper.tri(D_RR, diag = FALSE)]
geom_corr <- cor(v_EE, v_RR, method = "spearman")
```

5) Optional joint/composite metrics
- ROI‑level scatter: `ERA_acc` vs `geom_corr`
- Group‑level residuals or mixed models are better handled downstream; per‑ROI just store both.

6) Package results using rMVPA containers
- Although not a classic classifier, we can reuse `classification_result`:
  - Observed labels = item keys (`common`)
  - Predicted labels = `argmax` per row in `S_ER`
  - Probabilities = `S_ER` (optionally row‑normalized / softmax)

Sketch:
```
cres <- classification_result(
  observed     = factor(common, levels = common),
  predicted    = factor(common[apply(S_ER, 1, which.max)], levels = common),
  probs        = S_ER,
  testind      = seq_len(K),
  test_design  = dte[match(common, keys_te), , drop = FALSE]
)

base_perf <- performance(cres)
extra <- c(
  ERA_acc       = era_acc,
  ERA_diag_mean = mean(diag_vals),
  ERA_off_mean  = mean(off_vals),
  ERA_delta     = era_delta,
  ER_geom_corr  = geom_corr
)
perf_vec <- c(base_perf, extra)

out <- tibble::tibble(
  result       = list(cres),
  indices      = list(neuroim2::indices(roi$train_roi)),
  performance  = list(perf_vec),
  id           = rnum,
  error        = FALSE,
  error_message= "~",
  warning      = FALSE,
  warning_message = "~"
)
```

- This integrates cleanly with `run_regional(era_rsa_model(...), region_mask)` and `run_searchlight`.

## Enhancements and Extensions
These remain API‑friendly within rMVPA.

### 1) Domain‑adaptive ERA via REMAP‑RRR
- Use REMAP to transform encoding prototypes: `E' = E + EΔ`.
- Compute ERA similarity and geometry on `E'` vs `R`.
- Compare `ERA_naive` vs `ERA_remap`, `ER_geom_corr_naive` vs `ER_geom_corr_remap` to test whether adaptation improves reactivation and/or geometry alignment.
- Implementation options:
  - `era_remap_model` wrapping REMAP inside `process_roi`, or
  - `era_rsa_model(transform = c("naive","remap"))` to factor shared code.

### 2) Geometric remapping diagnostics
- Procrustes alignment between E and R (after mean‑centering); use Procrustes distance as another geometry metric.
- Decompose RDM similarity using contrast RSA:
  - `RDM_E` vs `RDM_R`
  - `RDM_E` vs category/semantic model RDMs
  - `RDM_R` vs category/semantic model RDMs
- Ask whether retrieval geometry becomes more categorical relative to encoding.

### 3) Vector RSA / trial‑level ERA
- Treat trial‑level E–R similarity profiles as RSA vectors.
- For each retrieval trial, vectorize `S_ER[i, ·]` and regress against model predictors (semantic similarity, behavior).
- Answers whether E–R similarity structure reflects semantic organization beyond identity reinstatement.
- Could piggy‑back on `vector_rsa_model` by constructing the design, or add `era_vector_rsa_model` to build vectors internally.

### 4) Joint first‑/second‑order analysis
- ROI‑wise scatterplots: `ERA_acc` vs `ER_geom_corr` (per subject)
- Group‑level modeling: e.g., LMM where geometry predicts ERA across subjects/ROIs.
- Use rMVPA’s tidy `prediction_table` / `performance_table` with dplyr/ggplot.

## Positioning within rMVPA
- `era_rsa_model` is a new analysis family that:
  - Expects a train/test dataset (encoding/retrieval)
  - Uses `link_by` to define item keys
  - Returns `classification_result` plus a rich performance vector with ERA and geometry metrics
- Works seamlessly with `run_regional` / `run_searchlight` and can be combined with REMAP‑RRR, classic RSA, contrast RSA, vector/feature RSA in the same pipeline.

## Notes on Distances/Similarities
- `simfun` for `S_ER`: correlation, cosine, crossnobis; convert distance objects to similarity if needed.
- `distfun` for RDMs: often `1 − correlation`; can swap in Euclidean or crossnobis distances.
- Consider centering/shrinkage options (`center`, `shrink`) consistent with existing rMVPA utilities.

## Practical Considerations
- Ensure only items present in both phases are used (intersect keys).
- If retrieval has multiple trials per item, average or choose a policy parameter.
- Define off‑diagonal sampling consistently (upper triangle vs both off‑diagonal halves) for `ERA_off_mean`.
- Optionally row‑normalize `S_ER` (e.g., softmax) before storing in `probs`.
- Keep ROI indices and design alignment intact for tidy outputs downstream.

## Summary
- First‑order ERA quantifies item‑wise reinstatement; second‑order ER geometry quantifies preservation/warping of representational structure.
- `era_rsa_model` provides both within the rMVPA framework, returning standard classification‑like outputs plus ERA/geometry metrics.
- Extensions (REMAP, Procrustes, vector RSA, contrast RSA) allow rich, publishable analyses without changing the API surface.

---

# Discussion 2 Addendum: Confounders via `rsa_design`

This addendum makes ERA‑RSA feel like “just another RSA flavor” that supports multiple model RDMs and confounds using the existing `rsa_design` → `rsa_model` idiom, while still delivering first‑order ERA metrics.

## Design Reuse: `rsa_design` for confounders
- Reuse `rsa_design` to specify target model RDM(s) plus confounder RDMs (block, temporal position, etc.).
- `era_rsa_model(dataset, design, link_by, ...)` then:
  - builds encoding/retrieval prototypes (`E`, `R`),
  - computes neural RDMs (`D_EE`, `D_RR`),
  - regresses each neural RDM vector on the design’s model RDM vectors (multi‑RDM, optional semipartial),
  - computes first‑order ERA similarity `S_ER` and derived metrics.

### Example design with confounds
```r
# K items
K <- length(unique(item_ids))

# Model RDMs (K x K, symmetric)
D_sem   <- as.dist(semantic_dist)   # theory RDM
D_block <- as.dist(block_dist)      # confound: same vs different block
D_time  <- as.dist(time_dist)       # confound: temporal lag

era_design <- rsa_design(
  ~ D_sem + D_block + D_time,
  data_list = list(
    D_sem   = D_sem,
    D_block = D_block,
    D_time  = D_time,
    block   = dataset$design$block_var  # optional for perms
  ),
  block_var = "block"
)

era_spec <- era_rsa_model(
  dataset     = dataset,     # mvpa_dataset(train=encoding, test=retrieval)
  design      = era_design,  # rsa_design with multiple Dmats
  link_by     = "image_id",  # align E/R trials to items
  regtype     = "lm",
  semipartial = TRUE,        # semipartial RSA-style control
  distfun     = cordist()    # neural distances (often 1 - r)
)
```

## Per‑ROI Steps (with confounds)
1) Prototypes and key alignment
```r
Xenc <- as.matrix(neuroim2::values(roi$train_roi))
Xret <- as.matrix(neuroim2::values(roi$test_roi))

dtr <- mod_spec$design$train_design
dte <- mod_spec$design$test_design

key_tr <- factor(dtr[[mod_spec$link_by]])
key_te <- factor(dte[[mod_spec$link_by]])

E <- rowsum(Xenc, key_tr, reorder = TRUE) / as.vector(table(key_tr))
R <- rowsum(Xret, key_te, reorder = TRUE) / pmax(1, as.vector(table(key_te)))

common <- intersect(rownames(E), rownames(R))
E <- E[common, , drop = FALSE]
R <- R[common, , drop = FALSE]
K <- nrow(E)
```

2) First‑order ERA (E×R)
```r
S_ER <- cor(t(E), t(R))  # K x K

diag_vals <- diag(S_ER)
off_vals  <- S_ER[upper.tri(S_ER, diag = FALSE) | lower.tri(S_ER, diag = FALSE)]
era_acc   <- mean(apply(S_ER, 1, which.max) == seq_len(K))
era_delta <- mean(diag_vals) - mean(off_vals)
```

3) Second‑order geometries and ER correlation
```r
D_EE <- distfun(E)
D_RR <- distfun(R)

v_EE <- D_EE[upper.tri(D_EE, diag = FALSE)]
v_RR <- D_RR[upper.tri(D_RR, diag = FALSE)]
geom_corr <- cor(v_EE, v_RR, method = "spearman")
```

4) Fold in model RDMs + confounds (RSA‑style)
```r
# Vectorized model RDMs from rsa_design (same item order)
mm_list <- rsa_model_mat(mod_spec$design)  # list of vectors
X_model <- as.data.frame(mm_list)          # Npair x P

# Regress neural RDMs on model predictors (or semipartial as in rsa_model)
fit_E <- lm(v_EE ~ ., data = X_model)
fit_R <- lm(v_RR ~ ., data = X_model)

beta_E <- coef(fit_E)[-1]  # per predictor RDM
beta_R <- coef(fit_R)[-1]
```
- With `semipartial = TRUE`, compute semipartial correlations as in `rsa_model` (i.e., regress out confounds and correlate residuals), mirroring existing implementation.

5) Package outputs
```r
obs  <- factor(common, levels = common)
pred <- factor(common[apply(S_ER, 1, which.max)], levels = common)
probs <- S_ER  # optionally row-softmax

cres <- classification_result(
  observed    = obs,
  predicted   = pred,
  probs       = probs,
  testind     = seq_len(K),
  test_design = dte[match(common, key_te), , drop = FALSE]
)

base_perf <- performance(cres)
extra <- c(
  ERA_acc       = era_acc,
  ERA_diag_mean = mean(diag_vals),
  ERA_off_mean  = mean(off_vals),
  ERA_delta     = era_delta,
  ER_geom_corr  = geom_corr,
  setNames(beta_E, paste0(names(beta_E), "_E")),
  setNames(beta_R, paste0(names(beta_R), "_R"))
)
perf_vec <- c(base_perf, extra)
```

## Enhancements with confounds in place
- Domain‑adaptive ERA (REMAP): transform `E` → `E'` via REMAP; recompute `S_ER`, `D_EE'`, `geom_corr`, and RSA coefficients; compare deltas (first‑ and second‑order).
- Geometric remapping diagnostics: compare `beta_*_E` vs `beta_*_R` (e.g., semantics vs block/time) and add a Procrustes distance metric if desired.
- ERA‑vector RSA: use rows of `S_ER` as vectors in a vector‑RSA with the same `rsa_design` confounds.

## TL;DR
- Define model + confound RDMs in `rsa_design` (e.g., `~ D_sem + D_block + D_time`).
- `era_rsa_model(...)` computes:
  - First‑order ERA metrics (accuracy, delta, etc.).
  - Encoding and retrieval neural RDMs and their correlation.
  - RSA‑style coefficients for each predictor RDM, separately for encoding and retrieval, with semipartial control identical to `rsa_model`.
- Returns a single per‑ROI performance vector compatible with `run_regional` / `run_searchlight`.

---

# Discussion 3: Concrete Implementation Outline

This section captures the concrete implementation sketched for adding ERA‑RSA directly into rMVPA as:
- `era_rsa_model()` — model spec constructor (mirrors `rsa_model()`/`vector_rsa_model()` with `create_model_spec`).
- `process_roi.era_rsa_model()` — per‑ROI processor compatible with `run_regional()` / `run_searchlight()`.

It implements:
- First‑order ERA (encoding–retrieval similarity / matching)
- Second‑order ER geometry (RDM of encoding vs RDM of retrieval)
- Optional confound RDMs (block, temporal, etc.) with RSA‑style regression and semipartial correlations.

## 1) Constructor: `era_rsa_model()`
Expects a standard `mvpa_dataset` and `mvpa_design`, plus key/phase variables and options.

Arguments (key ones):
- `dataset`: `mvpa_dataset` (train_data + mask)
- `design`: `mvpa_design` (trial table)
- `key_var`: column/formula giving the item/key that links encoding and retrieval
- `phase_var`: column/formula indicating phase (encoding vs retrieval)
- `encoding_level`, `retrieval_level`: which levels of `phase_var` map to E/R (default: first two levels)
- `distfun`: distance function for within‑phase RDMs (object or shortcut handled by `create_dist()`)
- `rsa_simfun`: similarity for comparing `D_EE` vs `D_RR` ("pearson"|"spearman")
- `confound_rdms`: optional named list of K×K matrices or `dist` objects at item level
- `include_diag`, `return_item_scores`, `...`

Key constructor steps:
- Parse `phase_var`, `key_var` from `design$train_design` via a `parse_variable` helper.
- Validate phase levels and default `encoding_level` / `retrieval_level` if not supplied.
- Normalize `distfun` to a `distfun` object.
- Normalize `confound_rdms` to square matrices with row/col names matching item keys where possible.
- Call `create_model_spec("era_rsa_model", ...)` storing: dataset, design, parsed key/phase factors, E/R level labels, `distfun`, `rsa_simfun`, `confound_rdms`, flags, and standard spec fields (`compute_performance = TRUE`, `return_predictions = FALSE`).

## 2) Per‑ROI Processing: `process_roi.era_rsa_model()`
For each ROI/searchlight, compute ERA metrics, geometry correlation, and optional RSA regression with confounds.

High‑level steps:
1) Extract ROI data: `X <- values(roi$train_roi)` (trials × voxels); get `indices`.
2) Validate sizes; ensure `phase` and `key` vectors (from `mod_spec`) match `nrow(X)`.
3) Build item prototypes per phase via `group_means`:
   - `E_full <- group_means(X[enc_idx,], group = key[enc_idx])`
   - `R_full <- group_means(X[ret_idx,], group = key[ret_idx])`
   - Intersect and sort keys; subset to `E`, `R`; require `K ≥ 2`.
   - Drop voxels with zero variance across both phases.
4) First‑order ERA (E×R):
   - `S <- cor(t(E), t(R))` (rows = encoding items; cols = retrieval items; size K×K).
   - `era_diag_mean <- mean(diag(S))`.
   - Off‑diagonal mean: set `diag(S) <- NA`; `era_off_mean <- mean(S, na.rm = TRUE)`.
   - `era_diag_minus_off <- era_diag_mean − era_off_mean`.
   - Top‑1 matching accuracy: for each retrieval column, `which.max` over encoding rows and compare to identity index (`seq_len(K)`). Note: orientation is E(rows)×R(cols).
5) Second‑order geometries:
   - `D_enc <- pairwise_dist(distfun, E)`; `D_ret <- pairwise_dist(distfun, R)`.
   - Vectorize lower triangles: `dE`, `dR` and correlate with `method = rsa_simfun` to get `geom_cor`.
6) Optional geometry regression with confounds:
   - Build vectorized predictors aligned to `common_keys`:
     - Base predictor: `enc_geom = dE`.
     - For each confound RDM, subset to `common_keys` and vectorize lower triangle.
   - Optionally call `check_collinearity(mm)` if available.
   - Fit `lm(dR ~ ., data = as.data.frame(mm))`; extract non‑intercept betas as `beta_<name>`.
   - If `run_lm_semipartial()` exists, compute semipartial correlations; store as `sp_<name>`.

Return value: tibble with columns
- `result = NULL` (no per‑trial predictions for now)
- `indices = neuroim2::indices(roi$train_roi)`
- `performance = list(c(n_items = K, era_top1_acc, era_diag_mean, era_diag_minus_off, geom_cor, beta_*?, sp_*?))`
- `id = rnum`; `error|warning` flags + messages

### R snippet (illustrative only)
```
# E x R similarity
S <- suppressWarnings(stats::cor(t(E), t(R), use = "pairwise.complete.obs"))
diag_sim <- diag(S)
era_diag_mean <- mean(diag_sim, na.rm = TRUE)
off <- S; diag(off) <- NA_real_
era_off_mean <- mean(off, na.rm = TRUE)
era_diag_minus_off <- era_diag_mean - era_off_mean
max_enc_idx <- apply(S, 2, function(col) if (all(is.na(col))) NA_integer_ else which.max(col))
era_top1_acc <- mean((!is.na(max_enc_idx)) & (max_enc_idx == seq_len(K)))

# Geometry
D_enc <- pairwise_dist(distfun, E)
D_ret <- pairwise_dist(distfun, R)
dE <- as.numeric(D_enc[lower.tri(D_enc)])
dR <- as.numeric(D_ret[lower.tri(D_ret)])
geom_cor <- stats::cor(dE, dR, method = mod_spec$rsa_simfun, use = "complete.obs")

# Confounds (optional)
mm <- list(enc_geom = dE)
for (nm in names(mod_spec$confound_rdms)) {
  M <- as.matrix(mod_spec$confound_rdms[[nm]])
  M <- M[common_keys, common_keys, drop = FALSE]
  mm[[nm]] <- as.numeric(M[lower.tri(M)])
}
fit <- stats::lm(dR ~ ., data = as.data.frame(mm))
beta_vec <- stats::coef(fit)[-1]
```

## 3) Integration and Example Usage
Because the spec uses `create_model_spec("era_rsa_model", ...)` and the method `process_roi.era_rsa_model()` follows the standard dispatch, it plugs in exactly like other models:

```
mspec <- era_rsa_model(
  dataset = my_dataset,
  design  = my_design,
  key_var   = ~ ImageID,
  phase_var = ~ Phase,
  encoding_level  = "enc",
  retrieval_level = "ret",
  distfun    = cordist(method = "pearson"),
  rsa_simfun = "pearson",
  confound_rdms = list(
    block = block_rdm,
    temp  = temp_rdm
  )
)

res_reg <- run_regional(mspec, region_mask = my_roi_mask)
res_sl  <- run_searchlight(mspec, sl_mask = my_sl_mask, radius = 3)
```

Performance table includes:
- `n_items`, `era_top1_acc`, `era_diag_mean`, `era_diag_minus_off`, `geom_cor`
- plus any `beta_*` and `sp_*` terms from geometry regression.

## 4) Next Steps (from the discussion)
- Refine first‑order metrics (e.g., block‑limited ERA, lag‑specific ERA).
- Consider a small `era_rsa_design()` helper to auto‑build key/phase/confound RDMs from a tidy design table, mirroring the ergonomics of `rsa_design()`.

---

# Discussion 4: Block‑Limited and Lag‑Specific ERA + `era_rsa_design` Helper

This section adds two practical first‑order analyses and a tiny helper that builds item‑level confound RDMs from your design, while keeping the main ERA‑RSA API clean.

## 1) Concepts

### 1.1 Block‑limited ERA
- Standard ERA uses the full off‑diagonal of `S_ER` as the “non‑match” baseline.
- Off‑diagonals can be inflated by block/run confounds (shared context). To diagnose this, compute ERA using only selected off‑diagonal pairs:
  - Within‑block ERA: use only pairs where items share the same block.
  - Across‑block ERA: use only pairs where items come from different blocks.
- Metrics:
  - `Δ_ERA_same = mean(diag(S_ER)) − mean(off‑diag within same block)`
  - `Δ_ERA_diff = mean(diag(S_ER)) − mean(off‑diag across blocks)`
- Interpretation: if ERA is large only within‑block, block structure may drive the effect; if it survives across‑block, it’s more robustly mnemonic.

### 1.2 Lag‑specific ERA
- For each item `i` with encoding time `t_enc[i]` and retrieval time `t_ret[i]`, define lag `ℓ[i] = t_ret[i] − t_enc[i]`.
- Quantify whether reactivation decays with lag using a per‑ROI correlation:
  - `ρ_ERA,lag = cor(diag(S_ER), ℓ, method = "spearman")` (or Pearson).
- Single scalar per ROI reporting lag sensitivity of item‑wise reinstatement.

## 2) Helper: `era_rsa_design()`
Builds item lists and optional item‑level confound summaries (block, time), plus ready‑to‑use item‑level confound RDMs.

```r
#' Build item-level ERA-RSA confound RDMs from an mvpa_design
#'
#' @param design mvpa_design (same object you pass to era_rsa_model).
#' @param key_var column name or formula giving the item ID (e.g. ~ ImageID).
#' @param phase_var column name or formula indicating phase (e.g. ~ Phase).
#' @param encoding_level(optional) encoding phase label; default = first level.
#' @param retrieval_level(optional) retrieval phase label; default = second level.
#' @param block_var(optional) column giving block/run membership.
#' @param time_var(optional) column giving trial index or onset time.
#'
#' @return A list with elements:
#'   - items: factor of item IDs used (only keys with both enc & ret trials)
#'   - item_block: optional factor of per-item block labels
#'   - item_time_enc, item_time_ret: optional numeric vectors
#'   - item_lag: optional numeric vector (ret - enc)
#'   - confound_rdms: named list of item-level RDMs (block/time-based)
#' @export
era_rsa_design <- function(design,
                           key_var,
                           phase_var,
                           encoding_level  = NULL,
                           retrieval_level = NULL,
                           block_var = NULL,
                           time_var  = NULL) {

  assertthat::assert_that(inherits(design, "mvpa_design"))
  d <- design$train_design

  key   <- parse_variable(key_var,   d)
  phase <- factor(parse_variable(phase_var, d))
  phase_lev <- levels(phase)
  if (is.null(encoding_level))  encoding_level  <- phase_lev[1L]
  if (is.null(retrieval_level)) retrieval_level <- phase_lev[2L]

  enc_mask <- phase == encoding_level
  ret_mask <- phase == retrieval_level

  block <- if (!is.null(block_var)) parse_variable(block_var, d) else NULL
  timev <- if (!is.null(time_var))  parse_variable(time_var,  d) else NULL

  # items with at least one enc and one ret trial
  key_enc <- unique(key[enc_mask])
  key_ret <- unique(key[ret_mask])
  common_items <- sort(intersect(key_enc, key_ret))
  K <- length(common_items)
  if (K < 2L) warning("era_rsa_design: fewer than 2 items with both encoding and retrieval trials.")

  Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }

  item_block     <- NULL
  item_time_enc  <- NULL
  item_time_ret  <- NULL
  item_lag       <- NULL
  confound_rdms  <- list()

  if (!is.null(block)) {
    item_block <- sapply(common_items, function(k) {
      b_enc <- block[enc_mask & key == k]
      if (length(b_enc)) Mode(b_enc) else NA
    })
    names(item_block) <- common_items

    # same/different-block RDM (0 same, 1 diff)
    B <- outer(item_block, item_block, FUN = function(a, b) as.numeric(a != b))
    rownames(B) <- colnames(B) <- common_items
    confound_rdms$block <- B
  }

  if (!is.null(timev)) {
    item_time_enc <- sapply(common_items, function(k) {
      t_enc <- timev[enc_mask & key == k]
      if (length(t_enc)) mean(t_enc) else NA_real_
    })
    item_time_ret <- sapply(common_items, function(k) {
      t_ret <- timev[ret_mask & key == k]
      if (length(t_ret)) mean(t_ret) else NA_real_
    })
    names(item_time_enc) <- common_items
    names(item_time_ret) <- common_items

    item_lag <- item_time_ret - item_time_enc

    # encoding-time distance RDM
    Tenc <- outer(item_time_enc, item_time_enc, FUN = function(a, b) abs(a - b))
    rownames(Tenc) <- colnames(Tenc) <- common_items
    confound_rdms$time_enc <- Tenc
  }

  list(
    items          = factor(common_items, levels = common_items),
    item_block     = if (!is.null(item_block)) factor(item_block) else NULL,
    item_time_enc  = item_time_enc,
    item_time_ret  = item_time_ret,
    item_lag       = item_lag,
    confound_rdms  = confound_rdms
  )
}
```

Usage with `era_rsa_model`:

```r
era_des <- era_rsa_design(
  design       = my_design,
  key_var      = ~ ImageID,
  phase_var    = ~ Phase,
  encoding_level  = "enc",
  retrieval_level = "ret",
  block_var    = ~ Block,
  time_var     = ~ Trial
)

mspec <- era_rsa_model(
  dataset       = my_dataset,
  design        = my_design,
  key_var       = ~ ImageID,
  phase_var     = ~ Phase,
  encoding_level  = "enc",
  retrieval_level = "ret",
  distfun       = cordist("pearson"),
  rsa_simfun    = "spearman",
  confound_rdms = era_des$confound_rdms,
  # optional extras for block/lag metrics:
  item_block    = era_des$item_block,
  item_lag      = era_des$item_lag
)
```

## 3) Adding block‑limited and lag‑specific ERA to `process_roi.era_rsa_model`
Patch snippets to add after computing `S` and `diag_sim`:

```r
# --- Block-limited ERA (if item_block is available) ------------------------
era_diag_minus_off_same_block <- NA_real_
era_diag_minus_off_diff_block <- NA_real_

if (!is.null(mod_spec$item_block)) {
  # align item_block to the common_keys subset
  ib <- mod_spec$item_block
  if (!is.null(names(ib))) ib <- ib[match(common_keys, names(ib))] else ib <- ib[seq_along(common_keys)]

  same_vals <- c(); diff_vals <- c()
  for (i in seq_len(K)) {
    same_idx <- setdiff(which(ib == ib[i]), i)
    diff_idx <- which(ib != ib[i])
    if (length(same_idx)) same_vals <- c(same_vals, S[i, same_idx])
    if (length(diff_idx)) diff_vals <- c(diff_vals, S[i, diff_idx])
  }

  if (length(same_vals)) era_diag_minus_off_same_block <- mean(diag_sim, na.rm = TRUE) - mean(same_vals, na.rm = TRUE)
  if (length(diff_vals)) era_diag_minus_off_diff_block <- mean(diag_sim, na.rm = TRUE) - mean(diff_vals, na.rm = TRUE)
}

# append to perf
perf <- c(perf,
  era_diag_minus_off_same_block = era_diag_minus_off_same_block,
  era_diag_minus_off_diff_block = era_diag_minus_off_diff_block
)

# --- Lag-specific ERA (diag ERA vs E->R lag) ------------------------------
era_lag_cor <- NA_real_
if (!is.null(mod_spec$item_lag)) {
  lag <- mod_spec$item_lag
  if (!is.null(names(lag))) lag <- lag[match(common_keys, names(lag))] else lag <- lag[seq_along(common_keys)]
  if (length(lag) == length(diag_sim) && any(!is.na(lag))) {
    era_lag_cor <- suppressWarnings(stats::cor(diag_sim, lag, method = "spearman", use = "complete.obs"))
  }
}
perf <- c(perf, era_lag_cor = era_lag_cor)
```

## 4) Recommended Workflow
1) Build design helpers:
```r
era_des <- era_rsa_design(
  design       = my_design,
  key_var      = ~ ImageID,
  phase_var    = ~ Phase,
  encoding_level  = "enc",
  retrieval_level = "ret",
  block_var    = ~ Block,
  time_var     = ~ Trial
)
```
2) Create model spec:
```r
mspec <- era_rsa_model(
  dataset  = my_dataset,
  design   = my_design,
  key_var      = ~ ImageID,
  phase_var    = ~ Phase,
  encoding_level  = "enc",
  retrieval_level = "ret",
  distfun       = cordist("pearson"),
  rsa_simfun    = "spearman",
  confound_rdms = era_des$confound_rdms,
  item_block    = era_des$item_block,
  item_lag      = era_des$item_lag
)
```
3) Run and inspect:
```r
res <- run_regional(mspec, region_mask)
# res$performance_table now contains:
# - era_top1_acc, era_diag_mean, era_diag_minus_off
# - geom_cor
# - beta_* / sp_* (if geometry regression/confounds enabled)
# - era_diag_minus_off_same_block, era_diag_minus_off_diff_block, era_lag_cor
```

---

# Discussion 5: Run Confounds and Temporal Lags — Reusing `R/temporal_rdms.R`

Second‑order ER geometry can be inflated by shared run/block structure, especially when encoding and retrieval occur in the same run. To keep ERA‑RSA coherent with the existing codebase and avoid redundancy, we reuse helpers in `R/temporal_rdms.R` wherever possible and add minimal metrics for robustness.

## 1) Conceptual picture
- First‑order ERA (often cross‑run by design) tends to be less affected by run artifacts.
- Second‑order ER geometry (`D_EE` vs `D_RR`) can be boosted by shared run structure within each phase, even if mnemonic geometry is modest.

## 2) Strategies

### A) Remove run at the pattern level (upstream)
- Regress out run intercepts or demean patterns within run before building prototypes. This benefits all analyses.

### B) Partial out run at the geometry level (preferred minimal change)
- Build item‑level run RDMs for encoding and retrieval (0/1 same‑run) and include them as confounds. Then compute a run‑partial ER geometry metric by residualizing both `dE` and `dR` on these run vectors before correlating.
- This mirrors standard RSA semipartial practice and integrates with our `confound_rdms` story.

### C) Cross‑run‑only geometry (mask same‑run pairs)
- Compute ER geometry using only item pairs that do not share a run in the relevant phase(s). Report alongside the raw geometry corr.

## 3) Reusing `temporal_rdms.R` for lags and temporal confounds
- For item‑level temporal confound RDMs, avoid custom `outer()` code and call `temporal_rdm()` directly on per‑item times:
  - `confound_rdms$time_enc <- temporal_rdm(index = item_time_enc, block = item_block, kernel = "linear" /*or exp/gauss*/, within_blocks_only = TRUE)`
  - Optionally build a retrieval‑time RDM: `confound_rdms$time_ret <- temporal_rdm(index = item_time_ret, block = item_block, ...)`.
- If you have onsets rather than indices, use `temporal_from_onsets()` convenience sugar.
- For richer temporal nuisance sets, `temporal_confounds()` can generate multiple trial‑level nuisances that you can then reduce to item level if needed.
- Keep the lag‑specific ERA metric as a simple scalar (correlation of `diag(S)` with `item_lag`); it isn’t an RDM and doesn’t duplicate `temporal_rdms` functionality.

## 4) era_rsa_design: minimal, non‑redundant additions
- In addition to `item_block` and `item_lag` noted earlier, store phase‑specific run labels per item (using the mode across trials if repeated):
  - `item_run_enc` and `item_run_ret`.
- Produce run‑confound RDMs:
  - `confound_rdms$run_enc[i,j] = 1(a same encoding run as b)`
  - `confound_rdms$run_ret[i,j] = 1(a same retrieval run as b)`
- Produce temporal confounds via `temporal_rdm()` (rather than manual `outer()`):
  - `confound_rdms$time_enc <- temporal_rdm(item_time_enc, block = item_block, kernel = "exp", within_blocks_only = TRUE, metric = "distance")`
  - Optionally `confound_rdms$time_ret` similarly, if retrieval times vary.

## 5) process_roi additions (lightweight, leveraging confounds)
After computing `D_enc`, `D_ret`, and vectors `dE`, `dR`:

```r
# Run-partial ER geometry (residualize dE and dR on run confounds)
geom_cor_run_partial <- NA_real_
if (!is.null(mod_spec$confound_rdms$run_enc) && !is.null(mod_spec$confound_rdms$run_ret)) {
  Renc <- as.numeric(mod_spec$confound_rdms$run_enc[common_keys, common_keys][lower.tri(D_enc)])
  Rret <- as.numeric(mod_spec$confound_rdms$run_ret[common_keys, common_keys][lower.tri(D_ret)])
  conf_df <- data.frame(enc_run = Renc, ret_run = Rret)
  dE_res <- stats::resid(stats::lm(dE ~ ., data = conf_df))
  dR_res <- stats::resid(stats::lm(dR ~ ., data = conf_df))
  geom_cor_run_partial <- suppressWarnings(stats::cor(dE_res, dR_res, method = mod_spec$rsa_simfun, use = "complete.obs"))
}

# Cross-run-only geometry (mask pairs sharing a run in either phase)
geom_cor_xrun <- NA_real_
if (!is.null(mod_spec$item_run_enc) && !is.null(mod_spec$item_run_ret)) {
  ire <- mod_spec$item_run_enc[match(common_keys, names(mod_spec$item_run_enc))]
  irr <- mod_spec$item_run_ret[match(common_keys, names(mod_spec$item_run_ret))]
  same_enc <- outer(ire, ire, "==")
  same_ret <- outer(irr, irr, "==")
  mask <- lower.tri(D_enc) & !(same_enc | same_ret)
  if (any(mask)) {
    geom_cor_xrun <- suppressWarnings(stats::cor(D_enc[mask], D_ret[mask], method = mod_spec$rsa_simfun, use = "complete.obs"))
  }
}

# append to perf
perf <- c(perf,
  geom_cor_run_partial = geom_cor_run_partial,
  geom_cor_xrun = geom_cor_xrun
)
```

Notes:
- If you already include `run_enc`/`run_ret` in `confound_rdms`, their beta/semipartial entries will appear in `perf` from the existing geometry regression. The `geom_cor_run_partial` scalar complements those by reporting ER geometry after jointly removing run structure from both phases.
- Prefer encoding cross‑run ERA for first‑order metrics when feasible; the second‑order controls above then provide symmetry.

## 6) Methods‑section friendly reporting
- Report ERA as cross‑run (train‑run vs test‑run) where applicable.
- Report three ER geometry numbers per ROI: `geom_cor` (raw), `geom_cor_run_partial` (run‑partial), and optionally `geom_cor_xrun` (cross‑run only).
- Use `temporal_rdm()`‑based confounds (encoding time, optionally retrieval time) to show lag/carry‑over is not driving second‑order effects.
