# Representational Network Analysis (ReNA) – Planning Notes (Part 1)

These notes capture the first installment of a proposal for a new **Representational Network Analysis (ReNA)** module in `rMVPA`. The goal is to extend the package beyond per‑ROI/searchlight analyses of *what* is represented, to analyses of **how representational geometries flow through networks of regions and models**.

ReNA treats ROIs and searchlights as **nodes in a representational network**, linked by:

- **Representational connectivity**: who shares what representational geometry with whom?
- **Representational mediation**: does this node mediate the effect of a model or another region on memory/behavior?
- **Representational mapping (optional)**: how does this node transform the representational code it receives?

All analyses are:

- linear,
- R‑native (no GPUs),
- built on existing RSA + REMAP machinery,
- designed to expose **crisp, per‑ROI metrics** that can be turned into maps via `run_regional()` / `run_searchlight()`.

Conceptually, this is “hyperalignment / information connectivity / mediation analysis, but all in RDM space”, integrated into the rMVPA ecosystem.

---

## 0. Context: What rMVPA Already Has

Existing capabilities:

- MVPA (classification/regression),
- vanilla RSA,
- contrast RSA / MS‑ReVE,
- vector / feature RSA,
- REMAP (domain‑adaptive mapping),
- ERA‑RSA (encoding–retrieval similarity with confounds).

Missing piece:

- A first‑class framework to ask:
  - **How do representations flow through the system?**
  - **Which regions carry, transform, or mediate representational content from A → B → behavior?**

ReNA is intended to fill this gap.

---

## 1. High‑Level Design

ReNA will be implemented in three layers:

1. **Seed‑based Representational Connectivity (ReNA‑RC)**  
   Basic but powerful: “where in the brain instantiates a representational geometry similar to a seed, above confounds?”
2. **Representational Mediation (ReNA‑RM)**  
   Mediation in RDM space: “does this ROI’s geometry mediate the influence of a model/region on another region or behavior?”
3. **(Optional) Representational Mapping Layer (ReNA‑MAP)**  
   Uses RRR/REMAP‑style mappings to characterize how regions transform representational codes.

All three are framed as model specs compatible with:

- `run_regional()`
- `run_searchlight()`

and share design helpers and confound handling patterns with `rsa_design()` and REMAP.

---

## 2. Seed‑Based Representational Connectivity (ReNA‑RC)

### 2.1 Concept

Given:

- A **seed RDM** `D^{seed}` (K × K), which may come from:
  - another brain region,
  - another phase (e.g., encoding vs retrieval),
  - a model RDM (e.g., CLIP, semantics, behavior).
- A **target searchlight/ROI**, for which rMVPA already computes a neural RDM `D^{ROI}`.

For each target ROI, we ask:

> “Does this ROI instantiate a representational geometry similar to the seed’s, above and beyond cheap confounds (run, time, category, etc.)?”

Procedure per ROI:

- Compute neural RDM `D^{ROI}`.
- Align item order to the seed RDM `D^{seed}`.
- Vectorize lower triangles:
  - `d_roi  = vec_lower(D^{ROI})`
  - `d_seed = vec_lower(D^{seed})`
- Regress:

```r
d_ROI ~ d_seed + D^{block} + D^{lag} + D^{behav} + ...
```

From the fitted model we extract:

- **Representational connectivity coefficient**: `beta_seed`
- **Optional semipartial correlation**: `sp_seed` (like in RSA)
- **Raw correlation**: `conn_raw = cor(d_roi, d_seed)`

These become maps:

- `conn_raw`: “where in the brain looks like the seed?”
- `conn_seed_partial` / `sp_seed`: “where in the brain shares the seed geometry above confounds?”

This is essentially “information/representational connectivity”, implemented in:

- a searchlight‑native,
- confound‑aware,
- rMVPA‑idiomatic way.

### 2.2 API Sketch – Design Helper

`repnet_design()` is an RSA‑design‑lite helper for connectivity:

```r
repnet_design <- function(design,
                          key_var,
                          seed_rdm,
                          confound_rdms = list()) {
  # key_var: item id
  # seed_rdm: K x K matrix or "dist" with row/colnames = item ids
  # confound_rdms: list of KxK matrices (block, lag, behavior, etc.)

  list(
    key           = factor(parse_variable(key_var, design$train_design)),
    seed_rdm      = as.matrix(seed_rdm),
    confound_rdms = lapply(confound_rdms, as.matrix)
  )
}
```

Key points:

- `key_var` links items in the design to the seed RDM.
- Confounds are supplied as RDMs and converted to matrices; lower triangles are used downstream.

### 2.3 API Sketch – Model Spec

`repnet_model()` constructs a model specification compatible with `run_searchlight()` / `run_regional()`:

```r
repnet_model <- function(dataset,
                         design,
                         repnet_des,
                         distfun = cordist("pearson"),
                         simfun  = c("pearson", "spearman"),
                         ...) {

  simfun <- match.arg(simfun)

  create_model_spec(
    "repnet_model",
    dataset       = dataset,
    design        = design,
    key           = repnet_des$key,
    seed_rdm      = repnet_des$seed_rdm,
    confound_rdms = repnet_des$confound_rdms,
    distfun       = if (is.character(distfun)) create_dist(distfun) else distfun,
    simfun        = simfun,
    compute_performance = TRUE,
    return_predictions  = FALSE,
    ...
  )
}
```

Notes:

- `distfun` follows existing RSA conventions (`cordist("pearson")`, etc.).
- `simfun` specifies how to compute simple similarity metrics if desired.
- `compute_performance = TRUE` ensures per‑ROI metrics (connectivity coefficients) are returned.

### 2.4 Per‑ROI Logic – `process_roi.repnet_model`

Core per‑ROI steps:

1. **Compute ROI RDM**
   - Extract ROI patterns and compute `D_roi` using `distfun`.
2. **Align item order**
   - Ensure `D_roi` and `seed_rdm` are aligned to the same item order specified by `key`.
3. **Vectorize RDMs**

   ```r
   d_roi  <- D_roi[lower.tri(D_roi)]
   seed_mat <- seed_rdm_aligned
   d_seed <- seed_mat[lower.tri(seed_mat)]
   ```

4. **Build model matrix with seed and confounds**

   ```r
   mm <- list(seed = d_seed)
   for (nm in names(mod_spec$confound_rdms)) {
     M <- mod_spec$confound_rdms[[nm]]
     # align to items, vectorize lower tri
     mm[[nm]] <- M[lower.tri(M)]
   }
   df  <- as.data.frame(mm)
   fit <- lm(d_roi ~ ., data = df)
   ```

5. **Extract connectivity metrics**

   - `conn_raw  = cor(d_roi, d_seed)`
   - `beta_seed = coef(fit)["seed"]`
   - `sp_seed` (optional): computed via `run_lm_semipartial()` or similar helper.

6. **Return performance vector per ROI**

   ```r
   perf <- c(
     conn_raw  = conn_raw,
     beta_seed = beta_seed,
     sp_seed   = sp_seed,    # if computed
     n_pairs   = length(d_roi)
   )
   ```

Running:

- `run_searchlight(repnet_model(...))` or
- `run_regional(repnet_model(...))`

produces whole‑brain/ROI maps of representational connectivity to the seed geometry, with confound control.

---

## 3. Representational Mediation (ReNA‑RM)

### 3.1 Concept

Goal: implement **classical mediation analysis in RDM space**, per ROI.

We want to answer questions like:

- “Does region M mediate the representational influence of a model or region A on retrieval/behavior?”
- “Does a candidate region M carry transformations from a predictor geometry to an outcome geometry?”

Setup:

- **Predictor RDM** `D^X` (e.g., model semantics, low‑level vision, encoding region A).
- **Mediator RDM** `D^M` (searchlight/ROI under test).
- **Outcome RDM** `D^Y` (retrieval region, behavior, ERA RDM, etc.).

All three RDMs are defined over the same set of items.

Vectorize lower triangles:

- `x = vec(D^X)`
- `m = vec(D^M)`
- `y = vec(D^Y)`

Standard mediation model (per ROI):

1. Path **a**: `m ~ x`
2. Path **b**: `y ~ m + x`
3. Path **c′**: direct `x → y` controlling for `m`
4. Indirect effect: `a × b` (optionally with Sobel or bootstrap approximations)

Per ROI, we obtain:

- `a_roi`, `b_roi`, `cprime_roi`, `indirect_roi`, plus simple significance approximations if desired.

Example:

- `D^X`: semantic RDM from a language model.
- `D^M`: ATL region candidate.
- `D^Y`: hippocampal retrieval RDM.

ReNA‑RM then identifies **where in cortex the geometry both follows semantics and propagates that semantics to hippocampal retrieval geometry**.

### 3.2 API Sketch – Design Helper

`repmed_design()` defines the predictor and outcome RDMs plus confounds:

```r
repmed_design <- function(items,
                          X_rdm,  # predictor RDM
                          Y_rdm,  # outcome RDM
                          confound_rdms = list()) {

  X <- as.matrix(X_rdm); rownames(X) <- colnames(X) <- items
  Y <- as.matrix(Y_rdm); rownames(Y) <- colnames(Y) <- items

  list(
    items         = items,
    X_rdm         = X,
    Y_rdm         = Y,
    confound_rdms = lapply(confound_rdms, as.matrix)
  )
}
```

Notes:

- `items` defines the item universe and is used to align all RDMs.
- Confound RDMs again follow the same pattern as in connectivity.

### 3.3 API Sketch – Model Spec

`repmed_model()` wraps the design into a model spec:

```r
repmed_model <- function(dataset,
                         design,
                         repmed_des,
                         distfun = cordist("pearson"),
                         ...) {

  create_model_spec(
    "repmed_model",
    dataset       = dataset,
    design        = design,
    items         = repmed_des$items,
    X_rdm         = repmed_des$X_rdm,
    Y_rdm         = repmed_des$Y_rdm,
    confound_rdms = repmed_des$confound_rdms,
    distfun       = if (is.character(distfun)) create_dist(distfun) else distfun,
    compute_performance = TRUE,
    return_predictions  = FALSE,
    ...
  )
}
```

### 3.4 Per‑ROI Logic – `process_roi.repmed_model`

Steps per ROI:

1. **Compute mediator RDM (`D_M`)**
   - Extract ROI patterns and compute `D_M` using `distfun`.
2. **Align RDMs**
   - Align `X`, `M`, and `Y` to the same item order (`items`).
3. **Vectorize**

   ```r
   x <- X[lower.tri(X)]
   m <- M[lower.tri(M)]
   y <- Y[lower.tri(Y)]
   ```

4. **Handle confounds**
   - Option A: residualize `x`, `m`, and `y` against confound vectors.
   - Option B: include confounds directly as predictors in the regression formulas.
5. **Fit mediation regressions**

   ```r
   # Path a: m ~ x (+ confounds)
   a_fit <- lm(m ~ x + conf_mat)

   # Paths b and c': y ~ m + x (+ confounds)
   b_fit <- lm(y ~ m + x + conf_mat)
   ```

6. **Extract coefficients**

   - `a       = coef(a_fit)["x"]`
   - `b       = coef(b_fit)["m"]`
   - `cprime  = coef(b_fit)["x"]`
   - `indirect = a * b`
   - Optionally compute simple Sobel statistics or bootstrap estimates.

7. **Store performance metrics**

   ```r
   perf <- c(
     med_a        = a,
     med_b        = b,
     med_cprime   = cprime,
     med_indirect = indirect,
     n_pairs      = length(x)
   )
   ```

Running `run_searchlight(repmed_model(...))` yields maps of:

- where representational geometry mediates `X → Y`,
- where `X` has strong direct effects beyond that ROI’s geometry.

This provides a causal‑ish, network‑level story about representations in a standard, RDM‑native searchlight framework.

---

## 4. Optional Mapping Mode with REMAP / RRR (ReNA‑MAP)

### 4.1 Concept

A third variant would integrate REMAP / reduced‑rank regression (RRR) to characterize **transformations** between seed and target patterns, beyond simple connectivity/mediation:

- Instead of comparing RDMs, learn a mapping between seed patterns and target patterns using RRR (as in REMAP).
- Summarize:
  - mapping rank,
  - explained variance,
  - principal mapping axes,
  - how strongly mapping axes express particular model RDMs/features.

Per ROI:

- Fit a reduced‑rank mapping:

  ```r
  Y ≈ X C
  ```

  where:

  - `X`: seed patterns,
  - `Y`: ROI patterns,
  - `C`: reduced‑rank mapping matrix.

- Summaries:
  - **rank** of the mapping,
  - **pattern‑level R²** (`map_r2`),
  - SVD of `C` to obtain principal mapping axes.
- Optionally correlate these axes against external model RDMs/features.

This is “90% written already” via existing REMAP code; the main work is exposing it as a **representational connectivity/mapping primitive**.

---

## 5. Integration into rMVPA Ecosystem

ReNA would form a small “Representational Network” module with:

- **Design helpers**
  - `repnet_design()` – seed + confounds.
  - `repmed_design()` – predictor, outcome, confounds.
  - (Optional) `repmap_design()` – mapping‑oriented helper if needed.
- **Model specs**
  - `repnet_model()` → representational connectivity searchlights/ROIs.
  - `repmed_model()` → representational mediation searchlights/ROIs.
  - (Optional) `repmap_model()` → mapping searchlights/ROIs.
- **Standard runners**
  - `run_regional(repnet_model, ...)`
  - `run_searchlight(repnet_model, ...)`
  - analogous `run_regional()` / `run_searchlight()` for `repmed_*`, `repmap_*`.

Outputs:

- Connectivity maps:
  - `conn_raw`, `beta_seed`, `sp_seed`.
- Mediation maps:
  - `med_a`, `med_b`, `med_cprime`, `med_indirect`.
- Mapping maps (optional):
  - `map_rank`, `map_r2`, etc.

Integration points:

- **`rsa_design`**: reuse contrast RDM machinery where helpful.
- **ERA‑RSA**: e.g., set `Y` to retrieval RDMs from ERA‑RSA.
- **REMAP**: e.g., seed patterns could be P→M adapted encodings.
- **Block/lag/behavior confounds**: handled via RDMs/residualization as above.

---

## 6. Relation to Python Ecosystem

Python currently offers:

- searchlight MVPA,
- RSA implementations,
- various information connectivity scripts,
- ad‑hoc mediation‑style analyses.

What it lacks is a **coherent, RDM‑native representational network framework** that:

- treats connectivity, mediation, and mapping as **first‑class per‑searchlight analyses**,
- plugs into a **clean MVPA/RSA ecosystem**,
- is built on **honest cross‑validation and confound control**.

ReNA would provide:

> “Not just what each voxel represents, but who talks to whom, how, and on which representational axes, in a principled and inspectable way.”

All implemented with:

- linear models,
- reduced‑rank regression,
- RDMs, and
- S3 methods,

in line with the rMVPA aesthetic.

---

## 7. Next Steps (Outside This Part)

- Draft a concrete implementation of:
  - `repnet_model()` and
  - `process_roi.repnet_model()`
  
as a first installment of the Representational Network module, wired into existing searchlight/regional infrastructure.

---

## 8. Part 2 – Concrete R Implementations (ReNA‑RC / ReNA‑RM / ReNA‑Map)

Part 2 provides **drop‑in R code** for three S3 model families:

- `repnet_model` – Representational Connectivity (ReNA‑RC)
- `repmed_model` – Representational Mediation (ReNA‑RM)
- `repmap_model` – Representational Mapping (ReNA‑Map)

All three:

- operate in **RDM space** or item‑level feature space,
- plug into `run_regional()`, `run_searchlight()`, and `mvpa_iterate`,
- follow existing `mvpa_design` idioms (keys, confound RDMs),
- avoid GPUs/transformers while doing things Python users typically hand‑roll.

Below is a structured summary of the provided code, with function signatures and behavior.

---

### 8.1 Representational Connectivity – `repnet_model` (ReNA‑RC)

**Goal:** per ROI/searchlight, quantify how similar its representational geometry is to a **seed RDM**, optionally controlling for **confound RDMs** (block, lag, behavior, etc.).

#### 8.1.1 `repnet_design()`: connectivity design helper

Signature:

```r
repnet_design <- function(design,
                          key_var,
                          seed_rdm,
                          confound_rdms = NULL)
```

Behavior:

- Asserts `design` inherits `"mvpa_design"`.
- Uses `design$train_design` and `parse_variable(key_var, d)` to extract item identities; coerces to a factor `key_fac`.
- Coerces `seed_rdm` to a square matrix `S`:
  - If `"dist"`: uses `Labels` attribute if present; otherwise falls back to levels of `key_fac` and the `"Size"` attribute.
  - Ensures `rownames(S)` and `colnames(S)` are set to item labels.
- Handles `confound_rdms`:
  - Must be a **named list**.
  - Each element can be `"dist"` or matrix; converted to matrices with row/colnames (using labels, or `rownames(S)`, or `key_fac` levels as defaults).
- Returns a list:

  - `key` – factor of item IDs (length = number of training observations),
  - `seed_rdm` – seed RDM matrix,
  - `confound_rdms` – named list of confound RDM matrices.

#### 8.1.2 `repnet_model()`: connectivity model spec

Signature:

```r
repnet_model <- function(dataset,
                         design,
                         repnet_des,
                         distfun = cordist("pearson"),
                         simfun = c("pearson", "spearman"),
                         ...)
```

Behavior:

- Asserts `dataset` inherits `"mvpa_dataset"` and `design` inherits `"mvpa_design"`.
- `simfun` is matched to `"pearson"` or `"spearman"`.
- `distfun`:
  - If a character, passed through `create_dist()`.
  - Must inherit `"distfun"`.
- Calls `create_model_spec("repnet_model", ...)` with:

  - `dataset`, `design`,
  - `key`           = `repnet_des$key`,
  - `seed_rdm`      = `repnet_des$seed_rdm`,
  - `confound_rdms` = `repnet_des$confound_rdms`,
  - `distfun`, `simfun`,
  - `compute_performance = FALSE`,
  - `return_predictions  = FALSE`,
  - plus any extra fields in `...`.

Note: **Per‑ROI performance metrics are computed inside `process_roi.repnet_model()`**, not by generic `compute_performance`.

#### 8.1.3 `process_roi.repnet_model()`: per‑ROI connectivity metrics

Signature:

```r
process_roi.repnet_model <- function(mod_spec,
                                     roi,
                                     rnum,
                                     center_global_id = NA,
                                     ...)
```

Core steps:

1. **Extract ROI data**
   - `X <- as.matrix(neuroim2::values(roi$train_roi))` (`n_obs × n_vox`).
   - `ind <- neuroim2::indices(roi$train_roi)`.
   - If `n_obs < 2` or `n_vox < 1`: return a single‑row tibble with `error = TRUE` and NULL performance.
2. **Check key length**
   - `key <- mod_spec$key`; must have `length(key) == n_obs`, else error tibble.
3. **Item‑level prototypes**
   - `E_full <- group_means(X, margin = 1, group = key)` → `K_roi × V` matrix.
   - `items_roi <- rownames(E_full)`.
4. **Align with seed RDM**
   - `seed_mat <- mod_spec$seed_rdm`; `items_seed <- rownames(seed_mat)`.
   - `common_items <- intersect(items_roi, items_seed)`; require at least 3 items.
   - Sort `common_items`; subset:
     - `E      <- E_full[common_items, , drop = FALSE]`
     - `seed_sub <- seed_mat[common_items, common_items, drop = FALSE]`.
5. **Drop constant voxels**
   - `sd_E <- apply(E, 2, sd, na.rm = TRUE)`; keep voxels with `sd > 0`.
   - If no voxels remain: error tibble.
6. **Compute ROI RDM and vectorize**
   - `D_roi <- pairwise_dist(mod_spec$distfun, E)`.
   - `d_roi  <- as.numeric(D_roi[lower.tri(D_roi)])`.
   - `d_seed <- as.numeric(seed_sub[lower.tri(seed_sub)])`.
7. **Raw connectivity**
   - `conn_raw <- cor(d_roi, d_seed, method = mod_spec$simfun, use = "complete.obs")` (with `suppressWarnings()`).
8. **Build regression predictors (seed + confounds)**
   - Start `mm <- list(seed = d_seed)`.
   - For each named confound RDM:
     - Align rows/cols to `common_items` if possible (matching rownames).
     - Otherwise use `[seq_len(K), seq_len(K)]` as a fallback.
     - Vectorize lower triangle and add as numeric vector to `mm`.
9. **Multiple regression + semipartials**
   - `df <- as.data.frame(mm)`.
   - Fit `fit <- lm(d_roi ~ ., data = df)` with `try(...)`:
     - Extract non‑intercept coefficient estimates and name them `beta_<predictor>`.
   - If a function `run_lm_semipartial()` exists:
     - Build `tmp_obj <- list(design = list(model_mat = mm))`.
     - Call `run_lm_semipartial(dvec = d_roi, obj = tmp_obj)`; if numeric, prefix names with `"sp_"`.
10. **Assemble performance vector**

    - Always include:
      - `n_items`   = `K`,
      - `n_pairs`   = `length(d_roi)`,
      - `conn_raw`  = raw correlation.
    - Optionally include:
      - `beta_seed`, `beta_<confound>` etc.,
      - `sp_seed`, `sp_<confound>` semipartials.

11. **Return tibble row**

    - `result = list(NULL)` (no predictions),
    - `indices = list(ind)`,
    - `performance = list(perf)`,
    - `id = rnum`,
    - `error = FALSE`, `warning = FALSE` (or TRUE with messages in failure cases).

---

### 8.2 Representational Mediation – `repmed_model` (ReNA‑RM)

**Goal:** for each ROI/searchlight, test whether its RDM **mediates** between a predictor RDM `X` and an outcome RDM `Y` in RDM space, optionally with confound RDMs.

Classical paths in RDM space:

- Path **a**: `m ~ x + conf`
- Paths **b/c′**: `y ~ m + x + conf`
- Indirect effect: `a * b`, with optional Sobel z/p.

#### 8.2.1 `repmed_design()`: mediation design helper

Signature:

```r
repmed_design <- function(items,
                          X_rdm,
                          Y_rdm,
                          confound_rdms = NULL)
```

Behavior:

- Coerces `items` to a character vector.
- Uses an internal helper `as_mat()` to turn `"dist"` or matrices into square matrices with row/colnames:
  - `"dist"`: uses `Labels` if present; else uses `default_items[seq_len(Size)]`.
  - Matrices with missing row/colnames get both set from `default_items`.
- Constructs:
  - `X <- as_mat(X_rdm, items)`,
  - `Y <- as_mat(Y_rdm, items)`.
- `confound_rdms`:
  - optional named list of RDMs, each passed through `as_mat()` with `default_items = items`.
- Returns list:

  - `items` – canonical item IDs,
  - `X_rdm`, `Y_rdm` – aligned matrices,
  - `confound_rdms` – named list of matrices.

#### 8.2.2 `repmed_model()`: mediation model spec

Signature:

```r
repmed_model <- function(dataset,
                         design,
                         repmed_des,
                         key_var,
                         distfun = cordist("pearson"),
                         ...)
```

Behavior:

- Asserts `dataset` and `design` classes.
- Extracts `d <- design$train_design`, then:
  - `key <- parse_variable(key_var, d)`,
  - `key_fac <- factor(key)`.
- `distfun`:
  - If character: `distfun <- create_dist(distfun)`.
  - Must inherit `"distfun"`.
- Calls `create_model_spec("repmed_model", ...)` with:

  - `dataset`, `design`,
  - `key`           = `key_fac`,
  - `items`         = `repmed_des$items`,
  - `X_rdm`         = `repmed_des$X_rdm`,
  - `Y_rdm`         = `repmed_des$Y_rdm`,
  - `confound_rdms` = `repmed_des$confound_rdms`,
  - `distfun`,
  - `compute_performance = FALSE`,
  - `return_predictions  = FALSE`,
  - plus `...`.

Again, ROI‑level mediation metrics are supplied by `process_roi.repmed_model()`.

#### 8.2.3 `process_roi.repmed_model()`: per‑ROI mediation

Signature:

```r
process_roi.repmed_model <- function(mod_spec,
                                     roi,
                                     rnum,
                                     center_global_id = NA,
                                     ...)
```

Core steps:

1. **Extract ROI data**
   - `Xmat <- as.matrix(neuroim2::values(roi$train_roi))`,
   - `ind  <- neuroim2::indices(roi$train_roi)`.
   - If `n_obs < 2` or `n_vox < 1`: return error tibble.
2. **Check key length**
   - `key <- mod_spec$key`; require `length(key) == n_obs`.
3. **Item‑level mediator prototypes**
   - `M_full <- group_means(Xmat, margin = 1, group = key)`.
   - `items_roi <- rownames(M_full)`.
4. **Align items across X, M, Y**
   - `items <- mod_spec$items`,
   - `items_X <- rownames(mod_spec$X_rdm)`,
   - `items_Y <- rownames(mod_spec$Y_rdm)`.
   - `common_items <- Reduce(intersect, list(items, items_roi, items_X, items_Y))`.
   - Require `K = length(common_items) >= 3`.
   - Sort and subset:
     - `M  <- M_full[common_items, , drop = FALSE]`,
     - `Xr <- mod_spec$X_rdm[common_items, common_items, drop = FALSE]`,
     - `Yr <- mod_spec$Y_rdm[common_items, common_items, drop = FALSE]`.
5. **Drop constant voxels**
   - `sd_M <- apply(M, 2, sd, na.rm = TRUE)`; keep `> 0`.
   - If none: error tibble.
6. **Compute mediator RDM and vectorize**
   - `D_M <- pairwise_dist(mod_spec$distfun, M)`.
   - Lower‑tri vectors:
     - `x <- as.numeric(Xr[lower.tri(Xr)])`,
     - `y <- as.numeric(Yr[lower.tri(Yr)])`,
     - `m <- as.numeric(D_M[lower.tri(D_M)])`.
7. **Confound vectors (optional)**
   - If `confound_rdms` present:
     - For each `Cr`:
       - `Csub <- Cr[common_items, common_items, drop = FALSE]`,
       - vectorize lower tri.
     - Bind into `conf_mat` with column names = names of `confound_rdms`.
8. **Build data frames for mediation paths**
   - Path a (`m ~ x (+ conf)`):
     - `df_a <- data.frame(m = m, x = x)` or with `conf_mat` columns.
   - Paths b/c′ (`y ~ m + x (+ conf)`):
     - `df_b <- data.frame(y = y, m = m, x = x)` or plus confounds.
9. **Fit regression models**
   - **Path a:** `a_fit <- lm(m ~ ., data = df_a)`; coefficient:
     - `a = coef(summary(a_fit))["x", "Estimate"]` if available.
   - **Paths b/c′:** `b_fit <- lm(y ~ ., data = df_b)`; coefficients:
     - `b      = coef(summary(b_fit))["m", "Estimate"]`,
     - `cprime = coef(summary(b_fit))["x", "Estimate"]`.
   - Any failure yields `NA_real_` values.
10. **Indirect and Sobel**
    - `indirect <- a * b`.
    - Sobel z/p (if both fits succeed and required coefficients exist):
      - `sa <- cf_a["x", "Std. Error"]`,
      - `sb <- cf_b["m", "Std. Error"]`,
      - `se_ind <- sqrt(b^2 * sa^2 + a^2 * sb^2)`,
      - `sobel_z <- indirect / se_ind`,
      - `sobel_p <- 2 * pnorm(-abs(sobel_z))`.
11. **Performance vector**

    ```r
    perf <- c(
      n_items      = K,
      n_pairs      = length(x),
      med_a        = a,
      med_b        = b,
      med_cprime   = cprime,
      med_indirect = indirect,
      med_sobel_z  = sobel_z,
      med_sobel_p  = sobel_p
    )
    ```

12. **Return tibble row**

    - `performance = list(perf)`,
    - plus ROI indices and standard error/warning fields.

---

### 8.3 Representational Mapping – `repmap_model` (ReNA‑Map)

**Goal:** for each ROI/searchlight, fit a **linear reduced‑rank mapping** from **seed features** (e.g., model embeddings or patterns from another region) to ROI patterns. Summarize each ROI by mapping rank, fit (R²), and mapping norm.

#### 8.3.1 `repmap_design()`: mapping design helper

Signature:

```r
repmap_design <- function(items,
                          seed_features)
```

Behavior:

- Coerces `items` to character.
- Coerces `seed_features` to matrix `F`.
- If `rownames(F)` is `NULL`, sets it to `items[seq_len(nrow(F))]`.
- Returns:

  - `items` – canonical item IDs,
  - `seed_features` – `K × P` matrix, rows labeled by items.

#### 8.3.2 `repmap_model()`: mapping model spec

Signature:

```r
repmap_model <- function(dataset,
                         design,
                         repmap_des,
                         key_var,
                         rank = "auto",
                         max_rank = 20,
                         ridge_lambda = NULL,
                         ...)
```

Behavior:

- Asserts `dataset` and `design` classes.
- Extracts `d <- design$train_design`, then:
  - `key <- parse_variable(key_var, d)`,
  - `key_fac <- factor(key)`.
- Calls `create_model_spec("repmap_model", ...)` with:

  - `dataset`, `design`,
  - `key`           = `key_fac`,
  - `items`         = `repmap_des$items`,
  - `seed_features` = `repmap_des$seed_features`,
  - `rank`, `max_rank`, `ridge_lambda`,
  - `compute_performance = FALSE`,
  - `return_predictions  = FALSE`,
  - plus `...`.

#### 8.3.3 `.repmap_fit_rrr()`: internal reduced‑rank regression helper

Signature:

```r
.repmap_fit_rrr <- function(X, Y,
                            rank = "auto",
                            max_rank = 20,
                            ridge_lambda = NULL)
```

Behavior:

- Requires `rrpack`; if not available, returns a **zero mapping**:

  - `C` = zero `p × q` matrix,
  - `rank = 0L`,
  - `singvals = numeric(0)`.

- Computes `n`, `p`, `q`; sets:

  - `max_rank <- max(1L, min(max_rank, p, q, n - 1L))`.

- Handles `rank` argument:
  - If numeric and `<= 0`: zero mapping (rank 0).
  - If `ridge_lambda` not `NULL`:
    - Uses `rrpack::rrs.fit()` with `nrank` set to either `max_rank` (for `"auto"`) or fixed integer rank.
  - Else if `rank == "auto"`:
    - Uses `rrpack::cv.rrr()` with:
      - `nfold = min(5, max(2, floor(n/3)))`,
      - `maxrank = max_rank`.
  - Else:
    - Uses `rrpack::rrr.fit()` with `nrank = rank`.

- Extracts mapping matrix and singular values:
  - `C_hat` from `fit$coef` if present, otherwise `rrpack::coef(fit)`; wrapped in `tryCatch` with zero fallback.
  - `singvals` from `fit$Ad` (`numeric(0)` on error).
  - `r_used` from `fit$rank` if available, else `length(singvals)`, else `NA`.

- Returns:

  - `C`       = `C_hat`,
  - `rank`    = `r_used`,
  - `singvals` = `singvals`.

#### 8.3.4 `process_roi.repmap_model()`: per‑ROI mapping metrics

Signature:

```r
process_roi.repmap_model <- function(mod_spec,
                                     roi,
                                     rnum,
                                     center_global_id = NA,
                                     ...)
```

Core steps:

1. **Extract ROI data**
   - `Ymat <- as.matrix(neuroim2::values(roi$train_roi))`,
   - `ind  <- neuroim2::indices(roi$train_roi)`.
   - If `n_obs < 2` or `n_vox < 1`: error tibble.
2. **Check key length**
   - `key <- mod_spec$key`; require `length(key) == n_obs`.
3. **Item‑level ROI patterns**
   - `Y_full <- group_means(Ymat, margin = 1, group = key)`.
   - `items_roi <- rownames(Y_full)`.
4. **Align with seed features**
   - `items <- mod_spec$items`,
   - `items_seed <- rownames(mod_spec$seed_features)`.
   - `common_items <- Reduce(intersect, list(items, items_roi, items_seed))`.
   - Require `K >= 3`.
   - Sort and subset:
     - `X <- mod_spec$seed_features[common_items, , drop = FALSE]`,
     - `Y <- Y_full[common_items, , drop = FALSE]`.
5. **Drop constant voxels/features**
   - In `Y`: keep voxels with `sd_Y > 0`; else error.
   - In `X`: keep features with `sd_X > 0`; else error.
6. **Center X and Y**
   - `Xc <- scale(X, center = TRUE, scale = FALSE)`,
   - `Yc <- scale(Y, center = TRUE, scale = FALSE)`.
7. **Fit reduced‑rank mapping**
   - `fit <- .repmap_fit_rrr(Xc, Yc, rank = mod_spec$rank, max_rank = mod_spec$max_rank, ridge_lambda = mod_spec$ridge_lambda)`.
   - Extract:
     - `C_hat`, `r_used`, `svals`.
8. **In‑sample mapping fit (per‑voxel R²)**
   - `Y_hat <- Xc %*% C_hat`,
   - `resid <- Yc - Y_hat`.
   - `ss_tot <- colSums(Yc^2)`; `ss_res <- colSums(resid^2)`.
   - Guard `ss_tot` against zero with `.Machine$double.eps`.
   - `r2_vox <- 1 - ss_res / ss_tot`.
   - Summaries:
     - `mapping_r2_mean = mean(r2_vox, na.rm = TRUE)`,
     - `mapping_r2_med  = median(r2_vox, na.rm = TRUE)`.
9. **Mapping norm and singular values**
   - `frob_norm = sqrt(sum(C_hat^2))`.
   - `map_sv1` = first singular value (or `NA`).
   - `map_sv_mean` = mean of singular values (or `NA`).
10. **Performance vector**

    ```r
    perf <- c(
      n_items        = K,
      n_seed_feats   = ncol(X),
      n_vox          = ncol(Y),
      map_rank       = as.numeric(r_used),
      map_sv1        = if (length(svals)) svals[1] else NA_real_,
      map_sv_mean    = if (length(svals)) mean(svals) else NA_real_,
      map_r2_mean    = mapping_r2_mean,
      map_r2_median  = mapping_r2_med,
      map_frob_norm  = frob_norm
    )
    ```

11. **Return tibble row**

    - `performance = list(perf)`,
    - `indices = list(ind)`,
    - `result = list(NULL)`,
    - standard `error`/`warning` fields.

---

### 8.4 Integration & Combiners

- All three models (`repnet_model`, `repmed_model`, `repmap_model`) are **S3 model specs** with `process_roi.*` methods.
- They integrate with existing infrastructure:
  - `mvpa_iterate()` / `run_searchlight()` / `run_regional()` can dispatch on their classes similarly to `rsa_model`, `vector_rsa_model`, `contrast_rsa_model`, etc.
- They do **not** rely on generic cross‑validation helpers:
  - Per‑ROI metrics are computed directly inside `process_roi.*` and returned in the `performance` list‑column.
- Combiners can:
  - Reuse `combine_rsa_standard()` where appropriate, or
  - Implement a simple `combine_repnet_standard()` pattern:

    ```r
    perf_mat <- do.call(rbind, good_results$performance)
    wrap_out(perf_mat, model_spec$dataset, unlist(good_results$id))
    ```

Once wired, rMVPA will expose:

- **`repnet_model` – ReNA‑RC connectivity maps**  
  “Where does the geometry look like the seed, with/without confounds?”
- **`repmed_model` – ReNA‑RM mediation maps**  
  “Where does this region’s geometry mediate model → retrieval/behavior?”
- **`repmap_model` – ReNA‑Map mapping maps**  
  “Where does this region implement a strong low‑rank mapping from seed features to neural patterns?”

All within the existing rMVPA framework, in R, and in an RDM‑native, linear‑model aesthetic that’s deliberately **slicker and more coherent** than what’s typically available in Python toolboxes today.

