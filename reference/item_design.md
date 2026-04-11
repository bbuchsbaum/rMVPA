# Construct an ITEM design for trial-wise decoding

Creates a design object for ITEM-style decoding where the trial-wise
design matrix (\`X_t\`) and supervised trial targets (\`T_target\`) are
decoupled from TR-level observations used in searchlight/ROI extraction.

## Usage

``` r
item_design(
  train_design = NULL,
  X_t,
  T_target,
  run_id,
  C_transform = NULL,
  V = NULL,
  v_type = c("cov", "precision"),
  trial_id = NULL,
  trial_hash = NULL,
  split_by = NULL,
  Z = NULL,
  Nuisance = NULL,
  meta = list(),
  ...
)
```

## Arguments

- train_design:

  Data frame/tibble describing TR-level observations. If \`NULL\`, a
  minimal design with one row per TR is created.

- X_t:

  Numeric trial-wise design matrix (\`n_time x n_trials\`).

- T_target:

  Trial-level supervised targets (matrix/data.frame/vector/factor).
  Length/rows must equal \`n_trials\`.

- run_id:

  Trial-level run/session identifiers (length \`n_trials\`).

- C_transform:

  Optional trial-to-condition transform matrix.

- V:

  Optional temporal covariance/precision passed to
  \`fmrilss::item_compute_u()\`.

- v_type:

  One of \`"cov"\` or \`"precision"\`.

- trial_id:

  Optional trial identifiers.

- trial_hash:

  Optional trial hash used by alignment guards.

- split_by:

  Optional split variable passed to \`mvpa_design()\`.

- Z:

  Optional nuisance matrix (\`n_time x q\`) used by LS-A
  (\`fmrilss::lsa\`).

- Nuisance:

  Alias nuisance matrix for LS-A; if both provided, \`Z\` is used.

- meta:

  Optional metadata list.

- ...:

  Reserved for forward compatibility.

## Value

An object of class \`c("item_design", "mvpa_design", "list")\`.
