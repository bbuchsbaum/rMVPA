# Construct an hrfdecoder design (explicit mvpa_design extension)

Creates a specialized design object for continuous-time decoding that
carries event metadata in addition to the fields required by
\`mvpa_design\`.

## Usage

``` r
hrfdecoder_design(event_model, events, block_var, split_by = NULL)
```

## Arguments

- event_model:

  An event model object (e.g., from fmridesign::event_model). This must
  include a sampling_frame attribute defining the TR structure.

- events:

  A data.frame/tibble of trial/event rows with event-level labels. Must
  contain columns: \`onset\` (event onset times in seconds),
  \`condition\` (factor or character labels).

- block_var:

  Integer/factor vector of length T (TRs) specifying run/block ids.

- split_by:

  Optional formula or vector to define split groups (passed to
  mvpa_design).

## Value

An object inheriting from mvpa_design with explicit class
\`hrfdecoder_design\` and additional fields \`\$event_model\` and
\`\$events\`.

## Details

This constructor avoids adding ad-hoc members to a plain \`mvpa_design\`
by returning an object with class \`c("hrfdecoder_design",
"mvpa_design", "list")\`.

The hrfdecoder approach operates on continuous TR-level data rather than
trial-level beta estimates. \`cv_labels\` is set to TR indices (1:T) and
is used only for fold construction in the cross-validation machinery.
Critically, fold assignment is determined by \`block_var\` (run IDs),
NOT by \`cv_labels\` values. The actual decoding targets come from the
\`event_model\` and \`events\` fields (passed as \`targets\`), and
\`train_model.hrfdecoder_model()\` ignores the \`cv_labels\` parameter
entirely.

## Examples

``` r
if (FALSE) { # \dontrun{
library(fmridesign)
library(fmrihrf)

# Create event table
events_df <- data.frame(
  onset = c(5, 15, 35, 45),
  condition = factor(c("A", "B", "A", "B")),
  run = c(1, 1, 2, 2)
)

# Define sampling frame (TR structure)
sframe <- sampling_frame(blocklens = c(60, 60), TR = 2)

# Build event model
evmod <- event_model(
  onset ~ hrf(condition, basis = "spmg1"),
  data = events_df,
  block = ~run,
  sampling_frame = sframe
)

# Create hrfdecoder design
block_var <- rep(1:2, each = 60)  # 60 TRs per run
design <- hrfdecoder_design(
  event_model = evmod,
  events = events_df,
  block_var = block_var
)
} # }
```
