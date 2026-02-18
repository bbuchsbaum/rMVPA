# ITEM vs hrfdecoder in rMVPA

This note clarifies how the two event-level decoding paths differ and when to use each.

## Shared goal

Both approaches decode event/condition information while accounting for overlap in BOLD responses.

## ITEM path (`item_design()` + `item_model()`)

- Pipeline: estimate trial-wise responses first, then decode.
- Stage 1: LS-A estimation of trial responses (`fmrilss::lsa`).
- Stage 2: covariance-aware ITEM decoding (`fmrilss::item_cv`).
- Main strengths:
  - explicit trial-response matrix (`Gamma`) for diagnostics;
  - direct control over trial covariance construction (`U`, `u_storage`, `ridge_u`);
  - explicit trial alignment controls (`trial_id`, `trial_hash`, `check_hash`).

## hrfdecoder path (`hrfdecoder_design()` + `hrfdecoder_model()`)

- Pipeline: fit directly in TR space, then aggregate to events.
- Stage 1: continuous-time decoder fit on TR-level observations.
- Stage 2: aggregate TR-level soft labels to event-level probabilities.
- Main strengths:
  - avoids separate trial-estimation step;
  - naturally expressed in TR-domain with temporal smoothness/HRF-conformity penalties.

## Practical selection guidance

- Choose ITEM when trial-level interpretability and trial-covariance control are primary.
- Choose hrfdecoder when you prefer a direct TR-level model with post-hoc event aggregation.

## Reliability contracts now covered by tests

- ITEM:
  - covariance contract checks (`cov` vs `precision`, dense vs sparse);
  - fold slicing consistency (`U` matrix vs by-run blocks);
  - deterministic fold/tie behavior;
  - null/signal simulation behavior;
  - ROI mismatch guards and prediction payload behavior.
- hrfdecoder:
  - design validation (required columns, run/block checks, event-range checks);
  - backend/runtime guards;
  - merge behavior for zero-probability events;
  - end-to-end small-data training smoke test.
