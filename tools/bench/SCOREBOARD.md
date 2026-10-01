# Benchmark scoreboard

Generated 2026-10-01 07:18 from `tools/bench/receipts/` by `tools/bench/scoreboard.R`.
Single-threaded medians. **Ratio = rMVPA / competitor (< 1 means rMVPA is faster).**
A ratio is shown only when the two output digests agree to 1e-4 relative;
otherwise the methods are not computing the same thing and the row is marked.
Losses are recorded as losses.

## sl_synth12_r3 (ms per centre)

| method | rMVPA side | engine | rMVPA | competitor | tool | ratio | digests |
|---|---|---|---|---|---|---|---|
| corclass | bench/competitor-harness @8a13e6e | legacy | 32.8 | 2.5 | nilearn | **13.14** | match |
| gaussian_nb | bench/competitor-harness @8a13e6e | legacy | 44.0 | 4.5 | nilearn | **9.83** | match |
| lda_shrinkage | bench/competitor-harness @8a13e6e | dual_lda_fast | 2.2 | 17.6 | nilearn | n/a | differ (0.237436 vs 0.239508) |
| sda_notune | bench/competitor-harness @8a13e6e | legacy | 62.3 | — | — | — | rMVPA only |
| corclass | perf/logger-gating @fe364b8 | legacy | 21.1 | 2.5 | nilearn | **8.45** | match |
| gaussian_nb | perf/logger-gating @fe364b8 | legacy | 32.1 | 4.5 | nilearn | **7.17** | match |
| lda_shrinkage | perf/logger-gating @fe364b8 | dual_lda_fast | 2.1 | 17.6 | nilearn | n/a | differ (0.237436 vs 0.239508) |
| sda_notune | perf/logger-gating @fe364b8 | legacy | 50.7 | — | — | — | rMVPA only |
| corclass | perf/regional-scaffold @fe364b8 | legacy | 13.0 | 2.5 | nilearn | **5.20** | match |
| gaussian_nb | perf/regional-scaffold @fe364b8 | legacy | 24.1 | 4.5 | nilearn | **5.38** | match |
| lda_shrinkage | perf/regional-scaffold @fe364b8 | dual_lda_fast | 2.1 | 17.6 | nilearn | n/a | differ (0.237442 vs 0.239508) |
| sda_notune | perf/regional-scaffold @fe364b8 | legacy | 42.7 | — | — | — | rMVPA only |

## regional_haxby_vt (ms per roi)

| method | rMVPA side | engine | rMVPA | competitor | tool | ratio | digests |
|---|---|---|---|---|---|---|---|
| corclass | bench/competitor-harness @8a13e6e |  | 590.0 | 7.0 | scikit-learn | **83.88** | match |
| gaussian_nb | bench/competitor-harness @8a13e6e |  | 885.5 | 14.0 | scikit-learn | **63.17** | match |
| lda_shrinkage | bench/competitor-harness @8a13e6e |  | 1061.0 | 568.9 | scikit-learn | n/a | differ (0.927083 vs 0.885417) |
| sda_notune | bench/competitor-harness @8a13e6e |  | 880.5 | — | — | — | rMVPA only |
| corclass | perf/logger-gating @fe364b8 |  | 556.0 | 7.0 | scikit-learn | **79.05** | match |
| gaussian_nb | perf/logger-gating @fe364b8 |  | 853.0 | 14.0 | scikit-learn | **60.85** | match |
| lda_shrinkage | perf/logger-gating @fe364b8 |  | 1032.5 | 568.9 | scikit-learn | n/a | differ (0.927083 vs 0.885417) |
| sda_notune | perf/logger-gating @fe364b8 |  | 850.0 | — | — | — | rMVPA only |
| corclass | perf/regional-scaffold @fe364b8 |  | 73.0 | 7.0 | scikit-learn | **10.38** | match |
| gaussian_nb | perf/regional-scaffold @fe364b8 |  | 365.5 | 14.0 | scikit-learn | **26.08** | match |
| lda_shrinkage | perf/regional-scaffold @fe364b8 |  | 553.0 | 568.9 | scikit-learn | n/a | differ (0.927083 vs 0.885417) |
| sda_notune | perf/regional-scaffold @fe364b8 |  | 372.0 | — | — | — | rMVPA only |

## rsa_haxby (ms per rdm)

| method | rMVPA side | engine | rMVPA | competitor | tool | ratio | digests |
|---|---|---|---|---|---|---|---|
| rdm_correlation_condmeans | bench/competitor-harness @8a13e6e |  | 0.333 | 0.185 | rsatoolbox | **1.80** | match |
| rdm_crossnobis_identity | bench/competitor-harness @8a13e6e |  | 1.8 | 3.4 | rsatoolbox | **0.53** | match |
| rdm_correlation_condmeans | perf/logger-gating @fe364b8 |  | 0.346 | 0.185 | rsatoolbox | **1.87** | match |
| rdm_crossnobis_identity | perf/logger-gating @fe364b8 |  | 1.8 | 3.4 | rsatoolbox | **0.54** | match |
| rdm_correlation_condmeans | perf/regional-scaffold @fe364b8 |  | 0.334 | 0.185 | rsatoolbox | **1.81** | match |
| rdm_crossnobis_identity | perf/regional-scaffold @fe364b8 |  | 1.8 | 3.4 | rsatoolbox | **0.55** | match |

## Environment

- rMVPA: libRblas.0.dylib
- scikit-learn: accelerate
- nilearn: accelerate
- rsatoolbox: accelerate
- R uses the reference BLAS and numpy uses Apple Accelerate, so BLAS-heavy rows favour the competitor.
- PyMVPA 2.6.5 cannot be installed with numpy >= 1.26 / Python 3.12 (needs numpy.distutils); not benchmarked.
- CoSMoMVPA needs MATLAB/Octave; not available.
