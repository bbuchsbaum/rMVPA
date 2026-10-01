# Audit receipts, 2026-09-30

Inputs behind the "measured" claims in `../../cran-release-plan.md` §0.
Base: master @ 8a13e6e. Machine: Apple M2 Pro, R 4.3.3, reference BLAS, sequential.
Every timing is a single run. Treat it as an order of magnitude, not a benchmark result.

| File | What it is | Author |
|---|---|---|
| `t3.R` | Legacy corclass searchlight on 20³ (572 s observed in the session console) | main session; later overwritten by an audit subagent, so the current content may differ |
| `t4.R` | Rprof of legacy corclass on 10³, plus sda_notune auto/legacy and svmLinear auto | main session |
| `profile_summary_latest_prof_out.txt` | `summaryRprof` of the last `prof.out` in scratch | regenerated; the last writer may have been a subagent run |
| `t1.R`, `t2.R`, `t5.R`, `t6.R`, `t2.out`, `rsa.out` | Timing and profile scripts and outputs | performance audit subagent |
| `R-CMD-check.log`, `testthat.Rout.fail.tail`, `cran_deps.R` | `--as-cran` check (no vignettes) and CRAN availability query | CRAN audit subagent |

The console output of the main-session runs (t3/t4) was not captured to file.
The figures quoted in the plan are transcribed from the session. Rerun the
scripts to reproduce them.
