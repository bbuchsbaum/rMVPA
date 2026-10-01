# Benchmark harness: rMVPA vs competitors

Head-to-head timings on identical inputs, with an output digest for each run
so we only compare runs that computed the same thing. Results are in
[`SCOREBOARD.md`](SCOREBOARD.md). The plan and its rules are in
`adocs/cran-release-plan.md` (Phase 0 harness, Phase 2 hill-climb).

## Layout

| Path | Purpose |
|---|---|
| `make_data.R` | Writes shared inputs to `data/` (gitignored): a synthetic 12³ 4D NIfTI volume with mask and design, and the Haxby VT block patterns as CSV |
| `run_rmvpa.R` | rMVPA side. Appends receipts to `receipts/<date>-rmvpa.jsonl` |
| `competitors/run_python.py` | nilearn / scikit-learn / rsatoolbox side. Appends to `receipts/<date>-python.jsonl` |
| `competitors/requirements.{in,lock}` | Pinned competitor environment |
| `scoreboard.R` | Collates the latest receipts into `SCOREBOARD.md` |
| `receipts/` | Append-only JSON-lines records, committed |

## Running

From the package root:

```sh
# one-time: competitor environment (isolated; touches neither R nor system Python)
uv venv --python 3.12 tools/bench/.venv
VIRTUAL_ENV=tools/bench/.venv uv pip install -r tools/bench/competitors/requirements.lock

Rscript tools/bench/make_data.R
Rscript tools/bench/run_rmvpa.R 3                 # reps; optional groups: sl,regional,rsa
OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  tools/bench/.venv/bin/python tools/bench/competitors/run_python.py 3
Rscript tools/bench/scoreboard.R
```

A full rMVPA run on master code takes about 25 minutes because the general
searchlight path is slow. For jobs that long, launch with
`nohup ... > log 2>&1 &` and check the process with `pgrep`.

## Rules

- **Same model, same data, same folds.** Every scenario uses
  leave-one-run-out folds and the same neighbourhoods (radius in mm, 1 mm
  voxels). The scoreboard shows a ratio only when the output digests agree.
  Comparable-but-different estimators, such as `dual_lda` vs sklearn
  Ledoit-Wolf LDA, are shown as `n/a` until parity is established.
- **Single thread, medians.** Every receipt records reps, median and IQR.
  Operations below the timer resolution loop internally (`inner`).
- **Record the environment.** Each receipt names its BLAS. Our R uses the
  reference BLAS while numpy uses Apple Accelerate, so BLAS-heavy rows favour
  the competitor. Report both R BLAS configurations before claiming a win on
  those rows.
- **Fixed seeds.** rMVPA breaks predicted-class ties at random (release-plan
  finding C2), so `run_rmvpa.R` reseeds before every rep. Receipts recorded
  before that change may differ in the 5th or 6th decimal of a digest for
  that reason alone.
- **Append, never edit.** New runs add receipts. Losses stay on the
  scoreboard until a later receipt beats them.

## Not covered

- **PyMVPA:** 2.6.5 needs `numpy.distutils`, so it cannot be installed with
  numpy ≥ 1.26 / Python 3.12.
- **CoSMoMVPA:** needs MATLAB or Octave, which are not available here.
- **Linear SVM:** rMVPA's `svmLinear` needs e1071, which is not installed
  here, and the planned LinearSVC-equivalent model does not exist yet.
- **Peak memory:** not yet recorded.
- **Byte compilation:** runners use `pkgload::load_all`, which does not
  byte-compile. Some per-call cost is JIT compilation an installed package
  would not pay. Relative comparisons between rMVPA branches are fair;
  before quoting absolute numbers, run against an installed build in a
  temporary library.
