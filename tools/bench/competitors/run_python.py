"""Competitor side of the rMVPA benchmark harness.

Run from the package root, after `Rscript tools/bench/make_data.R`:

    OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
      tools/bench/.venv/bin/python tools/bench/competitors/run_python.py [reps] [groups]

groups: comma-separated subset of sl,regional,rsa (default: all).

Appends one JSON line per (scenario, method) to
tools/bench/receipts/<date>-python.jsonl. Every run is single-threaded.
`digest` summarises the outputs, so a timing is only compared with an rMVPA
timing that computed the same thing.
"""

import json
import os
import platform
import statistics
import subprocess
import sys
import time
from datetime import datetime, timezone
from importlib.metadata import version

import numpy as np
import nibabel as nib
import pandas as pd
from threadpoolctl import threadpool_limits
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import LeaveOneGroupOut, cross_val_score
from sklearn.naive_bayes import GaussianNB
from nilearn.decoding import SearchLight
import rsatoolbox

ROOT = os.path.join("tools", "bench")
DATA = os.path.join(ROOT, "data")


class CorrelationClassifier(ClassifierMixin, BaseEstimator):
    """Pearson correlation to class-mean prototypes, argmax (first max wins).

    The same model as CoSMoMVPA's cosmo_classify_correlation and rMVPA's
    corclass (method = "pearson").
    """

    def fit(self, X, y):
        self.classes_ = np.unique(y)
        self.means_ = np.vstack([X[y == c].mean(axis=0) for c in self.classes_])
        return self

    def predict(self, X):
        def z(a):
            a = a - a.mean(axis=1, keepdims=True)
            n = np.linalg.norm(a, axis=1, keepdims=True)
            n[n == 0] = 1.0
            return a / n
        r = z(X) @ z(self.means_).T
        return self.classes_[np.argmax(r, axis=1)]


def estimators():
    return {
        "corclass": CorrelationClassifier(),
        "gaussian_nb": GaussianNB(),
        "lda_shrinkage": LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto"),
    }


def time_reps(fn, reps, inner=1):
    """Median-ready per-call times; `inner` repeats fn within each rep."""
    fn()  # warm-up
    times, out = [], None
    for _ in range(reps):
        t0 = time.perf_counter()
        for _ in range(inner):
            out = fn()
        times.append((time.perf_counter() - t0) / inner)
    return times, out


def receipt(scenario, method, tool, times, n_units, unit, digest, notes=""):
    q = np.percentile(times, [25, 75])
    return {
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tool": tool,
        "tool_versions": {p: version(p) for p in
                          ["numpy", "scipy", "scikit-learn", "nilearn", "rsatoolbox"]},
        "python": platform.python_version(),
        "blas": np.show_config(mode="dicts")["Build Dependencies"]["blas"]["name"],
        "threads": 1,
        "cpu": platform.processor() or platform.machine(),
        "scenario": scenario,
        "method": method,
        "reps": len(times),
        "median_s": statistics.median(times),
        "iqr_s": float(q[1] - q[0]),
        "n_units": n_units,
        "unit": unit,
        "ms_per_unit": 1000 * statistics.median(times) / n_units,
        "digest": digest,
        "notes": notes,
    }


def run_searchlight(reps):
    bold = nib.load(os.path.join(DATA, "synth12_bold.nii.gz"))
    mask = nib.load(os.path.join(DATA, "synth12_mask.nii.gz"))
    design = pd.read_csv(os.path.join(DATA, "synth12_design.csv"))
    y, runs = design["label"].to_numpy(), design["run"].to_numpy()
    n_centres = int(np.asarray(mask.dataobj).astype(bool).sum())
    rows = []
    for name, est in estimators().items():
        def fit():
            sl = SearchLight(mask_img=mask, process_mask_img=mask, radius=3.0,
                             estimator=est, cv=LeaveOneGroupOut(), scoring="accuracy",
                             n_jobs=1, verbose=0)
            sl.fit(bold, y, groups=runs)
            return sl.scores_
        times, scores = time_reps(fit, reps)
        vals = scores[np.asarray(mask.dataobj).astype(bool)]
        rows.append(receipt("sl_synth12_r3", name, "nilearn", times, n_centres, "centre",
                            {"mean_score": float(np.nanmean(vals))},
                            "nilearn SearchLight; mean of per-fold accuracies"))
    return rows


def run_regional(reps):
    X = pd.read_csv(os.path.join(DATA, "haxby_patterns.csv")).to_numpy(dtype=float)
    design = pd.read_csv(os.path.join(DATA, "haxby_design.csv"))
    y, runs = design["label"].to_numpy(), design["run"].to_numpy()
    rows = []
    for name, est in estimators().items():
        def fit():
            return cross_val_score(est, X, y, groups=runs, cv=LeaveOneGroupOut(),
                                   scoring="accuracy", n_jobs=1)
        times, scores = time_reps(fit, reps)
        rows.append(receipt("regional_haxby_vt", name, "scikit-learn", times, 1, "roi",
                            {"accuracy": float(np.mean(scores))},
                            "LeaveOneGroupOut over 12 runs; 577 voxels"))
    return rows


def run_rsa(reps):
    X = pd.read_csv(os.path.join(DATA, "haxby_patterns.csv")).to_numpy(dtype=float)
    design = pd.read_csv(os.path.join(DATA, "haxby_design.csv"))
    data = rsatoolbox.data.Dataset(
        X, obs_descriptors={"conds": design["label"].to_numpy(),
                            "runs": design["run"].to_numpy()})
    rows = []

    def corr():
        return rsatoolbox.rdm.calc_rdm(data, method="correlation", descriptor="conds")
    times, rdm = time_reps(corr, reps, inner=500)
    rows.append(receipt("rsa_haxby", "rdm_correlation_condmeans", "rsatoolbox", times, 1, "rdm",
                        {"sum": float(rdm.get_vectors().sum())},
                        "8 condition means over 577 voxels"))

    def xnobis():
        return rsatoolbox.rdm.calc_rdm(data, method="crossnobis", descriptor="conds",
                                       cv_descriptor="runs")
    times, rdm = time_reps(xnobis, reps, inner=100)
    rows.append(receipt("rsa_haxby", "rdm_crossnobis_identity", "rsatoolbox", times, 1, "rdm",
                        {"sum": float(rdm.get_vectors().sum())},
                        "identity noise; runs as cv folds"))
    return rows


def run_rsa_searchlight(reps):
    """rsatoolbox searchlight: same 20 condition patterns and model RDM as rMVPA.

    Spheres from get_volume_searchlight(radius=3, threshold=1.0), RDMs from
    get_searchlight_RDMs(method='correlation'), model fit with the vectorised
    rsatoolbox.rdm.compare(method='corr') over all searchlight RDMs (faster than
    evaluate_models_searchlight's per-centre loop).
    """
    from rsatoolbox.util.searchlight import get_volume_searchlight, get_searchlight_RDMs
    from rsatoolbox.rdm import RDMs, compare
    import io, contextlib
    bold = np.asarray(nib.load(os.path.join(DATA, "synth12_bold.nii.gz")).dataobj, dtype=float)
    mask = np.asarray(nib.load(os.path.join(DATA, "synth12_mask.nii.gz")).dataobj).astype(bool)
    data = bold[..., :20]
    data_2d = data.reshape(-1, 20).T  # C-order voxel index, as ravel_multi_index
    model = RDMs(pd.read_csv(os.path.join(DATA, "rsa_model_rdm.csv"))["model"].to_numpy()[None, :])
    events = np.arange(20)

    def run():
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            centers, neighbors = get_volume_searchlight(mask, radius=3, threshold=1.0)
            sl_rdms = get_searchlight_RDMs(data_2d, centers, neighbors, events, method="correlation", verbose=False)
            r = compare(model, sl_rdms, method="corr").ravel()
        return centers, r
    times, (centers, r) = time_reps(run, reps)
    xyz = np.array(np.unravel_index(centers, mask.shape)).T
    pd.DataFrame({"x": xyz[:, 0], "y": xyz[:, 1], "z": xyz[:, 2], "r": r}).to_csv(
        os.path.join(DATA, "rsa_sl_rsatoolbox_values.csv"), index=False)
    return [receipt("rsa_sl_synth12_r3", "rdm_corr_pearson_fit", "rsatoolbox", times, int(mask.sum()), "centre",
                    {"mean_r": float(np.nanmean(r))},
                    "get_volume_searchlight r=3 + get_searchlight_RDMs + vectorised compare")]


def main():
    reps = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    groups = sys.argv[2].split(",") if len(sys.argv) > 2 else ["sl", "regional", "rsa", "rsa_sl"]
    os.makedirs(os.path.join(ROOT, "receipts"), exist_ok=True)
    sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    path = os.path.join(ROOT, "receipts", f"{datetime.now():%Y-%m-%d}-python.jsonl")
    with threadpool_limits(limits=1):
        rows = []
        if "rsa" in groups:
            rows += run_rsa(max(reps, 20))
        if "regional" in groups:
            rows += run_regional(max(reps, 10))
        if "sl" in groups:
            rows += run_searchlight(reps)
        if "rsa_sl" in groups:
            rows += run_rsa_searchlight(reps)
    with open(path, "a") as fh:
        for r in rows:
            r["rmvpa_git_sha"] = sha
            fh.write(json.dumps(r) + "\n")
            print(f"{r['scenario']:<18} {r['method']:<28} {r['tool']:<12} "
                  f"median={r['median_s']:.4f}s  {r['ms_per_unit']:.3f} ms/{r['unit']}  {r['digest']}")


if __name__ == "__main__":
    main()
