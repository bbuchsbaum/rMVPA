#!/usr/bin/env python3
"""Benchmark vendored rsatoolbox crossnobis on deterministic fixtures."""

from __future__ import annotations

import argparse
import gc
import json
import sys
import time
import types
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in cast",
    category=RuntimeWarning,
)


def install_cengine_stub() -> None:
    """Allow balanced calc_rdm_crossnobis import when the Cython engine is unbuilt."""
    sim_mod = types.ModuleType("rsatoolbox.cengine.similarity")
    sim_mod.similarity = lambda *args, **kwargs: None
    sim_mod.calc_one = lambda *args, **kwargs: None
    sim_mod.calc = lambda *args, **kwargs: None

    cengine_mod = types.ModuleType("rsatoolbox.cengine")
    cengine_mod.similarity = sim_mod.similarity
    cengine_mod.calc_one = sim_mod.calc_one
    cengine_mod.calc = sim_mod.calc

    sys.modules.setdefault("rsatoolbox.cengine", cengine_mod)
    sys.modules.setdefault("rsatoolbox.cengine.similarity", sim_mod)


def load_rsatoolbox(vendor_src: str | None):
    if vendor_src:
        sys.path.insert(0, vendor_src)
    install_cengine_stub()
    from rsatoolbox.data import Dataset
    from rsatoolbox.rdm import calc_rdm_crossnobis

    return Dataset, calc_rdm_crossnobis


def make_measurements(k_cond: int, n_voxels: int, n_folds: int, n_reps: int):
    vox = np.arange(1, n_voxels + 1, dtype=float)
    n_obs = k_cond * n_folds * n_reps
    measurements = np.empty((n_obs, n_voxels), dtype=float)
    conds = []
    folds = []

    row = 0
    for fold in range(1, n_folds + 1):
        for rep in range(1, n_reps + 1):
            for cond in range(1, k_cond + 1):
                measurements[row, :] = (
                    np.sin(0.07 * cond * vox + 0.11 * fold)
                    + np.cos(0.03 * rep * vox + 0.17 * cond)
                    + 0.02 * fold
                    - 0.015 * cond
                )
                conds.append(f"C{cond:03d}")
                folds.append(fold)
                row += 1

    return measurements, np.asarray(conds), np.asarray(folds)


def make_fold_means(measurements: np.ndarray, k_cond: int, n_voxels: int,
                    n_folds: int, n_reps: int) -> np.ndarray:
    return measurements.reshape(n_folds, n_reps, k_cond, n_voxels).mean(axis=1).transpose(1, 2, 0)


def crossnobis_numpy_formula(u_folds: np.ndarray) -> np.ndarray:
    k_cond, n_voxels, n_folds = u_folds.shape
    tri = np.triu_indices(k_cond, k=1)
    u_sum = np.sum(u_folds, axis=2)
    rdms = np.empty((n_folds, int(k_cond * (k_cond - 1) / 2)), dtype=float)

    for fold in range(n_folds):
        test = u_folds[:, :, fold]
        train = (u_sum - test) / (n_folds - 1)
        kernel = train @ test.T
        dmat = np.expand_dims(np.diag(kernel), 0) + np.expand_dims(np.diag(kernel), 1) - kernel - kernel.T
        rdms[fold, :] = dmat[tri] / n_voxels

    return np.mean(rdms, axis=0)


def summarize(times: list[float]) -> dict[str, float]:
    arr = np.asarray(times, dtype=float) * 1000.0
    return {
        "median_ms": float(np.median(arr)),
        "min_ms": float(np.min(arr)),
        "max_ms": float(np.max(arr)),
        "mean_ms": float(np.mean(arr)),
    }


def benchmark_case(case: dict, dataset_cls, calc_rdm_crossnobis, times: int, warmup: int):
    k_cond = int(case["K"])
    n_voxels = int(case["V"])
    n_folds = int(case["M"])
    n_reps = int(case["reps"])
    case_id = case["case_id"]

    measurements, conds, folds = make_measurements(k_cond, n_voxels, n_folds, n_reps)
    dataset = dataset_cls(
        measurements=measurements,
        obs_descriptors={"cond": conds, "fold": folds},
    )
    u_folds = make_fold_means(measurements, k_cond, n_voxels, n_folds, n_reps)

    def run_public_once():
        rdm = calc_rdm_crossnobis(dataset, descriptor="cond", cv_descriptor="fold")
        return np.asarray(rdm.dissimilarities).reshape(-1)

    def run_formula_once():
        return crossnobis_numpy_formula(u_folds)

    distances_public = run_public_once()
    distances_formula = run_formula_once()
    for _ in range(warmup):
        run_public_once()
        run_formula_once()

    elapsed_public = []
    elapsed_formula = []
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        for _ in range(times):
            start = time.perf_counter()
            run_public_once()
            elapsed_public.append(time.perf_counter() - start)

            start = time.perf_counter()
            run_formula_once()
            elapsed_formula.append(time.perf_counter() - start)
    finally:
        if was_enabled:
            gc.enable()

    common = {
        "case_id": case_id,
        "K": k_cond,
        "V": n_voxels,
        "M": n_folds,
        "reps": n_reps,
        "iterations": times,
    }

    public = dict(common)
    public.update({
        "method": "rsatoolbox_public",
        "dissimilarities": distances_public.tolist(),
    })
    public.update(summarize(elapsed_public))

    formula = dict(common)
    formula.update({
        "method": "rsatoolbox_numpy_formula",
        "dissimilarities": distances_formula.tolist(),
    })
    formula.update(summarize(elapsed_formula))

    return [public, formula]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    with open(args.config, "r", encoding="utf-8") as handle:
        config = json.load(handle)

    vendor_src = config.get("vendor_src")
    if vendor_src is not None:
        vendor_src = str(Path(vendor_src).resolve())

    dataset_cls, calc_rdm_crossnobis = load_rsatoolbox(vendor_src)
    results = []
    for case in config["cases"]:
        results.extend(
            benchmark_case(case, dataset_cls, calc_rdm_crossnobis,
                           int(config["times"]), int(config["warmup"]))
        )

    payload = {
        "toolbox": "rsatoolbox",
        "vendor_src": vendor_src,
        "cengine_stubbed": True,
        "cases": results,
    }
    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(payload, handle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
