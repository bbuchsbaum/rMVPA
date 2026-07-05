#!/usr/bin/env python3
"""Benchmark vendored rsatoolbox RDM calculators on deterministic fixtures."""

from __future__ import annotations

import argparse
import gc
import json
import os
import sys
import time
import types
import warnings
from pathlib import Path
from tempfile import mkdtemp

import numpy as np

warnings.filterwarnings(
    "ignore",
    message="invalid value encountered in cast",
    category=RuntimeWarning,
)


def configure_writable_caches() -> None:
    """Keep optional plotting imports from writing under the user's home dir."""
    cache_root = Path(mkdtemp(prefix="rsatoolbox-rdm-cache-"))
    mpl_dir = cache_root / "matplotlib"
    xdg_dir = cache_root / "xdg"
    mpl_dir.mkdir(parents=True, exist_ok=True)
    xdg_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(xdg_dir))


def install_cengine_stub() -> None:
    """Allow rsatoolbox imports when the optional Cython engine is unbuilt."""
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
    configure_writable_caches()
    install_cengine_stub()
    from rsatoolbox.data import Dataset
    from rsatoolbox.rdm import calc_rdm

    return Dataset, calc_rdm


def make_measurements(k_cond: int, n_voxels: int, n_reps: int):
    vox = np.arange(1, n_voxels + 1, dtype=float)
    n_obs = k_cond * n_reps
    measurements = np.empty((n_obs, n_voxels), dtype=float)
    conds = []

    row = 0
    for rep in range(1, n_reps + 1):
        for cond in range(1, k_cond + 1):
            measurements[row, :] = (
                np.sin(0.07 * cond * vox + 0.11 * rep)
                + np.cos(0.03 * rep * vox + 0.17 * cond)
                + 0.02 * rep
                - 0.015 * cond
                + 0.001 * (cond % 3) * vox / n_voxels
            )
            conds.append(f"C{cond:03d}")
            row += 1

    return measurements, np.asarray(conds)


def make_condition_means(measurements: np.ndarray, k_cond: int, n_reps: int) -> np.ndarray:
    n_voxels = measurements.shape[1]
    return measurements.reshape(n_reps, k_cond, n_voxels).mean(axis=0)


def make_precision(n_voxels: int) -> np.ndarray:
    idx = np.arange(1, n_voxels + 1, dtype=float)
    diag_vals = 1.0 + (np.mod(idx, 7.0) / 10.0)
    u = np.sin(0.17 * idx)
    return np.diag(diag_vals) + 0.02 * np.outer(u, u) / n_voxels


def triu_vec(dmat: np.ndarray) -> np.ndarray:
    return dmat[np.triu_indices(dmat.shape[0], k=1)]


def rdm_correlation_formula(means: np.ndarray) -> np.ndarray:
    centered = means - means.mean(axis=1, keepdims=True)
    norms = np.sqrt(np.einsum("ij,ij->i", centered, centered))
    scaled = centered / norms[:, None]
    return triu_vec(1.0 - scaled @ scaled.T)


def rdm_euclidean_formula(means: np.ndarray) -> np.ndarray:
    sum_sq = np.sum(means**2, axis=1, keepdims=True)
    dmat = sum_sq + sum_sq.T - 2.0 * means @ means.T
    return triu_vec(dmat) / means.shape[1]


def rdm_mahalanobis_formula(means: np.ndarray, precision: np.ndarray) -> np.ndarray:
    kernel = means @ precision @ means.T
    diag = np.diag(kernel)
    dmat = diag[None, :] + diag[:, None] - 2.0 * kernel
    return triu_vec(dmat) / means.shape[1]


def summarize(times: list[float]) -> dict[str, float]:
    arr = np.asarray(times, dtype=float) * 1000.0
    return {
        "median_ms": float(np.median(arr)),
        "min_ms": float(np.min(arr)),
        "max_ms": float(np.max(arr)),
        "mean_ms": float(np.mean(arr)),
    }


def benchmark_case(case: dict, dataset_cls, calc_rdm, times: int, warmup: int):
    k_cond = int(case["K"])
    n_voxels = int(case["V"])
    n_reps = int(case["reps"])
    case_id = case["case_id"]

    measurements, conds = make_measurements(k_cond, n_voxels, n_reps)
    means = make_condition_means(measurements, k_cond, n_reps)
    precision = make_precision(n_voxels)
    dataset = dataset_cls(
        measurements=measurements,
        obs_descriptors={"cond": conds},
    )

    def public_correlation():
        rdm = calc_rdm(dataset, method="correlation", descriptor="cond")
        return np.asarray(rdm.dissimilarities).reshape(-1)

    def public_euclidean():
        rdm = calc_rdm(dataset, method="euclidean", descriptor="cond")
        return np.asarray(rdm.dissimilarities).reshape(-1)

    def public_mahalanobis():
        rdm = calc_rdm(
            dataset,
            method="mahalanobis",
            descriptor="cond",
            noise=precision,
        )
        return np.asarray(rdm.dissimilarities).reshape(-1)

    runners = [
        ("rsatoolbox_public", "correlation", public_correlation),
        ("rsatoolbox_public", "euclidean", public_euclidean),
        ("rsatoolbox_public", "mahalanobis", public_mahalanobis),
        ("rsatoolbox_numpy_formula", "correlation", lambda: rdm_correlation_formula(means)),
        ("rsatoolbox_numpy_formula", "euclidean", lambda: rdm_euclidean_formula(means)),
        (
            "rsatoolbox_numpy_formula",
            "mahalanobis",
            lambda: rdm_mahalanobis_formula(means, precision),
        ),
    ]

    distances = {
        (method, distance): runner()
        for method, distance, runner in runners
    }

    for _ in range(warmup):
        for _, _, runner in runners:
            runner()

    elapsed = {(method, distance): [] for method, distance, _ in runners}
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        for _ in range(times):
            for method, distance, runner in runners:
                start = time.perf_counter()
                runner()
                elapsed[(method, distance)].append(time.perf_counter() - start)
    finally:
        if was_enabled:
            gc.enable()

    rows = []
    common = {
        "case_id": case_id,
        "K": k_cond,
        "V": n_voxels,
        "reps": n_reps,
        "iterations": times,
    }
    for method, distance, _ in runners:
        row = dict(common)
        row.update({
            "method": method,
            "distance": distance,
            "dissimilarities": distances[(method, distance)].tolist(),
        })
        row.update(summarize(elapsed[(method, distance)]))
        rows.append(row)

    return rows


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

    dataset_cls, calc_rdm = load_rsatoolbox(vendor_src)
    results = []
    for case in config["cases"]:
        results.extend(benchmark_case(case, dataset_cls, calc_rdm,
                                      int(config["times"]), int(config["warmup"])))

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
