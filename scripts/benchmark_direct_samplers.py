#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import io
import json
import pickle
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path
from typing import Iterable

import numpy as np
import scipy.sparse as sp


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from metacooc.format import format_data  # noqa: E402
from metacooc.null_models import make_null_sampler  # noqa: E402


MODELS = ("FE", "EF", "EE")
EE_VARIANTS = ("auto", "legacy", "oversample", "chunked", "numpy-choice", "floyd")
SCENARIOS = (
    "fixture",
    "ingredients",
    "ee-wide",
    "fe-k1",
    "fe-k2",
    "fe-low",
    "fe-medium",
    "fe-near-full",
    "fe-mixed",
    "all-synthetic",
    "all",
)


def _positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def _positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def _density(X: sp.spmatrix) -> float:
    n_rows, n_cols = map(int, X.shape)
    cells = n_rows * n_cols
    return 0.0 if cells == 0 else float(X.nnz) / float(cells)


def _shape_label(shape: tuple[int, int]) -> str:
    return f"{int(shape[0])}x{int(shape[1])}"


def _load_ingredients_matrix(path: Path) -> sp.csr_matrix:
    with path.open("rb") as handle:
        ingredients = pickle.load(handle)
    return ingredients.presence_matrix.tocsr()


def _fixture_matrix() -> sp.csr_matrix:
    fixture_dir = REPO_ROOT / "tests" / "data"
    with tempfile.TemporaryDirectory(prefix="metacooc-benchmark-fixture-") as tmp:
        with contextlib.redirect_stdout(io.StringIO()):
            format_data(
                tax_profile=str(fixture_dir / "tax_profile.tsv"),
                output_dir=tmp,
                sample_to_biome_file=str(fixture_dir / "sample_to_biome.tsv"),
                aggregated=True,
                tag="benchmark",
                data_version="benchmark",
            )
        return _load_ingredients_matrix(Path(tmp) / "ingredients_raw_benchmark.pkl")


def _matrix_from_row_degrees(degrees: np.ndarray, n_cols: int, seed: int) -> sp.csr_matrix:
    rng = np.random.default_rng(seed)
    degrees = np.asarray(degrees, dtype=np.int64)
    n_rows = int(degrees.size)
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []

    for row, k_raw in enumerate(degrees):
        k = int(k_raw)
        if k < 0 or k > int(n_cols):
            raise ValueError(f"invalid degree {k} for n_cols={n_cols}")
        if k == 0:
            continue
        if k == int(n_cols):
            draw = np.arange(int(n_cols), dtype=np.int64)
        else:
            draw = rng.choice(int(n_cols), size=k, replace=False).astype(np.int64, copy=False)
        rows.append(np.full(draw.size, row, dtype=np.int64))
        cols.append(draw)

    if not rows:
        return sp.csr_matrix((n_rows, int(n_cols)), dtype=np.int8)

    row_idx = np.concatenate(rows)
    col_idx = np.concatenate(cols)
    data = np.ones(row_idx.size, dtype=np.int8)
    X = sp.csr_matrix((data, (row_idx, col_idx)), shape=(n_rows, int(n_cols)))
    X.sum_duplicates()
    X.sort_indices()
    return X


def _matrix_from_col_degrees(degrees: np.ndarray, n_rows: int, seed: int) -> sp.csr_matrix:
    return _matrix_from_row_degrees(degrees, n_rows, seed).T.tocsr()


def _mixed_degrees(n_entities: int, n_items: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    raw = rng.lognormal(mean=1.35, sigma=1.15, size=int(n_entities))
    degrees = np.rint(raw).astype(np.int64)
    zeros = rng.random(int(n_entities)) < 0.04
    degrees[zeros] = 0
    dense = rng.random(int(n_entities)) < 0.015
    degrees[dense] = rng.integers(max(1, n_items - 20), n_items + 1, size=int(np.count_nonzero(dense)))
    return np.clip(degrees, 0, int(n_items))


def _synthetic_fe_cases(seed: int) -> list[dict]:
    cases: list[dict] = []
    specs = [
        ("fe-k1", np.ones(4_000, dtype=np.int64), 100_000),
        ("fe-k2", np.full(4_000, 2, dtype=np.int64), 100_000),
        ("fe-low", np.full(3_000, 8, dtype=np.int64), 100_000),
        ("fe-medium", np.full(800, 500, dtype=np.int64), 10_000),
        ("fe-near-full", np.full(500, 1_022, dtype=np.int64), 1_024),
        ("fe-mixed", _mixed_degrees(3_000, 2_000, seed + 99), 2_000),
    ]
    for offset, (name, degrees, n_items) in enumerate(specs):
        cases.append(
            {
                "name": f"{name}:rows",
                "scenario": name,
                "X": _matrix_from_row_degrees(degrees, n_items, seed + offset),
                "models": ("FE",),
            }
        )
        cases.append(
            {
                "name": f"{name}:cols",
                "scenario": name,
                "X": _matrix_from_col_degrees(degrees, n_items, seed + 100 + offset),
                "models": ("EF",),
            }
        )
    return cases


def _synthetic_ee_wide(seed: int) -> dict:
    rng = np.random.default_rng(seed)
    n_rows = 1_000
    n_cols = 2_000_000
    nnz = 10_000
    linear = rng.choice(n_rows * n_cols, size=nnz, replace=False)
    rows = linear // n_cols
    cols = linear - rows * n_cols
    data = np.ones(nnz, dtype=np.int8)
    X = sp.csr_matrix((data, (rows, cols)), shape=(n_rows, n_cols))
    X.sum_duplicates()
    X.sort_indices()
    return {"name": "ee-wide:sparse", "scenario": "ee-wide", "X": X, "models": ("EE",)}


def _expand_scenarios(requested: Iterable[str]) -> set[str]:
    out = set(requested)
    if "all" in out:
        return set(SCENARIOS) - {"all", "all-synthetic"}
    if "all-synthetic" in out:
        out.remove("all-synthetic")
        out.update({"ee-wide", "fe-k1", "fe-k2", "fe-low", "fe-medium", "fe-near-full", "fe-mixed"})
    return out


def _build_cases(args: argparse.Namespace) -> list[dict]:
    requested = _expand_scenarios(args.scenario)
    cases: list[dict] = []

    if "fixture" in requested:
        cases.append({"name": "fixture:tests-data", "scenario": "fixture", "X": _fixture_matrix(), "models": MODELS})

    if "ingredients" in requested:
        for path in args.ingredients:
            p = Path(path)
            cases.append({"name": f"ingredients:{p.name}", "scenario": "ingredients", "X": _load_ingredients_matrix(p), "models": MODELS})

    synthetic_requested = requested & {"fe-k1", "fe-k2", "fe-low", "fe-medium", "fe-near-full", "fe-mixed"}
    if synthetic_requested:
        for case in _synthetic_fe_cases(args.seed):
            if case["scenario"] in synthetic_requested:
                cases.append(case)

    if "ee-wide" in requested:
        cases.append(_synthetic_ee_wide(args.seed + 500))

    if not cases:
        raise SystemExit("No benchmark cases selected. Use --scenario fixture, --ingredients PATH, or synthetic scenarios.")
    return cases


def _validate_output(X: sp.csr_matrix, Y: sp.csr_matrix, model: str, sort_indices: bool) -> dict:
    coo = Y.tocoo()
    support = np.column_stack((coo.row, coo.col)) if Y.nnz else np.empty((0, 2), dtype=np.int64)
    duplicate_free = np.unique(support, axis=0).shape[0] == int(Y.nnz)
    checks = {
        "shape": tuple(Y.shape) == tuple(X.shape),
        "csr": Y.getformat() == "csr",
        "duplicate_free": bool(duplicate_free),
        "int8_data": Y.data.dtype == np.dtype(np.int8),
    }
    if sort_indices:
        checks["sorted_indices"] = bool(Y.has_sorted_indices)
    if model == "FE":
        checks["row_totals"] = bool(np.array_equal(Y.getnnz(axis=1), X.getnnz(axis=1)))
    elif model == "EF":
        checks["col_totals"] = bool(np.array_equal(Y.getnnz(axis=0), X.getnnz(axis=0)))
    elif model == "EE":
        checks["nnz"] = int(Y.nnz) == int(X.nnz)
    return checks


def _estimate_temp_mb(events: list[dict]) -> float:
    if not events:
        return 0.0
    max_entries = max(int(event.get("estimated_temp_entries", 0) or 0) for event in events)
    return (max_entries * np.dtype(np.int64).itemsize) / (1024 * 1024)


def _run_benchmark(case: dict, model: str, args: argparse.Namespace, *, ee_variant: str | None = None) -> dict:
    X = case["X"].tocsr()
    debug_events: list[dict] = []
    ee_variant = ee_variant if model == "EE" else None
    result = {
        "case": case["name"],
        "scenario": case["scenario"],
        "model": model,
        "ee_variant": ee_variant,
        "seed": int(args.seed),
        "reps": int(args.reps),
        "shape": [int(X.shape[0]), int(X.shape[1])],
        "nnz": int(X.nnz),
        "density": _density(X),
        "sort_indices": bool(args.sort_indices),
        "memory_mb": args.memory_mb,
        "timings_sec": [],
        "algorithm_counts": {},
        "estimated_temp_mb": 0.0,
        "invariants": {},
        "ok": False,
        "failure": None,
    }

    try:
        construct_start = time.perf_counter()
        sampler = make_null_sampler(
            X,
            model,
            random_state=args.seed,
            sort_indices=args.sort_indices,
            memory_mb=args.memory_mb,
            ee_strategy=ee_variant,
            debug_callback=debug_events,
        )
        result["construct_sec"] = time.perf_counter() - construct_start

        iterator = sampler.sample(args.reps, seed=args.seed)
        invariant_accumulator: dict[str, bool] = {}
        for _ in range(args.reps):
            rep_start = time.perf_counter()
            Y = next(iterator)
            result["timings_sec"].append(time.perf_counter() - rep_start)
            checks = _validate_output(X, Y, model, bool(args.sort_indices))
            for key, value in checks.items():
                invariant_accumulator[key] = bool(value) if key not in invariant_accumulator else invariant_accumulator[key] and bool(value)

        result["invariants"] = invariant_accumulator
        result["ok"] = bool(invariant_accumulator) and all(invariant_accumulator.values())

    except Exception as exc:  # benchmark output should capture failures and continue
        result["failure"] = f"{type(exc).__name__}: {exc}"

    algorithms = Counter(str(event.get("algorithm", "unknown")) for event in debug_events)
    result["algorithm_counts"] = dict(sorted(algorithms.items()))
    result["algorithm_events"] = debug_events
    result["estimated_temp_mb"] = _estimate_temp_mb(debug_events)
    timings = result["timings_sec"]
    if timings:
        result["mean_sec"] = float(np.mean(timings))
        result["min_sec"] = float(np.min(timings))
        result["max_sec"] = float(np.max(timings))
    else:
        result["mean_sec"] = None
        result["min_sec"] = None
        result["max_sec"] = None
    return result


def _markdown_summary(results: list[dict]) -> str:
    lines = [
        "# Direct Null Sampler Benchmark",
        "",
        "| case | model | EE variant | shape | nnz | density | reps | algorithms | mean s/rep | ok |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | --- | ---: | --- |",
    ]
    for row in results:
        algorithms = ", ".join(f"{name}={count}" for name, count in row["algorithm_counts"].items()) or "-"
        mean_sec = "-" if row["mean_sec"] is None else f"{row['mean_sec']:.6f}"
        ok = "yes" if row["ok"] else f"no ({row['failure'] or 'invariant failed'})"
        lines.append(
            "| {case} | {model} | {ee_variant} | {shape} | {nnz} | {density:.3e} | {reps} | {algorithms} | {mean_sec} | {ok} |".format(
                case=row["case"],
                model=row["model"],
                ee_variant=row.get("ee_variant") or "-",
                shape=_shape_label(tuple(row["shape"])),
                nnz=row["nnz"],
                density=row["density"],
                reps=row["reps"],
                algorithms=algorithms,
                mean_sec=mean_sec,
                ok=ok,
            )
        )
    return "\n".join(lines)


def _write_results(path: Path, results: list[dict], args: argparse.Namespace) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix.lower() == ".jsonl":
        with path.open("w", encoding="utf-8") as handle:
            for result in results:
                handle.write(json.dumps(result, sort_keys=True) + "\n")
        return

    payload = {
        "metadata": {
            "reps": int(args.reps),
            "seed": int(args.seed),
            "memory_mb": args.memory_mb,
            "sort_indices": bool(args.sort_indices),
            "scenarios": list(args.scenario),
            "ee_variants": list(args.ee_variant),
        },
        "results": results,
    }
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark FE/EF/EE direct null samplers.")
    parser.add_argument("--ingredients", action="append", default=[], help="Path to an Ingredients pickle. Repeat for multiple HPC inputs.")
    parser.add_argument(
        "--scenario",
        action="append",
        choices=SCENARIOS,
        help="Benchmark scenario to run. Defaults to ingredients when --ingredients is set, otherwise fixture.",
    )
    parser.add_argument("--model", action="append", choices=MODELS, help="Restrict to one or more models. Defaults to each case's models.")
    parser.add_argument(
        "--ee-variant",
        action="append",
        choices=EE_VARIANTS,
        default=None,
        help=(
            "EE implementation variant to benchmark. Repeat for multiple variants. "
            "Use with --model EE to compare EE variants directly."
        ),
    )
    parser.add_argument("--reps", type=_positive_int, default=3, help="Replicates per case/model.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    parser.add_argument("--memory-mb", type=_positive_float, default=None, help="Temporary-memory budget passed to direct samplers.")
    parser.add_argument("--sort-indices", action="store_true", help="Request sorted CSR output and validate sorted indices.")
    parser.add_argument("--out", type=Path, help="Detailed output path. Use .jsonl for JSON Lines, otherwise JSON.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.scenario is None:
        args.scenario = ["ingredients"] if args.ingredients else ["fixture"]
    if args.ee_variant is None:
        args.ee_variant = ["auto"]
    if "ingredients" in args.scenario and not args.ingredients:
        raise SystemExit("--scenario ingredients requires --ingredients PATH")

    model_filter = set(args.model or MODELS)
    results: list[dict] = []
    for case in _build_cases(args):
        for model in case["models"]:
            if model in model_filter:
                if model == "EE":
                    for ee_variant in args.ee_variant:
                        results.append(_run_benchmark(case, model, args, ee_variant=ee_variant))
                else:
                    results.append(_run_benchmark(case, model, args))

    print(_markdown_summary(results))
    if args.out is not None:
        _write_results(args.out, results, args)
    return 1 if any(not row["ok"] for row in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
