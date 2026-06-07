from __future__ import annotations

import os
from collections.abc import Mapping
from typing import Any

import pandas as pd


NULL_METADATA_KEYS = [
    "null_model",
    "null_replicates_requested",
    "null_replicates_completed",
    "null_replicates_ok",
    "null_replicates_error",
    "null_seed",
    "null_seed_source",
]


def result_metadata_path(result_path: str) -> str:
    stem, _ = os.path.splitext(result_path)
    return f"{stem}_metadata.tsv"


def null_metadata_from_reduction(reduction: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "null_model": reduction["null_model"],
        "null_replicates_requested": int(reduction["n_target"]),
        "null_replicates_completed": int(reduction["n_done"]),
        "null_replicates_ok": int(reduction["n_ok"]),
        "null_replicates_error": int(reduction["n_err"]),
        "null_seed": int(reduction["null_seed"]),
        "null_seed_source": reduction["null_seed_source"],
    }


def write_result_metadata_sidecar(
    result_path: str,
    metadata: Mapping[str, Any] | None,
) -> str | None:
    if not metadata:
        return None

    metadata_path = result_metadata_path(result_path)
    ordered_keys = [key for key in NULL_METADATA_KEYS if key in metadata]
    ordered_keys.extend(key for key in metadata if key not in ordered_keys)
    rows = [{"key": key, "value": metadata[key]} for key in ordered_keys]
    pd.DataFrame(rows, columns=["key", "value"]).to_csv(metadata_path, sep="\t", index=False)
    return metadata_path
