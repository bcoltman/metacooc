from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd


COMPACT_NULL_METADATA_COLUMNS = [
    "null_model",
    "null_replicates",
    "null_replicates_failed",
    "null_seed",
]


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


def with_compact_null_metadata(
    result: pd.DataFrame,
    metadata: Mapping[str, Any] | None,
    *,
    fallback_null_model: str | None = None,
) -> pd.DataFrame:
    """Return an export copy with compact run-level null metadata columns."""
    if metadata:
        null_model = metadata.get("null_model")
        null_replicates = metadata.get("null_replicates_ok")
        null_replicates_failed = metadata.get("null_replicates_error")
        null_seed = metadata.get("null_seed")
    else:
        null_model = fallback_null_model
        null_replicates = None
        null_replicates_failed = None
        null_seed = None

    out = result.copy()
    out["null_model"] = pd.array([null_model] * len(out), dtype="string")
    out["null_replicates"] = pd.array(
        [null_replicates] * len(out),
        dtype="Int64",
    )
    out["null_replicates_failed"] = pd.array(
        [null_replicates_failed] * len(out),
        dtype="Int64",
    )
    out["null_seed"] = pd.array([null_seed] * len(out), dtype="Int64")
    return out
