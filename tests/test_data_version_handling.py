from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest
import scipy.sparse as sp

from metacooc._data_config import (
    DataVersionError,
    format_local_data_versions,
    local_data_versions,
)
from metacooc.pantry import Ingredients, load_ingredients, save_ingredients
from metacooc.search import search_data_obj


def _write_minimal_ingredients(path, data_version="1.1.0_globdb"):
    matrix = sp.csr_matrix(np.array([[1]], dtype=np.int8))
    ingredients = Ingredients(
        samples=["S001"],
        taxa=["d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__x"],
        presence_matrix=matrix,
        coverage_matrix=matrix.copy(),
        data_version=data_version,
    )
    save_ingredients(ingredients, str(path.parent), data_version=data_version)


def test_local_data_versions_reports_full_sources(tmp_path):
    (tmp_path / "ingredients_raw_1.1.0_globdb").mkdir()
    (tmp_path / "ingredients_aggregated_1.1.0_gtdb").mkdir()
    (tmp_path / "sra_metadata_1.1.0.tsv").touch()

    assert local_data_versions(tmp_path) == ["1.1.0_globdb", "1.1.0_gtdb"]
    display = format_local_data_versions(tmp_path)
    assert "1.1.0_globdb (GlobDB r226)" in display
    assert "1.1.0_gtdb (GTDB)" in display


def test_load_ingredients_missing_default_lists_local_variant(tmp_path):
    (tmp_path / "ingredients_raw_1.1.0_globdb").mkdir()

    with pytest.raises(DataVersionError) as excinfo:
        load_ingredients(data_dir=str(tmp_path), aggregated=False, data_version=None)

    message = str(excinfo.value)
    assert "No --data_version was specified" in message
    assert "1.1.0_gtdb (GTDB)" in message
    assert "1.1.0_globdb (GlobDB r226)" in message
    assert "--data_version 1.1.0_globdb" in message
    assert "Available: 1.1.0" not in message


def test_load_ingredients_explicit_available_variant_succeeds(tmp_path):
    path = tmp_path / "ingredients_raw_1.1.0_globdb"
    _write_minimal_ingredients(path)

    ingredients = load_ingredients(
        data_dir=str(tmp_path),
        aggregated=False,
        data_version="1.1.0_globdb",
    )

    assert ingredients.data_version == "1.1.0_globdb"


def test_metadata_missing_default_uses_data_version_error(tmp_path):
    (tmp_path / "ingredients_raw_1.1.0_globdb").mkdir()

    with pytest.raises(DataVersionError) as excinfo:
        search_data_obj(
            search_mode="metadata",
            search_string="soil",
            data_dir=str(tmp_path),
            data_version=None,
        )

    message = str(excinfo.value)
    assert "Missing required metadata file" in message
    assert "1.1.0_globdb (GlobDB r226)" in message


def test_cli_missing_default_data_version_has_no_traceback(tmp_path):
    (tmp_path / "ingredients_raw_1.1.0_globdb").mkdir()
    out = tmp_path / "out"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "metacooc",
            "search",
            "--search_mode",
            "taxa_context",
            "--search_string",
            "g__G",
            "--data_dir",
            str(tmp_path),
            "--output_dir",
            str(out),
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "Traceback" not in result.stderr
    assert "No --data_version was specified" in result.stderr
    assert "1.1.0_globdb (GlobDB r226)" in result.stderr
    assert "Available: 1.1.0" not in result.stderr
