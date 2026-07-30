from __future__ import annotations

import json

import numpy as np
import pytest
import scipy.sparse as sp

from metacooc._data_config import (
    DataReleaseError,
    IngredientsFormatError,
    ReleaseSpec,
    available_releases,
    describe_data_release,
    get_file_info,
    local_data_releases,
    missing_data_release_message,
)
from metacooc.pantry import Ingredients, load_ingredients, save_ingredients


def _published_registry(monkeypatch):
    import metacooc._data_config as config

    monkeypatch.setitem(config.RELEASES, "R226", ReleaseSpec("test-record"))


def _write_minimal_ingredients(path, data_release="R226_globdb"):
    matrix = sp.csr_matrix(np.array([[1]], dtype=np.uint8))
    ingredients = Ingredients(
        samples=["S001"],
        taxa=["d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__x"],
        presence_matrix=matrix,
        coverage_matrix=matrix.astype(float),
        data_release=data_release,
    )
    save_ingredients(ingredients, str(path.parent), data_release=data_release)


def test_release_parser_and_display():
    assert available_releases() == ["R226_gtdb", "R226_globdb"]
    assert describe_data_release("R226_gtdb") == "R226_gtdb (GTDB, database release R226)"

    with pytest.raises(DataReleaseError, match="must match"):
        get_file_info("2.0.0_gtdb")


def test_current_download_info_is_explicit_and_variant_specific(monkeypatch):
    _published_registry(monkeypatch)

    filenames, download_urls = get_file_info("R226_gtdb")

    assert available_releases() == ["R226_gtdb", "R226_globdb"]
    assert filenames == {
        "ingredients_raw": "ingredients_raw_R226_gtdb",
        "ingredients_aggregated": "ingredients_aggregated_R226_gtdb",
        "sra_metadata": "sra_metadata_R226.tsv",
    }
    assert all("test-record" in url for url in download_urls.values())


def test_missing_release_message_uses_real_data_dir_flag(tmp_path):
    message = missing_data_release_message(
        data_dir=tmp_path,
        data_release="R226_gtdb",
        missing_path=str(tmp_path / "ingredients_raw_R226_gtdb"),
        file_kind="Ingredients",
        defaulted=False,
    )
    assert "--data_dir" in message
    assert "--data-dir" not in message


def test_missing_metadata_message_requests_optional_download(tmp_path):
    message = missing_data_release_message(
        data_dir=tmp_path,
        data_release="R226_gtdb",
        missing_path=str(tmp_path / "sra_metadata_R226.tsv"),
        file_kind="metadata",
        defaulted=False,
    )

    assert message.endswith("--data-release R226_gtdb --include-metadata")


def test_local_data_releases_reports_canonical_sources(tmp_path):
    (tmp_path / "ingredients_raw_R226_globdb").mkdir()
    (tmp_path / "ingredients_aggregated_R226_gtdb").mkdir()
    (tmp_path / "ingredients_raw_2.0.0_gtdb").mkdir()

    assert local_data_releases(tmp_path) == ["R226_globdb", "R226_gtdb"]


def test_load_ingredients_explicit_release_succeeds(monkeypatch, tmp_path):
    _published_registry(monkeypatch)
    path = tmp_path / "ingredients_raw_R226_globdb"
    _write_minimal_ingredients(path)

    ingredients = load_ingredients(
        data_dir=str(tmp_path),
        data_release="R226_globdb",
    )

    assert ingredients.data_release == "R226_globdb"


def test_load_ingredients_rejects_manifest_release_mismatch(monkeypatch, tmp_path):
    _published_registry(monkeypatch)
    path = tmp_path / "ingredients_raw_R226_globdb"
    _write_minimal_ingredients(path)

    with pytest.raises(DataReleaseError, match="expected 'R226_gtdb'"):
        load_ingredients(custom_ingredients=str(path), data_release="R226_gtdb")


def test_load_ingredients_rejects_in_memory_release_mismatch():
    matrix = sp.csr_matrix(np.array([[1]], dtype=np.uint8))
    ingredients = Ingredients(
        samples=["S001"],
        taxa=["d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__x"],
        presence_matrix=matrix,
        coverage_matrix=matrix.astype(float),
        data_release="R226_globdb",
    )

    with pytest.raises(DataReleaseError, match="expected 'R226_gtdb'"):
        load_ingredients(
            custom_ingredients=ingredients,
            data_release="R226_gtdb",
        )


@pytest.mark.parametrize("format_version", [None, 99])
def test_load_ingredients_rejects_unsupported_format(tmp_path, format_version):
    path = tmp_path / "ingredients_raw_R226_globdb"
    _write_minimal_ingredients(path)

    manifest_path = path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if format_version is None:
        manifest.pop("format_version")
    else:
        manifest["format_version"] = format_version
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(IngredientsFormatError, match="format version"):
        load_ingredients(custom_ingredients=str(path))
