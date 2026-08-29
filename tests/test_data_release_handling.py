from __future__ import annotations

import json

import numpy as np
import pytest
import scipy.sparse as sp

from metacooc._data_config import (
    DataReleaseError,
    IngredientsFormatError,
    available_releases,
    describe_data_release,
    get_file_info,
    local_data_releases,
    parse_data_release,
    resolve_release,
    validate_registry,
)
from metacooc.pantry import Ingredients, load_ingredients, save_ingredients


def _write_minimal_ingredients(path, data_release="R226_globdb_rev1"):
    matrix = sp.csr_matrix(np.array([[1]], dtype=np.uint8))
    ingredients = Ingredients(
        samples=["S001"],
        taxa=["d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__x"],
        presence_matrix=matrix,
        coverage_matrix=matrix.astype(float),
        data_release=data_release,
    )
    return save_ingredients(
        ingredients,
        str(path.parent),
        data_release=data_release,
    )


def test_release_parser_and_global_revision_listing(registry_factory):
    registry = validate_registry(registry_factory(revisions=(1, 2), current=2))

    parsed = parse_data_release("R226_gtdb_rev2")
    assert (parsed.base, parsed.variant, parsed.revision) == ("R226", "gtdb", 2)
    assert available_releases(registry) == [
        "R226_gtdb_rev1",
        "R226_globdb_rev1",
        "R226_gtdb_rev2",
        "R226_globdb_rev2",
    ]
    assert describe_data_release("R226_gtdb_rev2") == (
        "R226_gtdb_rev2 (GTDB, database release R226, revision 2)"
    )

    with pytest.raises(ValueError, match="must match"):
        parse_data_release("R226_gtdb")


def test_registry_enforces_complete_deterministic_snapshots(registry_factory):
    registry = registry_factory()
    del registry["releases"]["R226"]["revisions"]["1"]["ingredients_formats"]["1"]["globdb"]
    with pytest.raises(ValueError, match="lacks globdb"):
        validate_registry(registry)

    registry = registry_factory()
    raw = registry["releases"]["R226"]["revisions"]["1"]["ingredients_formats"]["1"]["gtdb"]["raw"]
    raw["filename"] = "wrong.tar.gz"
    with pytest.raises(ValueError, match="filename must be"):
        validate_registry(registry)


def test_resolve_release_selects_client_format_and_variant(registry_factory):
    registry = validate_registry(registry_factory(formats=(1, 2)))
    resolved = resolve_release("R226_gtdb_rev1", registry=registry)

    assert resolved.current is True
    assert resolved.format_version == 1
    assert resolved.filenames == {
        "ingredients_raw": "ingredients_raw_R226_gtdb_rev1_format1",
        "ingredients_aggregated": "ingredients_aggregated_R226_gtdb_rev1_format1",
        "sra_metadata": "sra_metadata_R226_rev1.tsv",
        "sample_to_biome": "sample_to_biome_R226_rev1.tsv",
    }
    assert all(artifact.zenodo_record.isdigit() for artifact in resolved.artifacts.values())


def test_get_file_info_uses_immutable_zenodo_record(registry_factory, monkeypatch):
    import metacooc._data_config as config

    registry = validate_registry(registry_factory())
    monkeypatch.setattr(config, "load_registry", lambda: registry)
    filenames, urls = get_file_info("R226_globdb_rev1")

    assert filenames["ingredients_raw"].endswith("globdb_rev1_format1")
    assert "sra_metadata_R226_rev1.tsv.gz" in urls
    assert all("zenodo.org/records/" in url for url in urls.values())


def test_resolve_release_rejects_unsupported_format(registry_factory):
    registry = validate_registry(registry_factory(formats=(2,)))
    with pytest.raises(DataReleaseError, match="Update MetaCoOc"):
        resolve_release("R226_gtdb_rev1", registry=registry)


def test_local_data_releases_require_revision_and_supported_format(tmp_path):
    (tmp_path / "ingredients_raw_R226_globdb_rev1_format1").mkdir()
    (tmp_path / "ingredients_aggregated_R226_gtdb_rev1_format1").mkdir()
    (tmp_path / "ingredients_raw_R226_gtdb_rev1_format2").mkdir()
    (tmp_path / "ingredients_raw_R226_gtdb").mkdir()

    assert local_data_releases(tmp_path) == [
        "R226_globdb_rev1",
        "R226_gtdb_rev1",
    ]


def test_official_save_adds_format_suffix_and_custom_label_does_not(tmp_path):
    official = _write_minimal_ingredients(
        tmp_path / "ingredients_raw_R226_gtdb_rev1_format1",
        "R226_gtdb_rev1",
    )
    custom = _write_minimal_ingredients(
        tmp_path / "ingredients_raw_R226_gtdb",
        "R226_gtdb",
    )

    assert official.endswith("ingredients_raw_R226_gtdb_rev1_format1")
    assert custom.endswith("ingredients_raw_R226_gtdb")


def test_load_published_ingredients_requires_exact_release(monkeypatch, tmp_path):
    path = _write_minimal_ingredients(
        tmp_path / "ingredients_raw_R226_globdb_rev1_format1",
        "R226_globdb_rev1",
    )
    monkeypatch.setattr(
        "metacooc.pantry.get_file_info",
        lambda release: (
            {"ingredients_raw": path.rsplit("/", 1)[-1]},
            {},
        ),
    )

    ingredients = load_ingredients(
        data_dir=str(tmp_path),
        data_release="R226_globdb_rev1",
    )
    assert ingredients.data_release == "R226_globdb_rev1"

    with pytest.raises(DataReleaseError, match="exact --data_release"):
        load_ingredients(data_dir=str(tmp_path))


def test_load_ingredients_rejects_manifest_release_mismatch(tmp_path):
    path = _write_minimal_ingredients(
        tmp_path / "ingredients_raw_R226_globdb_rev1_format1",
        "R226_globdb_rev1",
    )
    with pytest.raises(DataReleaseError, match="expected 'R226_gtdb_rev1'"):
        load_ingredients(
            custom_ingredients=path,
            data_release="R226_gtdb_rev1",
        )


def test_load_ingredients_rejects_in_memory_release_mismatch():
    matrix = sp.csr_matrix(np.array([[1]], dtype=np.uint8))
    ingredients = Ingredients(
        samples=["S001"],
        taxa=["d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__x"],
        presence_matrix=matrix,
        coverage_matrix=matrix.astype(float),
        data_release="R226_globdb_rev1",
    )

    with pytest.raises(DataReleaseError, match="expected 'R226_gtdb_rev1'"):
        load_ingredients(
            custom_ingredients=ingredients,
            data_release="R226_gtdb_rev1",
        )


@pytest.mark.parametrize("format_version", [None, 99])
def test_load_ingredients_rejects_unsupported_format(tmp_path, format_version):
    path = _write_minimal_ingredients(
        tmp_path / "ingredients_raw_R226_globdb_rev1_format1",
        "R226_globdb_rev1",
    )
    manifest_path = tmp_path / path.rsplit("/", 1)[-1] / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if format_version is None:
        manifest.pop("format_version")
    else:
        manifest["format_version"] = format_version
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(IngredientsFormatError, match="format version"):
        load_ingredients(custom_ingredients=str(manifest_path.parent))
