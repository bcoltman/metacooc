from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import numpy as np
import pytest
import scipy.sparse as sp

from metacooc.format import format_data
from metacooc.pantry import (
    Ingredients,
    _read_sample_to_biome,
    load_ingredients,
    presence_for_counts,
    save_ingredients_directory,
    threshold_ingredients_presence,
)
from metacooc.analysis import _cooccur_core


def _assert_date_generated(value):
    assert value.endswith("Z")
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    assert parsed.tzinfo == timezone.utc
    assert parsed.microsecond == 0


def test_format_data_writes_raw_and_aggregated(raw_ingredients_path, aggregated_ingredients_path):
    assert raw_ingredients_path.exists()
    assert aggregated_ingredients_path.exists()
    for path in (raw_ingredients_path, aggregated_ingredients_path):
        assert (path / "manifest.json").exists()
        assert (path / "samples.tsv").exists()
        assert (path / "taxa.tsv").exists()
        assert (path / "sample_to_biome.tsv").exists()
        assert (path / "presence.npz").exists()
        assert (path / "coverage.npz").exists()
        assert (path / "total_counts.npy").exists()
        manifest = json.loads((path / "manifest.json").read_text())
        _assert_date_generated(manifest["date_generated"])
        assert manifest["data_release"] == "test"


def test_load_ingredients_accepts_manifest_without_date_generated(tmp_path):
    ingredients = Ingredients(
        samples=["S1"],
        taxa=["d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__a"],
        presence_matrix=sp.csr_matrix([[1]], dtype=np.uint8),
        coverage_matrix=sp.csr_matrix([[1.0]], dtype=float),
    )
    out = tmp_path / "ingredients_legacy_manifest"
    save_ingredients_directory(ingredients, str(out))

    manifest_path = out / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.pop("date_generated")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    loaded = load_ingredients(custom_ingredients=str(out))
    assert loaded.samples == ["S1"]
    assert loaded.taxa == ["d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__a"]


def test_raw_ingredients_contents(raw_ingredients):
    assert raw_ingredients.data_release == "test"
    assert raw_ingredients.samples[0] == "S001"
    assert len(raw_ingredients.samples) == 100
    assert set(raw_ingredients.samples) == {f"S{i:03d}" for i in range(1, 101)}
    assert len(raw_ingredients.taxa) == 300
    assert raw_ingredients.presence_matrix.shape == (300, 100)
    assert raw_ingredients.coverage_matrix.shape == (300, 100)
    assert raw_ingredients.presence_matrix.getformat() == "csr"
    assert raw_ingredients.presence_matrix.dtype == np.uint8
    assert set(raw_ingredients.presence_matrix.data.tolist()) == {1}
    assert raw_ingredients.total_counts.shape == (300,)
    assert raw_ingredients.total_counts[0] == 60
    assert raw_ingredients.total_counts[50] == 51
    assert raw_ingredients.sample_to_biome["S001"] == ("terrestrial", "soil")
    assert raw_ingredients.sample_to_biome["S100"] == ("aquatic", "marine")
    assert raw_ingredients.biomes_order["level_1"] == ["terrestrial", "aquatic"]
    assert raw_ingredients.biomes_order["level_2"] == ["soil", "marine"]


def test_available_biomes_reports_counts(raw_ingredients):
    available = raw_ingredients.available_biomes()
    assert available == {
        "level_1": [
            {"biome": "terrestrial", "n_samples": 50},
            {"biome": "aquatic", "n_samples": 50},
        ],
        "level_2": [
            {"biome": "soil", "n_samples": 50},
            {"biome": "marine", "n_samples": 50},
        ],
    }


def test_available_biomes_handles_missing_mapping(raw_ingredients):
    ingredients = raw_ingredients.copy(deep=False)
    ingredients.sample_to_biome = {}
    if hasattr(ingredients, "biomes_order"):
        delattr(ingredients, "biomes_order")
    if hasattr(ingredients, "sample_biome_indices"):
        delattr(ingredients, "sample_biome_indices")

    assert ingredients.available_biomes() == {"level_1": [], "level_2": []}


def test_aggregated_ingredients_contains_ancestors(aggregated_ingredients):
    assert aggregated_ingredients.data_release == "test"
    assert aggregated_ingredients.presence_matrix.shape[1] == 100
    assert len(aggregated_ingredients.taxa) > 300
    assert any(t.endswith("AGGREGATED") for t in aggregated_ingredients.taxa)
    assert any("g__Rhizo AGGREGATED" in t for t in aggregated_ingredients.taxa)


def test_ingredients_directory_lazy_loads_matrices(raw_ingredients_path, monkeypatch):
    calls = []
    real_load_npz = sp.load_npz

    def tracked_load_npz(path):
        calls.append(os.path.basename(path))
        return real_load_npz(path)

    monkeypatch.setattr(sp, "load_npz", tracked_load_npz)
    loaded = load_ingredients(custom_ingredients=str(raw_ingredients_path))

    assert loaded.samples[0] == "S001"
    assert loaded.taxa[0].endswith("s__rhizo_000")
    assert loaded.sample_to_biome["S001"] == ("terrestrial", "soil")
    assert loaded.data_release == "test"
    assert np.array_equal(loaded.total_counts[:3], np.array([60, 13, 14], dtype=np.int32))
    assert calls == []

    assert "presence: (300, 100)" in repr(loaded)
    assert "coverage: (300, 100)" in repr(loaded)
    assert calls == []

    assert loaded.presence_matrix.dtype == np.uint8
    assert calls == ["presence.npz"]
    assert loaded.coverage_matrix.shape == (300, 100)
    assert calls == ["presence.npz", "coverage.npz"]


def test_biome_distribution_does_not_load_coverage(raw_ingredients_path, monkeypatch):
    calls = []
    real_load_npz = sp.load_npz

    def tracked_load_npz(path):
        calls.append(os.path.basename(path))
        return real_load_npz(path)

    monkeypatch.setattr(sp, "load_npz", tracked_load_npz)
    loaded = load_ingredients(custom_ingredients=str(raw_ingredients_path))
    assert loaded._coverage_matrix is None

    biomes, presence, n_dropped = loaded.biome_distribution()

    assert biomes == ["terrestrial", "aquatic"]
    assert presence.shape == (2, 300)
    assert n_dropped == 0
    assert calls == ["presence.npz"]
    assert loaded._coverage_matrix is None


def test_read_sample_to_biome_vectorized_missing_values(tmp_path):
    path = tmp_path / "sample_to_biome.tsv"
    path.write_text(
        "accession\tlevel_1\tlevel_2\textra\n"
        "S001\tterrestrial\tsoil\tignored\n"
        "S002\t\tmarine\tignored\n"
        "S003\taquatic\t\tignored\n"
    )

    assert _read_sample_to_biome(str(path)) == {
        "S001": ("terrestrial", "soil"),
        "S002": (None, "marine"),
        "S003": ("aquatic", None),
    }


def test_format_data_archive_option_writes_tarballs(tmp_path, fixture_dir):
    out = tmp_path / "formatted"
    format_data(
        tax_profile=str(fixture_dir / "tax_profile.tsv"),
        output_dir=str(out),
        sample_to_biome_file=str(fixture_dir / "sample_to_biome.tsv"),
        aggregated=True,
        tag="archive",
        data_release="test",
        archive_ingredients=True,
    )

    assert (out / "ingredients_raw_archive.tar.gz").exists()
    assert (out / "ingredients_aggregated_archive.tar.gz").exists()


def test_format_official_release_uses_revision_and_schema_in_names(
    tmp_path, fixture_dir
):
    out = tmp_path / "official"
    format_data(
        tax_profile=str(fixture_dir / "tax_profile.tsv"),
        output_dir=str(out),
        sample_to_biome_file=str(fixture_dir / "sample_to_biome.tsv"),
        aggregated=True,
        data_release="R226_gtdb_rev1",
        archive_ingredients=True,
    )

    raw = out / "ingredients_raw_R226_gtdb_rev1_format1"
    aggregated = out / "ingredients_aggregated_R226_gtdb_rev1_format1"
    assert raw.is_dir()
    assert aggregated.is_dir()
    assert (out / f"{raw.name}.tar.gz").exists()
    assert (out / f"{aggregated.name}.tar.gz").exists()
    manifest = json.loads((raw / "manifest.json").read_text())
    assert manifest["data_release"] == "R226_gtdb_rev1"
    assert manifest["format_version"] == 1


def test_format_rejects_tag_for_official_release_before_reading_input(tmp_path):
    from metacooc._data_config import DataReleaseError

    with pytest.raises(DataReleaseError, match="--tag cannot be combined"):
        format_data(
            tax_profile=str(tmp_path / "not-read.tsv"),
            output_dir=str(tmp_path / "out"),
            tag="custom",
            data_release="R226_gtdb_rev1",
        )


def test_count_operations_do_not_overflow_uint8_presence():
    n_samples = 301
    matrix = sp.lil_matrix((2, n_samples), dtype=np.uint8)
    matrix[0, :300] = 1
    matrix[1, 1:] = 1
    matrix = matrix.tocsr()
    ingredients = Ingredients(
        samples=[f"S{i}" for i in range(n_samples)],
        taxa=[
            "d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__a",
            "d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__b",
        ],
        presence_matrix=matrix,
        coverage_matrix=matrix.astype(float),
    )

    edge_arrays, _ = _cooccur_core(
        ingredients,
        ingredients.taxa,
        min_conditional_probability=0.0,
    )
    assert edge_arrays.cols["shared_sample_count"][0] == 299


def test_presence_for_counts_prevents_uint8_sparse_product_overflow():
    matrix = sp.csr_matrix(np.ones((1, 300), dtype=np.uint8))

    direct = matrix @ matrix.T
    safe = presence_for_counts(matrix) @ presence_for_counts(matrix).T

    assert direct.dtype == np.uint8
    assert direct[0, 0] == 44
    assert safe.dtype == np.int32
    assert safe[0, 0] == 300


def test_coverage_threshold_updates_presence_and_masks_coverage():
    ingredients = Ingredients(
        samples=["S1", "S2"],
        taxa=[
            "d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__a",
            "d__Bacteria; p__P; c__C; o__O; f__F; g__G",
        ],
        presence_matrix=sp.csr_matrix([[1, 1], [1, 1]], dtype=np.uint8),
        coverage_matrix=sp.csr_matrix([[0.5, 1.5], [0.1, 0.8]], dtype=float),
    )

    thresholded = threshold_ingredients_presence(
        ingredients,
        min_coverage=1.0,
        min_coverage_by_rank={"genus": 0.5},
    )

    assert thresholded.presence_matrix.toarray().tolist() == [[0, 1], [0, 1]]
    assert thresholded.coverage_matrix.toarray().tolist() == [[0.0, 1.5], [0.0, 0.8]]
    assert thresholded.total_counts.tolist() == [1, 1]
    assert ingredients.presence_matrix.toarray().tolist() == [[1, 1], [1, 1]]


def test_relative_abundance_threshold_uses_sample_fractions():
    ingredients = Ingredients(
        samples=["S1", "S2"],
        taxa=[
            "d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__a",
            "d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__b",
        ],
        presence_matrix=sp.csr_matrix([[1, 1], [1, 1]], dtype=np.uint8),
        coverage_matrix=sp.csr_matrix([[3.0, 1.0], [1.0, 9.0]], dtype=float),
    )

    thresholded = threshold_ingredients_presence(
        ingredients,
        min_relative_abundance=0.5,
    )

    assert thresholded.presence_matrix.toarray().tolist() == [[1, 0], [0, 1]]
    assert thresholded.coverage_matrix.toarray().tolist() == [[3.0, 0.0], [0.0, 9.0]]
    assert thresholded.presence_thresholds["min_relative_abundance"] == 0.5


def test_combined_coverage_and_relative_threshold_warns():
    ingredients = Ingredients(
        samples=["S1", "S2"],
        taxa=[
            "d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__a",
            "d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__b",
        ],
        presence_matrix=sp.csr_matrix([[1, 1], [1, 1]], dtype=np.uint8),
        coverage_matrix=sp.csr_matrix([[3.0, 1.0], [1.0, 9.0]], dtype=float),
    )

    with pytest.warns(UserWarning, match="Both coverage and relative-abundance"):
        thresholded = threshold_ingredients_presence(
            ingredients,
            min_coverage=2.0,
            min_relative_abundance=0.5,
        )

    assert thresholded.presence_matrix.toarray().tolist() == [[1, 0], [0, 1]]


def test_aggregated_relative_abundance_recomputed_from_descendants():
    ingredients = Ingredients(
        samples=["S1", "S2"],
        taxa=[
            "Root; d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__a",
            "Root; d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__b",
            "Root; d__Bacteria; p__P; c__C; o__O; f__F; g__G AGGREGATED",
        ],
        presence_matrix=sp.csr_matrix([[1, 1], [1, 1], [1, 1]], dtype=np.uint8),
        coverage_matrix=sp.csr_matrix(
            [[3.0, 1.0], [1.0, 9.0], [4.0, 10.0]],
            dtype=float,
        ),
    )

    thresholded = threshold_ingredients_presence(
        ingredients,
        min_relative_abundance_by_rank={
            "species": 0.5,
            "genus": 0.9,
        },
    )

    assert thresholded.presence_matrix.toarray().tolist() == [
        [1, 0],
        [0, 1],
        [1, 1],
    ]
    assert thresholded.coverage_matrix.toarray().tolist() == [
        [3.0, 0.0],
        [0.0, 9.0],
        [4.0, 10.0],
    ]


def test_thresholded_ingredients_manifest_provenance_roundtrip(tmp_path):
    ingredients = Ingredients(
        samples=["S1", "S2"],
        taxa=["d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__a"],
        presence_matrix=sp.csr_matrix([[1, 1]], dtype=np.uint8),
        coverage_matrix=sp.csr_matrix([[0.5, 2.0]], dtype=float),
    )
    thresholded = threshold_ingredients_presence(ingredients, min_coverage=1.0)

    out = tmp_path / "ingredients_thresholded"
    save_ingredients_directory(thresholded, str(out))

    loaded = load_ingredients(custom_ingredients=str(out))
    assert loaded.presence_matrix.toarray().tolist() == [[0, 1]]
    assert loaded.coverage_matrix.toarray().tolist() == [[0.0, 2.0]]
    assert loaded.presence_thresholds["min_coverage"] == 1.0
    assert loaded.presence_thresholds["comparison"] == ">="


def test_ingredients_filter_and_copy_consistency(raw_ingredients):
    shallow = raw_ingredients.copy(deep=False)
    deep = raw_ingredients.copy(deep=True)
    assert shallow.samples == raw_ingredients.samples
    assert deep.samples == raw_ingredients.samples

    mask = [i < 50 for i in range(100)]
    filtered_samples = raw_ingredients.filtered_samples(mask)
    assert filtered_samples.samples[0] == "S001"
    assert filtered_samples.samples[-1] == "S050"
    assert filtered_samples.presence_matrix.shape == (300, 50)
    assert filtered_samples.total_counts[0] == 50

    filtered_taxa = raw_ingredients.filtered_taxa([0, 1, 2])
    assert len(filtered_taxa.taxa) == 3
    assert filtered_taxa.presence_matrix.shape == (3, 100)
    assert filtered_taxa.total_counts[0] == 60
