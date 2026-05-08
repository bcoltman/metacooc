from __future__ import annotations

import numpy as np
import os
import scipy.sparse as sp

from metacooc.format import format_data
from metacooc.pantry import Ingredients, _read_sample_to_biome, load_ingredients
from metacooc.analysis import _cooccur_core

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


def test_raw_ingredients_contents(raw_ingredients):
    assert raw_ingredients.data_version == "test"
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
    assert aggregated_ingredients.data_version == "test"
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
    assert loaded.data_version == "test"
    assert np.array_equal(loaded.total_counts[:3], np.array([60, 13, 14], dtype=np.int32))
    assert calls == []

    assert "presence: (300, 100)" in repr(loaded)
    assert "coverage: (300, 100)" in repr(loaded)
    assert calls == []

    assert loaded.presence_matrix.dtype == np.uint8
    assert calls == ["presence.npz"]
    assert loaded.coverage_matrix.shape == (300, 100)
    assert calls == ["presence.npz", "coverage.npz"]


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
        data_version="test",
        archive_ingredients=True,
    )

    assert (out / "ingredients_raw_archive.tar.gz").exists()
    assert (out / "ingredients_aggregated_archive.tar.gz").exists()


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

    edge_arrays, _ = _cooccur_core(ingredients, ingredients.taxa, threshold=0.0)
    assert edge_arrays.cols["inter"][0] == 299


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
