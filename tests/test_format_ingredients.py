from __future__ import annotations

import pickle

import numpy as np


def test_format_data_writes_raw_and_aggregated(raw_ingredients_path, aggregated_ingredients_path):
    assert raw_ingredients_path.exists()
    assert aggregated_ingredients_path.exists()


def test_raw_ingredients_contents(raw_ingredients):
    assert raw_ingredients.data_version == "test"
    assert raw_ingredients.samples[0] == "S001"
    assert len(raw_ingredients.samples) == 100
    assert set(raw_ingredients.samples) == {f"S{i:03d}" for i in range(1, 101)}
    assert len(raw_ingredients.taxa) == 300
    assert raw_ingredients.presence_matrix.shape == (300, 100)
    assert raw_ingredients.coverage_matrix.shape == (300, 100)
    assert raw_ingredients.presence_matrix.getformat() == "csr"
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


def test_ingredients_pickle_roundtrip(raw_ingredients_path):
    with raw_ingredients_path.open("rb") as f:
        loaded = pickle.load(f)
    dumped = pickle.dumps(loaded)
    restored = pickle.loads(dumped)
    assert restored.samples == loaded.samples
    assert restored.taxa == loaded.taxa
    assert restored.data_version == "test"
    assert np.array_equal(restored.total_counts, loaded.total_counts)
    assert (restored.presence_matrix != loaded.presence_matrix).nnz == 0


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
