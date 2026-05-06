from __future__ import annotations

import pickle

import numpy as np


def test_format_data_writes_raw_and_aggregated(raw_ingredients_path, aggregated_ingredients_path):
    assert raw_ingredients_path.exists()
    assert aggregated_ingredients_path.exists()


def test_raw_ingredients_contents(raw_ingredients):
    assert raw_ingredients.data_version == "test"
    assert raw_ingredients.samples == ["S1", "S2", "S3", "S4", "S5", "S6"]
    assert len(raw_ingredients.taxa) == 6
    assert raw_ingredients.presence_matrix.shape == (6, 6)
    assert raw_ingredients.coverage_matrix.shape == (6, 6)
    assert raw_ingredients.presence_matrix.getformat() == "csr"
    assert set(raw_ingredients.presence_matrix.data.tolist()) == {1}
    assert raw_ingredients.total_counts.tolist() == [4, 3, 3, 3, 3, 3]
    assert raw_ingredients.sample_to_biome["S1"] == ("terrestrial", "soil")
    assert raw_ingredients.biomes_order["level_1"] == ["terrestrial", "aquatic"]
    assert raw_ingredients.biomes_order["level_2"] == ["soil", "marine"]


def test_aggregated_ingredients_contains_ancestors(aggregated_ingredients):
    assert aggregated_ingredients.data_version == "test"
    assert aggregated_ingredients.presence_matrix.shape[1] == 6
    assert len(aggregated_ingredients.taxa) > 6
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

    filtered_samples = raw_ingredients.filtered_samples([True, True, True, False, False, False])
    assert filtered_samples.samples == ["S1", "S2", "S3"]
    assert filtered_samples.presence_matrix.shape == (6, 3)
    assert filtered_samples.total_counts.tolist() == [3, 3, 2, 0, 1, 1]

    filtered_taxa = raw_ingredients.filtered_taxa([0, 1, 2])
    assert len(filtered_taxa.taxa) == 3
    assert filtered_taxa.presence_matrix.shape == (3, 6)
    assert filtered_taxa.total_counts.tolist() == [4, 3, 3]
