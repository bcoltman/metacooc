from __future__ import annotations

import pickle

from metacooc.filter import filter_data, filter_data_obj
from metacooc.search import search_data_obj


def test_taxon_biome_and_metadata_search(raw_ingredients, metadata_file):
    taxon_hits = search_data_obj(
        search_mode="taxon",
        search_string="g__Rhizo",
        custom_ingredients=raw_ingredients,
    )
    assert taxon_hits == {"S1", "S2", "S3", "S4"}

    and_hits = search_data_obj(
        search_mode="taxon",
        search_string="g__Rhizo+g__Micro",
        custom_ingredients=raw_ingredients,
    )
    assert and_hits == {"S2", "S3"}

    inverse_hits = search_data_obj(
        search_mode="taxon",
        search_string="g__Rhizo",
        custom_ingredients=raw_ingredients,
        inverse=True,
    )
    assert inverse_hits == {"S5", "S6"}

    biome_hits = search_data_obj(
        search_mode="biome",
        search_string="soil",
        custom_ingredients=raw_ingredients,
    )
    assert biome_hits == {"S1", "S2", "S3"}

    metadata_hits = search_data_obj(
        search_mode="metadata",
        search_string="soil",
        metadata_file=str(metadata_file),
        column_names=["env_biome_sam"],
    )
    assert metadata_hits == {"S1", "S2", "S3"}


def test_filter_data_obj_by_counts_rank_and_accessions(raw_ingredients):
    filtered, ok = filter_data_obj(
        raw_ingredients,
        accession_set={"S1", "S2", "S3"},
        min_taxa_count=2,
        min_sample_count=2,
        filter_rank="species",
        taxa_count_rank="species",
    )
    assert ok
    assert filtered.samples == ["S1", "S2", "S3"]
    assert all("s__" in t for t in filtered.taxa)
    assert filtered.presence_matrix.shape[1] == 3
    assert filtered.presence_matrix.shape[0] >= 3


def test_file_filter_writes_null_and_filtered(tmp_path, raw_ingredients_path):
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("S1\nS2\nS3\n")

    out = tmp_path / "filter_out"
    filter_data(
        accessions_file=str(accessions),
        data_dir=None,
        output_dir=str(out),
        custom_ingredients=str(raw_ingredients_path),
        min_taxa_count=1,
        min_sample_count=1,
        filter_rank="species",
        taxa_count_rank="species",
        tag="test_",
    )

    null_path = out / "test_ingredients_null.pkl"
    filtered_path = out / "test_ingredients_filtered.pkl"
    assert null_path.exists()
    assert filtered_path.exists()

    with filtered_path.open("rb") as f:
        filtered = pickle.load(f)
    assert filtered.samples == ["S1", "S2", "S3"]
