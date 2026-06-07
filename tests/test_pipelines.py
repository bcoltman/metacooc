from __future__ import annotations

import pandas as pd
import pytest

from conftest import pipeline_args
from metacooc import pipelines


def test_run_association_pipeline(tmp_path, raw_ingredients_path, monkeypatch):
    monkeypatch.setattr(pipelines, "plot_analysis_obj", lambda *args, **kwargs: None)
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        search_mode="taxa_context",
        search_string="g__Rhizo",
        null_model="FE",
    )
    pipelines.run_association(args)
    out = tmp_path / "test_global_association.tsv"
    assert out.exists()
    assert not pd.read_csv(out, sep="\t").empty


def test_run_association_pipeline_writes_null_metadata_sidecar(
    tmp_path,
    raw_ingredients_path,
    monkeypatch,
):
    monkeypatch.setattr(pipelines, "plot_analysis_obj", lambda *args, **kwargs: None)
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        search_mode="taxa_context",
        search_string="g__Rhizo",
        null_model="EE",
        nm_n_reps=1,
        nm_seed=7,
    )
    pipelines.run_association(args)

    out = tmp_path / "test_global_association.tsv"
    metadata = tmp_path / "test_global_association_metadata.tsv"
    association_df = pd.read_csv(out, sep="\t")
    metadata_df = pd.read_csv(metadata, sep="\t")

    assert "jaccard_null_mean" in association_df.columns
    assert {"null_seed", "null_seed_source", "null_model"}.isdisjoint(association_df.columns)
    assert dict(zip(metadata_df["key"], metadata_df["value"].astype(str))) == {
        "null_model": "EE",
        "null_replicates_requested": "1",
        "null_replicates_completed": "1",
        "null_replicates_ok": "1",
        "null_replicates_error": "0",
        "null_seed": "7",
        "null_seed_source": "user",
    }


def test_run_cooccurrence_pipeline(tmp_path, raw_ingredients_path):
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        search_mode="taxa_context",
        search_string="g__Rhizo",
        null_model="FE",
        min_conditional_probability=0.0,
    )
    pipelines.run_cooccurrence(args)
    assert (tmp_path / "test_global_nodes.tsv").exists()
    assert (tmp_path / "test_global_edges.tsv").exists()


def test_run_structure_pipeline(tmp_path, raw_ingredients_path):
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        search_mode="taxa_context",
        search_string="g__Rhizo",
        null_model="FE",
        nm_n_reps=1,
    )
    pipelines.run_structure(args)
    out = tmp_path / "test_global_structure.tsv"
    assert out.exists()
    assert pd.read_csv(out, sep="\t")["metric"].tolist() == ["c_score", "mean_jaccard", "nodf"]


def test_run_biome_distribution_pipeline(tmp_path, raw_ingredients_path):
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        return_all_taxa=True,
        min_taxa_count=None,
        min_sample_count=None,
        filter_rank=None,
    )
    pipelines.run_biome_distribution(args)
    assert (tmp_path / "test_taxa_biome_distribution.tsv").exists()


def test_run_biome_distribution_skips_filtering_by_default(
    tmp_path,
    raw_ingredients_path,
    monkeypatch,
):
    def fail_filter_data_obj(*args, **kwargs):
        raise AssertionError("biome_distribution should not filter unless requested")

    monkeypatch.setattr(pipelines, "filter_data_obj", fail_filter_data_obj)
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        return_all_taxa=True,
        min_taxa_count=None,
        min_sample_count=None,
        filter_rank=None,
    )

    pipelines.run_biome_distribution(args)
    assert (tmp_path / "test_taxa_biome_distribution.tsv").exists()


def test_run_biome_distribution_taxa_query_and_biome_level(tmp_path, raw_ingredients_path):
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        return_all_taxa=False,
        taxa_query="s__rhizo_000,s__micro_000",
        biome_level="level_2",
        min_taxa_count=None,
        min_sample_count=None,
        filter_rank=None,
    )
    pipelines.run_biome_distribution(args)

    out = pd.read_csv(tmp_path / "test_taxa_biome_distribution.tsv", sep="\t", index_col=0)
    assert out.columns.tolist() == ["soil", "marine"]
    assert len(out) == 2
    assert any(row.endswith("g__Rhizo; s__rhizo_000") for row in out.index)
    assert any(row.endswith("g__Micro; s__micro_000") for row in out.index)


def test_run_biome_distribution_applies_count_filter_before_distribution(
    tmp_path,
    raw_ingredients_path,
):
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        return_all_taxa=False,
        taxa_query="s__rhizo_000",
        min_taxa_count=50,
        min_sample_count=None,
        filter_rank=None,
        taxa_count_rank="species",
    )
    pipelines.run_biome_distribution(args)

    out = pd.read_csv(tmp_path / "test_taxa_biome_distribution.tsv", sep="\t", index_col=0)
    row = out.iloc[0]
    assert row["terrestrial"] == 7
    assert row["aquatic"] == 10


def test_run_biome_distribution_rejects_invalid_taxa_query(tmp_path, raw_ingredients_path):
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        return_all_taxa=False,
        taxa_query="g__Rhizo|g__Micro",
        min_taxa_count=None,
        min_sample_count=None,
        filter_rank=None,
    )

    with pytest.raises(ValueError, match="focal_taxa mode does not support"):
        pipelines.run_biome_distribution(args)
