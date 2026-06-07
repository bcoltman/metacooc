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


def test_run_cooccurrence_pipeline_writes_null_metadata_sidecar(tmp_path, raw_ingredients_path):
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        search_mode="taxa_context",
        search_string="g__Rhizo",
        null_model="EE",
        nm_n_reps=1,
        nm_seed=7,
        min_conditional_probability=0.0,
    )
    pipelines.run_cooccurrence(args)

    edges = pd.read_csv(tmp_path / "test_global_edges.tsv", sep="\t")
    metadata = pd.read_csv(tmp_path / "test_global_edges_metadata.tsv", sep="\t")
    nodes = pd.read_csv(tmp_path / "test_global_nodes.tsv", sep="\t")

    assert "source_taxon" in edges.columns
    assert "target_taxon" in edges.columns
    assert "jaccard_null_mean" in edges.columns
    assert {"null_seed", "null_seed_source", "null_model"}.isdisjoint(edges.columns)
    assert "taxon_sample_count" in nodes.columns
    assert "out_degree_p_target_given_source_gt_0.0" in nodes.columns
    assert dict(zip(metadata["key"], metadata["value"].astype(str))) == {
        "null_model": "EE",
        "null_replicates_requested": "1",
        "null_replicates_completed": "1",
        "null_replicates_ok": "1",
        "null_replicates_error": "0",
        "null_seed": "7",
        "null_seed_source": "user",
    }


def test_run_cooccurrence_pipeline_compute_fisher_outputs_readable_columns(
    tmp_path,
    raw_ingredients_path,
):
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        search_mode="taxa_context",
        search_string="g__Rhizo",
        null_model="FE",
        min_conditional_probability=0.0,
        compute_fisher=True,
    )
    pipelines.run_cooccurrence(args)

    edges = pd.read_csv(tmp_path / "test_global_edges.tsv", sep="\t")
    assert {
        "fisher_odds_ratio",
        "fisher_p_value",
        "fisher_log_p_value",
    }.issubset(edges.columns)
    assert {"fisher_odds", "fisher_p", "log_fisher_p"}.isdisjoint(edges.columns)


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
    metadata = tmp_path / "test_global_structure_metadata.tsv"
    assert out.exists()
    structure_df = pd.read_csv(out, sep="\t")
    metadata_df = pd.read_csv(metadata, sep="\t")
    assert structure_df["metric"].tolist() == ["c_score", "mean_jaccard", "nodf"]
    assert "observed_value" in structure_df.columns
    assert "null_mean" in structure_df.columns
    assert {"null_seed", "null_seed_source", "null_model"}.isdisjoint(structure_df.columns)
    assert dict(zip(metadata_df["key"], metadata_df["value"].astype(str))) == {
        "null_model": "FE",
        "null_replicates_requested": "1",
        "null_replicates_completed": "1",
        "null_replicates_ok": "1",
        "null_replicates_error": "0",
        "null_seed": "7",
        "null_seed_source": "user",
    }


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
    out_path = tmp_path / "test_taxa_biome_distribution.tsv"
    out = pd.read_csv(out_path, sep="\t")
    assert out.columns.tolist() == ["taxon", "terrestrial", "aquatic"]
    assert len(out) == 300
    assert not (tmp_path / "test_taxa_biome_distribution_metadata.tsv").exists()


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

    out = pd.read_csv(tmp_path / "test_taxa_biome_distribution.tsv", sep="\t")
    assert out.columns.tolist() == ["taxon", "soil", "marine"]
    assert len(out) == 2
    assert out["taxon"].str.endswith("g__Rhizo; s__rhizo_000").any()
    assert out["taxon"].str.endswith("g__Micro; s__micro_000").any()


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

    out = pd.read_csv(tmp_path / "test_taxa_biome_distribution.tsv", sep="\t")
    row = out.iloc[0]
    assert row["taxon"].endswith("g__Rhizo; s__rhizo_000")
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
