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


def test_run_cooccurrence_pipeline(tmp_path, raw_ingredients_path):
    args = pipeline_args(
        custom_ingredients=str(raw_ingredients_path),
        output_dir=str(tmp_path),
        search_mode="taxa_context",
        search_string="g__Rhizo",
        null_model="FE",
        threshold=0.0,
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
    )

    with pytest.raises(ValueError, match="focal_taxa mode does not support"):
        pipelines.run_biome_distribution(args)
