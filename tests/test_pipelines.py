from __future__ import annotations

import pandas as pd

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
