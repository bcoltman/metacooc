from __future__ import annotations

import subprocess
import sys

import pandas as pd
import pytest

SOIL_SAMPLE_LINES = [f"S{i:03d}" for i in range(1, 51)]
RHIZO_SAMPLE_LINES = [f"S{i:03d}" for i in range(1, 61)]


def run_cli(*args, cwd=None):
    return subprocess.run(
        [sys.executable, "-m", "metacooc", *map(str, args)],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )


@pytest.fixture
def cli_formatted_dir(tmp_path, fixture_dir):
    out = tmp_path / "cli_formatted"
    run_cli(
        "format",
        "--tax_profile",
        fixture_dir / "tax_profile.tsv",
        "--sample_to_biome_file",
        fixture_dir / "sample_to_biome.tsv",
        "--output_dir",
        out,
        "--aggregated",
        "--tag",
        "cli",
        "--data_version",
        "test",
    )
    assert (out / "ingredients_raw_cli.pkl").exists()
    assert (out / "ingredients_aggregated_cli.pkl").exists()
    return out


@pytest.mark.cli
def test_cli_search_taxon_biome_and_metadata(tmp_path, cli_formatted_dir, metadata_file):
    raw = cli_formatted_dir / "ingredients_raw_cli.pkl"

    taxon_out = tmp_path / "search_taxon"
    run_cli(
        "search",
        "--search_mode",
        "taxa_context",
        "--search_string",
        "g__Rhizo",
        "--custom_ingredients",
        raw,
        "--output_dir",
        taxon_out,
    )
    assert (taxon_out / "search_results.txt").read_text().splitlines() == RHIZO_SAMPLE_LINES

    biome_out = tmp_path / "search_biome"
    run_cli(
        "search",
        "--search_mode",
        "biome",
        "--search_string",
        "soil",
        "--custom_ingredients",
        raw,
        "--output_dir",
        biome_out,
    )
    assert (biome_out / "search_results.txt").read_text().splitlines() == SOIL_SAMPLE_LINES

    metadata_out = tmp_path / "search_metadata"
    run_cli(
        "search",
        "--search_mode",
        "metadata",
        "--search_string",
        "soil",
        "--metadata_file",
        metadata_file,
        "--column_names",
        "env_biome_sam",
        "--output_dir",
        metadata_out,
    )
    assert (metadata_out / "search_results.txt").read_text().splitlines() == SOIL_SAMPLE_LINES


@pytest.mark.cli
def test_cli_filter_and_analysis_commands(tmp_path, cli_formatted_dir):
    raw = cli_formatted_dir / "ingredients_raw_cli.pkl"
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("".join(f"{sample}\n" for sample in SOIL_SAMPLE_LINES))

    filter_out = tmp_path / "filter"
    run_cli(
        "filter",
        "--custom_ingredients",
        raw,
        "--accessions_file",
        accessions,
        "--output_dir",
        filter_out,
        "--filter_rank",
        "species",
        "--tag",
        "cli",
    )
    null_file = filter_out / "cli_ingredients_null.pkl"
    filtered_file = filter_out / "cli_ingredients_filtered.pkl"
    assert null_file.exists()
    assert filtered_file.exists()

    assoc_out = tmp_path / "analysis_assoc"
    run_cli(
        "analysis",
        "--analysis_type",
        "association",
        "--null_file",
        null_file,
        "--filtered_file",
        filtered_file,
        "--output_dir",
        assoc_out,
        "--null_model",
        "FE",
        "--nm_n_reps",
        "1",
    )
    assert not pd.read_csv(assoc_out / "association.tsv", sep="\t").empty

    cooc_out = tmp_path / "analysis_cooc"
    run_cli(
        "analysis",
        "--analysis_type",
        "cooccurrence",
        "--null_file",
        null_file,
        "--filtered_file",
        filtered_file,
        "--output_dir",
        cooc_out,
        "--large",
        "--threshold",
        "0",
        "--null_model",
        "FE",
        "--nm_n_reps",
        "1",
    )
    assert (cooc_out / "taxon_nodes.tsv").exists()

    structure_out = tmp_path / "analysis_structure"
    run_cli(
        "analysis",
        "--analysis_type",
        "structure",
        "--filtered_file",
        filtered_file,
        "--output_dir",
        structure_out,
        "--null_model",
        "FE",
        "--nm_n_reps",
        "1",
    )
    assert (structure_out / "structure.tsv").exists()


@pytest.mark.cli
def test_cli_full_workflow_commands(tmp_path, cli_formatted_dir):
    raw = cli_formatted_dir / "ingredients_raw_cli.pkl"

    assoc_out = tmp_path / "workflow_assoc"
    run_cli(
        "association",
        "--search_mode",
        "taxa_context",
        "--search_string",
        "g__Rhizo",
        "--custom_ingredients",
        raw,
        "--output_dir",
        assoc_out,
        "--null_model",
        "FE",
        "--nm_n_reps",
        "1",
        "--tag",
        "cli",
    )
    assert (assoc_out / "cli_global_association.tsv").exists()

    cooc_out = tmp_path / "workflow_cooc"
    run_cli(
        "cooccurrence",
        "--search_mode",
        "taxa_context",
        "--search_string",
        "g__Rhizo",
        "--custom_ingredients",
        raw,
        "--output_dir",
        cooc_out,
        "--large",
        "--threshold",
        "0",
        "--null_model",
        "FE",
        "--nm_n_reps",
        "1",
        "--tag",
        "cli",
    )
    assert (cooc_out / "cli_global_nodes.tsv").exists()

    structure_out = tmp_path / "workflow_structure"
    run_cli(
        "structure",
        "--search_mode",
        "taxa_context",
        "--search_string",
        "g__Rhizo",
        "--custom_ingredients",
        raw,
        "--output_dir",
        structure_out,
        "--null_model",
        "FE",
        "--nm_n_reps",
        "1",
        "--tag",
        "cli",
    )
    assert (structure_out / "cli_global_structure.tsv").exists()

    biome_out = tmp_path / "workflow_biome"
    run_cli(
        "biome_distribution",
        "--custom_ingredients",
        raw,
        "--output_dir",
        biome_out,
        "--tag",
        "cli",
        "--return_all_taxa",
    )
    assert (biome_out / "cli_taxa_biome_distribution.tsv").exists()


@pytest.mark.cli
def test_cli_full_workflow_accepts_metadata_file(tmp_path, cli_formatted_dir, metadata_file):
    raw = cli_formatted_dir / "ingredients_raw_cli.pkl"
    out = tmp_path / "workflow_metadata"
    run_cli(
        "association",
        "--search_mode",
        "metadata",
        "--search_string",
        "soil",
        "--metadata_file",
        metadata_file,
        "--column_names",
        "env_biome_sam",
        "--custom_ingredients",
        raw,
        "--output_dir",
        out,
        "--null_model",
        "FE",
        "--nm_n_reps",
        "1",
        "--tag",
        "cli",
    )
    assert (out / "cli_global_association.tsv").exists()


@pytest.mark.cli
def test_cli_focal_lhs_rhs_cooccurrence_outputs_metrics(tmp_path, cli_formatted_dir):
    raw = cli_formatted_dir / "ingredients_raw_cli.pkl"
    out = tmp_path / "workflow_focal_rhs"

    run_cli(
        "cooccurrence",
        "--search_mode",
        "focal_taxa",
        "--search_string",
        "s__rhizo_000 -> g__Micro",
        "--custom_ingredients",
        raw,
        "--output_dir",
        out,
        "--large",
        "--threshold",
        "0",
        "--null_model",
        "FE",
        "--nm_n_reps",
        "1",
        "--tag",
        "cli",
    )

    edges = pd.read_csv(out / "cli_global_edges.tsv", sep="\t")
    assert not edges.empty
    assert {"focal_query", "focal_taxon"}.issubset(edges.columns)
    assert set(edges["focal_query"]) == {"s__rhizo_000"}
    assert edges["focal_taxon"].str.endswith("g__Rhizo; s__rhizo_000").all()
    assert edges["A_taxon"].str.endswith("g__Rhizo; s__rhizo_000").all()
    assert edges["B_taxon"].str.contains("g__Micro; s__micro_", regex=False).all()

    micro_000 = edges.loc[edges["B_taxon"].str.endswith("g__Micro; s__micro_000")].iloc[0]
    assert micro_000["A_B_intersection_count"] == 36
    assert micro_000["A_total_count"] == 60
    assert micro_000["B_total_count"] == 51
    assert micro_000["a"] == 36
    assert micro_000["b"] == 24
    assert micro_000["c"] == 15
    assert micro_000["d"] == 25
    assert micro_000["N"] == 100
    assert micro_000["jaccard"] == pytest.approx(36 / 75)
    assert micro_000["P_B_given_A"] == pytest.approx(36 / 60)
    assert micro_000["P_A_given_B"] == pytest.approx(36 / 51)
    assert 0 <= micro_000["q_bh"] <= 1
