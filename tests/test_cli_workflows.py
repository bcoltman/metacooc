from __future__ import annotations

import importlib
import subprocess
import sys
import multiprocessing as mp

import pandas as pd
import pytest

from metacooc.output import COMPACT_NULL_METADATA_COLUMNS

SOIL_SAMPLE_LINES = [f"S{i:03d}" for i in range(1, 51)]
RHIZO_SAMPLE_LINES = [f"S{i:03d}" for i in range(1, 61)]


def _assert_compact_null_metadata(result, *, model, replicates, failed, seed):
    assert result.columns.tolist()[-4:] == COMPACT_NULL_METADATA_COLUMNS
    assert set(result["null_model"]) == {model}
    assert set(result["null_replicates"]) == {replicates}
    assert set(result["null_replicates_failed"]) == {failed}
    assert set(result["null_seed"]) == {seed}


def run_cli(*args, cwd=None):
    return subprocess.run(
        [sys.executable, "-m", "metacooc", *map(str, args)],
        cwd=cwd,
        check=True,
        capture_output=True,
        text=True,
    )


def run_cli_no_check(*args, cwd=None):
    return subprocess.run(
        [sys.executable, "-m", "metacooc", *map(str, args)],
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
    )


@pytest.mark.cli
def test_cli_import_does_not_create_package_data_dir(monkeypatch):
    import pathlib
    import metacooc.cli as cli

    def fail_mkdir(*args, **kwargs):
        raise AssertionError("metacooc.cli should not create directories at import time")

    monkeypatch.setattr(pathlib.Path, "mkdir", fail_mkdir)
    importlib.reload(cli)


@pytest.mark.cli
def test_cli_data_dir_default_uses_env_var(monkeypatch, tmp_path):
    from metacooc.cli import build_parser

    expected = tmp_path / "metacooc-env-data"
    monkeypatch.setenv("METACOOC_DATA_DIR", str(expected))

    parser = build_parser()
    args = parser.parse_args(["download", "--list_data_releases"])

    assert args.data_dir == str(expected)


@pytest.mark.cli
def test_cli_data_dir_default_uses_platformdirs(monkeypatch, tmp_path):
    from metacooc.cli import build_parser
    import metacooc._paths as paths

    monkeypatch.delenv("METACOOC_DATA_DIR", raising=False)
    monkeypatch.setattr(
        paths,
        "user_data_dir",
        lambda *args, **kwargs: str(tmp_path / "share" / "metacooc"),
    )

    parser = build_parser()
    args = parser.parse_args(["download", "--list_data_releases"])

    assert args.data_dir == str(tmp_path / "share" / "metacooc" / "data")


@pytest.mark.cli
def test_cli_uses_underscore_option_names_and_rejects_legacy_flags():
    from metacooc.cli import build_parser

    parser = build_parser()
    args = parser.parse_args(["download", "--data_release", "R226_gtdb_rev1"])
    assert args.data_release == "R226_gtdb_rev1"
    assert args.include_metadata is False

    args = parser.parse_args(
        ["download", "--data_release", "R226_gtdb_rev1", "--include_metadata"]
    )
    assert args.include_metadata is True

    with pytest.raises(SystemExit):
        parser.parse_args(["download", "--data_version", "2.0.0_gtdb"])

    with pytest.raises(SystemExit):
        parser.parse_args(["download", "--list_data_versions"])

    for argv in (
        ["download", "--data-release", "R226_gtdb_rev1"],
        ["download", "--list-data-releases"],
        ["download", "--include-metadata"],
    ):
        with pytest.raises(SystemExit):
            parser.parse_args(argv)

    invalid_release = run_cli_no_check(
        "download",
        "--data_release",
        "R226_gtdb",
    )
    assert invalid_release.returncode != 0
    assert "Traceback" not in invalid_release.stderr
    assert "must match R<number>_<variant>_rev<number>" in invalid_release.stderr


@pytest.mark.cli
@pytest.mark.parametrize(
    "argv",
    [
        ["download"],
        [
            "search",
            "--search_mode",
            "biome",
            "--search_string",
            "soil",
            "--output_dir",
            "results",
        ],
        ["filter", "--output_dir", "results", "--min_taxa_count", "2"],
        [
            "cooccurrence",
            "--search_mode",
            "taxa_context",
            "--search_string",
            "g__Rhizo",
            "--output_dir",
            "results",
        ],
        [
            "association",
            "--search_mode",
            "biome",
            "--search_string",
            "soil",
            "--output_dir",
            "results",
        ],
        [
            "structure",
            "--search_mode",
            "biome",
            "--search_string",
            "soil",
            "--output_dir",
            "results",
        ],
        ["biome_distribution", "--output_dir", "results"],
    ],
)
def test_cli_applies_registry_default_to_published_data_commands(
    monkeypatch,
    capsys,
    registry_factory,
    argv,
):
    import metacooc._data_config as config
    from metacooc.cli import build_parser, prepare_cli_args

    registry = registry_factory()
    monkeypatch.setattr(config, "load_registry", lambda: registry)
    args = build_parser().parse_args(argv)

    prepare_cli_args(args)

    assert args.data_release == "R226_globdb_rev1"
    assert "Using default data release: R226_globdb_rev1" in capsys.readouterr().out


@pytest.mark.cli
def test_cli_custom_and_explicit_metadata_sources_suppress_default(
    tmp_path,
    capsys,
    monkeypatch,
    registry_factory,
):
    import metacooc._data_config as config
    from metacooc.cli import build_parser, prepare_cli_args

    custom = tmp_path / "ingredients"
    custom.mkdir()
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("acc\tterm\n", encoding="utf-8")
    parser = build_parser()

    custom_args = parser.parse_args(
        [
            "search",
            "--search_mode",
            "biome",
            "--search_string",
            "soil",
            "--output_dir",
            str(tmp_path / "biome"),
            "--custom_ingredients",
            str(custom),
        ]
    )
    prepare_cli_args(custom_args)
    assert custom_args.data_release is None

    metadata_args = parser.parse_args(
        [
            "search",
            "--search_mode",
            "metadata",
            "--search_string",
            "soil",
            "--output_dir",
            str(tmp_path / "metadata"),
            "--metadata_file",
            str(metadata),
        ]
    )
    prepare_cli_args(metadata_args)
    assert metadata_args.data_release is None
    assert "Using default data release" not in capsys.readouterr().out

    registry = registry_factory()
    monkeypatch.setattr(config, "load_registry", lambda: registry)
    workflow_args = parser.parse_args(
        [
            "association",
            "--search_mode",
            "metadata",
            "--search_string",
            "soil",
            "--output_dir",
            str(tmp_path / "association"),
            "--metadata_file",
            str(metadata),
        ]
    )
    prepare_cli_args(workflow_args)
    assert workflow_args.data_release == "R226_globdb_rev1"

    tax_profile = tmp_path / "profile.tsv"
    tax_profile.write_text("sample\ttaxonomy\tcoverage\n", encoding="utf-8")
    format_args = parser.parse_args(
        [
            "format",
            "--tax_profile",
            str(tax_profile),
            "--output_dir",
            str(tmp_path / "formatted"),
        ]
    )
    prepare_cli_args(format_args)
    assert format_args.data_release is None

    list_args = parser.parse_args(["download", "--list_data_releases"])
    prepare_cli_args(list_args)
    assert list_args.data_release is None


@pytest.mark.cli
def test_cli_contract_rejects_incompatible_and_stale_options(tmp_path):
    from metacooc.cli import build_parser, prepare_cli_args

    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "plot",
                "--analysis_file",
                str(tmp_path / "analysis.tsv"),
                "--output_dir",
                str(tmp_path),
                "--analysis_type",
                "association",
                "--aggregated",
            ]
        )
    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "analysis",
                "--analysis_type",
                "structure",
                "--filtered_file",
                str(tmp_path / "ingredients"),
                "--output_dir",
                str(tmp_path),
                "--aggregated",
            ]
        )

    custom = tmp_path / "ingredients"
    custom.mkdir()
    args = parser.parse_args(
        [
            "cooccurrence",
            "--search_mode",
            "taxa_context",
            "--search_string",
            "g__Rhizo",
            "--output_dir",
            str(tmp_path / "cooc"),
            "--custom_ingredients",
            str(custom),
            "--null_scope",
            "taxa",
        ]
    )
    with pytest.raises(SystemExit):
        prepare_cli_args(args)

    args = parser.parse_args(
        [
            "search",
            "--search_mode",
            "focal_taxa",
            "--search_string",
            "g__Rhizo -> g__Micro",
            "--output_dir",
            str(tmp_path / "search"),
        ]
    )
    with pytest.raises(SystemExit):
        prepare_cli_args(args)

    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("acc\torganism\n", encoding="utf-8")
    args = parser.parse_args(
        [
            "search",
            "--search_mode",
            "metadata",
            "--search_string",
            "soil",
            "--output_dir",
            str(tmp_path / "metadata-search"),
            "--metadata_file",
            str(metadata),
            "--strict",
            "--column_names",
            "organism",
        ]
    )
    with pytest.raises(SystemExit):
        prepare_cli_args(args)


@pytest.mark.cli
def test_cli_missing_custom_path_and_unknown_command_do_not_traceback(tmp_path):
    missing = run_cli_no_check(
        "search",
        "--search_mode",
        "biome",
        "--search_string",
        "soil",
        "--output_dir",
        tmp_path / "out",
        "--custom_ingredients",
        tmp_path / "missing",
    )
    assert missing.returncode != 0
    assert "must be an existing directory" in missing.stderr
    assert "Traceback" not in missing.stderr

    unknown = run_cli_no_check("not-a-command")
    assert unknown.returncode == 2
    assert "Unknown command" in unknown.stderr
    assert "Traceback" not in unknown.stderr

@pytest.mark.cli
def test_cli_conditional_probability_flags_parse_and_reject_old_names(tmp_path):
    from metacooc.cli import build_parser

    parser = build_parser()

    assoc_args = parser.parse_args(
        [
            "association",
            "--search_mode",
            "taxa_context",
            "--search_string",
            "g__Rhizo",
            "--output_dir",
            str(tmp_path / "assoc"),
            "--min_conditional_probability",
            "0.2",
        ]
    )
    assert assoc_args.min_conditional_probability == 0.2

    cooc_args = parser.parse_args(
        [
            "cooccurrence",
            "--search_mode",
            "taxa_context",
            "--search_string",
            "g__Rhizo",
            "--output_dir",
            str(tmp_path / "cooc"),
            "--min_conditional_probability",
            "0.3",
        ]
    )
    assert cooc_args.min_conditional_probability == 0.3

    plot_args = parser.parse_args(
        [
            "plot",
            "--analysis_file",
            str(tmp_path / "analysis.tsv"),
            "--output_dir",
            str(tmp_path / "plot"),
            "--analysis_type",
            "cooccurrence",
            "--q_threshold",
            "0.4",
            "--q_metric",
            "chi2_q_value_bh",
            "--plot_all",
            "--label_top_n",
            "4",
            "--x_metric",
            "p_target_given_source",
            "--y_metric",
            "phi_coefficient",
        ]
    )
    assert plot_args.q_threshold == 0.4
    assert plot_args.analysis_type == "cooccurrence"
    assert plot_args.q_metric == "chi2_q_value_bh"
    assert plot_args.plot_all is True
    assert plot_args.label_top_n == 4

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "plot",
                "--analysis_file",
                str(tmp_path / "analysis.tsv"),
                "--output_dir",
                str(tmp_path / "missing_type"),
            ]
        )

    unpaired_plot_args = parser.parse_args(
        [
            "plot",
            "--analysis_file",
            str(tmp_path / "analysis.tsv"),
            "--output_dir",
            str(tmp_path / "unpaired"),
            "--analysis_type",
            "association",
            "--x_metric",
            "phi_coefficient",
        ]
    )
    with pytest.raises(SystemExit):
        unpaired_plot_args.func(unpaired_plot_args)

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "association",
                "--search_mode",
                "taxa_context",
                "--search_string",
                "g__Rhizo",
                "--output_dir",
                str(tmp_path / "old"),
                "--threshold",
                "0.1",
            ]
        )

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "cooccurrence",
                "--search_mode",
                "taxa_context",
                "--search_string",
                "g__Rhizo",
                "--output_dir",
                str(tmp_path / "old_cooc"),
                "--threshold",
                "0.1",
            ]
        )

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "plot",
                "--analysis_file",
                str(tmp_path / "analysis.tsv"),
                "--output_dir",
                str(tmp_path / "old_plot"),
                "--analysis_type",
                "association",
                "--threshold",
                "0.1",
            ]
        )


@pytest.mark.cli
def test_cli_presence_threshold_flags_parse_and_validate(tmp_path):
    from metacooc.cli import build_parser

    parser = build_parser()

    args = parser.parse_args(
        [
            "association",
            "--search_mode",
            "taxa_context",
            "--search_string",
            "g__Rhizo",
            "--output_dir",
            str(tmp_path / "assoc"),
            "--min_coverage",
            "1.5",
            "--min_coverage_by_rank",
            "species=0.2",
            "genus=0.1",
            "--min_relative_abundance",
            "0.02",
            "--min_relative_abundance_by_rank",
            "species=0.03",
        ]
    )

    assert args.min_coverage == 1.5
    assert args.min_coverage_by_rank == [("species", 0.2), ("genus", 0.1)]
    assert args.min_relative_abundance == 0.02
    assert args.min_relative_abundance_by_rank == [("species", 0.03)]

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "association",
                "--search_mode",
                "taxa_context",
                "--search_string",
                "g__Rhizo",
                "--output_dir",
                str(tmp_path / "bad_rank"),
                "--min_coverage_by_rank",
                "strain=1",
            ]
        )

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "association",
                "--search_mode",
                "taxa_context",
                "--search_string",
                "g__Rhizo",
                "--output_dir",
                str(tmp_path / "bad_cov"),
                "--min_coverage",
                "-1",
            ]
        )

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "association",
                "--search_mode",
                "taxa_context",
                "--search_string",
                "g__Rhizo",
                "--output_dir",
                str(tmp_path / "bad_rel"),
                "--min_relative_abundance",
                "1.1",
            ]
        )


@pytest.mark.cli
def test_cli_filter_accepts_presence_threshold_as_filter(monkeypatch, tmp_path):
    from metacooc.cli import build_parser
    import metacooc.filter as filter_module

    captured = {}

    def fake_filter_data(**kwargs):
        captured.update(kwargs)

    monkeypatch.setattr(filter_module, "filter_data", fake_filter_data)

    parser = build_parser()
    args = parser.parse_args(
        [
            "filter",
            "--custom_ingredients",
            str(tmp_path / "ingredients"),
            "--output_dir",
            str(tmp_path / "filter"),
            "--min_coverage",
            "1.0",
        ]
    )
    args.func(args)

    assert captured["min_coverage"] == 1.0
    assert captured["min_relative_abundance"] is None


@pytest.mark.cli
def test_cli_analysis_cutoff_flag_parses_and_is_rejected_for_structure(tmp_path):
    from metacooc.cli import build_parser

    parser = build_parser()

    args = parser.parse_args(
        [
            "analysis",
            "--analysis_type",
            "association",
            "--null_file",
            str(tmp_path / "null"),
            "--filtered_file",
            str(tmp_path / "filtered"),
            "--output_dir",
            str(tmp_path / "assoc"),
            "--min_conditional_probability",
            "0.1",
        ]
    )
    assert args.min_conditional_probability == 0.1

    structure_args = parser.parse_args(
        [
            "analysis",
            "--analysis_type",
            "structure",
            "--filtered_file",
            str(tmp_path / "filtered"),
            "--output_dir",
            str(tmp_path / "structure"),
            "--min_conditional_probability",
            "0.1",
        ]
    )
    with pytest.raises(SystemExit):
        structure_args.func(structure_args)


@pytest.mark.cli
def test_cli_null_hpc_options_parse_and_forward(monkeypatch, tmp_path):
    from metacooc.cli import build_parser
    import metacooc.pipelines as pipelines

    captured = {}

    def fake_run_structure(args):
        captured.update(vars(args))

    monkeypatch.setattr(pipelines, "run_structure", fake_run_structure)
    mp_start = mp.get_all_start_methods()[0]

    parser = build_parser()
    args = parser.parse_args(
        [
            "structure",
            "--search_mode",
            "taxa_context",
            "--search_string",
            "g__Rhizo",
            "--output_dir",
            str(tmp_path),
            "--nm_seed",
            "123",
            "--nm_n_workers",
            "2",
            "--nm_mp_start",
            mp_start,
            "--nm_burn_in_steps",
            "3",
            "--nm_steps_per_rep",
            "4",
            "--nm_progress_every",
            "5",
        ]
    )
    args.func(args)

    assert captured["nm_seed"] == 123
    assert captured["nm_n_workers"] == 2
    assert captured["nm_mp_start"] == mp_start
    assert captured["nm_burn_in_steps"] == 3
    assert captured["nm_steps_per_rep"] == 4
    assert captured["nm_progress_every"] == 5


@pytest.mark.cli
def test_cli_biome_distribution_filter_defaults_are_command_specific(tmp_path):
    from metacooc.cli import build_parser

    parser = build_parser()
    biome_args = parser.parse_args(
        [
            "biome_distribution",
            "--output_dir",
            str(tmp_path / "biome"),
        ]
    )
    assert biome_args.min_taxa_count is None
    assert biome_args.min_sample_count is None
    assert biome_args.taxa_count_rank == "species"

    structure_args = parser.parse_args(
        [
            "structure",
            "--search_mode",
            "taxa_context",
            "--search_string",
            "g__Rhizo",
            "--output_dir",
            str(tmp_path / "structure"),
        ]
    )
    assert structure_args.min_taxa_count == 1
    assert structure_args.min_sample_count == 1
    assert structure_args.taxa_count_rank == "species"


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
        "--data_release",
        "test",
    )
    assert (out / "ingredients_raw_cli").exists()
    assert (out / "ingredients_aggregated_cli").exists()
    return out


@pytest.mark.cli
def test_cli_rejects_custom_ingredients_with_data_release(
    tmp_path, cli_formatted_dir
):
    result = run_cli_no_check(
        "search",
        "--search_mode",
        "taxa_context",
        "--search_string",
        "g__Rhizo",
        "--custom_ingredients",
        cli_formatted_dir / "ingredients_raw_cli",
        "--data_release",
        "R226_gtdb_rev1",
        "--output_dir",
        tmp_path / "invalid",
    )
    assert result.returncode != 0
    assert "mutually exclusive" in result.stderr


@pytest.mark.cli
def test_cli_search_taxon_biome_and_metadata(tmp_path, cli_formatted_dir, metadata_file):
    raw = cli_formatted_dir / "ingredients_raw_cli"

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
def test_cli_search_lists_available_biomes(cli_formatted_dir):
    raw = cli_formatted_dir / "ingredients_raw_cli"

    result = run_cli(
        "search",
        "--list_biomes",
        "--custom_ingredients",
        raw,
    )

    assert "level_1:" in result.stdout
    assert "  terrestrial (50 samples)" in result.stdout
    assert "  aquatic (50 samples)" in result.stdout
    assert "level_2:" in result.stdout
    assert "  soil (50 samples)" in result.stdout
    assert "  marine (50 samples)" in result.stdout


@pytest.mark.cli
def test_cli_filter_and_analysis_commands(tmp_path, cli_formatted_dir):
    raw = cli_formatted_dir / "ingredients_raw_cli"
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
    null_file = filter_out / "cli_ingredients_null"
    filtered_file = filter_out / "cli_ingredients_filtered"
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
        "--min_conditional_probability",
        "0",
    )
    assert not pd.read_csv(assoc_out / "association.tsv", sep="\t").empty
    assert (assoc_out / "association_summary.tsv").exists()

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
        "--min_conditional_probability",
        "0",
        "--null_model",
        "EE",
        "--nm_n_reps",
        "1",
    )
    assert (cooc_out / "taxon_nodes.tsv").exists()
    assert (cooc_out / "taxon_edges_summary.tsv").exists()

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
        "--nm_seed",
        "7",
    )
    structure = pd.read_csv(structure_out / "structure.tsv", sep="\t")
    assert "null_mean" in structure.columns
    _assert_compact_null_metadata(
        structure,
        model="FE",
        replicates=1,
        failed=0,
        seed=7,
    )
    assert not (structure_out / "structure_metadata.tsv").exists()


@pytest.mark.cli
def test_cli_full_workflow_commands(tmp_path, cli_formatted_dir):
    raw = cli_formatted_dir / "ingredients_raw_cli"

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
        "EE",
        "--nm_n_reps",
        "1",
        "--nm_seed",
        "7",
        "--min_conditional_probability",
        "0",
        "--tag",
        "cli",
    )
    association = pd.read_csv(assoc_out / "cli_global_association.tsv", sep="\t")
    association_summary = pd.read_csv(
        assoc_out / "cli_global_association_summary.tsv",
        sep="\t",
    )
    assert "jaccard_null_mean" in association.columns
    assert "jaccard_null_mean" not in association_summary.columns
    assert "jaccard_null_ses" in association_summary.columns
    _assert_compact_null_metadata(
        association,
        model="EE",
        replicates=1,
        failed=0,
        seed=7,
    )
    _assert_compact_null_metadata(
        association_summary,
        model="EE",
        replicates=1,
        failed=0,
        seed=7,
    )
    assert (assoc_out / "cli_global_association_plot.png").exists()
    assert not (assoc_out / "cli_global_association_metadata.tsv").exists()

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
        "--min_conditional_probability",
        "0",
        "--null_model",
        "EE",
        "--nm_n_reps",
        "1",
        "--nm_seed",
        "7",
        "--tag",
        "cli",
    )
    assert (cooc_out / "cli_global_nodes.tsv").exists()
    assert (cooc_out / "cli_global_edges_summary.tsv").exists()
    assert (cooc_out / "cli_global_cooccurrence_plot.png").exists()

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
        "--nm_seed",
        "7",
        "--tag",
        "cli",
    )
    structure = pd.read_csv(structure_out / "cli_global_structure.tsv", sep="\t")
    assert "null_mean" in structure.columns
    _assert_compact_null_metadata(
        structure,
        model="FE",
        replicates=1,
        failed=0,
        seed=7,
    )
    assert not (structure_out / "cli_global_structure_metadata.tsv").exists()

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
    biome = pd.read_csv(biome_out / "cli_taxa_biome_distribution.tsv", sep="\t")
    assert biome.columns.tolist() == ["taxon", "terrestrial", "aquatic"]
    assert not (biome_out / "cli_taxa_biome_distribution_metadata.tsv").exists()

    filtered_biome_out = tmp_path / "workflow_biome_filtered"
    run_cli(
        "biome_distribution",
        "--custom_ingredients",
        raw,
        "--output_dir",
        filtered_biome_out,
        "--tag",
        "cli",
        "--taxa_query",
        "s__rhizo_000,s__micro_000",
        "--biome_level",
        "level_2",
        "--min_taxa_count",
        "2",
        "--taxa_count_rank",
        "species",
    )
    filtered_biome = pd.read_csv(
        filtered_biome_out / "cli_taxa_biome_distribution.tsv",
        sep="\t",
    )
    assert filtered_biome.columns.tolist() == ["taxon", "soil", "marine"]
    assert len(filtered_biome) == 2

    invalid_biome = run_cli_no_check(
        "biome_distribution",
        "--custom_ingredients",
        raw,
        "--output_dir",
        tmp_path / "workflow_biome_invalid",
        "--taxa_query",
        "g__Rhizo|g__Micro",
    )
    assert invalid_biome.returncode != 0
    assert "focal_taxa mode does not support" in invalid_biome.stderr


@pytest.mark.cli
def test_cli_full_workflow_accepts_metadata_file(tmp_path, cli_formatted_dir, metadata_file):
    raw = cli_formatted_dir / "ingredients_raw_cli"
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
    raw = cli_formatted_dir / "ingredients_raw_cli"
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
        "--min_conditional_probability",
        "0",
        "--null_model",
        "EE",
        "--nm_n_reps",
        "1",
        "--nm_seed",
        "7",
        "--tag",
        "cli",
    )

    edges = pd.read_csv(out / "cli_global_edges.tsv", sep="\t")
    summary = pd.read_csv(out / "cli_global_edges_summary.tsv", sep="\t")
    assert not edges.empty
    assert summary.columns.tolist()[:4] == [
        "focal_query",
        "focal_taxon",
        "source_taxon",
        "target_taxon",
    ]
    assert len(summary) == len(edges)
    assert {"focal_query", "focal_taxon"}.issubset(edges.columns)
    assert set(edges["focal_query"]) == {"s__rhizo_000"}
    assert edges["focal_taxon"].str.endswith("g__Rhizo; s__rhizo_000").all()
    assert edges["source_taxon"].str.endswith("g__Rhizo; s__rhizo_000").all()
    assert edges["target_taxon"].str.contains("g__Micro; s__micro_", regex=False).all()

    micro_000 = edges.loc[edges["target_taxon"].str.endswith("g__Micro; s__micro_000")].iloc[0]
    assert micro_000["shared_sample_count"] == 36
    assert micro_000["source_taxon_sample_count"] == 60
    assert micro_000["target_taxon_sample_count"] == 51
    assert micro_000["source_only_sample_count"] == 24
    assert micro_000["target_only_sample_count"] == 15
    assert micro_000["neither_source_nor_target_sample_count"] == 25
    assert micro_000["background_sample_count"] == 100
    assert micro_000["jaccard_taxon_pair"] == pytest.approx(36 / 75)
    assert micro_000["p_target_given_source"] == pytest.approx(36 / 60)
    assert micro_000["p_source_given_target"] == pytest.approx(36 / 51)
    assert 0 <= micro_000["chi2_q_value_bh"] <= 1
    assert "jaccard_null_mean" in edges.columns
    assert "jaccard_null_p_empirical" in edges.columns
    _assert_compact_null_metadata(
        edges,
        model="EE",
        replicates=1,
        failed=0,
        seed=7,
    )
    _assert_compact_null_metadata(
        summary,
        model="EE",
        replicates=1,
        failed=0,
        seed=7,
    )
    assert not (out / "cli_global_edges_metadata.tsv").exists()


@pytest.mark.cli
def test_cli_cooccurrence_warns_for_large_null_edge_request(
    monkeypatch,
    capsys,
    tmp_path,
    cli_formatted_dir,
):
    from metacooc.cli import build_parser
    import metacooc.analysis as analysis

    raw = cli_formatted_dir / "ingredients_raw_cli"
    out = tmp_path / "workflow_warn_large_null"
    mp_start = mp.get_all_start_methods()[0]
    monkeypatch.setattr(analysis, "_COOCCURRENCE_NULL_EDGE_WARNING_THRESHOLD", 1)

    parser = build_parser()
    args = parser.parse_args(
        [
            "cooccurrence",
            "--search_mode",
            "focal_taxa",
            "--search_string",
            "s__rhizo_000 -> g__Micro",
            "--custom_ingredients",
            str(raw),
            "--output_dir",
            str(out),
            "--large",
            "--min_conditional_probability",
            "0",
            "--null_model",
            "EE",
            "--nm_n_reps",
            "1",
            "--nm_n_workers",
            "1",
            "--nm_mp_start",
            mp_start,
            "--tag",
            "cli",
        ]
    )
    args.func(args)

    captured = capsys.readouterr()
    assert "WARNING: Cooccurrence null simulation will evaluate" in captured.out
    assert "bounded dot-product batching" in captured.out
    assert (out / "cli_global_edges.tsv").exists()


@pytest.mark.cli
def test_cli_invalid_query_grammar_fails_clearly(tmp_path, cli_formatted_dir):
    raw = cli_formatted_dir / "ingredients_raw_cli"

    taxon_mode = run_cli_no_check(
        "search",
        "--search_mode",
        "taxon",
        "--search_string",
        "g__Rhizo",
        "--custom_ingredients",
        raw,
        "--output_dir",
        tmp_path / "bad_taxon",
    )
    assert taxon_mode.returncode != 0
    assert "invalid choice: 'taxon'" in taxon_mode.stderr

    invalid_focal = run_cli_no_check(
        "cooccurrence",
        "--search_mode",
        "focal_taxa",
        "--search_string",
        "g__Rhizo+g__Micro",
        "--custom_ingredients",
        raw,
        "--output_dir",
        tmp_path / "bad_focal",
        "--large",
    )
    assert invalid_focal.returncode != 0
    assert "focal_taxa mode does not support" in invalid_focal.stderr

    missing_metadata_source = run_cli_no_check(
        "search",
        "--search_mode",
        "metadata",
        "--search_string",
        "soil",
        "--output_dir",
        tmp_path / "bad_metadata",
    )
    assert missing_metadata_source.returncode != 0
    assert "Traceback" not in missing_metadata_source.stderr
    assert "R226_globdb_rev1" in missing_metadata_source.stderr
    assert "--include_metadata" in missing_metadata_source.stderr
