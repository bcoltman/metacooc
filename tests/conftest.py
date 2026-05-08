from __future__ import annotations

import pickle
from pathlib import Path
from types import SimpleNamespace

import pytest

from metacooc.format import format_data


@pytest.fixture
def fixture_dir() -> Path:
    return Path(__file__).parent / "data"


@pytest.fixture
def formatted_data_dir(tmp_path: Path, fixture_dir: Path) -> Path:
    out = tmp_path / "formatted"
    format_data(
        tax_profile=str(fixture_dir / "tax_profile.tsv"),
        output_dir=str(out),
        sample_to_biome_file=str(fixture_dir / "sample_to_biome.tsv"),
        aggregated=True,
        tag="test",
        data_version="test",
    )
    return out


@pytest.fixture
def raw_ingredients_path(formatted_data_dir: Path) -> Path:
    return formatted_data_dir / "ingredients_raw_test.pkl"


@pytest.fixture
def aggregated_ingredients_path(formatted_data_dir: Path) -> Path:
    return formatted_data_dir / "ingredients_aggregated_test.pkl"


@pytest.fixture
def raw_ingredients(raw_ingredients_path: Path):
    with raw_ingredients_path.open("rb") as f:
        return pickle.load(f)


@pytest.fixture
def aggregated_ingredients(aggregated_ingredients_path: Path):
    with aggregated_ingredients_path.open("rb") as f:
        return pickle.load(f)


@pytest.fixture
def metadata_file(fixture_dir: Path) -> Path:
    return fixture_dir / "metadata.tsv"


def pipeline_args(**overrides):
    base = dict(
        data_dir=None,
        data_version=None,
        aggregated=False,
        custom_ingredients=None,
        metadata_file=None,
        search_mode="taxa_context",
        search_string="g__Rhizo",
        ranks_for_search_inclusion=None,
        strict=False,
        column_names=None,
        inverse=False,
        min_taxa_count=1,
        min_sample_count=1,
        filter_rank="species",
        taxa_count_rank="species",
        remove_null_threshold=False,
        null_scope=None,
        null_biome_query=None,
        null_taxa_query=None,
        null_metadata_query=None,
        taxa_degree=1,
        min_shared_samples_between_taxa=1,
        output_dir=None,
        tag="test_",
        null_model="FE",
        nm_n_reps=1,
        nm_seed=7,
        nm_n_workers=None,
        nm_mp_start=None,
        nm_sort_indices=False,
        nm_burn_in_steps=None,
        nm_steps_per_rep=None,
        nm_progress_every=25,
        compute_fisher=False,
        threshold=0.0,
        large=True,
        max_pairs=100_000,
        return_all_taxa=True,
    )
    base.update(overrides)
    return SimpleNamespace(**base)
