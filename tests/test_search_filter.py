from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from metacooc.filter import filter_data, filter_data_obj
from metacooc.pantry import Ingredients, load_ingredients
from metacooc.search import search_data_obj, search_in_metadata


SOIL_SAMPLES = {f"S{i:03d}" for i in range(1, 51)}
RHIZO_HITS = {f"S{i:03d}" for i in range(1, 61)}
MICRO_HITS = {f"S{i:03d}" for i in range(25, 76)}
RHIZO_MICRO_HITS = RHIZO_HITS & MICRO_HITS
MARINE_SAMPLES = {f"S{i:03d}" for i in range(51, 101)}


def test_taxon_biome_and_metadata_search(raw_ingredients, metadata_file):
    taxon_hits = search_data_obj(
        search_mode="taxa_context",
        search_string="g__Rhizo",
        custom_ingredients=raw_ingredients,
    )
    assert taxon_hits == RHIZO_HITS

    and_hits = search_data_obj(
        search_mode="taxa_context",
        search_string="g__Rhizo+g__Micro",
        custom_ingredients=raw_ingredients,
    )
    assert and_hits == RHIZO_MICRO_HITS

    inverse_hits = search_data_obj(
        search_mode="taxa_context",
        search_string="g__Rhizo",
        custom_ingredients=raw_ingredients,
        inverse=True,
    )
    assert inverse_hits == set(raw_ingredients.samples) - RHIZO_HITS

    biome_hits = search_data_obj(
        search_mode="biome",
        search_string="soil",
        custom_ingredients=raw_ingredients,
    )
    assert biome_hits == SOIL_SAMPLES

    metadata_hits = search_data_obj(
        search_mode="metadata",
        search_string="soil",
        metadata_file=str(metadata_file),
        column_names=["env_biome_sam"],
    )
    assert metadata_hits == SOIL_SAMPLES


def test_metadata_search_python_backend(monkeypatch, metadata_file):
    monkeypatch.setenv("METACOOC_METADATA_SEARCH_BACKEND", "python")

    assert search_in_metadata(
        str(metadata_file),
        "SOIL",
        column_names=["env_biome_sam"],
    ) == SOIL_SAMPLES
    assert search_in_metadata(
        str(metadata_file),
        "soil",
        column_names=["env_biome_sam"],
        inverse=True,
    ) == MARINE_SAMPLES


def test_metadata_search_auto_falls_back_to_python(monkeypatch, metadata_file):
    import metacooc.search as search

    monkeypatch.setenv("METACOOC_METADATA_SEARCH_BACKEND", "auto")
    monkeypatch.setattr(search.shutil, "which", lambda tool: None)

    hits = search_data_obj(
        search_mode="metadata",
        search_string="soil",
        metadata_file=str(metadata_file),
        column_names=["env_biome_sam"],
    )

    assert hits == SOIL_SAMPLES


def test_metadata_search_auto_falls_back_when_external_probe_fails(
    monkeypatch, metadata_file
):
    import metacooc.search as search

    class FailedProbe:
        returncode = 2
        stdout = ""
        stderr = "probe failed"

    monkeypatch.setenv("METACOOC_METADATA_SEARCH_BACKEND", "auto")
    monkeypatch.setattr(search.shutil, "which", lambda tool: f"/usr/bin/{tool}")
    monkeypatch.setattr(search.subprocess, "run", lambda *args, **kwargs: FailedProbe())
    monkeypatch.setattr(search, "_EXTERNAL_METADATA_SEARCH_PROBE", {})

    hits = search_data_obj(
        search_mode="metadata",
        search_string="soil",
        metadata_file=str(metadata_file),
        column_names=["env_biome_sam"],
    )

    assert hits == SOIL_SAMPLES


def test_metadata_search_forced_external_missing_tools_raises(monkeypatch, metadata_file):
    import metacooc.search as search

    monkeypatch.setenv("METACOOC_METADATA_SEARCH_BACKEND", "external")
    monkeypatch.setattr(search.shutil, "which", lambda tool: None)

    with pytest.raises(RuntimeError, match="required tool"):
        search_in_metadata(
            str(metadata_file),
            "soil",
            column_names=["env_biome_sam"],
        )


def test_metadata_search_forced_external_probe_failure_raises(
    monkeypatch, metadata_file
):
    import metacooc.search as search

    class FailedProbe:
        returncode = 2
        stdout = ""
        stderr = "probe failed"

    monkeypatch.setenv("METACOOC_METADATA_SEARCH_BACKEND", "external")
    monkeypatch.setattr(search.shutil, "which", lambda tool: f"/usr/bin/{tool}")
    monkeypatch.setattr(search.subprocess, "run", lambda *args, **kwargs: FailedProbe())
    monkeypatch.setattr(search, "_EXTERNAL_METADATA_SEARCH_PROBE", {})

    with pytest.raises(RuntimeError, match="probe failed"):
        search_in_metadata(
            str(metadata_file),
            "soil",
            column_names=["env_biome_sam"],
        )


def test_metadata_search_invalid_backend_raises(monkeypatch, metadata_file):
    monkeypatch.setenv("METACOOC_METADATA_SEARCH_BACKEND", "bad")

    with pytest.raises(ValueError, match="METACOOC_METADATA_SEARCH_BACKEND"):
        search_in_metadata(
            str(metadata_file),
            "soil",
            column_names=["env_biome_sam"],
        )


def test_taxa_search_uses_coverage_threshold():
    ingredients = Ingredients(
        samples=["S1", "S2"],
        taxa=["d__Bacteria; p__P; c__C; o__O; f__F; g__G; s__a"],
        presence_matrix=sp.csr_matrix(np.ones((1, 2), dtype=np.uint8)),
        coverage_matrix=sp.csr_matrix([[0.5, 2.0]], dtype=float),
    )

    hits = search_data_obj(
        search_mode="taxa_context",
        search_string="s__a",
        custom_ingredients=ingredients,
        min_coverage=1.0,
    )

    assert hits == {"S2"}


def test_focal_lhs_rhs_query_resolution(raw_ingredients):
    resolved = search_data_obj(
        search_mode="focal_taxa",
        search_string="s__rhizo_000 -> g__Micro",
        custom_ingredients=raw_ingredients,
        return_details=True,
    )

    assert set(resolved.focal_rows_by_query_lhs) == {"s__rhizo_000"}
    assert set(resolved.focal_rows_by_query_rhs) == {"g__Micro"}
    assert resolved.focal_sample_union == RHIZO_HITS

    lhs_taxa = [raw_ingredients.taxa[i] for i in resolved.focal_rows_by_query_lhs["s__rhizo_000"]]
    rhs_taxa = [raw_ingredients.taxa[i] for i in resolved.focal_rows_by_query_rhs["g__Micro"]]

    assert len(lhs_taxa) == 1
    assert lhs_taxa[0].endswith("g__Rhizo; s__rhizo_000")
    assert len(rhs_taxa) == 50
    assert all("g__Micro; s__micro_" in taxon for taxon in rhs_taxa)


def test_focal_multi_lhs_query_resolution(raw_ingredients):
    resolved = search_data_obj(
        search_mode="focal_taxa",
        search_string="s__rhizo_000,s__micro_000",
        custom_ingredients=raw_ingredients,
        return_details=True,
    )

    assert set(resolved.focal_rows_by_query_lhs) == {"s__rhizo_000", "s__micro_000"}
    assert resolved.focal_sample_union == RHIZO_HITS | MICRO_HITS


def test_focal_endpoint_and_aggregated_resolution(raw_ingredients, aggregated_ingredients):
    species = search_data_obj(
        search_mode="focal_taxa",
        search_string="s__rhizo_000",
        custom_ingredients=raw_ingredients,
        return_details=True,
    )
    species_taxa = [
        raw_ingredients.taxa[i]
        for i in species.focal_rows_by_query_lhs["s__rhizo_000"]
    ]
    assert len(species_taxa) == 1
    assert species_taxa[0].endswith("g__Rhizo; s__rhizo_000")

    exact_aggregated = search_data_obj(
        search_mode="focal_taxa",
        search_string="g__Rhizo AGGREGATED",
        custom_ingredients=aggregated_ingredients,
        return_details=True,
    )
    aggregated_taxa = [
        aggregated_ingredients.taxa[i]
        for i in exact_aggregated.focal_rows_by_query_lhs["g__Rhizo AGGREGATED"]
    ]
    assert aggregated_taxa == [
        "Root; d__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; "
        "o__Rhizobiales; f__Rhizobiaceae; g__Rhizo AGGREGATED"
    ]
    assert exact_aggregated.focal_sample_union == RHIZO_HITS

    expanded = search_data_obj(
        search_mode="focal_taxa",
        search_string="g__Rhizo",
        custom_ingredients=aggregated_ingredients,
        return_details=True,
    )
    expanded_taxa = [
        aggregated_ingredients.taxa[i]
        for i in expanded.focal_rows_by_query_lhs["g__Rhizo"]
    ]
    assert len(expanded_taxa) == 51
    assert any(t.endswith("g__Rhizo AGGREGATED") for t in expanded_taxa)
    assert expanded.focal_sample_union == RHIZO_HITS


@pytest.mark.parametrize(
    ("search_mode", "search_string", "message"),
    [
        ("taxon", "g__Rhizo", "search_mode='taxon' has been removed"),
        ("metadata", "soil -> reef", "'->' syntax is not supported in metadata mode"),
        ("biome", "soil -> marine", "'->' syntax is not supported in biome mode"),
        ("taxa_context", "g__Rhizo -> g__Micro", "'->' syntax is not supported in taxa_context mode"),
        ("focal_taxa", "g__Rhizo+g__Micro", "focal_taxa mode does not support"),
        ("focal_taxa", "g__Rhizo|g__Micro", "focal_taxa mode does not support"),
        ("focal_taxa", "g__Rhizo -> ", "Right-hand side of '->' cannot be empty"),
        ("focal_taxa", " -> g__Micro", "Left-hand side of '->' cannot be empty"),
    ],
)
def test_invalid_query_grammar_raises(raw_ingredients, metadata_file, search_mode, search_string, message):
    kwargs = {"custom_ingredients": raw_ingredients}
    if search_mode == "metadata":
        kwargs = {"metadata_file": str(metadata_file)}

    with pytest.raises(ValueError, match=message):
        search_data_obj(
            search_mode=search_mode,
            search_string=search_string,
            **kwargs,
        )


def test_metadata_without_file_or_data_dir_raises():
    with pytest.raises(ValueError, match="data_dir must be provided"):
        search_data_obj(
            search_mode="metadata",
            search_string="soil",
        )


def test_filter_data_obj_by_counts_rank_and_accessions(raw_ingredients):
    filtered, ok = filter_data_obj(
        raw_ingredients,
        accession_set=SOIL_SAMPLES,
        min_taxa_count=2,
        min_sample_count=2,
        filter_rank="species",
        taxa_count_rank="species",
    )
    assert ok
    assert filtered.samples == [f"S{i:03d}" for i in range(1, 51)]
    assert all("s__" in t for t in filtered.taxa)
    assert filtered.presence_matrix.shape[1] == 50
    assert filtered.presence_matrix.shape[0] >= 100


def test_file_filter_writes_null_and_filtered(tmp_path, raw_ingredients_path):
    accessions = tmp_path / "accessions.txt"
    accessions.write_text("".join(f"S{i:03d}\n" for i in range(1, 51)))

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

    null_path = out / "test_ingredients_null"
    filtered_path = out / "test_ingredients_filtered"
    assert null_path.exists()
    assert filtered_path.exists()

    filtered = load_ingredients(custom_ingredients=str(filtered_path))
    assert filtered.samples == [f"S{i:03d}" for i in range(1, 51)]
