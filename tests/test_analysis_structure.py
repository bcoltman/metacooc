from __future__ import annotations

import numpy as np

from metacooc.analysis import association_obj, cooccurrence_obj
from metacooc.structure import _structure_core, structure_obj


def test_association_obj_fe_and_ee(raw_ingredients):
    filtered = raw_ingredients.filtered_samples([True, True, True, False, False, False])

    fe = association_obj(
        raw_ingredients,
        filtered,
        threshold=0.0,
        null_model="FE",
        nm_n_reps=1,
        nm_random_state=3,
    )
    assert not fe.empty
    assert {"taxon", "a", "b", "c", "d", "p_T_given_X", "jaccard", "q_bh"}.issubset(fe.columns)
    assert not any(c.startswith("jaccard_null_mean_FE") for c in fe.columns)

    ee = association_obj(
        raw_ingredients,
        filtered,
        threshold=0.0,
        null_model="EE",
        nm_n_reps=1,
        nm_random_state=3,
    )
    assert "jaccard_null_mean_EE" in ee.columns
    assert "jaccard_p_EE" in ee.columns
    assert np.isfinite(ee["jaccard"].to_numpy()).all()


def test_cooccurrence_obj_fe_and_ee(raw_ingredients):
    taxa_universe = list(raw_ingredients.taxa)

    edge_arrays, nodes_df = cooccurrence_obj(
        raw_ingredients,
        taxa_universe,
        large=True,
        threshold=0.0,
        null_model="FE",
        nm_n_reps=1,
        nm_random_state=4,
    )
    assert nodes_df is not None
    assert not nodes_df.empty
    assert edge_arrays is not None
    assert edge_arrays.n_rows > 0
    assert "jaccard" in edge_arrays.cols

    edge_arrays_ee, _ = cooccurrence_obj(
        raw_ingredients,
        taxa_universe,
        large=True,
        threshold=0.0,
        null_model="EE",
        nm_n_reps=1,
        nm_random_state=4,
    )
    assert edge_arrays_ee is not None
    assert "jaccard_null_mean_EE" in edge_arrays_ee.cols
    assert "jaccard_p_EE" in edge_arrays_ee.cols


def test_structure_obj_observed_only(raw_ingredients):
    out = structure_obj(raw_ingredients, compute_null=False)
    assert out["metric"].tolist() == ["c_score", "mean_jaccard", "nodf"]
    assert "obs" in out.columns


def test_structure_core_all_null_models(raw_ingredients):
    for model in ("FE", "EF", "EE", "FF"):
        out = _structure_core(
            raw_ingredients,
            null_model=model,
            nm_n_reps=1,
            nm_random_state=5,
            compute_null=True,
            nm_burn_in_steps=2,
            nm_steps_per_rep=2,
            nm_progress_every=1,
        )
        assert out["metric"].tolist() == ["c_score", "mean_jaccard", "nodf"]
        assert f"null_mean_{model}" in out.columns


def test_biome_distribution(raw_ingredients):
    biomes, presence, coverage, n_dropped = raw_ingredients.biome_distribution()
    assert biomes == ["terrestrial", "aquatic"]
    assert presence.shape == (2, 6)
    assert coverage.shape == (2, 6)
    assert n_dropped == 0
