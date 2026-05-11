from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from metacooc.analysis import association_obj, cooccurrence_obj
from metacooc.structure import (
    _structure_core,
    compute_c_score,
    compute_nodf_streamed,
    mean_jaccard_dot,
    structure_obj,
)


def _reference_c_score(X):
    X = X.tocsr()
    R = np.diff(X.indptr).astype(float)
    n_taxa = R.size
    if n_taxa < 2:
        return np.nan

    S = (X @ X.T).toarray()
    vals = []
    for i in range(n_taxa):
        for j in range(i + 1, n_taxa):
            vals.append((R[i] - S[i, j]) * (R[j] - S[i, j]))
    return float(np.mean(vals))


def _reference_mean_jaccard(X):
    X = X.tocsr()
    deg_all = np.diff(X.indptr)
    keep = deg_all > 0
    if int(keep.sum()) < 2:
        return np.nan

    X = X[keep, :].tocsr()
    deg = np.diff(X.indptr).astype(float)
    S = (X @ X.T).toarray()
    vals = []
    for i in range(X.shape[0]):
        for j in range(i + 1, X.shape[0]):
            union = deg[i] + deg[j] - S[i, j]
            vals.append(0.0 if union <= 0 else S[i, j] / union)
    return float(np.mean(vals))


def _reference_nodf_for_orientation(X):
    X = X.tocsr()
    deg = np.diff(X.indptr).astype(float)
    S = (X @ X.T).toarray()
    total = 0.0
    pairs = 0
    for i in range(X.shape[0]):
        for j in range(i + 1, X.shape[0]):
            if S[i, j] > 0 and deg[i] > deg[j] > 0:
                total += (S[i, j] / deg[j]) * 100.0
                pairs += 1
    return total, pairs


def _reference_nodf(X):
    X = X.tocsr()
    if X.shape[0] < 2 or X.shape[1] < 2:
        return np.nan
    row_total, row_pairs = _reference_nodf_for_orientation(X)
    col_total, col_pairs = _reference_nodf_for_orientation(X.T.tocsr())
    pairs = row_pairs + col_pairs
    if pairs == 0:
        return np.nan
    return float((row_total + col_total) / pairs)


def test_association_obj_fe_and_ee(raw_ingredients):
    filtered = raw_ingredients.filtered_samples([s <= "S050" for s in raw_ingredients.samples])

    fe = association_obj(
        raw_ingredients,
        filtered,
        min_conditional_probability=0.0,
        null_model="FE",
        nm_n_reps=1,
        nm_seed=3,
    )
    assert not fe.empty
    assert {"taxon", "a", "b", "c", "d", "p_T_given_X", "jaccard", "q_bh"}.issubset(fe.columns)
    assert not any(c.startswith("jaccard_null_mean_FE") for c in fe.columns)

    ee = association_obj(
        raw_ingredients,
        filtered,
        min_conditional_probability=0.0,
        null_model="EE",
        nm_n_reps=1,
        nm_seed=3,
    )
    assert "jaccard_null_mean_EE" in ee.columns
    assert "jaccard_p_EE" in ee.columns
    assert {"null_seed", "null_seed_source", "null_model"}.issubset(ee.columns)
    assert np.isfinite(ee["jaccard"].to_numpy()).all()


def test_cooccurrence_obj_fe_and_ee(raw_ingredients):
    taxa_universe = list(raw_ingredients.taxa)

    edge_arrays, nodes_df = cooccurrence_obj(
        raw_ingredients,
        taxa_universe,
        large=True,
        min_conditional_probability=0.0,
        null_model="FE",
        nm_n_reps=1,
        nm_seed=4,
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
        min_conditional_probability=0.0,
        null_model="EE",
        nm_n_reps=1,
        nm_seed=4,
    )
    assert edge_arrays_ee is not None
    assert "jaccard_null_mean_EE" in edge_arrays_ee.cols
    assert "jaccard_p_EE" in edge_arrays_ee.cols
    assert edge_arrays_ee.meta["null_seed"] == 4
    assert edge_arrays_ee.meta["null_seed_source"] == "user"
    assert edge_arrays_ee.meta["null_model"] == "EE"


def test_focal_rhs_cooccurrence_edges_and_metrics(raw_ingredients):
    taxa_universe = list(raw_ingredients.taxa)
    taxa_by_name = {taxon: i for i, taxon in enumerate(taxa_universe)}
    taxa_arr = np.asarray(taxa_universe, dtype=object)

    focal_taxon = next(t for t in taxa_universe if t.endswith("g__Rhizo; s__rhizo_000"))
    rhs_taxa = [t for t in taxa_universe if "g__Micro; s__micro_" in t]
    micro_000 = next(t for t in rhs_taxa if t.endswith("g__Micro; s__micro_000"))

    edge_arrays, _ = cooccurrence_obj(
        raw_ingredients,
        taxa_universe,
        focal_query_to_taxa={"s__rhizo_000": [focal_taxon]},
        rhs_query_to_taxa={"s__rhizo_000": rhs_taxa},
        large=True,
        min_conditional_probability=0.0,
        null_model="FE",
        nm_n_reps=1,
        nm_seed=4,
    )

    assert edge_arrays is not None
    assert edge_arrays.n_rows == 50
    assert {"focal_query", "focal_taxon"}.issubset(edge_arrays.cols)
    assert set(edge_arrays.cols["focal_query"]) == {"s__rhizo_000"}
    assert set(edge_arrays.cols["focal_taxon"]) == {focal_taxon}
    assert np.all(edge_arrays.cols["iA"] == taxa_by_name[focal_taxon])
    assert all("g__Micro; s__micro_" in taxon for taxon in taxa_arr[edge_arrays.cols["iB"]])

    focal_samples = set(raw_ingredients.presence_matrix[taxa_by_name[focal_taxon]].indices)
    for row_idx in range(min(5, edge_arrays.n_rows)):
        b_taxon = taxa_arr[edge_arrays.cols["iB"][row_idx]]
        b_samples = set(raw_ingredients.presence_matrix[taxa_by_name[b_taxon]].indices)
        inter = len(focal_samples & b_samples)
        union = len(focal_samples | b_samples)
        assert edge_arrays.cols["inter"][row_idx] == inter
        assert np.isclose(edge_arrays.cols["jaccard"][row_idx], inter / union)

    row_idx = np.flatnonzero(taxa_arr[edge_arrays.cols["iB"]] == micro_000)
    assert row_idx.size == 1
    row_idx = int(row_idx[0])

    assert edge_arrays.cols["inter"][row_idx] == 36
    assert np.isclose(edge_arrays.cols["jaccard"][row_idx], 36 / 75)
    assert edge_arrays.cols["log_p"][row_idx] <= 0
    assert np.isfinite(edge_arrays.cols["log_q_bh"][row_idx])


def test_focal_rhs_cooccurrence_can_resolve_to_no_surviving_edges(raw_ingredients):
    taxa_universe = list(raw_ingredients.taxa)
    focal_taxon = next(t for t in taxa_universe if t.endswith("g__Rhizo; s__rhizo_000"))
    rhs_taxa = [t for t in taxa_universe if "g__Micro; s__micro_" in t]

    edge_arrays, nodes_df = cooccurrence_obj(
        raw_ingredients,
        taxa_universe,
        focal_query_to_taxa={"s__rhizo_000": [focal_taxon]},
        rhs_query_to_taxa={"s__rhizo_000": rhs_taxa},
        large=True,
        min_conditional_probability=1.0,
        null_model="FE",
        nm_n_reps=1,
        nm_seed=4,
    )

    assert edge_arrays is None
    assert nodes_df is not None


def test_structure_obj_observed_only(raw_ingredients):
    out = structure_obj(raw_ingredients, compute_null=False)
    assert out["metric"].tolist() == ["c_score", "mean_jaccard", "nodf"]
    assert "obs" in out.columns


def test_blocked_structure_metrics_match_full_reference():
    X = sp.csr_matrix(
        np.array(
            [
                [1, 1, 1, 0, 0],
                [1, 1, 0, 0, 0],
                [0, 0, 1, 1, 0],
                [0, 0, 0, 0, 0],
                [1, 0, 0, 1, 1],
            ],
            dtype=np.int32,
        )
    )

    assert np.isclose(compute_c_score(X, chunk_rows=1), _reference_c_score(X))
    assert np.isclose(mean_jaccard_dot(X, chunk_rows=1), _reference_mean_jaccard(X))
    assert np.isclose(compute_nodf_streamed(X, chunk_rows=1), _reference_nodf(X))

    assert np.isclose(compute_c_score(X, chunk_rows=2), _reference_c_score(X))
    assert np.isclose(mean_jaccard_dot(X, chunk_rows=2), _reference_mean_jaccard(X))
    assert np.isclose(compute_nodf_streamed(X, chunk_rows=2), _reference_nodf(X))


def test_blocked_structure_metrics_degenerate_inputs():
    assert np.isnan(
        compute_c_score(sp.csr_matrix([[1, 0, 1]], dtype=np.int32), chunk_rows=1)
    )
    assert np.isnan(
        mean_jaccard_dot(
            sp.csr_matrix([[0, 0], [1, 0]], dtype=np.int32),
            chunk_rows=1,
        )
    )
    assert np.isnan(
        compute_nodf_streamed(
            sp.csr_matrix([[1], [0], [1]], dtype=np.int32),
            chunk_rows=1,
        )
    )
    assert np.isnan(
        compute_nodf_streamed(
            sp.csr_matrix((0, 3), dtype=np.int32),
            chunk_rows=1,
        )
    )


def test_structure_core_uses_blocked_metrics(raw_ingredients):
    out = _structure_core(raw_ingredients, compute_null=False, chunk_rows=7)
    X = raw_ingredients.presence_matrix
    expected = {
        "c_score": compute_c_score(X, chunk_rows=7),
        "mean_jaccard": mean_jaccard_dot(X, chunk_rows=7),
        "nodf": compute_nodf_streamed(X, chunk_rows=7),
    }

    assert out["metric"].tolist() == ["c_score", "mean_jaccard", "nodf"]
    for metric, value in expected.items():
        observed = out.loc[out["metric"] == metric, "obs"].item()
        assert np.isclose(observed, value)


def test_structure_core_all_null_models(raw_ingredients):
    for model in ("FE", "EF", "EE", "FF"):
        out = _structure_core(
            raw_ingredients,
            null_model=model,
            nm_n_reps=1,
            nm_seed=5,
            compute_null=True,
            nm_burn_in_steps=2,
            nm_steps_per_rep=2,
            nm_progress_every=1,
        )
        assert out["metric"].tolist() == ["c_score", "mean_jaccard", "nodf"]
        assert f"null_mean_{model}" in out.columns
        assert {"null_seed", "null_seed_source", "null_model"}.issubset(out.columns)


def test_biome_distribution(raw_ingredients):
    biomes, presence, n_dropped = raw_ingredients.biome_distribution()
    assert biomes == ["terrestrial", "aquatic"]
    assert presence.shape == (2, 300)
    assert n_dropped == 0


def test_association_metrics_have_expected_counts(raw_ingredients):
    filtered = raw_ingredients.filtered_samples([s <= "S050" for s in raw_ingredients.samples])
    focal_taxon = next(t for t in raw_ingredients.taxa if t.endswith("g__Rhizo; s__rhizo_000"))

    out = association_obj(
        raw_ingredients,
        filtered,
        min_conditional_probability=0.0,
        null_model="FE",
        nm_n_reps=1,
        nm_seed=3,
    )
    row = out.loc[out["taxon"].eq(focal_taxon)].iloc[0]

    assert row["a"] == 50
    assert row["b"] == 10
    assert row["c"] == 0
    assert row["d"] == 40
    assert row["N_T"] == 50
    assert row["N_notT"] == 50
    assert row["N_null"] == 100
    assert np.isclose(row["jaccard"], 50 / 60)
    assert np.isclose(row["phi"], 0.8164965809277261)
    assert 0 <= row["q_bh"] <= 1

    filtered_sample_set = set(filtered.samples)
    null_sample_set = set(raw_ingredients.samples)
    taxon_to_idx = {taxon: i for i, taxon in enumerate(raw_ingredients.taxa)}

    for taxon in out["taxon"].head(10):
        taxon_samples = {
            raw_ingredients.samples[i]
            for i in raw_ingredients.presence_matrix[taxon_to_idx[taxon]].indices
        }
        a = len(taxon_samples & filtered_sample_set)
        b = len(taxon_samples - filtered_sample_set)
        c = len(filtered_sample_set - taxon_samples)
        d = len(null_sample_set - filtered_sample_set - taxon_samples)
        expected_jaccard = a / (a + b + c)

        metric_row = out.loc[out["taxon"].eq(taxon)].iloc[0]
        assert metric_row["a"] == a
        assert metric_row["b"] == b
        assert metric_row["c"] == c
        assert metric_row["d"] == d
        assert np.isclose(metric_row["jaccard"], expected_jaccard)
