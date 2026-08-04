from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

import metacooc.analysis as analysis_module
from metacooc.analysis import (
    COOCCURRENCE_BASE_COLUMNS,
    COOCCURRENCE_FISHER_COLUMNS,
    CooccurrenceArrays,
    _build_cooccurrence_summary,
    _reconstruct_edge_frame,
    _write_full_edges_parquet_from_arrays,
    association,
    association_obj,
    bh_logq_from_logp,
    cooccurrence_obj,
    export_cooccurrence_outputs,
)
from metacooc.output import COMPACT_NULL_METADATA_COLUMNS
from metacooc.structure import (
    _structure_core,
    compute_c_score,
    compute_nodf_streamed,
    mean_jaccard_dot,
    structure,
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


def _synthetic_edge_arrays(*, compute_fisher=False, include_null=False):
    shared = np.array([5, 5, 3, 3], dtype=np.int32)
    source = np.array([0, 0, 1, 0], dtype=np.int32)
    target = np.array([1, 2, 3, 3], dtype=np.int32)
    totals = np.array([10, 8, 8, 5], dtype=np.int32)
    union = totals[source] + totals[target] - shared

    cols = {
        "source_taxon_index": source,
        "target_taxon_index": target,
        "shared_sample_count": shared,
        "jaccard_taxon_pair": (shared / union).astype(np.float32),
        "chi2_log_p_value": np.log(np.array([0.001, 0.8, 0.9, 0.7])),
        "chi2_log_q_value_bh": np.log(np.array([0.001, 0.8, 0.9, 0.7])),
    }
    if include_null:
        cols.update(
            {
                "jaccard_null_mean": np.full(4, 0.2, dtype=np.float32),
                "jaccard_null_sd": np.full(4, 0.1, dtype=np.float32),
                "jaccard_null_ses": np.array([0.2, 0.3, 0.4, 0.5]),
                "jaccard_null_p_empirical": np.array([0.1, 0.1, 0.01, np.nan]),
                "jaccard_null_log_q_value_bh": np.log(
                    np.array([0.1, 0.1, 0.01, np.nan])
                ),
            }
        )

    return CooccurrenceArrays(
        cols=cols,
        meta={
            "totals": totals,
            "N_total": 20,
            "compute_fisher": compute_fisher,
        },
    )


def _assert_compact_null_metadata(
    result,
    *,
    model,
    replicates,
    failed,
    seed,
):
    assert result.columns.tolist()[-4:] == COMPACT_NULL_METADATA_COLUMNS
    expected = {
        "null_model": model,
        "null_replicates": replicates,
        "null_replicates_failed": failed,
        "null_seed": seed,
    }
    for column, value in expected.items():
        if value is None:
            assert result[column].isna().all()
        else:
            assert set(result[column]) == {value}


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


def test_bh_logq_caps_adjusted_probabilities_at_one():
    log_q = bh_logq_from_logp(
        np.log(np.array([0.8, 0.9, np.nan])),
        m_total=10,
    )

    np.testing.assert_allclose(np.exp(log_q[:2]), [1.0, 1.0])
    assert np.isnan(log_q[2])


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
    assert {
        "taxon",
        "taxon_in_cohort_count",
        "taxon_in_background_not_cohort_count",
        "cohort_without_taxon_count",
        "neither_taxon_nor_cohort_count",
        "p_cohort_given_taxon",
        "p_taxon_given_cohort",
        "jaccard_taxon_cohort",
        "rr_cohort_given_taxon_vs_without_taxon",
        "rr_taxon_given_cohort_vs_without_cohort",
        "ln_rr_cohort_given_taxon_vs_without_taxon",
        "ln_rr_taxon_given_cohort_vs_without_cohort",
        "chi2_q_value_bh",
    }.issubset(fe.columns)
    assert {"RR_A_to_B", "RR_B_to_A", "logRR_A_to_B", "logRR_B_to_A"}.isdisjoint(fe.columns)

    row = fe.iloc[0]
    a = row["taxon_in_cohort_count"]
    b = row["taxon_in_background_not_cohort_count"]
    c = row["cohort_without_taxon_count"]
    d = row["neither_taxon_nor_cohort_count"]
    expected_rr_cohort = ((a + 0.5) / (a + b + 0.5)) / ((c + 0.5) / (c + d + 0.5))
    expected_rr_taxon = ((a + 0.5) / (a + c + 0.5)) / ((b + 0.5) / (b + d + 0.5))
    assert np.isclose(row["rr_cohort_given_taxon_vs_without_taxon"], expected_rr_cohort)
    assert np.isclose(row["rr_taxon_given_cohort_vs_without_cohort"], expected_rr_taxon)
    assert np.isclose(row["ln_rr_cohort_given_taxon_vs_without_taxon"], np.log(expected_rr_cohort))
    assert np.isclose(row["ln_rr_taxon_given_cohort_vs_without_cohort"], np.log(expected_rr_taxon))
    assert "jaccard_null_mean" not in fe.columns

    ee = association_obj(
        raw_ingredients,
        filtered,
        min_conditional_probability=0.0,
        null_model="EE",
        nm_n_reps=1,
        nm_seed=3,
    )
    assert "jaccard_null_mean" in ee.columns
    assert "jaccard_null_p_empirical" in ee.columns
    assert {"null_seed", "null_seed_source", "null_model"}.isdisjoint(ee.columns)
    assert ee.attrs["null_metadata"] == {
        "null_model": "EE",
        "null_replicates_requested": 1,
        "null_replicates_completed": 1,
        "null_replicates_ok": 1,
        "null_replicates_error": 0,
        "null_seed": 3,
        "null_seed_source": "user",
    }
    assert np.isfinite(ee["jaccard_taxon_cohort"].to_numpy()).all()


def test_association_obj_compute_fisher_uses_readable_columns(raw_ingredients):
    filtered = raw_ingredients.filtered_samples([s <= "S050" for s in raw_ingredients.samples])

    out = association_obj(
        raw_ingredients,
        filtered,
        min_conditional_probability=0.0,
        null_model="FE",
        nm_n_reps=1,
        nm_seed=3,
        compute_fisher=True,
    )

    assert {
        "fisher_odds_ratio",
        "fisher_p_value",
        "fisher_log_p_value",
    }.issubset(out.columns)
    assert {"fisher_odds", "fisher_p", "log_fisher_p"}.isdisjoint(out.columns)
    assert np.isfinite(out["fisher_p_value"].to_numpy()).all()


def test_association_wrapper_embeds_compact_null_metadata(tmp_path, raw_ingredients):
    filtered = raw_ingredients.filtered_samples([s <= "S050" for s in raw_ingredients.samples])

    association(
        raw_ingredients,
        filtered,
        output_dir=str(tmp_path),
        tag="direct_",
        min_conditional_probability=0.0,
        null_model="EE",
        nm_n_reps=1,
        nm_seed=3,
    )

    association_df = pd.read_csv(tmp_path / "direct_association.tsv", sep="\t")
    summary_df = pd.read_csv(
        tmp_path / "direct_association_summary.tsv",
        sep="\t",
    )
    assert "jaccard_null_mean" in association_df.columns
    assert summary_df.columns.tolist() == [
        "taxon",
        "phi_coefficient",
        "p_cohort_given_taxon",
        "p_taxon_given_cohort",
        "taxon_in_cohort_count",
        "cohort_sample_count",
        "taxon_in_background_not_cohort_count",
        "background_not_cohort_sample_count",
        "chi2_q_value_bh",
        "jaccard_null_ses",
        "jaccard_null_p_empirical",
        *COMPACT_NULL_METADATA_COLUMNS,
    ]
    assert len(summary_df) == len(association_df)
    assert summary_df["chi2_q_value_bh"].is_monotonic_increasing
    _assert_compact_null_metadata(
        association_df,
        model="EE",
        replicates=1,
        failed=0,
        seed=3,
    )
    _assert_compact_null_metadata(
        summary_df,
        model="EE",
        replicates=1,
        failed=0,
        seed=3,
    )
    assert not (tmp_path / "direct_association_metadata.tsv").exists()


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
    assert "jaccard_taxon_pair" in edge_arrays.cols

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
    assert "jaccard_null_mean" in edge_arrays_ee.cols
    assert "jaccard_null_p_empirical" in edge_arrays_ee.cols
    assert "jaccard_null_log_q_value_bh" in edge_arrays_ee.cols
    assert edge_arrays_ee.meta["null_metadata"] == {
        "null_model": "EE",
        "null_replicates_requested": 1,
        "null_replicates_completed": 1,
        "null_replicates_ok": 1,
        "null_replicates_error": 0,
        "null_seed": 4,
        "null_seed_source": "user",
    }
    assert edge_arrays_ee.meta["null_seed"] == 4
    assert edge_arrays_ee.meta["null_seed_source"] == "user"
    assert edge_arrays_ee.meta["null_model"] == "EE"


def test_cooccurrence_large_export_embeds_compact_null_metadata(
    tmp_path,
    raw_ingredients,
):
    taxa_universe = list(raw_ingredients.taxa)
    edge_arrays, nodes_df = cooccurrence_obj(
        raw_ingredients,
        taxa_universe,
        large=True,
        min_conditional_probability=0.0,
        null_model="EE",
        nm_n_reps=1,
        nm_seed=4,
    )
    assert edge_arrays is not None
    assert edge_arrays.n_rows > 1

    export_cooccurrence_outputs(
        edge_arrays=edge_arrays,
        nodes_df=nodes_df,
        taxa_universe=taxa_universe,
        output_dir=str(tmp_path),
        edges_base="large_edges",
        nodes_base="large_nodes",
        null_model="EE",
        summary_n=1,
    )

    assert (tmp_path / "large_nodes.tsv").exists()
    assert (tmp_path / "large_edges.parquet").exists()
    assert (tmp_path / "large_edges_taxa.parquet").exists()
    assert (tmp_path / "large_edges_summary.tsv").exists()
    assert not (tmp_path / "large_edges_metadata.tsv").exists()

    full_edges = pd.read_parquet(tmp_path / "large_edges.parquet")
    summary_edges = pd.read_csv(tmp_path / "large_edges_summary.tsv", sep="\t")
    assert {
        "source_taxon_index",
        "target_taxon_index",
        "shared_sample_count",
        "jaccard_taxon_pair",
        "jaccard_null_mean",
    }.issubset(full_edges.columns)
    expected_full_edges = _reconstruct_edge_frame(
        edge_arrays,
        taxa_universe,
        idx=slice(None),
        compute_fisher=False,
    )
    expected_full_edges.drop(
        columns=["source_taxon", "target_taxon"],
        inplace=True,
    )
    expected_full_edges.insert(
        0,
        "source_taxon_index",
        edge_arrays.cols["source_taxon_index"],
    )
    expected_full_edges.insert(
        1,
        "target_taxon_index",
        edge_arrays.cols["target_taxon_index"],
    )
    pd.testing.assert_frame_equal(
        full_edges.drop(columns=COMPACT_NULL_METADATA_COLUMNS),
        expected_full_edges,
        check_dtype=False,
    )
    _assert_compact_null_metadata(
        full_edges,
        model="EE",
        replicates=1,
        failed=0,
        seed=4,
    )
    assert summary_edges.columns.tolist() == [
        "source_taxon",
        "target_taxon",
        "phi_coefficient",
        "p_target_given_source",
        "p_source_given_target",
        "shared_sample_count",
        "source_taxon_sample_count",
        "target_taxon_sample_count",
        "chi2_q_value_bh",
        "jaccard_null_ses",
        "jaccard_null_p_empirical",
        "jaccard_null_q_value_bh",
        *COMPACT_NULL_METADATA_COLUMNS,
    ]
    assert len(summary_edges) == 1
    _assert_compact_null_metadata(
        summary_edges,
        model="EE",
        replicates=1,
        failed=0,
        seed=4,
    )


def test_cooccurrence_summary_uses_empirical_q_and_exact_tie_breaking():
    taxa_universe = ["z_source", "b_target", "a_target", "c_target"]
    edge_arrays = _synthetic_edge_arrays(include_null=True)

    capped = _build_cooccurrence_summary(
        edge_arrays,
        taxa_universe,
        summary_n=2,
    )
    assert list(zip(capped["source_taxon"], capped["target_taxon"])) == [
        ("b_target", "c_target"),
        ("z_source", "a_target"),
    ]
    np.testing.assert_allclose(
        capped["jaccard_null_q_value_bh"],
        [0.01, 0.1],
    )

    uncapped = _build_cooccurrence_summary(
        edge_arrays,
        taxa_universe,
        summary_n=4,
    )
    assert list(zip(uncapped["source_taxon"], uncapped["target_taxon"])) == [
        ("b_target", "c_target"),
        ("z_source", "a_target"),
        ("z_source", "b_target"),
        ("z_source", "c_target"),
    ]
    assert np.isnan(uncapped["jaccard_null_q_value_bh"].iloc[-1])

    analytical = _build_cooccurrence_summary(
        _synthetic_edge_arrays(),
        taxa_universe,
        summary_n=1,
    )
    assert list(zip(analytical["source_taxon"], analytical["target_taxon"])) == [
        ("z_source", "b_target")
    ]


def test_large_parquet_batches_store_all_fisher_metrics_without_summary_recompute(
    tmp_path,
    monkeypatch,
    capsys,
):
    taxa_universe = ["z_source", "b_target", "a_target", "c_target"]
    edge_arrays = _synthetic_edge_arrays(compute_fisher=True)
    edge_arrays.cols["focal_query"] = np.array(["query"] * 4, dtype=object)
    edge_arrays.cols["focal_taxon"] = np.array(["z_source"] * 4, dtype=object)
    fisher_calls = 0
    fisher_exact = analysis_module._fisher_exact

    def counting_fisher_exact(*args, **kwargs):
        nonlocal fisher_calls
        fisher_calls += 1
        return fisher_exact(*args, **kwargs)

    monkeypatch.setattr(analysis_module, "_fisher_exact", counting_fisher_exact)
    output_path = str(tmp_path / "batched_edges")
    _write_full_edges_parquet_from_arrays(
        edge_arrays,
        taxa_universe,
        output_path,
        batch_size=2,
    )

    full_edges = pd.read_parquet(f"{output_path}.parquet")
    expected_columns = [
        "focal_query",
        "focal_taxon",
        "source_taxon_index",
        "target_taxon_index",
        *COOCCURRENCE_BASE_COLUMNS[2:],
        *COOCCURRENCE_FISHER_COLUMNS,
        *COMPACT_NULL_METADATA_COLUMNS,
    ]
    assert full_edges.columns.tolist() == expected_columns
    assert len(full_edges) == edge_arrays.n_rows
    assert full_edges["focal_taxon"].eq("z_source").all()
    assert full_edges["background_sample_count"].eq(20).all()
    assert full_edges[COOCCURRENCE_FISHER_COLUMNS].notna().all().all()
    _assert_compact_null_metadata(
        full_edges,
        model="FE",
        replicates=None,
        failed=None,
        seed=None,
    )
    assert fisher_calls == edge_arrays.n_rows

    summary = _build_cooccurrence_summary(
        edge_arrays,
        taxa_universe,
        summary_n=1,
    )
    assert fisher_calls == edge_arrays.n_rows
    assert set(COOCCURRENCE_FISHER_COLUMNS).isdisjoint(summary.columns)

    output = capsys.readouterr().out
    assert "all 4 edges" in output
    assert "batch 1/2" in output
    assert "batch 2/2" in output


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
    assert np.all(edge_arrays.cols["source_taxon_index"] == taxa_by_name[focal_taxon])
    assert all(
        "g__Micro; s__micro_" in taxon
        for taxon in taxa_arr[edge_arrays.cols["target_taxon_index"]]
    )

    focal_samples = set(raw_ingredients.presence_matrix[taxa_by_name[focal_taxon]].indices)
    for row_idx in range(min(5, edge_arrays.n_rows)):
        target_taxon = taxa_arr[edge_arrays.cols["target_taxon_index"][row_idx]]
        target_samples = set(raw_ingredients.presence_matrix[taxa_by_name[target_taxon]].indices)
        shared_sample_count = len(focal_samples & target_samples)
        union = len(focal_samples | target_samples)
        assert edge_arrays.cols["shared_sample_count"][row_idx] == shared_sample_count
        assert np.isclose(
            edge_arrays.cols["jaccard_taxon_pair"][row_idx],
            shared_sample_count / union,
        )

    row_idx = np.flatnonzero(taxa_arr[edge_arrays.cols["target_taxon_index"]] == micro_000)
    assert row_idx.size == 1
    row_idx = int(row_idx[0])

    assert edge_arrays.cols["shared_sample_count"][row_idx] == 36
    assert np.isclose(edge_arrays.cols["jaccard_taxon_pair"][row_idx], 36 / 75)
    assert edge_arrays.cols["chi2_log_p_value"][row_idx] <= 0
    assert np.isfinite(edge_arrays.cols["chi2_log_q_value_bh"][row_idx])


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
    assert {"observed_value", "observed_error"}.issubset(out.columns)
    assert "obs" not in out.columns
    assert "null_mean" not in out.columns


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
        observed = out.loc[out["metric"] == metric, "observed_value"].item()
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
        assert {
            "null_mean",
            "null_sd",
            "null_standardized_effect_size",
            "null_p_empirical",
        }.issubset(out.columns)
        assert {"null_seed", "null_seed_source", "null_model"}.isdisjoint(out.columns)
        assert out.attrs["null_metadata"] == {
            "null_model": model,
            "null_replicates_requested": 1,
            "null_replicates_completed": 1,
            "null_replicates_ok": 1,
            "null_replicates_error": 0,
            "null_seed": 5,
            "null_seed_source": "user",
        }


def test_structure_wrapper_embeds_compact_null_metadata(tmp_path, raw_ingredients):
    structure(
        raw_ingredients,
        output_dir=str(tmp_path),
        tag="direct_",
        null_model="EE",
        nm_n_reps=1,
        nm_seed=5,
    )

    structure_df = pd.read_csv(tmp_path / "direct_structure.tsv", sep="\t")

    assert "null_mean" in structure_df.columns
    _assert_compact_null_metadata(
        structure_df,
        model="EE",
        replicates=1,
        failed=0,
        seed=5,
    )
    assert not (tmp_path / "direct_structure_metadata.tsv").exists()


def test_structure_wrapper_without_null_leaves_compact_metadata_empty(
    tmp_path,
    raw_ingredients,
):
    structure(
        raw_ingredients,
        output_dir=str(tmp_path),
        tag="observed_",
        compute_null=False,
    )

    structure_df = pd.read_csv(tmp_path / "observed_structure.tsv", sep="\t")
    _assert_compact_null_metadata(
        structure_df,
        model=None,
        replicates=None,
        failed=None,
        seed=None,
    )


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

    assert row["taxon_in_cohort_count"] == 50
    assert row["taxon_in_background_not_cohort_count"] == 10
    assert row["cohort_without_taxon_count"] == 0
    assert row["neither_taxon_nor_cohort_count"] == 40
    assert row["cohort_sample_count"] == 50
    assert row["background_not_cohort_sample_count"] == 50
    assert row["background_sample_count"] == 100
    assert np.isclose(row["jaccard_taxon_cohort"], 50 / 60)
    assert np.isclose(row["phi_coefficient"], 0.8164965809277261)
    assert 0 <= row["chi2_q_value_bh"] <= 1

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
        assert metric_row["taxon_in_cohort_count"] == a
        assert metric_row["taxon_in_background_not_cohort_count"] == b
        assert metric_row["cohort_without_taxon_count"] == c
        assert metric_row["neither_taxon_nor_cohort_count"] == d
        assert np.isclose(metric_row["jaccard_taxon_cohort"], expected_jaccard)
