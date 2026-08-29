from __future__ import annotations

import itertools
from collections import Counter

import numpy as np
import pytest
import scipy.sparse as sp

import metacooc.null_models as null_models
from metacooc.null_models import (
    _seed_seq_spawn,
    make_null_sampler,
    parallel_null_reduce_vector,
)

pytestmark = pytest.mark.nullmodels


def _stat_support_checksum(X):
    coo = X.tocoo()
    linear = coo.row.astype(np.float64) * X.shape[1] + coo.col.astype(np.float64)
    return np.array([linear.sum()], dtype=float)


def _sample_one(X, model, seed=11):
    sampler = make_null_sampler(
        X,
        model,
        seed=seed,
        burn_in_steps=2,
        steps_per_rep=2,
        sort_indices=True,
    )
    return next(iter(sampler.sample(1, seed=seed)))


def _reference_cooccurrence_jaccard(X, subset_idx, iA, iB):
    X_sub = X[subset_idx, :].astype(np.int32, copy=False).tocsr()
    totals = np.asarray(X_sub.sum(axis=1)).ravel().astype(np.float64, copy=False)
    co = (X_sub @ X_sub.T).tocsr()
    inter = co[iA, iB].A1.astype(np.float64, copy=False)
    union = totals[iA] + totals[iB] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.divide(
            inter,
            union,
            out=np.zeros_like(inter, dtype=float),
            where=union > 0,
        )


def _adaptive_cooccurrence_jaccard(monkeypatch, X, subset_idx, iA, iB):
    monkeypatch.setattr(null_models, "_G_subset_idx", subset_idx)
    monkeypatch.setattr(null_models, "_G_iA", iA)
    monkeypatch.setattr(null_models, "_G_iB", iB)
    return null_models.stat_fn_cooccurrence_jaccard(X)


def test_fe_preserves_row_sums(raw_ingredients):
    X = raw_ingredients.presence_matrix
    Y = _sample_one(X, "FE")
    assert Y.getformat() == "csr"
    assert Y.shape == X.shape
    assert np.array_equal(Y.getnnz(axis=1), X.getnnz(axis=1))
    assert set(Y.data.tolist()) == {1}


def test_ef_preserves_column_sums(raw_ingredients):
    X = raw_ingredients.presence_matrix
    Y = _sample_one(X, "EF")
    assert Y.getformat() == "csr"
    assert Y.shape == X.shape
    assert np.array_equal(Y.getnnz(axis=0), X.getnnz(axis=0))
    assert set(Y.data.tolist()) == {1}


def test_ee_preserves_shape_and_fill(raw_ingredients):
    X = raw_ingredients.presence_matrix
    Y = _sample_one(X, "EE")
    assert Y.getformat() == "csr"
    assert Y.shape == X.shape
    assert Y.nnz == X.nnz
    assert set(Y.data.tolist()) == {1}


def test_ff_preserves_row_and_column_sums(raw_ingredients):
    X = raw_ingredients.presence_matrix
    Y = _sample_one(X, "FF")
    assert Y.shape == X.shape
    assert np.array_equal(Y.getnnz(axis=1), X.getnnz(axis=1))
    assert np.array_equal(Y.getnnz(axis=0), X.getnnz(axis=0))


def test_direct_samplers_are_seed_reproducible(raw_ingredients):
    X = raw_ingredients.presence_matrix
    for model in ("FE", "EF", "EE"):
        A = _sample_one(X, model, seed=99)
        B = _sample_one(X, model, seed=99)
        assert (A != B).nnz == 0


def test_direct_samplers_are_sorted_csr_without_duplicate_columns(raw_ingredients):
    X = raw_ingredients.presence_matrix

    for model in ("FE", "EF", "EE"):
        Y = _sample_one(X, model, seed=101)
        assert Y.getformat() == "csr"
        assert Y.has_sorted_indices

        for row in range(Y.shape[0]):
            start, end = Y.indptr[row], Y.indptr[row + 1]
            row_indices = Y.indices[start:end]
            assert np.unique(row_indices).size == row_indices.size


def test_direct_samplers_have_no_duplicate_support_when_unsorted(raw_ingredients):
    X = raw_ingredients.presence_matrix

    for model in ("FE", "EF", "EE"):
        sampler = make_null_sampler(X, model, seed=123, sort_indices=False)
        Y = next(iter(sampler.sample(1, seed=123)))
        coo = Y.tocoo()
        support = np.column_stack((coo.row, coo.col))
        assert np.unique(support, axis=0).shape[0] == Y.nnz


def test_ee_large_sparse_shape_avoids_population_allocation():
    rows = np.array([0, 17, 999, 20_000, 39_999], dtype=np.int64)
    cols = np.array([2, 123, 4567, 29_999, 1], dtype=np.int64)
    X = sp.csr_matrix(
        (np.ones(rows.size, dtype=np.int8), (rows, cols)),
        shape=(40_000, 30_000),
    )

    Y = _sample_one(X, "EE", seed=321)

    assert X.shape[0] * X.shape[1] > 1_000_000_000
    assert Y.shape == X.shape
    assert Y.nnz == X.nnz


def test_ee_small_matrix_support_frequencies_are_uniform():
    X = sp.csr_matrix(
        (np.ones(2, dtype=np.int8), ([0, 1], [0, 2])),
        shape=(2, 3),
    )
    sampler = make_null_sampler(X, "EE", seed=5, sort_indices=True)
    n_reps = 12_000
    counts = Counter()

    for Y in sampler.sample(n_reps, seed=5):
        coo = Y.tocoo()
        linear = tuple(sorted((coo.row * X.shape[1] + coo.col).tolist()))
        counts[linear] += 1

    supports = list(itertools.combinations(range(X.shape[0] * X.shape[1]), X.nnz))
    expected = n_reps / len(supports)

    assert set(counts) == set(supports)
    for support in supports:
        assert abs(counts[support] - expected) < expected * 0.20


def test_cooccurrence_jaccard_matches_full_matrix_reference(monkeypatch):
    X = sp.csr_matrix(
        (
            np.ones(14, dtype=np.int8),
            (
                [0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4],
                [0, 1, 4, 1, 2, 0, 2, 5, 0, 3, 5, 1, 3, 4],
            ),
        ),
        shape=(6, 6),
    )
    subset_idx = np.array([3, 0, 4, 1, 5], dtype=np.int64)
    iA = np.array([0, 0, 1, 2, 3, 3, 4], dtype=np.int64)
    iB = np.array([1, 3, 2, 0, 0, 3, 4], dtype=np.int64)

    observed = _adaptive_cooccurrence_jaccard(monkeypatch, X, subset_idx, iA, iB)
    expected = _reference_cooccurrence_jaccard(X, subset_idx, iA, iB)

    assert np.allclose(observed, expected)


def test_cooccurrence_jaccard_matches_reference_when_batched(monkeypatch):
    X = sp.csr_matrix(
        (
            np.ones(18, dtype=np.int8),
            (
                [0, 0, 1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 7],
                [0, 2, 1, 2, 5, 0, 5, 1, 3, 4, 0, 4, 2, 3, 5, 1, 5, 4],
            ),
        ),
        shape=(8, 6),
    )
    subset_idx = np.array([7, 0, 5, 2, 6, 3, 1], dtype=np.int64)
    iA = np.array([0, 1, 1, 2, 3, 4, 5, 5, 6], dtype=np.int64)
    iB = np.array([1, 2, 5, 0, 4, 6, 1, 3, 0], dtype=np.int64)
    monkeypatch.setattr(null_models, "_COOCCURRENCE_JACCARD_MAX_SOURCE_ROWS", 1)
    monkeypatch.setattr(null_models, "_COOCCURRENCE_JACCARD_MAX_EDGE_BLOCK", 2)

    observed = _adaptive_cooccurrence_jaccard(monkeypatch, X, subset_idx, iA, iB)
    expected = _reference_cooccurrence_jaccard(X, subset_idx, iA, iB)

    assert np.allclose(observed, expected)


def test_generated_null_seed_reproduces_when_reused():
    X = sp.csr_matrix(
        (np.ones(5, dtype=np.int8), ([0, 1, 2, 2, 3], [1, 2, 0, 3, 1])),
        shape=(4, 5),
    )
    obs = np.array([0.0], dtype=float)

    first = parallel_null_reduce_vector(
        X=X,
        model="EE",
        n_reps=8,
        obs=obs,
        stat_fn=_stat_support_checksum,
        seed=None,
        n_workers=2,
        progress_every=1,
    )
    replay = parallel_null_reduce_vector(
        X=X,
        model="EE",
        n_reps=8,
        obs=obs,
        stat_fn=_stat_support_checksum,
        seed=first["null_seed"],
        n_workers=2,
        progress_every=1,
    )

    assert first["null_seed_source"] == "generated"
    assert replay["null_seed_source"] == "user"
    assert replay["null_seed"] == first["null_seed"]
    assert np.array_equal(replay["mean"], first["mean"])
    assert np.array_equal(replay["sd"], first["sd"])


def test_spawned_worker_seeds_are_distinct():
    seeds = _seed_seq_spawn(12345, 8)
    assert len(seeds) == 8
    assert len(set(seeds)) == len(seeds)


def test_ff_sampler_streams_mutable_snapshots(raw_ingredients):
    sampler = make_null_sampler(
        raw_ingredients.presence_matrix,
        "FF",
        seed=123,
        burn_in_steps=0,
        steps_per_rep=1,
        sort_indices=True,
    )
    stream = iter(sampler.sample(2))
    first = next(stream)
    second = next(stream)

    assert first is second


def test_direct_sampler_snapshots_do_not_alias(raw_ingredients):
    sampler = make_null_sampler(
        raw_ingredients.presence_matrix,
        "FE",
        seed=123,
        sort_indices=True,
    )
    stream = iter(sampler.sample(2))
    first = next(stream)
    first_copy = first.copy()
    second = next(stream)

    assert first is not second
    assert (first != first_copy).nnz == 0
