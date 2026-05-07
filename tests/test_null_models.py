from __future__ import annotations

import itertools
from collections import Counter

import numpy as np
import pytest
import scipy.sparse as sp

from metacooc.null_models import make_null_sampler

pytestmark = pytest.mark.nullmodels


def _sample_one(X, model, seed=11):
    sampler = make_null_sampler(
        X,
        model,
        random_state=seed,
        burn_in_steps=2,
        steps_per_rep=2,
        sort_indices=True,
    )
    return next(iter(sampler.sample(1, seed=seed)))


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
        sampler = make_null_sampler(X, model, random_state=123, sort_indices=False)
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
    sampler = make_null_sampler(X, "EE", random_state=5, sort_indices=True)
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
