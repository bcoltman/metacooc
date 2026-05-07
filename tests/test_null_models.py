from __future__ import annotations

import numpy as np
import pytest

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
