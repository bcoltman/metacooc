#!/usr/bin/env python3

import os
import json
import tarfile
import warnings
import re
from collections import defaultdict
from datetime import datetime, timezone

from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

from metacooc._data_config import (
    DataReleaseError,
    IngredientsFormatError,
    LATEST_DATA_RELEASE,
    get_file_info,
    missing_data_release_message,
)
from metacooc.utils import (
    _RANK_PREFIXES, 
    _PREFIX_TO_RANK, 
    _parse_tokens, 
    _token_rank, 
    _terminal_rank_prefix, 
    _deepest_rank_token
)

INGREDIENTS_FORMAT_VERSION = 1
AGGREGATED_SUFFIX = " AGGREGATED"


def _utc_timestamp_seconds() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def presence_for_counts(obj):
    X = obj.presence_matrix if hasattr(obj, "presence_matrix") else obj
    return X.astype(np.int32, copy=False)


def _presence_csr(matrix: sp.spmatrix) -> sp.csr_matrix:
    out = matrix.tocsr(copy=True).astype(np.uint8, copy=False)
    out.data[:] = 1
    out.eliminate_zeros()
    out.sum_duplicates()
    out.sort_indices()
    return out


def _csr(matrix: sp.spmatrix) -> sp.csr_matrix:
    out = matrix.tocsr(copy=False)
    out.eliminate_zeros()
    out.sum_duplicates()
    out.sort_indices()
    return out


def presence_for_counts(obj):
    X = obj.presence_matrix if hasattr(obj, "presence_matrix") else obj
    return X.astype(np.int32, copy=False)


def _rank_thresholds_dict(thresholds, name: str) -> dict[str, float]:
    if thresholds is None:
        return {}

    items = thresholds.items() if isinstance(thresholds, dict) else thresholds
    out: dict[str, float] = {}
    for rank, value in items:
        rank = str(rank).strip().lower()
        if rank not in _RANK_PREFIXES:
            raise ValueError(
                f"Unknown rank '{rank}' in {name}. Expected one of: "
                f"{', '.join(_RANK_PREFIXES.keys())}"
            )
        out[rank] = float(value)
    return out


def _presence_threshold_active(
    min_coverage=None,
    min_coverage_by_rank=None,
    min_relative_abundance=None,
    min_relative_abundance_by_rank=None,
) -> bool:
    return any(
        value is not None
        for value in (min_coverage, min_relative_abundance)
    ) or bool(min_coverage_by_rank) or bool(min_relative_abundance_by_rank)


def _taxon_for_rank(taxon: str) -> str:
    if taxon.endswith(AGGREGATED_SUFFIX):
        return taxon[: -len(AGGREGATED_SUFFIX)]
    return taxon


def _terminal_rank_name(taxon: str) -> Optional[str]:
    pref = _terminal_rank_prefix(_parse_tokens(_taxon_for_rank(taxon)))
    if pref is None:
        return None
    return _PREFIX_TO_RANK.get(pref)


def _row_cutoffs(
    taxa: list[str],
    *,
    global_threshold,
    by_rank,
) -> np.ndarray:
    default = 0.0 if global_threshold is None else float(global_threshold)
    cutoffs = np.full(len(taxa), default, dtype=float)
    if not by_rank:
        return cutoffs

    for i, taxon in enumerate(taxa):
        rank = _terminal_rank_name(taxon)
        if rank in by_rank:
            cutoffs[i] = float(by_rank[rank])
    return cutoffs


def _pass_matrix_from_values(matrix: sp.spmatrix, row_cutoffs: np.ndarray) -> sp.csr_matrix:
    mat = matrix.tocsr(copy=True)
    mat.eliminate_zeros()
    if mat.nnz == 0:
        return mat.astype(np.uint8, copy=False)

    row_nnz = np.diff(mat.indptr)
    thresholds = np.repeat(row_cutoffs.astype(float, copy=False), row_nnz)
    keep = mat.data >= thresholds
    mat.data = keep.astype(np.uint8, copy=False)
    mat.eliminate_zeros()
    mat.sum_duplicates()
    mat.sort_indices()
    return mat.astype(np.uint8, copy=False)


def _ancestor_map_for_base_taxa(base_taxa: list[str]) -> dict[str, list[int]]:
    ancestors: dict[str, list[int]] = defaultdict(list)
    for local_idx, taxon in enumerate(base_taxa):
        tokens = _parse_tokens(taxon)
        for i, token in enumerate(tokens):
            if _token_rank(token) in _RANK_PREFIXES:
                ancestors["; ".join(tokens[: i + 1])].append(local_idx)
    return ancestors


def _aggregated_descendant_matrix(
    taxa: list[str],
    base_rows: np.ndarray,
    aggregated_rows: np.ndarray,
) -> sp.csr_matrix:
    if aggregated_rows.size == 0 or base_rows.size == 0:
        return sp.csr_matrix((aggregated_rows.size, base_rows.size), dtype=np.uint8)

    base_taxa = [taxa[i] for i in base_rows]
    ancestor_map = _ancestor_map_for_base_taxa(base_taxa)

    rows = []
    cols = []
    for agg_local_idx, row_idx in enumerate(aggregated_rows):
        ancestor = _taxon_for_rank(taxa[row_idx])
        descendants = ancestor_map.get(ancestor, [])
        rows.extend([agg_local_idx] * len(descendants))
        cols.extend(descendants)

    data = np.ones(len(rows), dtype=np.uint8)
    return sp.csr_matrix(
        (data, (rows, cols)),
        shape=(aggregated_rows.size, base_rows.size),
        dtype=np.uint8,
    )


def _relative_abundance_base_matrix(coverage_base: sp.csr_matrix) -> sp.csr_matrix:
    denom = np.asarray(coverage_base.sum(axis=0), dtype=float).ravel()
    rel = coverage_base.tocsr(copy=True).astype(float, copy=False)
    rel.eliminate_zeros()
    if rel.nnz == 0:
        return rel

    denom_for_data = denom[rel.indices]
    rel.data = np.divide(
        rel.data,
        denom_for_data,
        out=np.zeros_like(rel.data, dtype=float),
        where=denom_for_data > 0,
    )
    rel.eliminate_zeros()
    rel.sum_duplicates()
    rel.sort_indices()
    return rel


def _relative_abundance_pass_matrix(
    ingredients: "Ingredients",
    coverage: sp.csr_matrix,
    row_cutoffs: np.ndarray,
) -> sp.csr_matrix:
    taxa = list(ingredients.taxa)
    base_mask = np.fromiter(
        (not taxon.endswith(AGGREGATED_SUFFIX) for taxon in taxa),
        dtype=bool,
        count=len(taxa),
    )
    base_rows = np.flatnonzero(base_mask).astype(np.int64, copy=False)
    aggregated_rows = np.flatnonzero(~base_mask).astype(np.int64, copy=False)

    if base_rows.size == 0:
        return sp.csr_matrix(coverage.shape, dtype=np.uint8)

    rel_base = _relative_abundance_base_matrix(coverage[base_rows, :].tocsr())
    base_pass = _pass_matrix_from_values(rel_base, row_cutoffs[base_rows])

    row_parts = []
    col_parts = []
    data_parts = []

    if base_pass.nnz:
        coo = base_pass.tocoo()
        row_parts.append(base_rows[coo.row])
        col_parts.append(coo.col)
        data_parts.append(np.ones(coo.nnz, dtype=np.uint8))

    if aggregated_rows.size:
        descendant_map = _aggregated_descendant_matrix(taxa, base_rows, aggregated_rows)
        rel_agg = descendant_map @ rel_base
        rel_agg.eliminate_zeros()
        agg_pass = _pass_matrix_from_values(rel_agg, row_cutoffs[aggregated_rows])
        if agg_pass.nnz:
            coo = agg_pass.tocoo()
            row_parts.append(aggregated_rows[coo.row])
            col_parts.append(coo.col)
            data_parts.append(np.ones(coo.nnz, dtype=np.uint8))

    if not data_parts:
        return sp.csr_matrix(coverage.shape, dtype=np.uint8)

    rows = np.concatenate(row_parts)
    cols = np.concatenate(col_parts)
    data = np.concatenate(data_parts)
    out = sp.csr_matrix((data, (rows, cols)), shape=coverage.shape, dtype=np.uint8)
    out.eliminate_zeros()
    out.sum_duplicates()
    out.data[:] = 1
    out.sort_indices()
    return out


def _presence_threshold_provenance(
    *,
    min_coverage,
    min_coverage_by_rank,
    min_relative_abundance,
    min_relative_abundance_by_rank,
) -> dict[str, object]:
    provenance = {
        "comparison": ">=",
        "masked_coverage": "original values retained for passing cells; failed cells zeroed",
    }
    if min_coverage is not None:
        provenance["min_coverage"] = float(min_coverage)
    if min_coverage_by_rank:
        provenance["min_coverage_by_rank"] = {
            rank: float(value) for rank, value in sorted(min_coverage_by_rank.items())
        }
    if min_relative_abundance is not None:
        provenance["min_relative_abundance"] = float(min_relative_abundance)
    if min_relative_abundance_by_rank:
        provenance["min_relative_abundance_by_rank"] = {
            rank: float(value)
            for rank, value in sorted(min_relative_abundance_by_rank.items())
        }
    if min_relative_abundance is not None or min_relative_abundance_by_rank:
        provenance["relative_abundance"] = {
            "units": "fraction",
            "base_denominator": "non-AGGREGATED rows per sample",
            "aggregated_rows": "sum of descendant non-AGGREGATED relative abundances",
        }
    return provenance


def threshold_ingredients_presence(
    ingredients: "Ingredients",
    *,
    min_coverage=None,
    min_coverage_by_rank=None,
    min_relative_abundance=None,
    min_relative_abundance_by_rank=None,
) -> "Ingredients":
    """
    Derive binary presence from coverage/relative-abundance thresholds.

    The returned Ingredients object keeps the original coverage values for
    passing cells and zeroes coverage where the derived presence is absent.
    """
    min_coverage_by_rank = _rank_thresholds_dict(
        min_coverage_by_rank,
        "min_coverage_by_rank",
    )
    min_relative_abundance_by_rank = _rank_thresholds_dict(
        min_relative_abundance_by_rank,
        "min_relative_abundance_by_rank",
    )

    if not _presence_threshold_active(
        min_coverage=min_coverage,
        min_coverage_by_rank=min_coverage_by_rank,
        min_relative_abundance=min_relative_abundance,
        min_relative_abundance_by_rank=min_relative_abundance_by_rank,
    ):
        return ingredients

    if (
        min_coverage is not None or min_coverage_by_rank
    ) and (
        min_relative_abundance is not None or min_relative_abundance_by_rank
    ):
        warnings.warn(
            "Both coverage and relative-abundance presence thresholds were supplied; "
            "a taxon-sample cell must pass both to count as present.",
            UserWarning,
        )

    coverage = ingredients.coverage_matrix.tocsr(copy=True)
    coverage.eliminate_zeros()
    coverage.sum_duplicates()
    coverage.sort_indices()

    pass_matrix = coverage.copy()
    pass_matrix.data = np.ones(pass_matrix.nnz, dtype=np.uint8)
    pass_matrix = pass_matrix.astype(np.uint8, copy=False)

    if min_coverage is not None or min_coverage_by_rank:
        coverage_cutoffs = _row_cutoffs(
            list(ingredients.taxa),
            global_threshold=min_coverage,
            by_rank=min_coverage_by_rank,
        )
        coverage_pass = _pass_matrix_from_values(coverage, coverage_cutoffs)
        pass_matrix = pass_matrix.multiply(coverage_pass).tocsr()

    if min_relative_abundance is not None or min_relative_abundance_by_rank:
        rel_cutoffs = _row_cutoffs(
            list(ingredients.taxa),
            global_threshold=min_relative_abundance,
            by_rank=min_relative_abundance_by_rank,
        )
        rel_pass = _relative_abundance_pass_matrix(ingredients, coverage, rel_cutoffs)
        pass_matrix = pass_matrix.multiply(rel_pass).tocsr()

    presence = _presence_csr(pass_matrix)
    masked_coverage = coverage.multiply(presence).tocsr()
    masked_coverage.eliminate_zeros()
    masked_coverage.sum_duplicates()
    masked_coverage.sort_indices()

    out = ingredients.copy_shallow()
    out.presence_matrix = presence
    out.coverage_matrix = masked_coverage
    out.presence_thresholds = _presence_threshold_provenance(
        min_coverage=min_coverage,
        min_coverage_by_rank=min_coverage_by_rank,
        min_relative_abundance=min_relative_abundance,
        min_relative_abundance_by_rank=min_relative_abundance_by_rank,
    )
    return out


class Ingredients:
    """
    Container for metagenomic data with sparse matrices.
    
    Attributes:
        samples (List[str]): Sample identifiers.
        taxa (List[str]): Taxonomic identifiers.
        presence_matrix (sp.csr_matrix): Binary presence/absence.
        coverage_matrix (sp.csr_matrix): Coverage values.
        total_counts (np.ndarray): Cached per-taxon counts.
        data_release (Optional[str]): Source database release label, when known.
        sample_to_biome (Dict[str, Tuple[str, str]]): Mapping sample to (biome_level_1, biome_level_2).
        biomes_order (Dict[str, List[str]]): Unique biomes for each level.
        sample_biome_indices (Dict[str, np.ndarray]): Per-sample biome index for each level (-1 if missing).
    """
    def __init__(
        self,
        samples: List[str],
        taxa: List[str],
        presence_matrix: sp.csr_matrix,
        coverage_matrix: sp.csr_matrix,
        sample_to_biome: Dict[str, str] = None,
        data_release: Optional[str] = None,
        ):
        self.taxa = taxa
        self.samples = samples
        
        object.__setattr__(self, "_presence_matrix", _presence_csr(presence_matrix))
        object.__setattr__(self, "_coverage_matrix", _csr(coverage_matrix))
        object.__setattr__(self, "_presence_matrix_path", None)
        object.__setattr__(self, "_coverage_matrix_path", None)
        object.__setattr__(self, "_presence_matrix_shape", self._presence_matrix.shape)
        object.__setattr__(self, "_coverage_matrix_shape", self._coverage_matrix.shape)
        
        self.total_counts = self._compute_total_counts()
        
        # Load and preallocate biome mapping
        self.sample_to_biome = sample_to_biome or {}
        if self.sample_to_biome:
            self._allocate_biomes()
            
        self._rank_lookups = None
        self._terminal_rank_prefixes = None
        self.data_release = data_release
        self.presence_thresholds = None
    
    def __getstate__(self):
        state = {
            "samples": self.samples,
            "taxa": self.taxa,
            "_presence_matrix": self.presence_matrix,
            "_coverage_matrix": self.coverage_matrix,
            "total_counts": self.total_counts,
            "sample_to_biome": self.sample_to_biome,
            "data_release": self.data_release,
            "presence_thresholds": getattr(self, "presence_thresholds", None),
        }
        
        if hasattr(self, "biomes_order"):
            state["biomes_order"] = self.biomes_order
            state["sample_biome_indices"] = self.sample_biome_indices
            
        return state
    
    def __setstate__(self, state):
        self.samples = state["samples"]
        self.taxa = state["taxa"]
        
        object.__setattr__(self, "_presence_matrix", state.get("_presence_matrix"))
        object.__setattr__(self, "_coverage_matrix", state.get("_coverage_matrix"))
        object.__setattr__(self, "_presence_matrix_path", None)
        object.__setattr__(self, "_coverage_matrix_path", None)
        object.__setattr__(self, "_presence_matrix_shape", self._presence_matrix.shape)
        object.__setattr__(self, "_coverage_matrix_shape", self._coverage_matrix.shape)
        
        self._rank_lookups = None
        self._terminal_rank_prefixes = None
        
        # restore or compute total_counts
        # restore or compute total_counts
        tc = state.get("total_counts", None)
        self.total_counts = tc if tc is not None else self._compute_total_counts()
        
        # restore or default biome mapping
        self.sample_to_biome = state.get("sample_to_biome", {})
        
        # restore allocation if present
        if "biomes_order" in state and "sample_biome_indices" in state:
            self.biomes_order = state["biomes_order"]
            self.sample_biome_indices = state["sample_biome_indices"]
        elif self.sample_to_biome:
            self._allocate_biomes()
        
        self.data_release = state.get("data_release", None)
        self.presence_thresholds = state.get("presence_thresholds", None)
    
    def _invalidate_taxa_caches(self):
        self._rank_lookups = None
        self._terminal_rank_prefixes = None
    
    def _build_taxa_lookups(self):
        """
        Build per-rank exact-token lookups and terminal rank prefixes.
        _rank_lookups[rank][token] -> set(taxon_idx)
        """
        lookups = {rank: defaultdict(set) for rank in _RANK_PREFIXES.keys()}
        term_prefixes = []
        
        for i, taxon in enumerate(self.taxa):
            tokens = _parse_tokens(taxon)
            term_prefixes.append(_terminal_rank_prefix(tokens))
            for tok in tokens:
                r = _token_rank(tok)
                if r is None or r == "root":
                    continue
                lookups[r][tok].add(i)
                    
        self._rank_lookups = lookups
        self._terminal_rank_prefixes = term_prefixes
        
    def _ensure_taxa_lookups(self):
        if self._rank_lookups is None or self._terminal_rank_prefixes is None:
            self._build_taxa_lookups()
    
    def build_taxa_lookup(self) -> None:
        """Build (or rebuild) the per-rank taxa lookup caches."""
        self._build_taxa_lookups()
    
    @property
    def presence_matrix(self) -> sp.csr_matrix:
        if self._presence_matrix is None:
            self._presence_matrix = sp.load_npz(self._presence_matrix_path).tocsr()
            self._presence_matrix_shape = self._presence_matrix.shape
        return self._presence_matrix
    
    @presence_matrix.setter
    def presence_matrix(self, mat: sp.csr_matrix):
        self._presence_matrix = _presence_csr(mat)
        self._presence_matrix_path = None
        self._presence_matrix_shape = self._presence_matrix.shape
        self.total_counts = self._compute_total_counts()
    
    @property
    def coverage_matrix(self) -> sp.csr_matrix:
        if self._coverage_matrix is None:
            self._coverage_matrix = sp.load_npz(self._coverage_matrix_path).tocsr()
            self._coverage_matrix_shape = self._coverage_matrix.shape
        return self._coverage_matrix
    
    @coverage_matrix.setter
    def coverage_matrix(self, mat: sp.csr_matrix):
        self._coverage_matrix = _csr(mat)
        self._coverage_matrix_path = None
        self._coverage_matrix_shape = self._coverage_matrix.shape
    
    def _compute_total_counts(self) -> np.ndarray:
        """
        Per-taxon presence counts across samples.
        
        Because presence_matrix is already binary / sparse-presence, use getnnz(axis=1)
        instead of (matrix > 0).sum(axis=1), which allocates a large temporary sparse matrix.
        """
        return np.asarray(self.presence_matrix.getnnz(axis=1), dtype=np.int32).ravel()

        # return np.array((self._presence_matrix > 0).sum(axis=1)).flatten()
    
    def __repr__(self):
        presence_shape = getattr(self, "_presence_matrix_shape", None)
        if presence_shape is None:
            presence_shape = self.presence_matrix.shape
        coverage_shape = getattr(self, "_coverage_matrix_shape", None)
        if coverage_shape is None:
            coverage_shape = self.coverage_matrix.shape
        return (
            f"<Ingredients: {len(self.taxa)} taxa, "
            f"{len(self.samples)} samples, "
            f"presence: {presence_shape}, "
            f"coverage: {coverage_shape}>"
        )
    
    @staticmethod
    def _normalise_indexer(mask, size: int) -> np.ndarray:
        """
        Normalise a boolean mask or integer indexer into a 1D integer index array.
        """
        arr = np.asarray(mask)
        if arr.ndim != 1:
            arr = arr.ravel()
            
        if arr.dtype == bool:
            if arr.size != size:
                raise ValueError(
                    f"Boolean mask length {arr.size} does not match expected size {size}"
                )
            return np.flatnonzero(arr).astype(np.int64, copy=False)
            
        return arr.astype(np.int64, copy=False)


    def copy_shallow(self) -> "Ingredients":
        """
        Cheap structural copy for pipeline/filtering use.
        
        Shares sparse matrices and read-mostly metadata until the new object is sliced.
        This avoids duplicating GB-scale matrices during chained filtering.
        """
        new = Ingredients.__new__(Ingredients)
        
        new.samples = self.samples
        new.taxa = self.taxa
        
        object.__setattr__(new, "_presence_matrix", self._presence_matrix)
        object.__setattr__(new, "_coverage_matrix", self._coverage_matrix)
        object.__setattr__(new, "_presence_matrix_path", self._presence_matrix_path)
        object.__setattr__(new, "_coverage_matrix_path", self._coverage_matrix_path)
        object.__setattr__(new, "_presence_matrix_shape", self._presence_matrix_shape)
        object.__setattr__(new, "_coverage_matrix_shape", self._coverage_matrix_shape)
        
        new.total_counts = self.total_counts
        new.sample_to_biome = self.sample_to_biome
        new.data_release = self.data_release
        new.presence_thresholds = getattr(self, "presence_thresholds", None)
        
        new._rank_lookups = self._rank_lookups
        new._terminal_rank_prefixes = self._terminal_rank_prefixes
        
        if hasattr(self, "biomes_order"):
            new.biomes_order = self.biomes_order
        if hasattr(self, "sample_biome_indices"):
            new.sample_biome_indices = self.sample_biome_indices
            
        return new
    
    def copy(self, deep: bool = True):
        """
        deep=True preserves the old behaviour.
        deep=False is the fast path used internally by filtering helpers.
        """
        if not deep:
            return self.copy_shallow()
            
        copied = Ingredients(
            samples=self.samples.copy(),
            taxa=self.taxa.copy(),
            presence_matrix=self.presence_matrix.copy(),
            coverage_matrix=self.coverage_matrix.copy(),
            sample_to_biome=self.sample_to_biome.copy() if self.sample_to_biome is not None else None,
            data_release=self.data_release,
        )
        copied.presence_thresholds = getattr(self, "presence_thresholds", None)
        return copied
    
    def filter_samples(self, mask) -> None:
        """
        In-place filter of samples, keeping matrices and biome indices in sync.
        
        Args:
            mask (List[bool] | List[int] | np.ndarray): Boolean mask or list of indices of samples to keep.
        """
        idx = self._normalise_indexer(mask, len(self.samples))
        
        if idx.size == len(self.samples):
            return
            
        samples_arr = np.asarray(self.samples, dtype=object)
        self.samples = samples_arr[idx].tolist()
        
        self._presence_matrix = self.presence_matrix[:, idx]
        self._coverage_matrix = self.coverage_matrix[:, idx]
        self._presence_matrix_path = None
        self._coverage_matrix_path = None
        self._presence_matrix_shape = self._presence_matrix.shape
        self._coverage_matrix_shape = self._coverage_matrix.shape
        
        # sample filtering changes row totals
        self.total_counts = self._compute_total_counts()
        
        if getattr(self, "sample_to_biome", None):
            self._allocate_biomes()
    
    def filtered_samples(self, mask) -> 'Ingredients':
        """
        Return a new Ingredients instance filtered by samples, without deep-copying
        the backing sparse matrices first.
        """
        new = self.copy_shallow()
        new.filter_samples(mask)
        return new
    
    def filter_taxa(self, mask) -> None:
        """
        In-place filter of taxa (columns), keeping matrices and counts in sync.
        
        Args:
            mask (List[bool] | List[int] | np.ndarray): Boolean mask or list of indices of taxa to keep.
        """
        idx = self._normalise_indexer(mask, len(self.taxa))
        
        if idx.size == len(self.taxa):
            return
        
        taxa_arr = np.asarray(self.taxa, dtype=object)
        self.taxa = taxa_arr[idx].tolist()
        
        self._presence_matrix = self.presence_matrix[idx, :]
        self._coverage_matrix = self.coverage_matrix[idx, :]
        self._presence_matrix_path = None
        self._coverage_matrix_path = None
        self._presence_matrix_shape = self._presence_matrix.shape
        self._coverage_matrix_shape = self._coverage_matrix.shape
        
        # taxa filtering does not change per-row counts except subsetting them
        self.total_counts = self.total_counts[idx].astype(np.int32, copy=False)
        
        self._invalidate_taxa_caches()
        
    def filtered_taxa(self, mask) -> 'Ingredients':
        """
        Return a new Ingredients instance filtered by taxa, without deep-copying
        the backing sparse matrices first.
        """
        new = self.copy_shallow()
        new.filter_taxa(mask)
        return new
    
    def _allocate_biomes(self):
        """
        Precompute biome order and per-sample biome indices for both levels.
        """
        biomes_level_1: list[str] = []
        biomes_level_2: list[str] = []
        idxs_level_1: list[int] = []
        idxs_level_2: list[int] = []
        
        # maps biome -> index
        idx_map1: dict[str, int] = {}
        idx_map2: dict[str, int] = {}
        
        for s in self.samples:
            b1, b2 = self.sample_to_biome.get(s, (None, None))
            
            if b1 is None:
                idxs_level_1.append(-1)
            else:
                idx = idx_map1.get(b1)
                if idx is None:
                    idx = len(biomes_level_1)
                    biomes_level_1.append(b1)
                    idx_map1[b1] = idx
                idxs_level_1.append(idx)
                
            if b2 is None:
                idxs_level_2.append(-1)
            else:
                idx = idx_map2.get(b2)
                if idx is None:
                    idx = len(biomes_level_2)
                    biomes_level_2.append(b2)
                    idx_map2[b2] = idx
                idxs_level_2.append(idx)
                
        self.biomes_order = {
            "level_1": biomes_level_1,
            "level_2": biomes_level_2,
        }
        self.sample_biome_indices = {
            "level_1": np.array(idxs_level_1, dtype=int),
            "level_2": np.array(idxs_level_2, dtype=int),
        }

    def available_biomes(self) -> dict[str, list[dict[str, object]]]:
        """
        Return available biome query terms grouped by biome level.

        Each item contains:
          - biome: biome name
          - n_samples: number of samples assigned to that biome
        """
        if not getattr(self, "sample_to_biome", None):
            return {"level_1": [], "level_2": []}

        if not hasattr(self, "biomes_order") or not hasattr(self, "sample_biome_indices"):
            self._allocate_biomes()

        out: dict[str, list[dict[str, object]]] = {}
        for level in ("level_1", "level_2"):
            biomes = self.biomes_order.get(level, [])
            idxs = self.sample_biome_indices.get(level, np.array([], dtype=int))
            rows = []
            for idx, biome in enumerate(biomes):
                rows.append(
                    {
                        "biome": biome,
                        "n_samples": int(np.count_nonzero(idxs == idx)),
                    }
                )
            out[level] = rows
        return out
        
    def biome_distribution(self, level: str = "level_1"):
        """
        Compute biome distribution for the specified level.
        
        Args:
            level (str): Biome level to use ("level_1" or "level_2").
            
        Returns:
            Tuple[List[str], sp.csr_matrix, int]:
                - Unique biomes for the level.
                - Presence matrix for the level.
                - Number of samples with missing biome assignments.
        """
        if level not in ("level_1", "level_2"):
            raise ValueError("level must be 'level_1' or 'level_2'")
        
        biomes = self.biomes_order[level]
        idxs = self.sample_biome_indices[level]
        n_biomes = len(biomes)
        n_samples = len(idxs)
        
        # 1) Build biome-assignment matrix
        assigned = idxs >= 0
        rows = idxs[assigned]
        cols = np.nonzero(assigned)[0]
        data = np.ones_like(rows, dtype=int)
        B = sp.csr_matrix((data, (rows, cols)), shape=(n_biomes, n_samples))
        
        # 2) Presence counts
        Pbin = presence_for_counts(self)
        presence = B @ Pbin.T

        n_dropped = int((idxs < 0).sum())
        return biomes, presence, n_dropped

def _read_one_column_tsv(path: str, column: str) -> list[str]:
    df = pd.read_csv(path, sep="\t", dtype=str, usecols=lambda c: c == column)
    if column not in df.columns:
        raise ValueError(f"{path} is missing required column '{column}'")
    return df[column].fillna("").tolist()


def _write_one_column_tsv(path: str, column: str, values: list[str]) -> None:
    pd.DataFrame({column: values}).to_csv(path, sep="\t", index=False)


def _read_sample_to_biome(path: str) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    required = {"accession", "level_1", "level_2"}
    df = pd.read_csv(path, sep="\t", dtype=str, usecols=lambda c: c in required)
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")

    accessions = df["accession"].fillna("").to_numpy(dtype=object)
    level_1_series = df["level_1"].astype(object)
    level_2_series = df["level_2"].astype(object)
    level_1 = level_1_series.where(level_1_series.notna() & (level_1_series != ""), None).to_numpy(dtype=object)
    level_2 = level_2_series.where(level_2_series.notna() & (level_2_series != ""), None).to_numpy(dtype=object)

    return {
        sample: (b1, b2)
        for sample, b1, b2 in zip(accessions, level_1, level_2)
    }


def _write_sample_to_biome(path: str, sample_to_biome: Dict[str, Tuple[Optional[str], Optional[str]]]) -> None:
    rows = [
        {
            "accession": sample,
            "level_1": values[0] if values and values[0] is not None else "",
            "level_2": values[1] if values and values[1] is not None else "",
        }
        for sample, values in sample_to_biome.items()
    ]
    pd.DataFrame(rows, columns=["accession", "level_1", "level_2"]).to_csv(path, sep="\t", index=False)


def _ingredients_from_directory(directory: str) -> Ingredients:
    manifest_path = os.path.join(directory, "manifest.json")
    if not os.path.isdir(directory):
        raise ValueError(f"Ingredients path must be a directory: {directory}")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(f"Ingredients manifest not found: {manifest_path}")

    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    format_version = manifest.get("format_version")
    if (
        type(format_version) is not int
        or format_version != INGREDIENTS_FORMAT_VERSION
    ):
        raise IngredientsFormatError(
            f"Unsupported Ingredients format version {format_version!r} in "
            f"{manifest_path}; this metacooc version supports "
            f"{INGREDIENTS_FORMAT_VERSION}."
        )

    components = manifest.get("components", {})
    matrix_shapes = manifest.get("matrix_shapes", {})
    samples = _read_one_column_tsv(os.path.join(directory, components.get("samples", "samples.tsv")), "sample")
    taxa = _read_one_column_tsv(os.path.join(directory, components.get("taxa", "taxa.tsv")), "taxon")
    sample_to_biome = _read_sample_to_biome(
        os.path.join(directory, components.get("sample_to_biome", "sample_to_biome.tsv"))
    )
    total_counts = np.load(os.path.join(directory, components.get("total_counts", "total_counts.npy")))

    new = Ingredients.__new__(Ingredients)
    new.samples = samples
    new.taxa = taxa
    object.__setattr__(new, "_presence_matrix", None)
    object.__setattr__(
        new,
        "_presence_matrix_path",
        os.path.join(directory, components.get("presence_matrix", "presence.npz")),
    )
    object.__setattr__(new, "_coverage_matrix", None)
    object.__setattr__(
        new,
        "_coverage_matrix_path",
        os.path.join(directory, components.get("coverage_matrix", "coverage.npz")),
    )
    object.__setattr__(new, "_presence_matrix_shape", tuple(matrix_shapes.get("presence", (len(taxa), len(samples)))))
    object.__setattr__(new, "_coverage_matrix_shape", tuple(matrix_shapes.get("coverage", (len(taxa), len(samples)))))
    new.total_counts = total_counts.astype(np.int32, copy=False)
    new.sample_to_biome = sample_to_biome
    if new.sample_to_biome:
        new._allocate_biomes()
    new._rank_lookups = None
    new._terminal_rank_prefixes = None
    new.data_release = manifest.get("data_release")
    new.presence_thresholds = manifest.get("presence_thresholds")
    return new


def _write_ingredients_directory(
    ingredients: "Ingredients",
    directory: str,
    *,
    aggregated: bool = False,
) -> str:
    os.makedirs(directory, exist_ok=True)

    presence = _presence_csr(ingredients.presence_matrix)
    coverage = _csr(ingredients.coverage_matrix)
    total_counts = np.asarray(ingredients.total_counts, dtype=np.int32)

    components = {
        "samples": "samples.tsv",
        "taxa": "taxa.tsv",
        "sample_to_biome": "sample_to_biome.tsv",
        "presence_matrix": "presence.npz",
        "coverage_matrix": "coverage.npz",
        "total_counts": "total_counts.npy",
    }

    _write_one_column_tsv(os.path.join(directory, components["samples"]), "sample", list(ingredients.samples))
    _write_one_column_tsv(os.path.join(directory, components["taxa"]), "taxon", list(ingredients.taxa))
    _write_sample_to_biome(
        os.path.join(directory, components["sample_to_biome"]),
        getattr(ingredients, "sample_to_biome", {}) or {},
    )
    sp.save_npz(os.path.join(directory, components["presence_matrix"]), presence)
    sp.save_npz(os.path.join(directory, components["coverage_matrix"]), coverage)
    np.save(os.path.join(directory, components["total_counts"]), total_counts)

    manifest = {
        "format_version": INGREDIENTS_FORMAT_VERSION,
        "date_generated": _utc_timestamp_seconds(),
        "data_release": getattr(ingredients, "data_release", None),
        "aggregated": bool(aggregated),
        "matrix_shapes": {
            "presence": list(presence.shape),
            "coverage": list(coverage.shape),
        },
        "matrix_dtypes": {
            "presence": str(presence.dtype),
            "coverage": str(coverage.dtype),
        },
        "components": components,
    }
    presence_thresholds = getattr(ingredients, "presence_thresholds", None)
    if presence_thresholds:
        manifest["presence_thresholds"] = presence_thresholds

    with open(os.path.join(directory, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")
    return directory


def archive_ingredients_directory(directory: str, archive_path: Optional[str] = None) -> str:
    archive_path = archive_path or f"{directory}.tar.gz"
    base_name = os.path.basename(directory.rstrip(os.sep))
    with tarfile.open(archive_path, "w:gz") as tar:
        tar.add(directory, arcname=base_name)
    return archive_path


def save_ingredients_directory(
    ingredients: "Ingredients",
    directory: str,
    *,
    aggregated: bool = False,
    archive: bool = False,
) -> str:
    path = _write_ingredients_directory(ingredients, directory, aggregated=aggregated)
    if archive:
        archive_ingredients_directory(path)
    print(f"Saved Ingredients directory -> {path}")
    return path


def load_ingredients(
    data_dir: Optional[str] = None,
    aggregated: bool = False,
    custom_ingredients=None,
    data_release: Optional[str] = None,
    sample_to_biome_file=None) -> Ingredients:
    """Load an Ingredients object and associated biome mapping."""
    
    # determine ingredients file path
    if not custom_ingredients:
        defaulted_data_release = data_release is None
        data_release = data_release or LATEST_DATA_RELEASE
        filenames, _ = get_file_info(data_release)
        if not data_dir:
            raise ValueError(
                "data_dir must be provided when not using custom_ingredients"
            )
        key = f"ingredients_aggregated" if aggregated else "ingredients_raw"
        filepath = os.path.join(data_dir, filenames[key])
    else:
        if isinstance(custom_ingredients, Ingredients):
            return custom_ingredients
        filepath = custom_ingredients
    
    
    if not os.path.exists(filepath):
        if custom_ingredients:
            raise FileNotFoundError(
            f"{custom_ingredients} is either not found or isn't an Ingredients object"
            )
        
        raise DataReleaseError(
            missing_data_release_message(
                data_dir=data_dir,
                data_release=data_release,
                missing_path=filepath,
                file_kind="Ingredients",
                defaulted=defaulted_data_release,
            )
        )
    ingredients = _ingredients_from_directory(filepath)
    
    # A requested release is an assertion, including for an explicit custom
    # directory.  Without one, the manifest remains authoritative.
    if data_release is not None:
        ev = getattr(ingredients, "data_release", None)
        if ev != data_release:
            raise DataReleaseError(
                f"Loaded data release {ev!r}, expected {data_release!r}."
            )
    if not custom_ingredients:
        print(f"Using {filepath}")
    
    return ingredients

def save_ingredients(ingredients: "Ingredients", 
                     output_dir: str, 
                     *, 
                     aggregated: bool = False,
                     tag: Optional[str] = None, 
                     data_release: Optional[str] = None,
                     archive: bool = False) -> str:
    os.makedirs(output_dir, exist_ok=True)
    
    if data_release is not None:
        ingredients.data_release = data_release
        
    kind = "ingredients_aggregated" if aggregated else "ingredients_raw"
    label = tag or data_release or getattr(ingredients, "data_release", None)
    suffix = f"_{label}" if label else ""
    dirpath = os.path.join(output_dir, f"{kind}{suffix}")
    path = _write_ingredients_directory(ingredients, dirpath, aggregated=aggregated)
    if archive:
        archive_path = archive_ingredients_directory(path)
        print(f"Archived Ingredients directory -> {archive_path}")
    
    print(f"Saved Ingredients directory -> {path}")
    return path
