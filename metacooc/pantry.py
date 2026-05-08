#!/usr/bin/env python3

import os
import json
import tarfile
import warnings
import re
from collections import defaultdict

from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

from metacooc._data_config import *
from metacooc.utils import (
    _RANK_PREFIXES, 
    _PREFIX_TO_RANK, 
    _parse_tokens, 
    _token_rank, 
    _terminal_rank_prefix, 
    _deepest_rank_token
)

INGREDIENTS_FORMAT_VERSION = 1


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


class Ingredients:
    """
    Container for metagenomic data with sparse matrices.
    
    Attributes:
        samples (List[str]): Sample identifiers.
        taxa (List[str]): Taxonomic identifiers.
        presence_matrix (sp.csr_matrix): Binary presence/absence.
        coverage_matrix (sp.csr_matrix): Coverage values.
        total_counts (np.ndarray): Cached per-taxon counts.
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
        data_version: Optional[str] = None,
        ):
        self.taxa = taxa
        self.samples = samples
        
        object.__setattr__(self, "_presence_matrix", _presence_csr(presence_matrix))
        object.__setattr__(self, "_coverage_matrix", _csr(coverage_matrix))
        object.__setattr__(self, "_presence_matrix_path", None)
        object.__setattr__(self, "_coverage_matrix_path", None)
        
        self.total_counts = self._compute_total_counts()
        
        # Load and preallocate biome mapping
        self.sample_to_biome = sample_to_biome or {}
        if self.sample_to_biome:
            self._allocate_biomes()
            
        self._rank_lookups = None
        self._terminal_rank_prefixes = None
        self.data_version = data_version
    
    def __getstate__(self):
        state = {
            "samples": self.samples,
            "taxa": self.taxa,
            "_presence_matrix": self.presence_matrix,
            "_coverage_matrix": self.coverage_matrix,
            "total_counts": self.total_counts,
            "sample_to_biome": self.sample_to_biome,
            "data_version": self.data_version,
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
        
        self.data_version = state.get("data_version", None)
    
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
        return self._presence_matrix
    
    @presence_matrix.setter
    def presence_matrix(self, mat: sp.csr_matrix):
        self._presence_matrix = _presence_csr(mat)
        self._presence_matrix_path = None
        self.total_counts = self._compute_total_counts()
    
    @property
    def coverage_matrix(self) -> sp.csr_matrix:
        if self._coverage_matrix is None:
            self._coverage_matrix = sp.load_npz(self._coverage_matrix_path).tocsr()
        return self._coverage_matrix
    
    @coverage_matrix.setter
    def coverage_matrix(self, mat: sp.csr_matrix):
        self._coverage_matrix = _csr(mat)
        self._coverage_matrix_path = None
    
    def _compute_total_counts(self) -> np.ndarray:
        """
        Per-taxon presence counts across samples.
        
        Because presence_matrix is already binary / sparse-presence, use getnnz(axis=1)
        instead of (matrix > 0).sum(axis=1), which allocates a large temporary sparse matrix.
        """
        return np.asarray(self.presence_matrix.getnnz(axis=1), dtype=np.int32).ravel()

        # return np.array((self._presence_matrix > 0).sum(axis=1)).flatten()
    
    def __repr__(self):
        return (
            f"<Ingredients: {len(self.taxa)} taxa, "
            f"{len(self.samples)} samples, "
            f"presence: {self.presence_matrix.shape}, "
            f"coverage: {self.coverage_matrix.shape}>"
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
        
        new.total_counts = self.total_counts
        new.sample_to_biome = self.sample_to_biome
        new.data_version = self.data_version
        
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
            
        return Ingredients(
            samples=self.samples.copy(),
            taxa=self.taxa.copy(),
            presence_matrix=self.presence_matrix.copy(),
            coverage_matrix=self.coverage_matrix.copy(),
            sample_to_biome=self.sample_to_biome.copy() if self.sample_to_biome is not None else None,
            data_version=self.data_version,
        )
    
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
            Tuple[List[str], sp.csr_matrix, sp.csr_matrix, int]:
                - Unique biomes for the level.
                - Presence matrix for the level.
                - Coverage matrix for the level.
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
        
        # 3) Coverage means
        coverage_sums = B @ self.coverage_matrix.T
        counts = np.array(B.sum(axis=1)).ravel()
        # Make a copy for means
        coverage = coverage_sums.tolil()
        for b in range(n_biomes):
            if counts[b] > 0:
                coverage[b, :] = coverage_sums[b, :].toarray() / counts[b]
        coverage = coverage.tocsr()
        
        n_dropped = int((idxs < 0).sum())
        return biomes, presence, coverage, n_dropped

def _samples_path(directory: str) -> str:
    return os.path.join(directory, "samples.tsv")


def _taxa_path(directory: str) -> str:
    return os.path.join(directory, "taxa.tsv")


def _read_one_column_tsv(path: str, column: str) -> list[str]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    if column not in df.columns:
        raise ValueError(f"{path} is missing required column '{column}'")
    return df[column].fillna("").tolist()


def _write_one_column_tsv(path: str, column: str, values: list[str]) -> None:
    pd.DataFrame({column: values}).to_csv(path, sep="\t", index=False)


def _read_sample_to_biome(path: str) -> Dict[str, Tuple[Optional[str], Optional[str]]]:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"accession", "level_1", "level_2"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
    out = {}
    for _, row in df.iterrows():
        b1 = row["level_1"] or None
        b2 = row["level_2"] or None
        out[row["accession"]] = (b1, b2)
    return out


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

    components = manifest.get("components", {})
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
    new.total_counts = total_counts.astype(np.int32, copy=False)
    new.sample_to_biome = sample_to_biome
    if new.sample_to_biome:
        new._allocate_biomes()
    new._rank_lookups = None
    new._terminal_rank_prefixes = None
    new.data_version = manifest.get("data_version")
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
        "data_version": getattr(ingredients, "data_version", None),
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
    data_version: Optional[str] = None,
    sample_to_biome_file=None) -> Ingredients:
    """Load an Ingredients object and associated biome mapping."""
    
    # determine ingredients file path
    if not custom_ingredients:
        defaulted_data_version = data_version is None
        data_version = data_version or LATEST_VERSION
        filenames, _ = get_file_info(data_version)
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
        
        raise DataVersionError(
            missing_data_version_message(
                data_dir=data_dir,
                data_version=data_version,
                missing_path=filepath,
                file_kind="Ingredients",
                defaulted=defaulted_data_version,
            )
        )
    ingredients = _ingredients_from_directory(filepath)
    
    # version mismatch warning
    if not custom_ingredients:
        ev = getattr(ingredients, "data_version", None)
        if ev and ev != data_version:
            warnings.warn(
                f"Loaded version {ev}, expected {data_version}.", UserWarning
            )
        print(f"Using {filepath}")
    
    return ingredients

def save_ingredients(ingredients: "Ingredients", 
                     output_dir: str, 
                     *, 
                     aggregated: bool = False,
                     tag: Optional[str] = None, 
                     data_version: Optional[str] = None,
                     archive: bool = False) -> str:
    os.makedirs(output_dir, exist_ok=True)
    
    if data_version is not None:
        ingredients.data_version = data_version
        
    kind = "ingredients_aggregated" if aggregated else "ingredients_raw"
    label = tag or data_version
    suffix = f"_{label}" if label else ""
    dirpath = os.path.join(output_dir, f"{kind}{suffix}")
    path = _write_ingredients_directory(ingredients, dirpath, aggregated=aggregated)
    if archive:
        archive_path = archive_ingredients_directory(path)
        print(f"Archived Ingredients directory -> {archive_path}")
    
    print(f"Saved Ingredients directory -> {path}")
    return path
