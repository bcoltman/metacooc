#!/usr/bin/env python3

import os
import pickle
import warnings
import re
from collections import defaultdict

from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

from metacooc._data_config import *
from metacooc.utils import (
    _RANK_PREFIXES, 
    _PREFIX_TO_RANK, 
    _parse_tokens, 
    _token_rank, 
    _terminal_rank_prefix, 
    _deepest_rank_token
)


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
        
        object.__setattr__(self, "_presence_matrix", presence_matrix)
        object.__setattr__(self, "_coverage_matrix", coverage_matrix)
        
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
            "_presence_matrix": self._presence_matrix,
            "_coverage_matrix": self._coverage_matrix,
            "total_counts": self.total_counts,
            "sample_to_biome": self.sample_to_biome,
            "data_version": self.data_version,
        }
        
        if hasattr(self, "biomes_order"):
            state["biomes_order"] = self.biomes_order
            state["sample_biome_indices"] = self.sample_biome_indices
            
        if hasattr(self, "_rank_lookups") and self._rank_lookups is not None:
            state["_rank_lookups"] = self._rank_lookups
            state["_terminal_rank_prefixes"] = self._terminal_rank_prefixes
            
        return state
    
    def __setstate__(self, state):
        self.samples = state["samples"]
        self.taxa = state["taxa"]
        
        object.__setattr__(self, "_presence_matrix", state.get("_presence_matrix"))
        object.__setattr__(self, "_coverage_matrix", state.get("_coverage_matrix"))
        
        self._rank_lookups = state.get("_rank_lookups", None)
        self._terminal_rank_prefixes = state.get("_terminal_rank_prefixes", None)
        
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
        return self._presence_matrix
    
    @presence_matrix.setter
    def presence_matrix(self, mat: sp.csr_matrix):
        self._presence_matrix = mat
        self.total_counts = self._compute_total_counts()
    
    @property
    def coverage_matrix(self) -> sp.csr_matrix:
        return self._coverage_matrix
    
    @coverage_matrix.setter
    def coverage_matrix(self, mat: sp.csr_matrix):
        self._coverage_matrix = mat
    
    def _compute_total_counts(self) -> np.ndarray:
        return np.array((self._presence_matrix > 0).sum(axis=1)).flatten()
    
    def __repr__(self):
        return (
            f"<Ingredients: {len(self.taxa)} taxa, "
            f"{len(self.samples)} samples, "
            f"presence: {self.presence_matrix.shape}, "
            f"coverage: {self.coverage_matrix.shape}>"
        )
    
    def copy(self):
        return Ingredients(
            samples=self.samples.copy(),
            taxa=self.taxa.copy(),
            presence_matrix=self.presence_matrix.copy(),
            coverage_matrix=self.coverage_matrix.copy(),
            sample_to_biome=self.sample_to_biome.copy(),
        )
    
    def filter_samples(self, mask) -> None:
        """
        In-place filter of samples, keeping matrices and biome indices in sync.
        
        Args:
            mask (List[bool] | List[int] | np.ndarray): Boolean mask or list of indices of samples to keep.
        """
        import numpy as _np
        if isinstance(mask, (_np.ndarray, list)):
            arr = _np.array(mask)
            if arr.dtype == bool:
                idxs = _np.nonzero(arr)[0].tolist()
            else:
                idxs = arr.astype(int).tolist()
        else:
            raise ValueError("mask must be a list or numpy array of bools or ints")
        # apply
        self.samples = [self.samples[i] for i in idxs]
        self._presence_matrix = self._presence_matrix[:, idxs]
        self._coverage_matrix = self._coverage_matrix[:, idxs]
        self.total_counts = self._compute_total_counts()
        if getattr(self, 'sample_to_biome', None):
            self._allocate_biomes()
    
    def filtered_samples(self, mask) -> 'Ingredients':
        """
        Return a new Ingredients instance filtered by samples.
        """
        new = self.copy()
        new.filter_samples(mask)
        return new
    
    def filter_taxa(self, mask) -> None:
        """
        In-place filter of taxa (columns), keeping matrices and counts in sync.
        
        Args:
            mask (List[bool] | List[int] | np.ndarray): Boolean mask or list of indices of taxa to keep.
        """
        import numpy as _np
        if isinstance(mask, (_np.ndarray, list)):
            arr = _np.array(mask)
            if arr.dtype == bool:
                idxs = _np.nonzero(arr)[0].tolist()
            else:
                idxs = arr.astype(int).tolist()
        else:
            raise ValueError("mask must be a list or numpy array of bools or ints")
        # apply
        self.taxa = [self.taxa[i] for i in idxs]
        self._presence_matrix = self._presence_matrix[idxs, :]
        self._coverage_matrix = self._coverage_matrix[idxs, :]
        # update counts and caches
        self.total_counts = self._compute_total_counts()
        self._invalidate_taxa_caches()
        
    def filtered_taxa(self, mask) -> 'Ingredients':
        """
        Return a new Ingredients instance filtered by taxa.
        """
        new = self.copy()
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
        Pbin = (self._presence_matrix > 0).astype(int)
        presence = B @ Pbin.T
        
        # 3) Coverage means
        coverage_sums = B @ self._coverage_matrix.T
        counts = np.array(B.sum(axis=1)).ravel()
        # Make a copy for means
        coverage = coverage_sums.tolil()
        for b in range(n_biomes):
            if counts[b] > 0:
                coverage[b, :] = coverage_sums[b, :].toarray() / counts[b]
        coverage = coverage.tocsr()
        
        n_dropped = int((idxs < 0).sum())
        return biomes, presence, coverage, n_dropped



def load_ingredients(
    data_dir: Optional[str] = None,
    aggregated: bool = False,
    custom_ingredients=None,
    data_version: Optional[str] = None,
    sample_to_biome_file=None) -> Ingredients:
    """Load an Ingredients object and associated biome mapping."""
    
    # determine ingredients file path
    if not custom_ingredients:
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
        
        avail = ", ".join(sorted(RELEASES.keys()))
        raise FileNotFoundError(
            f"Ingredients file '{filepath}' not found.\n"
            f"Missing version {data_version}. Available: {avail}"
        )
    
    # load ingredients object
    with open(filepath, "rb") as f:
        ingredients = pickle.load(f)
    
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
                     data_version: Optional[str] = None) -> str:
    os.makedirs(output_dir, exist_ok=True)
    
    if data_version is not None:
        ingredients.data_version = data_version
        
    kind = "ingredients_aggregated" if aggregated else "ingredients_raw"
    suffix = f"_{tag}" if tag else ""
    filepath = os.path.join(output_dir, f"{kind}{suffix}.pkl")
    
    with open(filepath, "wb") as f:
        pickle.dump(ingredients, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    print(f"Saved Ingredients → {filepath}")
    return filepath
