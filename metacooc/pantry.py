#!/usr/bin/env python3

import os
import pickle
import warnings
from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
# from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from metacooc._data_config import *


class Ingredients:
    """
    Container for metagenomic data with sparse matrices.
    
    Attributes:
        samples (List[str]): Sample identifiers.
        taxa (List[str]): Taxonomic identifiers.
        presence_matrix (sp.csr_matrix): Binary presence/absence.
        coverage_matrix (sp.csr_matrix): Coverage values.
        total_counts (np.ndarray): Cached per-taxon counts.
        sample_to_biome (Dict[str,str]): Mapping sample→biome.
        biomes_order (List[str]): Unique biomes.
        sample_biome_indices (np.ndarray): Per-sample biome index (-1 if missing).
    """
    def __init__(
        self,
        samples: List[str],
        taxa: List[str],
        presence_matrix: sp.csr_matrix,
        coverage_matrix: sp.csr_matrix,
        sample_to_biome: Dict[str, str] = None,
    ):
        self.samples = samples
        self.taxa = taxa
        object.__setattr__(self, "_presence_matrix", presence_matrix)
        object.__setattr__(self, "_coverage_matrix", coverage_matrix)
        self.total_counts = self._compute_total_counts()
        
        # Load and preallocate biome mapping if provided
        self.sample_to_biome = sample_to_biome or {}
        if self.sample_to_biome:
            self._allocate_biomes()
    
    def __getstate__(self):
        state = {
            "samples": self.samples,
            "taxa": self.taxa,
            "_presence_matrix": self._presence_matrix,
            "_coverage_matrix": self._coverage_matrix,
            "total_counts": self.total_counts,
            "sample_to_biome": self.sample_to_biome,
        }
        # include precomputed if exists
        if hasattr(self, "biomes_order"):
            state["biomes_order"] = self.biomes_order
            state["sample_biome_indices"] = self.sample_biome_indices
        return state
    
    def __setstate__(self, state):
        self.samples = state["samples"]
        self.taxa = state["taxa"]
        object.__setattr__(self, "_presence_matrix", state.get("_presence_matrix"))
        object.__setattr__(self, "_coverage_matrix", state.get("_coverage_matrix"))
        
        # restore or compute total_counts
        if "total_counts" in state and state["total_counts"] is not None:
            self.total_counts = state["total_counts"]
        else:
            self.total_counts = self._compute_total_counts()
        
        # restore or default biome mapping
        self.sample_to_biome = state.get("sample_to_biome", {})
        # restore allocation if present
        if "biomes_order" in state and "sample_biome_indices" in state:
            self.biomes_order = state["biomes_order"]
            self.sample_biome_indices = state["sample_biome_indices"]
        elif self.sample_to_biome:
            self._allocate_biomes()
    
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
        return np.array((self._presence_matrix > 0).sum(axis=0)).flatten()
    
    def __repr__(self):
        return (
            f"<Ingredients: {len(self.samples)} samples, "
            f"{len(self.taxa)} taxa, "
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
        self._presence_matrix = self._presence_matrix[idxs, :]
        self._coverage_matrix = self._coverage_matrix[idxs, :]
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
        self._presence_matrix = self._presence_matrix[:, idxs]
        self._coverage_matrix = self._coverage_matrix[:, idxs]
        # update counts and caches
        self.total_counts = self._compute_total_counts()
        
    def filtered_taxa(self, mask) -> 'Ingredients':
        """
        Return a new Ingredients instance filtered by taxa.
        """
        new = self.copy()
        new.filter_taxa(mask)
        return new
    
    def _allocate_biomes(self):
        """
        Precompute biome order & per-sample biome indices.
        """
        biomes: List[str] = []
        idxs: List[int] = []
        for s in self.samples:
            b = self.sample_to_biome.get(s)
            if b is None:
                idxs.append(-1)
            else:
                if b not in biomes:
                    biomes.append(b)
                idxs.append(biomes.index(b))
        self.biomes_order = biomes
        self.sample_biome_indices = np.array(idxs, dtype=int)
    
    def biome_distribution(self):
        biomes = self.biomes_order
        idxs   = self.sample_biome_indices
        n_biomes = len(biomes)
        n_samples = len(idxs)
        
        # 1) build biome‐assignment matrix
        assigned = idxs >= 0
        rows = idxs[assigned]
        cols = np.nonzero(assigned)[0]
        data = np.ones_like(rows, dtype=int)
        B = sp.csr_matrix((data, (rows, cols)), shape=(n_biomes, n_samples))
        
        # 2) presence counts
        Pbin = (self._presence_matrix > 0).astype(int)
        presence = B @ Pbin
        
        # 3) coverage means
        coverage_sums = B @ self._coverage_matrix
        counts = np.array(B.sum(axis=1)).ravel()
        # make a copy for means
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
    sandpiper_version: Optional[str] = None,
    sample_to_biome_file=None) -> Ingredients:
    """Load an Ingredients object and associated biome mapping."""
    
    # determine ingredients file path
    if not custom_ingredients:
        version = sandpiper_version or LATEST_VERSION
        filenames, _ = get_file_info(version)
        if not data_dir:
            raise ValueError(
                "data_dir must be provided when not using custom_ingredients"
            )
        # key = f"ingredients_aggregated_{aggregation_level}" if aggregated else "ingredients_raw"
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
            f"Missing version {version}. Available: {avail}"
        )
    
    # load ingredients object
    with open(filepath, "rb") as f:
        ingredients = pickle.load(f)
    
    # version mismatch warning
    if not custom_ingredients:
        ev = getattr(ingredients, "version", None)
        if ev and ev != version:
            warnings.warn(
                f"Loaded version {ev}, expected {version}.", UserWarning
            )
        print(f"Using {filepath}")
    
    return ingredients
