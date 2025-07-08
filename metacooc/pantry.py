#!/usr/bin/env python3

import os
import warnings
import pickle

import scipy.sparse as sp
from typing import List
import numpy as np

# from metacooc._data_config import DATA_VERSION, FILENAMES
from metacooc._data_config import *

class Ingredients:
    """
    Container for metagenomic data with sparse matrices.
    
    Attributes:
        samples (List[str]): List of sample identifiers.
        taxa (List[str]): List of taxonomic identifiers.
        presence_matrix (sp.csr_matrix): Sparse binary matrix indicating presence/absence.
        coverage_matrix (sp.csr_matrix): Sparse matrix containing coverage values.
        total_counts (np.ndarray): Cached counts of each taxon (number of samples with nonzero presence).
    """
    def __init__(self, samples: List[str], taxa: List[str],
                 presence_matrix: sp.csr_matrix, coverage_matrix: sp.csr_matrix):
        self.samples = samples
        self.taxa = taxa
        self.presence_matrix = presence_matrix
        self.coverage_matrix = coverage_matrix
        self.total_counts = self._compute_total_counts()

    def _compute_total_counts(self) -> np.ndarray:
        """
        Compute the total counts for each taxon (number of samples with presence).
        """
        counts = np.array((self.presence_matrix > 0).sum(axis=0)).flatten()
        return counts

    def __repr__(self):
        return (f"<Ingredients: {len(self.samples)} samples, "
                f"{len(self.taxa)} taxa, "
                f"presence_matrix shape: {self.presence_matrix.shape}, "
                f"coverage_matrix shape: {self.coverage_matrix.shape}>")
    
    def copy(self):
        """
        Create and return a deep copy of the current Ingredients instance.
        """
        return Ingredients(
            samples=self.samples.copy(),
            taxa=self.taxa.copy(),
            presence_matrix=self.presence_matrix.copy(),
            coverage_matrix=self.coverage_matrix.copy()
        )



def load_ingredients(data_dir, aggregated=False, custom_ingredients=None, sandpiper_version=None):
    """Load an Ingredients object, checking version and download status."""
    
    version = sandpiper_version or LATEST_VERSION
    filenames, _ = get_file_info(version)
    
    if not custom_ingredients:
        filename = filenames["ingredients_aggregated"] if aggregated else filenames["ingredients_raw"]
        filepath = os.path.join(data_dir, filename)
    else:
        if isinstance(custom_ingredients, Ingredients):
            return custom_ingredients
        filepath = custom_ingredients
    
    if not os.path.exists(filepath):
        version_list = ", ".join(sorted(RELEASES.keys()))
        raise FileNotFoundError(
            f"Ingredients file '{filepath}' not found.\n"
            f"You may be missing version {version}. Available versions: {version_list}\n"
            f"Either download the correct data or specify a version using --sandpiper_version."
        )
    
    with open(filepath, "rb") as f:
        ingredients = pickle.load(f)
    
    # Skip version check if custom ingredient provided
    if not custom_ingredients:
        embedded_version = getattr(ingredients, "version", None)
        if embedded_version and embedded_version != version:
            warnings.warn(
                f"Loaded Ingredients object is version {embedded_version}, "
                f"but expected {version}. This might indicate a mismatch.",
                UserWarning,
            )
    
    return ingredients
    
    # """Load an Ingredients object and check its version."""
    # if not custom_ingredients:
        # if aggregated:
            # filename = FILENAMES["ingredients_aggregated"]
        # else:
            # filename = FILENAMES["ingredients_raw"]
        # filepath = os.path.join(data_dir, filename)
    # else:
        # if isinstance(custom_ingredients, Ingredients):
            # return custom_ingredients
        
        # filepath = custom_ingredients
    
    # if not os.path.exists(filepath):
        # raise FileNotFoundError(f"Ingredients file '{filepath}' not found.")
        
    # with open(filepath, "rb") as f:
        # ingredients = pickle.load(f)
        
    # if custom_ingredients:
        # return ingredients
        
    # # Check embedded version if present
    # embedded_version = getattr(ingredients, "version", None)
    # if embedded_version and embedded_version != DATA_VERSION:
        # warnings.warn(
            # f"Loaded Ingredients object is version {embedded_version}, "
            # f"but expected {DATA_VERSION}. You may be using outdated or mismatched data.",
            # UserWarning
        # )
        
    # return ingredients