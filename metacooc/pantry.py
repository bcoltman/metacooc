#!/usr/bin/env python3

import os
import warnings
import pickle

import scipy.sparse as sp
from typing import List
import numpy as np

from metacooc._data_config import DATA_VERSION, FILENAMES

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



def load_ingredients(data_dir, aggregated=False):
    """Load an Ingredients object and check its version."""
    filename = FILENAMES["ingredients_aggregated"] if aggregated else FILENAMES["ingredients_raw"]
    filepath = os.path.join(data_dir, filename)
    
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Ingredients file '{filepath}' not found.")
        
    with open(filepath, "rb") as f:
        ingredients = pickle.load(f)
        
    # Check embedded version if present
    embedded_version = getattr(ingredients, "version", None)
    if embedded_version and embedded_version != DATA_VERSION:
        warnings.warn(
            f"Loaded Ingredients object is version {embedded_version}, "
            f"but expected {DATA_VERSION}. You may be using outdated or mismatched data.",
            UserWarning
        )
        
    return ingredients