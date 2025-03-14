#!/usr/bin/env python3
"""
format.py

Converts a raw sandpiper taxonomic profile (TSV file) into two sparse matrix representations:
    1. A raw Ingredients object (species-level only)
    2. An aggregated Ingredients object (with additional taxonomic levels added)

Also provides functions for building metadata indices:
    - build_broad_metadata_index: Creates a tokenized index over specified metadata columns.
    - build_exact_metadata_index: Creates an index for exact matching on a given column.

Usage (for formatting):
  metacooc format --data_dir /path/to/data --output_dir /path/to/out [--pattern g__]

The script expects a file "sandpiper.tsv" in the data directory.
"""

import os
import re
import argparse
import pickle
import copy

import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix, hstack

from metacooc.pantry import Ingredients  # Simple container class for your data

# --- Functions for creating Ingredients objects ---

def build_indices(tax_profile: str):
    """
    Build indices for samples and taxa from the input TSV file.
    
    Returns:
        samples (list): Ordered list of unique sample identifiers.
        taxa (list): Ordered list of unique taxon identifiers.
        sample_to_index (dict): Mapping from sample to its index.
        taxon_to_index (dict): Mapping from taxon to its index.
    """
    chunk_iterator = pd.read_csv(
        tax_profile,
        delimiter='\t',
        dtype=str,
        usecols=['sample', 'taxonomy'],
        chunksize=100000,
        engine="c"
    )
    
    sample_to_index = {}
    taxon_to_index = {}
    
    for chunk in chunk_iterator:
        for sample in chunk['sample'].unique():
            if sample not in sample_to_index:
                sample_to_index[sample] = len(sample_to_index)
        for taxon in chunk['taxonomy'].unique():
            if taxon not in taxon_to_index:
                taxon_to_index[taxon] = len(taxon_to_index)
    
    samples = [None] * len(sample_to_index)
    for sample, idx in sample_to_index.items():
        samples[idx] = sample
        
    taxa = [None] * len(taxon_to_index)
    for taxon, idx in taxon_to_index.items():
        taxa[idx] = taxon
        
    return samples, taxa, sample_to_index, taxon_to_index

def create_sparse_matrices(tax_profile: str, sample_to_index: dict, taxon_to_index: dict, samples: list, taxa: list):
    """
    Create presence and coverage sparse matrices from the TSV file.
    
    Returns:
        presence_matrix (csr_matrix): Sparse binary matrix.
        coverage_matrix (csr_matrix): Sparse matrix with coverage values.
    """
    data_presence = []
    data_coverage = []
    row_indices = []
    col_indices = []
    
    chunk_iterator = pd.read_csv(
        tax_profile,
        delimiter='\t',
        dtype=str,
        chunksize=100000,
        engine="c"
    )
    
    for chunk in chunk_iterator:
        sample_indices = chunk['sample'].map(sample_to_index)
        taxon_indices = chunk['taxonomy'].map(taxon_to_index)
        
        coverage_values = pd.to_numeric(chunk['coverage'], errors='coerce').fillna(0)
        
        data_presence.extend([1] * len(chunk))
        data_coverage.extend(coverage_values)
        row_indices.extend(sample_indices)
        col_indices.extend(taxon_indices)
    
    num_samples = len(samples)
    num_taxa = len(taxa)
    
    presence_matrix = csr_matrix(
        (data_presence, (row_indices, col_indices)),
        shape=(num_samples, num_taxa),
        dtype=int
    )
    coverage_matrix = csr_matrix(
        (data_coverage, (row_indices, col_indices)),
        shape=(num_samples, num_taxa),
        dtype=float
    )
    
    return presence_matrix, coverage_matrix

def format_ingredients(tax_profile: str):
    """
    Process the sandpiper TSV file and return a raw Ingredients instance.
    """
    samples, taxa, sample_to_index, taxon_to_index = build_indices(tax_profile)
    presence_matrix, coverage_matrix = create_sparse_matrices(tax_profile, sample_to_index, taxon_to_index, samples, taxa)
    
    # Generate the raw Ingredients object.
    raw_ingredients = Ingredients(samples, taxa, presence_matrix, coverage_matrix)
    
    return raw_ingredients
    

def format_data(tax_profile: str, output_dir: str, aggregated=False, aggregation_pattern=None):
    """
    Process the sandpiper TSV file and return a raw Ingredients instance.
    """
    # samples, taxa, sample_to_index, taxon_to_index = build_indices(tax_profile)
    # presence_matrix, coverage_matrix = create_sparse_matrices(tax_profile, sample_to_index, taxon_to_index, samples, taxa)
    
    # Generate the raw Ingredients object.
    raw_ingredients = format_ingredients(tax_profile)
    
    # Save both objects.
    output_raw = os.path.join(output_dir, "ingredients_raw.pkl")
    
    with open(output_raw, "wb") as f:
        pickle.dump(raw_ingredients, f)
    
    print(f"Raw ingredients saved to {output_raw}")
    
    if aggregated & (aggregation_pattern is not None):
        # Generate the aggregated Ingredients object.
        aggregated_ingredients = copy.deepcopy(raw_ingredients)
        aggregated_ingredients = add_taxa_levels_to_ingredients(aggregated_ingredients, aggregation_pattern=aggregation_pattern)
        
        output_agg = os.path.join(output_dir, "ingredients_aggregated.pkl")
        
        with open(output_agg, "wb") as f:
            pickle.dump(aggregated_ingredients, f)
        print(f"Aggregated ingredients saved to {output_agg}")

def add_taxa_levels_to_ingredients(ingredients: Ingredients, aggregation_pattern="g__") -> Ingredients:
    """
    Add aggregated taxonomic levels to an Ingredients instance.
    This function splits the raw taxa strings, extracts the level specified by the regex aggregation_pattern,
    and computes a new aggregated presence matrix which is concatenated to the original.
    
    Args:
        ingredients (Ingredients): The Ingredients instance to process.
        aggregation_pattern (str): Regex aggregation_pattern for taxonomic level aggregation (default: "g__").
    
    Returns:
        Ingredients: Updated Ingredients with extra aggregated columns.
    """
    level_names = ["Root", "d__", "p__", "c__", "o__", "f__", "g__", "s__"]
    
    taxa_series = pd.Series(ingredients.taxa)
    taxa_df = taxa_series.str.split("; ", expand=True)
    taxa_df.columns = level_names[:taxa_df.shape[1]]
    
    target_level = next((level for level in level_names if re.fullmatch(aggregation_pattern, level)), None)
    if not target_level or target_level not in taxa_df.columns:
        print(f"No matching level found for aggregation_pattern: {aggregation_pattern}")
        return ingredients
    
    higher_level_taxa = taxa_df[target_level]
    unique_levels = higher_level_taxa.dropna().unique()
    
    num_taxa = len(ingredients.taxa)
    num_unique_levels = len(unique_levels)
    row_indices = []
    col_indices = []
    data = []
    
    level_to_index = {level: idx for idx, level in enumerate(unique_levels)}
    
    for taxon_idx, level in enumerate(higher_level_taxa):
        if pd.isnull(level):
            continue
        level_idx = level_to_index[level]
        row_indices.append(taxon_idx)
        col_indices.append(level_idx)
        data.append(1)
    
    T = csr_matrix((data, (row_indices, col_indices)), shape=(num_taxa, num_unique_levels), dtype=int)
    P = ingredients.presence_matrix
    L = P @ T
    L.data = np.ones_like(L.data) # This binarizes the aggregated presence values.
    
    aggregated_taxa = list(unique_levels)
    new_taxa = ingredients.taxa + aggregated_taxa
    
    new_presence_matrix = hstack([P, L], format='csr')
    
    C = ingredients.coverage_matrix
    M = C @ T
        
    new_coverage_matrix = hstack([C, M], format='csr')
    
    return Ingredients(ingredients.samples, new_taxa, new_presence_matrix, new_coverage_matrix)
