#!/usr/bin/env python3
"""
format.py

Converts a raw sandpiper taxonomic profile (TSV file) into two sparse matrix representations:
    1. A raw Ingredients object (species-level only)
    2. An aggregated Ingredients object (with additional taxonomic levels added)

"""

import numpy as np
import os
import pandas as pd
import pickle
import re
from scipy.sparse import csr_matrix, hstack
import warnings
from typing import Optional

from metacooc.pantry import Ingredients


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

def format_ingredients(tax_profile: str, sample_to_biome=None):
    """
    Process the sandpiper TSV file and return a raw Ingredients instance.
    """
    samples, taxa, sample_to_index, taxon_to_index = build_indices(tax_profile)
    presence_matrix, coverage_matrix = create_sparse_matrices(tax_profile, sample_to_index, taxon_to_index, samples, taxa)
    
    # Generate the raw Ingredients object.
    raw_ingredients = Ingredients(samples, taxa, presence_matrix, coverage_matrix, sample_to_biome)
    
    return raw_ingredients

def format_data(
    tax_profile: str,
    output_dir: str,
    sample_to_biome_file: Optional[str] = None,
    aggregated: bool = False,
    tag: str = "",
):
    """
    Process the sandpiper TSV file and return a raw Ingredients instance.
    
    Args:
        tax_profile (str): Path to the tax profile TSV file.
        output_dir (str): Directory to save the output.
        sample_to_biome_file (str, optional): Path to the biome mapping file.
        aggregated (bool, optional): Whether to generate aggregated ingredients.
        tag (str, optional): Tag for output filenames.
    """
    # Load the biome mapping file if specified
    sample_to_biome = {}
    if sample_to_biome_file:
        if os.path.exists(sample_to_biome_file):
            df = pd.read_csv(sample_to_biome_file, dtype=str, sep="\t")
            # CSV has columns: accession, level_1, level_2
            sample_to_biome = {
                row["accession"]: (row["level_1"], row["level_2"])
                for _, row in df.iterrows()
            }
        else:
            warnings.warn(
                f"Biome mapping file '{sample_to_biome_file}' not found",
                UserWarning,
            )
    
    # Generate the raw Ingredients object
    raw_ingredients = format_ingredients(tax_profile, sample_to_biome)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Save the raw Ingredients object
    output_raw = os.path.join(output_dir, f"ingredients_raw{f'_{tag}' if tag else ""}.pkl")
    with open(output_raw, "wb") as f:
        pickle.dump(raw_ingredients, f)
    print(f"Raw ingredients saved to {output_raw}")
    
    if aggregated:
        # Generate the aggregated Ingredients object
        aggregated_ingredients = raw_ingredients.copy()
        aggregated_ingredients = add_taxa_levels_to_ingredients(aggregated_ingredients)
        output_agg = os.path.join(output_dir, f"ingredients_aggregated{f'_{tag}' if tag else ""}.pkl")
        with open(output_agg, "wb") as f:
            pickle.dump(aggregated_ingredients, f)
        print(f"Aggregated ingredients saved to {output_agg}")


def add_taxa_levels_to_ingredients(ingredients: Ingredients) -> Ingredients:
    """
    Flat direct mapping taxonomy aggregation while preserving 'Root; d__…' as own entry:
    - original taxa (finest resolution) kept as-is
    - for each higher rank (genus → family → order → class → phylum → domain),
      build a single mapping from every raw taxon to its ancestor at that rank
      and aggregate in one matmul.
    """
    # 1) Grab raw taxa and split so first column is 'Root; d__…'
    raw_taxa = list(ingredients.taxa)
    split_regex = r"; (?=(?:p__|c__|o__|f__|g__|s__))"
    taxa_df = pd.Series(raw_taxa).str.split(split_regex, expand=True)
    
    # Assign column prefixes and true codes
    col_prefixes = ["Root; d__", "p__", "c__", "o__", "f__", "g__", "s__"]
    rank_codes   = ["d__",      "p__",  "c__",  "o__",  "f__",  "g__",  "s__"]
    taxa_df.columns = col_prefixes[: taxa_df.shape[1]]
    
    # 2) Build lineage_map from raw strings for full-label lookup
    lineage_map = {}
    for raw in raw_taxa:
        parts = raw.split("; ")
        for j, val in enumerate(parts):
            code = val[:3]
            lineage_map[(code, val)] = lineage_map.get(
                (code, val),
                "; ".join(parts[: j + 1])
            )
    
    # 3) Original matrices
    raw_P = ingredients.presence_matrix  # (n_samples, n_raw)
    raw_C = ingredients.coverage_matrix  # (n_samples, n_raw)
    
    # Initialize blocks with raw data
    P_blocks = [raw_P]
    C_blocks = [raw_C]
    labels   = list(raw_taxa)
    n_raw    = len(raw_taxa)
    
    # 4) Flat mapping per higher rank
    # iterate through each prefix+code pair
    for col_prefix, code in zip(col_prefixes, rank_codes):
        # skip species level (already raw)
        if code == "s__":
            continue
        if col_prefix not in taxa_df.columns:
            print(col_prefix)
            continue
        
        # series of labels at this rank (e.g. 'Root; d__Bacteria', 'p__Proteobacteria', ...)
        series = taxa_df[col_prefix]
        mask   = series.notna()
        if not mask.any():
            continue
            
        rows = series.index[mask].tolist()
        labels_at_rank = sorted(series[mask].unique())
        
        # build sparse mapping T: raw -> labels_at_rank
        cols = [labels_at_rank.index(series[i]) for i in rows]
        data = np.ones(len(rows), dtype=int)
        T = csr_matrix((data, (rows, cols)), shape=(n_raw, len(labels_at_rank)), dtype=int)
        
        # aggregate presence & coverage
        P_cur = raw_P @ T
        P_cur.data[:] = 1
        C_cur = raw_C @ T
        
        # label columns: use lineage_map with correct lookup for root case
        agg_labels = []
        for lab in labels_at_rank:
            # for domain: lab == 'Root; d__X'
            if code == "d__" and lab.startswith("Root; "):
                # strip 'Root; ' to get 'd__X'
                lookup = lab.split("; ", 1)[1]
            else:
                lookup = lab
            full_lin = lineage_map[(code, lookup)]
            agg_labels.append(f"{full_lin} AGGREGATED")
            
        # collect
        P_blocks.append(P_cur)
        C_blocks.append(C_cur)
        labels.extend(agg_labels)
    
    # 5) Stitch everything side by side
    full_P = hstack(P_blocks, format="csr")
    full_C = hstack(C_blocks, format="csr")
    
    return Ingredients(
        samples=ingredients.samples,
        taxa=labels,
        presence_matrix=full_P,
        coverage_matrix=full_C,
        sample_to_biome=ingredients.sample_to_biome,
    )

