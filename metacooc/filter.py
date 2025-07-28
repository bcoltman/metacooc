#!/usr/bin/env python3
"""
filter.py

This module provides two interfaces for filtering:

1. File-based filtering (filter_data):
   Loads an Ingredients object from disk, applies filters, and saves output.
2. Object-based filtering (filter_data_obj):
   Accepts an Ingredients object and returns a filtered Ingredients object.

Usage (file-based):
    metacooc filter --accessions_file path/to/accessions.txt --data_dir /path/to/data --output_dir /path/to/out --min_taxa_count 5 --min_sample_count 10
"""

import os
import pickle
import numpy as np
from metacooc.pantry import *

def filter_by_accessions(ingredients, accession_set):
    # build a boolean mask of samples to keep
    mask = [s in accession_set for s in ingredients.samples]
    if not any(mask):
        print("Warning: No samples match the given accessions.")
        return None
    return ingredients.filtered_samples(mask)

def filter_samples_by_taxa_count(ingredients, min_taxa_count):
    # keep samples that have at least min_taxa_count taxa
    sample_counts = np.array(ingredients.presence_matrix.sum(axis=1)).flatten()
    mask = sample_counts >= min_taxa_count
    if not mask.any():
        print("Warning: No samples meet the taxa count threshold.")
        return None
    return ingredients.filtered_samples(mask)

def filter_taxa_by_sample_count(ingredients, min_sample_count):
    # keep taxa present in at least min_sample_count samples
    taxa_counts = np.array((ingredients.presence_matrix > 0).sum(axis=0)).flatten()
    mask = taxa_counts >= min_sample_count
    if not mask.any():
        print("Warning: No taxa meet the sample count threshold.")
        return None
    return ingredients.filtered_taxa(mask)

def filter_taxa_by_rank(ingredients, filter_rank):
    rank_prefixes = {
        "domain": "d__", "phylum": "p__", "class": "c__",
        "order": "o__",   "family": "f__", "genus": "g__", "species": "s__",
    }
    prefix = rank_prefixes.get(filter_rank.lower()) if filter_rank else None
    mask = [bool(prefix and prefix in t) for t in ingredients.taxa]
    if not any(mask):
        print(f"Warning: Filtering on {filter_rank} resulted in no taxa.")
        return None
    return ingredients.filtered_taxa(mask)

def filter_data_obj(ingredients, accession_set=None, min_taxa_count=None, min_sample_count=None, filter_rank=None):
    
    filtered = ingredients.copy()
    
    if filter_rank is not None:
        filtered = filter_taxa_by_rank(filtered, filter_rank)
        if filtered is None:
            raise ValueError(f"Filtering on {filter_rank} resulted in no taxa.")
    
    if min_taxa_count is not None:
        filtered = filter_samples_by_taxa_count(filtered, min_taxa_count)
        if filtered is None:
            raise ValueError(f"Filtering by mnimum taxa count of {min_taxa_count} resulted in no samples. {'Could be affected by rank filtering' if filter_rank is not None else ''}")
    
    if min_sample_count is not None:
        filtered = filter_taxa_by_sample_count(filtered, min_sample_count)
        if filtered is None:
            raise ValueError(f"Filtering by minimum sample count of {min_sample_count} resulted in no taxa.")
    
    if accession_set is not None:
        filtered = filter_by_accessions(filtered, accession_set)
        if filtered is None:
            raise ValueError("Filtering by accessions resulted in no samples.")
    
    
    return filtered

def filter_data(accessions_file, 
                data_dir, 
                output_dir, 
                aggregated=False, 
                min_taxa_count=None, 
                min_sample_count=None, 
                filter_rank=None, 
                tag=None, 
                custom_ingredients=None,
                sandpiper_version=None):
    
    os.makedirs(output_dir, exist_ok=True)

    # Load Ingredients object from disk.
    ingredients = load_ingredients(data_dir, 
                                   aggregated, 
                                   custom_ingredients, 
                                   sandpiper_version)
    
    # Apply count-based filters.
    filtered = filter_data_obj(ingredients, 
                               accession_set=None, 
                               min_taxa_count=min_taxa_count, 
                               min_sample_count=min_sample_count, 
                               filter_rank=filter_rank)
    intermediate_path = os.path.join(output_dir, f"ingredients_counts_filtered{tag}.pkl")
    with open(intermediate_path, "wb") as f:
        pickle.dump(filtered, f)
    print(f"Intermediate filtered Ingredients saved to {intermediate_path}")
    
    # If an accessions file is provided, load it and filter.
    if accessions_file:
        with open(accessions_file, "r") as f:
            acc_list = [line.strip() for line in f if line.strip()]
        accession_set = set(acc_list)
        filtered = filter_data_obj(filtered, accession_set, None, None, None)
    
    final_path = os.path.join(output_dir, f"ingredients_all_filtered{tag}.pkl")
    with open(final_path, "wb") as f:
        pickle.dump(filtered, f)
    print(f"Final filtered Ingredients saved to {final_path}")
