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
    indices = [i for i, acc in enumerate(ingredients.samples) if acc in accession_set]
    if not indices:
        print("Warning: No samples match the given accessions.")
        return None
    new_samples = [ingredients.samples[i] for i in indices]
    new_presence = ingredients.presence_matrix[indices, :]
    new_coverage = ingredients.coverage_matrix[indices, :]
    return Ingredients(new_samples, ingredients.taxa, new_presence, new_coverage)

def filter_taxa_by_sample_count(ingredients, min_sample_count):
    taxa_counts = np.array((ingredients.presence_matrix > 0).sum(axis=0)).flatten()
    indices = np.where(taxa_counts >= min_sample_count)[0]
    if indices.size == 0:
        print("Warning: No taxa meet the sample count threshold.")
        return None
    new_taxa = [ingredients.taxa[i] for i in indices]
    new_presence = ingredients.presence_matrix[:, indices]
    new_coverage = ingredients.coverage_matrix[:, indices]
    return Ingredients(ingredients.samples, new_taxa, new_presence, new_coverage)

def filter_samples_by_taxa_count(ingredients, min_taxa_count):
    sample_counts = np.array(ingredients.presence_matrix.sum(axis=1)).flatten()
    indices = np.where(sample_counts >= min_taxa_count)[0]
    if indices.size == 0:
        print("Warning: No samples meet the taxa count threshold.")
        return None
    new_samples = [ingredients.samples[i] for i in indices]
    new_presence = ingredients.presence_matrix[indices, :]
    new_coverage = ingredients.coverage_matrix[indices, :]
    return Ingredients(new_samples, ingredients.taxa, new_presence, new_coverage)

def filter_taxa_by_rank(ingredients, filter_rank):
    rank_prefixes = {"domain": "d__", 
                     "phylum": "p__", 
                     "class": "c__", 
                     "order": "o__",
                     "family": "f__",  
                     "genus": "g__", 
                     "species": "s__"}
    prefix = rank_prefixes.get(filter_rank.lower()) if filter_rank else None
    
    indices = [i for i, taxa in enumerate(ingredients.taxa) if prefix in taxa]
    
    # indices = []
    # for i, taxon in enumerate(ingredients.taxa):
        # if prefix and prefix in taxon:
            # indices.append(i)
    # if not indices:
        # return set()
    
    if not indices:
        print(f"Warning: Filtering on {filter_rank} resulted in no taxa.")
        return None
        
    new_taxa = [ingredients.taxa[i] for i in indices]
    new_presence = ingredients.presence_matrix[:, indices]
    new_coverage = ingredients.coverage_matrix[:, indices]
    return Ingredients(ingredients.samples, new_taxa, new_presence, new_coverage)

def filter_data_obj(ingredients, accession_set=None, min_taxa_count=None, min_sample_count=None, filter_rank=None):
    
    filtered = ingredients.copy()
    
    if min_taxa_count is not None:
        filtered = filter_samples_by_taxa_count(filtered, min_taxa_count)
        if filtered is None:
            raise ValueError("Filtering by taxa count resulted in no samples.")
    
    if min_sample_count is not None:
        filtered = filter_taxa_by_sample_count(filtered, min_sample_count)
        if filtered is None:
            raise ValueError("Filtering by sample count resulted in no taxa.")
    
    if accession_set is not None:
        filtered = filter_by_accessions(filtered, accession_set)
        if filtered is None:
            raise ValueError("Filtering by accessions resulted in no samples.")
    
    if filter_rank is not None:
        filtered = filter_taxa_by_rank(filtered, filter_rank)
        if filtered is None:
            raise ValueError(f"Filtering on {filter_rank} resulted in no taxa.")
    
    return filtered

def filter_data(accessions_file, 
                data_dir, 
                output_dir, 
                aggregated=False, 
                min_taxa_count=None, 
                min_sample_count=None, 
                filter_rank=None, 
                tag=None, 
                custom_ingredients=None):
    
    # Load Ingredients object from disk.
    ingredients = load_ingredients(data_dir, aggregated, custom_ingredients)
    
    # Apply count-based filters.
    filtered = filter_data_obj(ingredients, None, min_taxa_count, min_sample_count, filter_rank)
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
