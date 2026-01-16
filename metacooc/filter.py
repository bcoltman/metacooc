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
from metacooc.utils import _RANK_PREFIXES
from metacooc.clustering import determine_taxa_context

def filter_by_accessions(ingredients, accession_set):
    # build a boolean mask of samples to keep
    mask = [s in accession_set for s in ingredients.samples]
    not_present = set(accession_set).difference(set(ingredients.samples))
    if not_present:
        print(f"{len(not_present)} accessions were not present in ingredients.samples. This may be due "
              f"to accessions being removed by other parameters (e.g. --min_taxa_count, --min_sample_count)")
    if not any(mask):
        return None
    return ingredients.filtered_samples(mask)

def filter_samples_by_taxa_count(ingredients, min_taxa_count, taxa_count_rank):
    """
    Keep samples that have at least `min_taxa_count` taxa whose *terminal* token
    is at the requested rank.
    
    By default, this means: keep samples with at least `min_taxa_count` species.
    """
    if min_taxa_count is None or min_taxa_count <= 0:
        # Nothing to filter on: return as-is
        return ingredients
    
    rp = taxa_count_rank.strip().lower()
    if rp not in _RANK_PREFIXES:
        raise ValueError(
            f"Unknown rank '{taxa_count_rank}'. Expected one of: {', '.join(_RANK_PREFIXES.keys())}"
        )
    req_pref = _RANK_PREFIXES[rp]
    
    # Ensure taxa lookup caches are built
    ingredients._ensure_taxa_lookups()
    term = ingredients._terminal_rank_prefixes 
    
    # taxa at the requested terminal rank (e.g. species)
    taxa_mask = np.array([p == req_pref for p in term], dtype=bool)
    if not taxa_mask.any():
        warnings.warn(
            f"No taxa with terminal rank '{taxa_count_rank}' found; cannot "
            "filter samples by taxa count at this rank."
        )
        return None
        
    # Count presence of those taxa per sample
    sample_counts = np.array(
        ingredients.presence_matrix[taxa_mask, :].sum(axis=0)
    ).ravel()
    
    # Keep samples with at least min_taxa_count taxa at this rank
    sample_mask = sample_counts >= min_taxa_count
    if not sample_mask.any():
        return None
        
    return ingredients.filtered_samples(sample_mask)

def filter_taxa_by_sample_count(ingredients, min_sample_count):
    # keep taxa present in at least min_sample_count samples
    taxa_counts = np.array((ingredients.presence_matrix > 0).sum(axis=1)).flatten()
    mask = taxa_counts >= min_sample_count
    if not mask.any():
        return None
    return ingredients.filtered_taxa(mask)

def filter_taxa_by_rank(ingredients, filter_rank):
    """
    Keep taxa whose *terminal* token is at the requested rank.
    e.g., filter_rank='species' keeps only s__... features.
    """
    if not filter_rank:
        return ingredients  # no-op
    
    rp = filter_rank.strip().lower()
    if rp not in _RANK_PREFIXES:
        raise ValueError(
            f"Unknown rank '{filter_rank}'. Expected one of: {', '.join(_RANK_PREFIXES.keys())}"
        )
    req_pref = _RANK_PREFIXES[rp]
    
    # ensure the cache (built lazily, not pickled)
    ingredients._ensure_taxa_lookups()
    term = ingredients._terminal_rank_prefixes  # list parallel to ingredients.taxa
    
    mask = [p == req_pref for p in term]
    if not any(mask):
        return None
    return ingredients.filtered_taxa(mask)


class FilteringError(Exception):
    pass

def filter_data_obj(ingredients, 
                    accession_set=None, 
                    min_taxa_count=None, 
                    min_sample_count=None, 
                    filter_rank=None, 
                    taxa_count_rank=None):
    
    filtered = ingredients.copy()
    
    try:        
        if min_taxa_count is not None:
            filtered = filter_samples_by_taxa_count(filtered, min_taxa_count, taxa_count_rank)
            if filtered is None:
                raise FilteringError(f"Warning: Filtering by minimum taxa count of {min_taxa_count} resulted in no samples. {'Could be affected by rank filtering' if filter_rank is not None else ''}")
        
        if min_sample_count is not None:
            filtered = filter_taxa_by_sample_count(filtered, min_sample_count)
            if filtered is None:
                raise FilteringError(f"Warning: Filtering by minimum sample count of {min_sample_count} resulted in no taxa.")
        
        if filter_rank is not None:
            filtered = filter_taxa_by_rank(filtered, filter_rank)
            if filtered is None:
                raise FilteringError(f"Warning: Filtering on {filter_rank} resulted in no taxa.")
        
        if accession_set is not None:
            filtered = filter_by_accessions(filtered, accession_set)
            if filtered is None:
                raise FilteringError("Warning: Filtering by accessions resulted in no samples.")
                
        return filtered, True
        
    except FilteringError as e:
        print(e)
        return None, False

def filter_data(accessions_file, 
                data_dir, 
                output_dir, 
                aggregated=False, 
                min_taxa_count=None, 
                min_sample_count=None, 
                filter_rank=None, 
                taxa_count_rank=None,
                tag=None,
                null_scope=None,
                null_taxa_query=None,
                null_biome_query=None,
                null_metadata_query=None,
                threshold_null=False,
                local_degree=1,
                min_shared_samples_between_taxa=1,
                custom_ingredients=None,
                sandpiper_version=None):
                
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load Ingredients object from disk.
    ingredients = load_ingredients(data_dir, 
                                   aggregated, 
                                   custom_ingredients, 
                                   sandpiper_version)
    
    if null_scope is None:
        
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=None, 
                                              min_taxa_count=min_taxa_count if threshold_null else 0, 
                                              min_sample_count=min_sample_count if threshold_null else 0, 
                                              filter_rank=filter_rank,
                                              taxa_count_rank=taxa_count_rank)
                                              
        if not is_successful:
            return
        
    elif null_scope == "local":
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=None, 
                                              min_taxa_count=min_taxa_count if threshold_null else 0, 
                                              min_sample_count=min_sample_count if threshold_null else 0, 
                                              filter_rank=filter_rank,
                                              taxa_count_rank=taxa_count_rank)
        
        if not is_successful:
            return
                                              
        null_ingredients = determine_taxa_context(null_ingredients,
                                      focal_taxa=null_taxa_query,
                                      degree=local_degree,
                                      min_shared_samples_between_taxa=min_shared_samples_between_taxa)
    
    elif null_scope == "biome" or null_scope == "metadata":
        search_string = null_biome_query if null_scope == "biome" else null_metadata_query
        
        null_matching_accessions = search_data_obj(search_mode=null_scope,
                                              search_string=search_string,
                                              custom_ingredients=ingredients)
        
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=null_matching_accessions, 
                                              min_taxa_count=min_taxa_count if threshold_null else 0, 
                                              min_sample_count=min_sample_count if threshold_null else 0, 
                                              filter_rank=filter_rank,
                                              taxa_count_rank=taxa_count_rank)
        
        if not is_successful:
            return None, None, None, None
    
    elif null_scope == "biome_taxa" or null_scope == "metadata_taxa" :
        
        search_mode = "biome" if null_scope == "biome_taxa" else "metadata"
        search_string = null_biome_query if search_mode == "biome" else null_metadata_query
        
        null_matching_accessions = search_data_obj(search_mode=search_mode,
                                              search_string=search_string,
                                              custom_ingredients=ingredients)
        
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=null_matching_accessions, 
                                              min_taxa_count=min_taxa_count if threshold_null else 0, 
                                              min_sample_count=min_sample_count if threshold_null else 0, 
                                              filter_rank=filter_rank,
                                              taxa_count_rank=taxa_count_rank)
        
        if not is_successful:
            return None, None, None, None
            
        null_ingredients = determine_taxa_context(null_ingredients,
                                                  focal_taxa=null_taxa_query,
                                                  degree=taxa_degree,
                                                  min_shared_samples_between_taxa=min_shared_samples_between_taxa)
    
    intermediate_path = os.path.join(output_dir, f"{tag}ingredients_null.pkl")
    with open(intermediate_path, "wb") as f:
        pickle.dump(null_ingredients, f)
    print(f"Null Ingredients saved to {intermediate_path}")
    
    # If an accessions file is provided, load it and filter.
    if accessions_file:
        with open(accessions_file, "r") as f:
            acc_list = [line.strip() for line in f if line.strip()]
        accession_set = set(acc_list)
        
        filtered, is_successful = filter_data_obj(null_ingredients, 
                                                  accession_set=accession_set)
        
        if is_successful:
            final_path = os.path.join(output_dir, f"{tag}ingredients_filtered.pkl")
            with open(final_path, "wb") as f:
                pickle.dump(filtered, f)
            print(f"Final filtered Ingredients saved to {final_path}")
