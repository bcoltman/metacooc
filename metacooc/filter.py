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
import numpy as np
from metacooc.pantry import *
from metacooc.utils import _RANK_PREFIXES
from metacooc.clustering import determine_taxa_context
from metacooc.search import search_data_obj
import warnings

def _sample_mask_from_accessions(ingredients, accession_set):
    if accession_set is None:
        return None
    
    accession_set = accession_set if isinstance(accession_set, set) else set(accession_set)
    sample_arr = np.asarray(ingredients.samples, dtype=object)
    
    mask = np.fromiter((s in accession_set for s in sample_arr),
                       dtype=bool,
                       count=sample_arr.size)
    
    sample_set = set(sample_arr.tolist())
    not_present = accession_set.difference(sample_set)
    
    if not_present:
        print(
            f"{len(not_present)} accessions were not present in ingredients.samples. "
            "This may be due to accessions being removed by other parameters "
            "(e.g. --min_taxa_count, --min_sample_count)"
        )
        
    return mask if mask.any() else None

def filter_by_accessions(ingredients, accession_set):
    mask = _sample_mask_from_accessions(ingredients, accession_set)
    if mask is None:
        return None
    return ingredients.filtered_samples(mask)


def _sample_mask_by_taxa_count(ingredients, min_taxa_count, taxa_count_rank):
    """
    Boolean sample mask: keep samples with at least `min_taxa_count` taxa whose
    terminal rank is `taxa_count_rank`.
    """
    if min_taxa_count is None or min_taxa_count <= 0:
        return np.ones(len(ingredients.samples), dtype=bool)
    
    rp = taxa_count_rank.strip().lower()
    if rp not in _RANK_PREFIXES:
        raise ValueError(
            f"Unknown rank '{taxa_count_rank}'. Expected one of: {', '.join(_RANK_PREFIXES.keys())}"
        )
    req_pref = _RANK_PREFIXES[rp]
    
    ingredients._ensure_taxa_lookups()
    term = ingredients._terminal_rank_prefixes
    
    taxa_mask = np.fromiter((p == req_pref for p in term), dtype=bool, count=len(term))
    if not taxa_mask.any():
        warnings.warn(
            f"No taxa with terminal rank '{taxa_count_rank}' found; cannot "
            "filter samples by taxa count at this rank."
        )
        return None
    
    P = ingredients.presence_matrix
    if taxa_mask.all():
        sample_counts = np.asarray(P.getnnz(axis=0)).ravel()
    else:
        row_idx = np.flatnonzero(taxa_mask)
        sample_counts = np.asarray(P[row_idx, :].getnnz(axis=0)).ravel()
        
    sample_mask = sample_counts >= min_taxa_count
    if not sample_mask.any():
        return None
    return sample_mask


def filter_samples_by_taxa_count(ingredients, min_taxa_count, taxa_count_rank):
    mask = _sample_mask_by_taxa_count(ingredients, min_taxa_count, taxa_count_rank)
    if mask is None:
        return None
    return ingredients.filtered_samples(mask)


def _taxa_mask_by_sample_count(ingredients, min_sample_count):
    """
    Boolean taxa mask: keep taxa present in at least `min_sample_count` samples.
    
    Uses cached total_counts instead of recomputing from sparse comparisons.
    """
    if min_sample_count is None or min_sample_count <= 0:
        return np.ones(len(ingredients.taxa), dtype=bool)
    
    mask = ingredients.total_counts >= int(min_sample_count)
    if not mask.any():
        return None
    return mask


def filter_taxa_by_sample_count(ingredients, min_sample_count):
    mask = _taxa_mask_by_sample_count(ingredients, min_sample_count)
    if mask is None:
        return None
    return ingredients.filtered_taxa(mask)


def _taxa_mask_by_rank(ingredients, filter_rank):
    """
    Boolean taxa mask: keep taxa whose terminal token is at the requested rank.
    """
    if not filter_rank:
        return np.ones(len(ingredients.taxa), dtype=bool)
    
    rp = filter_rank.strip().lower()
    if rp not in _RANK_PREFIXES:
        raise ValueError(
            f"Unknown rank '{filter_rank}'. Expected one of: {', '.join(_RANK_PREFIXES.keys())}"
        )
    req_pref = _RANK_PREFIXES[rp]
    
    ingredients._ensure_taxa_lookups()
    term = ingredients._terminal_rank_prefixes
    
    mask = np.fromiter((p == req_pref for p in term), dtype=bool, count=len(term))
    if not mask.any():
        return None
    return mask


def filter_taxa_by_rank(ingredients, filter_rank):
    mask = _taxa_mask_by_rank(ingredients, filter_rank)
    if mask is None:
        return None
    return ingredients.filtered_taxa(mask)


class FilteringError(Exception):
    pass


def filter_data_obj(
    ingredients,
    accession_set=None,
    min_taxa_count=None,
    min_sample_count=None,
    filter_rank=None,
    taxa_count_rank=None,
):
    """
    Apply filters to a single shallow working copy, in place.
    
    This avoids the previous pattern of:
      copy -> filtered_samples() -> copy -> filtered_taxa() -> copy -> ...
    """
    filtered = ingredients.copy_shallow()
    
    try:
        if min_taxa_count is not None and min_taxa_count > 0:
            sample_mask = _sample_mask_by_taxa_count(filtered, min_taxa_count, taxa_count_rank)
            if sample_mask is None:
                raise FilteringError(
                    f"Warning: Filtering by minimum taxa count of {min_taxa_count} resulted in no samples. "
                    f"{'Could be affected by rank filtering' if filter_rank is not None else ''}"
                )
            filtered.filter_samples(sample_mask)
        
        if min_sample_count is not None and min_sample_count > 0:
            taxa_mask = _taxa_mask_by_sample_count(filtered, min_sample_count)
            if taxa_mask is None:
                raise FilteringError(
                    f"Warning: Filtering by minimum sample count of {min_sample_count} resulted in no taxa."
                )
            filtered.filter_taxa(taxa_mask)
        
        if filter_rank is not None:
            taxa_mask = _taxa_mask_by_rank(filtered, filter_rank)
            if taxa_mask is None:
                raise FilteringError(f"Warning: Filtering on {filter_rank} resulted in no taxa.")
            filtered.filter_taxa(taxa_mask)
        
        if accession_set is not None:
            sample_mask = _sample_mask_from_accessions(filtered, accession_set)
            if sample_mask is None:
                raise FilteringError("Warning: Filtering by accessions resulted in no samples.")
            filtered.filter_samples(sample_mask)
        
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
                remove_null_threshold=False,
                taxa_degree=1,
                min_shared_samples_between_taxa=1,
                custom_ingredients=None,
                data_version=None,
                metadata_file=None,
                min_coverage=None,
                min_coverage_by_rank=None,
                min_relative_abundance=None,
                min_relative_abundance_by_rank=None):
                
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load Ingredients object from disk.
    ingredients = load_ingredients(data_dir, 
                                   aggregated, 
                                   custom_ingredients, 
                                   data_version)
    ingredients = threshold_ingredients_presence(
        ingredients,
        min_coverage=min_coverage,
        min_coverage_by_rank=min_coverage_by_rank,
        min_relative_abundance=min_relative_abundance,
        min_relative_abundance_by_rank=min_relative_abundance_by_rank,
    )
    
    if null_scope is None:
        
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=None, 
                                              min_taxa_count=0 if remove_null_threshold else min_taxa_count, 
                                              min_sample_count=0 if remove_null_threshold else min_sample_count, 
                                              filter_rank=filter_rank,
                                              taxa_count_rank=taxa_count_rank)
                                              
        if not is_successful:
            return
        
    elif null_scope == "taxa":
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=None,
                                              min_taxa_count=0 if remove_null_threshold else min_taxa_count, 
                                              min_sample_count=0 if remove_null_threshold else min_sample_count,
                                              filter_rank=filter_rank,
                                              taxa_count_rank=taxa_count_rank)
        
        if not is_successful:
            return
                                              
        null_ingredients = determine_taxa_context(null_ingredients,
                                      focal_taxa=null_taxa_query,
                                      degree=taxa_degree,
                                      min_shared_samples_between_taxa=min_shared_samples_between_taxa)
    
    elif null_scope == "biome" or null_scope == "metadata":
        search_string = null_biome_query if null_scope == "biome" else null_metadata_query
        
        null_matching_accessions = search_data_obj(search_mode=null_scope,
                                              search_string=search_string,
                                              data_dir=data_dir,
                                              custom_ingredients=ingredients,
                                              data_version=data_version,
                                              aggregated=aggregated,
                                              metadata_file=metadata_file)
        
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=null_matching_accessions, 
                                              min_taxa_count=0 if remove_null_threshold else min_taxa_count, 
                                              min_sample_count=0 if remove_null_threshold else min_sample_count,
                                              filter_rank=filter_rank,
                                              taxa_count_rank=taxa_count_rank)
        
        if not is_successful:
            return None, None, None, None
    
    elif null_scope == "biome_taxa" or null_scope == "metadata_taxa" :
        
        search_mode = "biome" if null_scope == "biome_taxa" else "metadata"
        search_string = null_biome_query if search_mode == "biome" else null_metadata_query
        
        null_matching_accessions = search_data_obj(search_mode=search_mode,
                                              search_string=search_string,
                                              data_dir=data_dir,
                                              custom_ingredients=ingredients,
                                              data_version=data_version,
                                              aggregated=aggregated,
                                              metadata_file=metadata_file)
        
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=null_matching_accessions, 
                                              min_taxa_count=0 if remove_null_threshold else min_taxa_count, 
                                              min_sample_count=0 if remove_null_threshold else min_sample_count, 
                                              filter_rank=filter_rank,
                                              taxa_count_rank=taxa_count_rank)
        
        if not is_successful:
            return None, None, None, None
            
        null_ingredients = determine_taxa_context(null_ingredients,
                                                  focal_taxa=null_taxa_query,
                                                  degree=taxa_degree,
                                                  min_shared_samples_between_taxa=min_shared_samples_between_taxa)
    
    intermediate_path = os.path.join(output_dir, f"{tag}ingredients_null")
    save_ingredients_directory(null_ingredients, intermediate_path, aggregated=aggregated)
    print(f"Null Ingredients saved to {intermediate_path}")
    
    # If an accessions file is provided, load it and filter.
    if accessions_file:
        with open(accessions_file, "r") as f:
            acc_list = [line.strip() for line in f if line.strip()]
        accession_set = set(acc_list)
        
        filtered, is_successful = filter_data_obj(null_ingredients, 
                                                  accession_set=accession_set)
        
        if is_successful:
            final_path = os.path.join(output_dir, f"{tag}ingredients_filtered")
            save_ingredients_directory(filtered, final_path, aggregated=aggregated)
            print(f"Final filtered Ingredients saved to {final_path}")
    else:
        print(
            "No accession file provided. Count/rank filters were applied only when constructing "
            "the null Ingredients object, which was saved as "
            f"{intermediate_path}. No separate filtered Ingredients file was written."
        )
        
