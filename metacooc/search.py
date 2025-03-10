#!/usr/bin/env python3
"""
search.py

Combined search functionality for metacooc.

Exposes:
  1. search_data(mode, data_dir, output_dir, search_string, rank=None, index_type='broad',
                 column=None, strict=False, index_file=None)
     - File-based interface: loads files from disk, writes the matching accessions to
       output_dir/search_results.txt, and returns the set of matching accessions.

  2. search_data_obj(mode, data, search_string, rank=None)
     - Object-based interface:
         * For mode "taxon", data should be an in-memory Ingredients object.
         * For mode "metadata", data should be a pre-built metadata index (a dict).
       It returns the set of matching accessions.
       
Modes:
  - Taxon:
      Loads an Ingredients pickle (e.g., "ingredients_aggregated.pkl") from data_dir and
      searches its taxa list for the given search_string (optionally restricting by rank).
  
  - Metadata:
      Either loads a pre-built metadata index (if index_file is provided) or loads an index
      from pre-generated files (e.g., "broad_metadata_index.pkl" or "broad_strict_metadata_index.pkl"
      for broad mode, or "{column}_metadata_index.pkl" for exact mode) from data_dir. The search
      is then performed on that index.
"""

import os
import pickle
import pandas as pd

def search_by_taxon(ingredients, search_string, rank=None):
    """
    Search the Ingredients object's taxa list.
    
    Returns:
        set: Accession IDs from samples where a matching taxon is present.
    """
    rank_prefixes = {"genus": "g__", "family": "f__", "order": "o__",
                     "class": "c__", "phylum": "p__", "species": "s__", "domain": "d__"}
    prefix = rank_prefixes.get(rank.lower()) if rank else None
    matching_indices = []
    for i, taxon in enumerate(ingredients.taxa):
        if prefix and not taxon.startswith(prefix):
            continue
        if search_string.lower() in taxon.lower():
            matching_indices.append(i)
    if not matching_indices:
        return set()
    submatrix = ingredients.presence_matrix[:, matching_indices]
    rows, _ = submatrix.nonzero()
    return {ingredients.samples[r] for r in rows}

def search_data_obj(mode, data, search_string, rank=None):
    """
    Object-based search function.
    
    For mode "taxon": data is an Ingredients object.
    For mode "metadata": data is a pre-built metadata index (dict).
    
    Returns:
        set: Matching accession IDs.
    """
    if mode.lower() == "taxon":
        return search_by_taxon(data, search_string, rank)
    elif mode.lower() == "metadata":
        # Here, data is assumed to be a dict mapping tokens to sets of accessions.
        return data.get(search_string.strip().lower(), set())
    else:
        raise ValueError("Invalid mode. Must be 'taxon' or 'metadata'.")

def search_data(mode, data_dir, output_dir, search_string, rank=None, 
                index_type="broad", column=None, strict=False, index_file=None, tag=None):
    """
    File-based search function for metacooc.
    
    For 'taxon' mode, loads an Ingredients object (expected as "ingredients_aggregated.pkl")
    from data_dir and uses search_by_taxon.
    
    For 'metadata' mode, either loads a pre-built index (if index_file is provided) or uses
    default index filenames from data_dir, then performs a lookup.
    
    Writes matching accessions (one per line) to output_dir/search_results.txt.
    
    Returns:
        set: Matching accession IDs.
    """
    matching_accessions = set()
    if mode.lower() == "taxon":
        ingredients_file = os.path.join(data_dir, "ingredients_aggregated.pkl")
        if not os.path.exists(ingredients_file):
            raise FileNotFoundError(f"Ingredients file '{ingredients_file}' not found.")
        with open(ingredients_file, "rb") as f:
            ingredients = pickle.load(f)
        matching_accessions = search_by_taxon(ingredients, search_string, rank)
    elif mode.lower() == "metadata":
        if index_file:
            with open(index_file, "rb") as f:
                metadata_index = pickle.load(f)
            print(f"Loaded metadata index from {index_file}")
        else:
            if index_type.lower() == "broad":
                fname = "broad_metadata_index.pkl" if not strict else "broad_strict_metadata_index.pkl"
            else:
                if not column:
                    raise ValueError("Exact index requires a column name.")
                fname = f"{column}_metadata_index.pkl"
            index_path = os.path.join(data_dir, fname)
            if not os.path.exists(index_path):
                raise FileNotFoundError(f"Metadata index file '{index_path}' not found.")
            with open(index_path, "rb") as f:
                metadata_index = pickle.load(f)
            print(f"Used metadata index from {index_path}")
        matching_accessions = metadata_index.get(search_string.strip().lower(), set())
    else:
        raise ValueError("Invalid mode specified. Must be 'taxon' or 'metadata'.")
    
    output_file = os.path.join(output_dir, f"search_results{tag}.txt")
    with open(output_file, "w") as f:
        for acc in sorted(matching_accessions):
            f.write(f"{acc}\n")
    print(f"Found {len(matching_accessions)} matching accessions. Results saved to {output_file}")

