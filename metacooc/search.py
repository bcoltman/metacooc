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
import subprocess

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
        if prefix and not prefix in taxon:
            continue
        if search_string.lower() in taxon.lower():
            matching_indices.append(i)
    if not matching_indices:
        return set()
    submatrix = ingredients.presence_matrix[:, matching_indices]
    rows, _ = submatrix.nonzero()
    return {ingredients.samples[r] for r in rows}

def search_data_obj(mode, data_dir, search_string, rank=None, strict=False, column_names=None):
    """
    Object-based search function.
    
    For mode "taxon": data is an Ingredients object.
    For mode "metadata": data is a pre-built metadata index (dict).
    
    Returns:
        set: Matching accession IDs.
    """
    matching_accessions = set()
    if mode.lower() == "taxon":
        ingredients_file = os.path.join(data_dir, "ingredients_raw.pkl")
        if not os.path.exists(ingredients_file):
            raise FileNotFoundError(f"Ingredients file '{ingredients_file}' not found.")
        with open(ingredients_file, "rb") as f:
            ingredients = pickle.load(f)
        
        matching_accessions = search_by_taxon(ingredients, search_string, rank)
    
    elif mode.lower() == "metadata":
        metadata = os.path.join(data_dir, "broad_metadata.tsv")
        
        matching_accessions = search_in_metadata(metadata, search_string, strict, column_names)
    
    else:
        raise ValueError("Invalid mode specified. Must be 'taxon' or 'metadata'.")
        
    return matching_accessions

def search_data(mode, data_dir, output_dir, search_string, rank=None,
                column_names=None, strict=False, tag=None):
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
    
    matching_accessions = search_data_obj(mode, data_dir, search_string, rank, strict, column_names)
    
    output_file = os.path.join(output_dir, f"search_results{tag}.txt")
    with open(output_file, "w") as f:
        for acc in sorted(matching_accessions):
            f.write(f"{acc}\n")
    print(f"Found {len(matching_accessions)} matching accessions. Results saved to {output_file}")

def get_column_indices(metadata_file, column_names, delimiter="\t"):
    """
    Find the indices of given column names in the metadata file.
    
    Args:
        metadata_file (str): Path to the metadata file.
        column_names (list of str): The column names to locate.
        delimiter (str): Column separator (default: tab).
        
    Returns:
        list: The 1-based indices of the columns for AWK (since AWK uses 1-based indexing).
    """
    with open(metadata_file, "r") as f:
        headers = f.readline().strip().split(delimiter)
    
    indices = []
    for column_name in column_names:
        if column_name not in headers:
            raise ValueError(f"Column '{column_name}' not found in metadata file.")
        indices.append(headers.index(column_name) + 1)  # Convert 0-based to 1-based index for AWK
        
    return indices

def grep_metadata(search_string, metadata_file, column_names=None, delimiter="\t"):
    """
    Search for a token in a large metadata file using optimized grep & awk for multiple columns.
    
    Args:
        search_string (str): The token to search for.
        metadata_file (str): Path to the metadata file.
        column_names (list of str or None): List of column names to restrict search to. If None, search whole file.
        delimiter (str): Column separator (default: tab).
        
    Returns:
        set: A set of matching accession numbers.
    """
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file '{metadata_file}' not found.")
        
    search_string = search_string.strip().lower()
    
    # Find column indices if column names are specified
    column_indices = None
    if column_names:
        column_indices = get_column_indices(metadata_file, column_names, delimiter)
        
    # Extract accession + specified columns
    if column_indices:
        column_indices_str = ",".join(f"${idx}" for idx in column_indices)  # e.g., "$2,$3" for awk
        awk_command = f"awk -F'{delimiter}' '{{print $1, {column_indices_str}}}' {metadata_file}"
        
        # Grep should search only in the extracted columns (not the accession column)
        grep_command = f"awk -F' ' '"
        grep_command += " || ".join(f"$i ~ /{search_string}/" for i in range(2, len(column_indices) + 2))
        grep_command += " {print $1}'"
        command = f"{awk_command} | {grep_command}"
        
    else:
        # If no specific columns, grep the whole file and extract the first column
        grep_command = f"grep -i '{search_string}' {metadata_file} | cut -f1"
        command = grep_command
        
    # Run the command
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    # Return a set of unique accession numbers
    return set(result.stdout.splitlines())

def search_in_metadata(metadata, search_string, strict=False, column_names=None):
    """
    Search for a token in metadata files using grep, with optional column restriction to multiple columns.
    
    Args:
        search_string (str): The token to search for.
        metadata (str): Metadata file.
        strict (bool): Whether to use the strict metadata file.
        column_names (list of str or None): List of column names to restrict search to.
        
    Returns:
        set: A set of matching accession numbers.
    """
    
    
    if strict:
        column_names = ["acc", "organism", "env_biome_sam", "env_feature_sam", "env_material_sam", "biosamplemodel_sam"]
        
    # Perform the search
    return grep_metadata(search_string, metadata, column_names)