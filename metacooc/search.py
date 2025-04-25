#!/usr/bin/env python3
"""
search.py

Combined search functionality for metacooc.

Exposes:
  1. search_data(mode, data_dir, output_dir, search_string, rank=None, column_names=None,
                 strict=False, tag="", inverse=False)
     - File-based interface: loads files from disk, writes the matching accessions to
       output_dir/search_results.txt, and returns the set of matching accessions.

  2. search_data_obj(mode, data_dir, search_string, rank=None, strict=False,
                     column_names=None, inverse=False)
     - Object-based interface:
         * For mode "taxon", data is an Ingredients object.
         * For mode "metadata", data is a metadata file.
       It returns the set of matching accessions.
       
Modes:
  - Taxon:
      Loads an Ingredients pickle (e.g., "ingredients_raw.pkl") from data_dir and
      searches its taxa list for the given search_string (optionally restricting by rank).
  
  - Metadata:
      Searches a metadata file (e.g., "sra_metadata.tsv") for the search token,
      optionally restricting the search to specified columns. Inverse searches are done
      by activating grep’s -v flag (or by inverting awk conditions).
"""

import os
import subprocess

from metacooc.pantry import load_ingredients 
from metacooc._data_config import FILENAMES

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
        if prefix and prefix not in taxon:
            continue
        if search_string.lower() in taxon.lower():
            matching_indices.append(i)
    if not matching_indices:
        return set()
    submatrix = ingredients.presence_matrix[:, matching_indices]
    rows, _ = submatrix.nonzero()
    return {ingredients.samples[r] for r in rows}

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

def grep_metadata(search_string, metadata_file, column_names=None, delimiter="\t", inverse=False):
    """
    Search for a token in a large metadata file using optimized grep & awk for multiple columns,
    with support for inverse search via grep's -v flag or by inverting awk conditions.
    
    Args:
        search_string (str): The token to search for.
        metadata_file (str): Path to the metadata file.
        column_names (list of str or None): List of column names to restrict search to.
            If None, the whole file is searched.
        delimiter (str): Column separator (default: tab).
        inverse (bool): If True, perform an inverse search (lines that do NOT match).
        
    Returns:
        set: A set of matching accession numbers.
    """
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file '{metadata_file}' not found.")
        
    search_string = search_string.strip().lower()
    
    # If no specific columns are provided, use grep directly with the -v flag if inverse is True.
    if not column_names:
        flag = "-iv" if inverse else "-i"
        command = f"grep {flag} '{search_string}' {metadata_file} | cut -f1"
    else:
        # Find column indices for the specified column names.
        column_indices = get_column_indices(metadata_file, column_names, delimiter)
        column_indices_str = ",".join(f"${idx}" for idx in column_indices)  # e.g., "$2,$3" for awk
        awk_command = f"awk -F'{delimiter}' '{{print $1, {column_indices_str}}}' {metadata_file}"
        
        if inverse:
            # Inverse: build condition so that all specified columns do NOT match search_string.
            conditions = []
            for i in range(2, len(column_indices) + 2):
                conditions.append(f"$${i} !~ /{search_string}/")
            condition_str = " && ".join(conditions)
        else:
            # Normal: any column matches search_string.
            conditions = []
            for i in range(2, len(column_indices) + 2):
                conditions.append(f"$${i} ~ /{search_string}/")
            condition_str = " || ".join(conditions)
        grep_command = f"awk -F' ' '{{if({condition_str}) print $1}}'"
        command = f"{awk_command} | {grep_command}"
    
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return set(result.stdout.splitlines())

def search_in_metadata(metadata, search_string, strict=False, column_names=None, inverse=False):
    """
    Search for a token in a metadata file using grep, with optional column restriction and
    inverse search capability.
    
    Args:
        search_string (str): The token to search for.
        metadata (str): Metadata file.
        strict (bool): If True, use a predefined list of column names.
        column_names (list of str or None): List of column names to restrict search to.
        inverse (bool): If True, perform an inverse search.
        
    Returns:
        set: A set of matching accession numbers.
    """
    if strict:
        column_names = ["acc", "organism", "env_biome_sam", "env_feature_sam",
                        "env_material_sam", "biosamplemodel_sam"]
    return grep_metadata(search_string, metadata, column_names, inverse=inverse)

def search_data_obj(mode, data_dir, search_string, rank=None, strict=False, column_names=None, inverse=False):
    """
    Object-based search function.
    
    For mode "taxon": data is an Ingredients object.
    For mode "metadata": data is a metadata file.
    
    Returns:
        set: Matching accession IDs.
    """
    matching_accessions = set()
    if mode.lower() == "taxon":
        ingredients = load_ingredients(data_dir)
        
        matching_accessions = search_by_taxon(ingredients, search_string, rank)
        if inverse:
            # For inverse search, return samples that do NOT match.
            matching_accessions = set(ingredients.samples) - matching_accessions
    
    elif mode.lower() == "metadata":
        metadata_file = os.path.join(data_dir, FILENAMES["sra_metadata"])
        if not os.path.exists(metadata_file):
            raise FileNotFoundError(f"Expected metadata file '{metadata_file}' not found.")
        matching_accessions = search_in_metadata(metadata_file, search_string, strict, column_names, inverse)
    # metadata_file = os.path.join(data_dir, "sra_metadata.tsv")
        # matching_accessions = search_in_metadata(metadata_file, search_string, strict, column_names, inverse)
    
    else:
        raise ValueError("Invalid mode specified. Must be 'taxon' or 'metadata'.")
        
    return matching_accessions

def search_data(mode, data_dir, output_dir, search_string, rank=None,
                column_names=None, strict=False, tag="", inverse=False):
    """
    File-based search function for metacooc.
    
    For 'taxon' mode, loads an Ingredients object (expected as "ingredients_raw.pkl")
    from data_dir and uses search_by_taxon.
    
    For 'metadata' mode, performs a grep-based lookup on the metadata file.
    
    Writes matching accessions (one per line) to output_dir/search_results.txt.
    
    Returns:
        set: Matching accession IDs.
    """
    matching_accessions = search_data_obj(mode, data_dir, search_string, rank, strict, column_names, inverse)
    
    output_file = os.path.join(output_dir, f"search_results{tag if tag else ''}.txt")
    with open(output_file, "w") as f:
        for acc in sorted(matching_accessions):
            f.write(f"{acc}\n")
    print(f"Found {len(matching_accessions)} matching accessions. Results saved to {output_file}")

