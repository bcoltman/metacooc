#!/usr/bin/env python3
"""
search.py

Combined search functionality for metacooc.

Exposes:
  1. search_data(mode, data_dir, output_dir, search_string, ranks_for_search_inclusion=None, column_names=None,
                 strict=False, tag="", inverse=False)
     - File-based interface: loads files from disk, writes the matching accessions to
       output_dir/search_results.txt, and returns the set of matching accessions.

  2. search_data_obj(mode, data_dir, search_string, ranks_for_search_inclusion=None, strict=False,
                     column_names=None, inverse=False)
     - Object-based interface:
         * For mode "taxon", data is an Ingredients object.
         * For mode "metadata", data is a metadata file.
       It returns the set of matching accessions.
       
Modes:
  - Taxon:
      Loads an Ingredients pickle (e.g., "ingredients_raw.pkl") from data_dir and
      searches its taxa list for the given search_string (optionally restricting by ranks_for_search_inclusion).
  
  - Metadata:
      Searches a metadata file (e.g., "sra_metadata.tsv") for the search token,
      optionally restricting the search to specified columns. Inverse searches are done
      by activating grep’s -v flag (or by inverting awk conditions).
"""

import os
import subprocess
from typing import List, Set

from metacooc.pantry import load_ingredients 
from metacooc._data_config import *

def search_by_taxon(ingredients, search_string, ranks_for_search_inclusion=None):
    """
    Search the Ingredients object's taxa list.
    
    Returns:
        set: Accession IDs from samples where a matching taxon is present.
    """
    rank_prefixes = {"domain": "d__", 
                     "phylum": "p__", 
                     "class": "c__", 
                     "order": "o__",
                     "family": "f__",  
                     "genus": "g__", 
                     "species": "s__"}
    
    prefix = rank_prefixes.get(ranks_for_search_inclusion.lower()) if ranks_for_search_inclusion else None
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

def search_by_biome(ingredients, biome_names):
    """
    Return all sample IDs whose sample_to_biome matches one of biome_names.
    
    Validates that each requested biome is in ingredients.biomes_order.
    """
    # ensure we have precomputed biomes
    if not hasattr(ingredients, "biomes_order"):
        ingredients._allocate_biomes()
    available = set(ingredients.biomes_order)
    
    # normalize input to list
    if isinstance(biome_names, str):
        requested = [b.strip() for b in biome_names.split(",")]
    else:
        requested = list(biome_names)
    
    # validate
    bad = [b for b in requested if b not in available]
    if bad:
        raise ValueError(
            f"Unknown biome(s): {bad}. "
            f"Available biomes are: {sorted(available)}"
        )
    
    # gather matching samples
    return {
        sample
        for sample, biome in ingredients.sample_to_biome.items()
        if biome in requested
    }

def _parse_query(q: str) -> List[List[str]]:
    groups = []
    # top‑level OR
    for or_part in q.split("|"):
        or_part = or_part.strip()
        if not or_part:
            continue
        # within each, split on AND (‘+’)
        terms = [t.strip() for t in or_part.split("+") if t.strip()]
        if terms:
            groups.append(terms)
    return groups

def search_data_obj(
    mode: str,
    search_string: str,
    data_dir: str = None,
    ranks_for_search_inclusion=None,
    strict=False,
    column_names=None,
    inverse=False,
    custom_ingredients=None,
    sandpiper_version=None
) -> Set:
    mode = mode.lower()
    
    # 1) metadata: raw regex search
    if mode == "metadata":
        version = sandpiper_version or LATEST_VERSION
        filenames, _ = get_file_info(version)
        if not data_dir:
            raise ValueError(
                "data_dir must be provided if searching metadata"
            )
        metadata_file = os.path.join(data_dir, filenames["sra_metadata"])
        if not os.path.exists(metadata_file):
            raise FileNotFoundError(f"Missing '{metadata_file}'")
        return search_in_metadata(
            metadata_file, search_string, strict, column_names, inverse
        )
    
    # 2) taxon or biome: boolean logic
    loader, search_fn = {
        "taxon": (load_ingredients, search_by_taxon),
        "biome": (load_ingredients, search_by_biome)
    }.get(mode, (None, None))
    if loader is None:
        raise ValueError("mode must be 'taxon', 'metadata' or 'biome'")
    
    ingredients = loader(
        data_dir,
        custom_ingredients=custom_ingredients,
        sandpiper_version=sandpiper_version
    )
    
    # parse into OR‑groups of AND‑terms
    groups = _parse_query(search_string)
    
    # biome mode doesn’t support AND‑chains
    if mode == "biome":
        for terms in groups:
            if len(terms) > 1:
                raise ValueError(
                    "AND queries (using '+') are not supported in biome mode; "
                    f"cannot process group: {terms!r}"
                )
    
    # now build your hits: AND within each group, OR across groups
    total_hits: Set = set()
    for terms in groups:
        # first term
        if mode == "taxon":
            hits = search_fn(ingredients, terms[0], ranks_for_search_inclusion)
        else:  # biome
            hits = search_fn(ingredients, terms[0])
        
        # AND‐chain further terms (only ever for taxon)
        for term in terms[1:]:
            hits &= search_fn(ingredients, term, ranks_for_search_inclusion)
        
        total_hits |= hits
    
    # inversion if requested
    if inverse:
        return set(ingredients.samples) - total_hits
    return total_hits

def search_data(mode, data_dir, output_dir, search_string, ranks_for_search_inclusion=None,
                column_names=None, strict=False, tag="", inverse=False, custom_ingredients=None, sandpiper_version=None):
    """
    File‑based search wrapper for metacooc.
    
    Supports three modes:
      * 'taxon':   loads an Ingredients object (ingredients_raw.pkl) from data_dir
                   and finds all samples whose taxa match `search_string`. You can
                   restrict to a particular taxonomic rank via `ranks_for_search_inclusion`.
      * 'metadata': parses the SRA metadata file in data_dir (as defined by your
                   data_config) and greps for `search_string` in specified `column_names`
                   (or a predefined set if `strict=True`).
      * 'biome':   loads the same Ingredients object and returns samples whose
                   `sample_to_biome` value matches one or more comma‑separated
                   names in `search_string`.
    
    You can also invert the match by passing `inverse=True` in any mode.
    
    The results (one accession per line) are written to
        {output_dir}/search_results{tag}.txt
    
    Parameters
    ----------
    mode : str
        One of 'taxon', 'metadata', or 'biome'.
    data_dir : str
        Directory containing the ingredients and/or metadata files.
    output_dir : str
        Directory to write the results file (created if necessary).
    search_string : str
        Token to search for (or comma‑separated list of biome names in 'biome' mode).
    ranks_for_search_inclusion : Optional[str]
        Taxonomic rank prefix (e.g. 'genus') to filter on in 'taxon' mode.
    column_names : Optional[list of str]
        Metadata columns to restrict the search to in 'metadata' mode.
    strict : bool
        If True in 'metadata' mode, use a predefined set of columns.
    tag : str
        Suffix to append to the results filename (before ".txt").
    inverse : bool
        If True, return the complement of the matching set.
    custom_ingredients : Ingredients or str
        Path to a custom pickled Ingredients object, or an Ingredients instance.
    sandpiper_version : Optional[str]
        Version string to select alternate data files via your data_config.
    
    Returns
    -------
    set
        The set of matching accession IDs (written to the output file).
    """
    matching_accessions = search_data_obj(mode, 
                                          data_dir, 
                                          search_string, 
                                          ranks_for_search_inclusion, 
                                          strict, 
                                          column_names, 
                                          inverse, 
                                          custom_ingredients, 
                                          sandpiper_version)
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        
    output_file = os.path.join(output_dir, f"search_results{tag}.txt")
    with open(output_file, "w") as f:
        for acc in sorted(matching_accessions):
            f.write(f"{acc}\n")
    print(f"Mode={mode!r}: Found {len(matching_accessions)} matching accessions. Results saved to {output_file}")
