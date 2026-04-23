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

import os, subprocess, shlex
from typing import Dict, List, Set, Optional

from metacooc.pantry import load_ingredients
from metacooc.utils import (
    _RANK_PREFIXES, 
    _PREFIX_TO_RANK, 
    _parse_tokens, 
    _token_rank, 
    _terminal_rank_prefix, 
    _deepest_rank_token
)

from metacooc._data_config import *


def _search_taxon_rows(
    ingredients,
    search_string: str,
    ranks_for_search_inclusion: Optional[str] = None,
) -> Set[int]:
    """
    Core row-resolution logic used by search_by_taxon and other utilities.
    
    Returns:
        Set[int]: rows indices in ingredients.taxa matching the search_string.
    """
    if not search_string or not search_string.strip():
        return set()
        
    ingredients._ensure_taxa_lookups()
    
    rank, token = _deepest_rank_token(search_string)
    if rank is None:
        return set()
    
    # If searching by 'Root', that means “everything”
    if rank == "root":
        num_taxa = ingredients.presence_matrix.shape[0]
        candidates: Set[int] = set(range(num_taxa))
    else:
        candidates = set(ingredients._rank_lookups[rank].get(token, ()))
    
    if not candidates:
        return set()
    
    if ranks_for_search_inclusion:
        rp = ranks_for_search_inclusion.strip().lower()
        if rp not in _RANK_PREFIXES:
            raise ValueError(
                f"Unknown rank '{ranks_for_search_inclusion}'. "
                f"Expected one of: {', '.join(_RANK_PREFIXES.keys())}"
            )
        
        lookups_for_rank = ingredients._rank_lookups.get(rp, {})
        if lookups_for_rank:
            cols_with_rank = set().union(*lookups_for_rank.values())
            candidates &= cols_with_rank
        else:
            # No taxa have this rank at all
            return set()
        
        if not candidates:
            return set()
    
    return candidates

def _normalise_search_mode(search_mode: str) -> str:
    search_mode = search_mode.lower().strip()
    if search_mode == "taxon":
        return "taxa_context"
    return search_mode

def _parse_focal_query(q: str) -> List[str]:
    if "|" in q or "+" in q:
        raise ValueError(
            "focal_taxa mode does not support '|' or '+'. "
            "Use commas to specify multiple focal taxa queries."
        )
    queries = [part.strip() for part in q.split(",") if part.strip()]
    if not queries:
        raise ValueError("No focal taxa queries were provided.")
    return queries

def _is_endpoint_focal_query(search_string: str) -> bool:
    rank, token = _deepest_rank_token(search_string)
    if rank is None or token is None:
        return False
    return (rank == "species") or token.endswith(" AGGREGATED")

def resolve_focal_taxa_queries(
    ingredients,
    search_string: str,
) -> Dict[str, Set[int]]:
    """
    Resolve focal taxa queries into explicit taxon-row sets.
    
    Endpoint queries (species or '* AGGREGATED') resolve to themselves only.
    Non-terminal ranked queries resolve to the row itself, descendants, and
    the aggregated counterpart if present.
    """
    ingredients._ensure_taxa_lookups()
    
    out: Dict[str, Set[int]] = {}
    for query in _parse_focal_query(search_string):
        rows = _search_taxon_rows(ingredients, query)
        if not rows:
            raise ValueError(f"No taxa matched focal query: {query!r}")
            
        if not _is_endpoint_focal_query(query):
            _, token = _deepest_rank_token(query)
            if token is not None:
                rows |= _search_taxon_rows(ingredients, f"{token} AGGREGATED")
                
        if not rows:
            raise ValueError(f"No taxa matched focal query after expansion: {query!r}")
            
        out[query] = rows
        
    return out

def search_by_focal_taxa(
    ingredients,
    search_string: str,
) -> set:
    """
    Resolve focal taxa queries and return the union of samples containing any
    resolved focal taxon row.
    """
    focal_rows = resolve_focal_taxa_queries(ingredients, search_string)
    
    all_rows: Set[int] = set()
    for rows in focal_rows.values():
        all_rows |= rows
        
    if not all_rows:
        return set()
        
    sub = ingredients.presence_matrix[sorted(all_rows), :]
    _, cols = sub.nonzero()
    return {ingredients.samples[c] for c in cols}

def search_by_taxon(
    ingredients,
    search_string: str,
    ranks_for_search_inclusion: Optional[str] = None,
) -> set:
    """
    Exact taxonomy search using ONLY the deepest ranked token from `search_string`,
    backed by a per-rank lookup cache on the Ingredients object.
    
    Returns a set of sample IDs that contain any of the matching taxa.
    """
    candidates = _search_taxon_rows(
        ingredients,
        search_string,
        ranks_for_search_inclusion=ranks_for_search_inclusion,
    )
    if not candidates:
        return set()
    
    sub = ingredients.presence_matrix[sorted(candidates), :]
    _, cols = sub.nonzero()
    return {ingredients.samples[c] for c in cols}

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
            print(f"Column '{column_name}' not found in metadata file.")
            raise ValueError(f"Column '{column_name}' not found in metadata file.")
        indices.append(headers.index(column_name) + 1)  # Convert 0-based to 1-based index for AWK
        
    return indices


def grep_metadata(search_string, metadata_file, column_names=None, delimiter="\t", inverse=False):
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file '{metadata_file}' not found.")

    needle = search_string.strip()

    if not column_names:
        # Fast path: whole-line search. grep -iF is usually faster than awk+tower/IGNORECASE.
        flag = "-ivF" if inverse else "-iF"
        cmd = f"LC_ALL=C grep {flag} {shlex.quote(needle)} {shlex.quote(metadata_file)} | cut -f1"
    else:
        col_idxs = get_column_indices(metadata_file, column_names, delimiter)
        if inverse:
            conds = [f'index(${i}, needle)==0' for i in col_idxs]
            cond = " && ".join(conds)
        else:
            conds = [f'index(${i}, needle)>0' for i in col_idxs]
            cond = " || ".join(conds)
        # IGNORECASE=1 makes index() case-insensitive without tolower()
        cmd = (
            f"LC_ALL=C awk -F'{delimiter}' -v IGNORECASE=1 -v needle={shlex.quote(needle)} "
            f"'NR==1{{next}} {cond} {{print $1}}' {shlex.quote(metadata_file)}"
        )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return set(line for line in result.stdout.splitlines() if line)

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
    Return sample IDs whose biome (level_1 or level_2) matches any of `biome_names`.
    """
    if not hasattr(ingredients, "biomes_order"):
        ingredients._allocate_biomes()
        
    # Accept "name1,name2" or list-like
    requested = (
        [b.strip() for b in biome_names.split(",")]
        if isinstance(biome_names, str)
        else list(biome_names)
    )
    
    available = set(ingredients.biomes_order.get("level_1", [])) | \
                set(ingredients.biomes_order.get("level_2", []))
    bad = [b for b in requested if b not in available]
    if bad:
        raise ValueError(f"Unknown biome(s): {bad}. Available: {sorted(available)}")
        
    out = set()
    for sample, (b1, b2) in ingredients.sample_to_biome.items():
        if (b1 in requested) or (b2 in requested):
            out.add(sample)
    return out

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

# def _parse_query(q: str) -> List[List[List[str]]]:
    # """
    # ',' → separate queries  (highest-level OR)
    # '|' → OR within a query
    # '+' → AND terms
    # """
    # queries = []
    
    # for comma_block in q.split(","):
        # or_groups = []
        # for or_part in comma_block.split("|"):
            # terms = [t.strip() for t in or_part.split("+") if t.strip()]
            # if terms:
                # or_groups.append(terms)
                
        # if or_groups:
            # queries.append(or_groups)
            
    # return queries

def search_data_obj(
    search_mode: str,
    search_string: str,
    data_dir: str = None,
    ranks_for_search_inclusion=None,
    strict=False,
    column_names=None,
    inverse=False,
    custom_ingredients=None,
    data_version=None,
    aggregated: bool = False,
) -> Set:
    search_mode = _normalise_search_mode(search_mode)
    
    if search_mode == "metadata":
        data_version = data_version or LATEST_VERSION
        filenames, _ = get_file_info(data_version)
        if not data_dir:
            raise ValueError("data_dir must be provided if searching metadata")
        metadata_file = os.path.join(data_dir, filenames["sra_metadata"])
        if not os.path.exists(metadata_file):
            raise FileNotFoundError(f"Missing '{metadata_file}'")
        return search_in_metadata(
            metadata_file, search_string, strict, column_names, inverse
        )
    
    loader, search_fn = {
        "taxa_context": (load_ingredients, search_by_taxon),
        "focal_taxa": (load_ingredients, search_by_focal_taxa),
        "biome": (load_ingredients, search_by_biome),
    }.get(search_mode, (None, None))
    if loader is None:
        raise ValueError(
            "search_mode must be one of: 'focal_taxa', 'taxa_context', 'metadata', 'biome'"
        )
    
    if search_mode == "focal_taxa" and inverse:
        raise ValueError("inverse searches are not supported in focal_taxa mode")
    
    ingredients = loader(
        data_dir,
        aggregated=aggregated,
        custom_ingredients=custom_ingredients,
        data_version=data_version,
    )
    
    if search_mode == "focal_taxa":
        return search_fn(ingredients, search_string)
    
    groups = _parse_query(search_string)
    
    if search_mode == "biome":
        for terms in groups:
            if len(terms) > 1:
                raise ValueError(
                    "AND queries (using '+') are not supported in biome mode; "
                    f"cannot process group: {terms!r}"
                )
    
    total_hits: Set = set()
    for terms in groups:
        if search_mode == "taxa_context":
            hits = search_fn(ingredients, terms[0], ranks_for_search_inclusion)
        else:
            hits = search_fn(ingredients, terms[0])
        
        for term in terms[1:]:
            hits &= search_fn(ingredients, term, ranks_for_search_inclusion)
        
        total_hits |= hits
    
    if inverse:
        return set(ingredients.samples) - total_hits
    return total_hits

def search_data(mode, data_dir, output_dir, search_string, ranks_for_search_inclusion=None,
                column_names=None, strict=False, tag="", inverse=False, custom_ingredients=None,
                data_version=None, list_column_names=False, aggregated=False):
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
        {output_dir}/{tag}search_results.txt
    
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
    data_version : Optional[str]
        Version string to select alternate data files via your data_config.
    
    Returns
    -------
    set
        The set of matching accession IDs (written to the output file).
    """
    if list_column_names:
            
        data_version = data_version or LATEST_VERSION
        filenames, _ = get_file_info(data_version)
        if not data_dir:
            raise ValueError(
                "data_dir must be provided if searching metadata"
            )
        metadata_file = os.path.join(data_dir, filenames["sra_metadata"])
        if not os.path.exists(metadata_file):
            raise FileNotFoundError(f"Missing '{metadata_file}'")
        with open(metadata_file, "r") as f:
            headers = f.readline().strip().split("\t")
            print(headers)
        return
    
    matching_accessions = search_data_obj(mode,
                                          search_string, 
                                          data_dir, 
                                          ranks_for_search_inclusion, 
                                          strict, 
                                          column_names, 
                                          inverse, 
                                          custom_ingredients, 
                                          data_version,
                                          aggregated,
                                          )
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        
    output_file = os.path.join(output_dir, f"{tag}search_results.txt")
    with open(output_file, "w") as f:
        for acc in sorted(matching_accessions):
            f.write(f"{acc}\n")
    print(f"Mode={mode!r}: Found {len(matching_accessions)} matching accessions. Results saved to {output_file}")
