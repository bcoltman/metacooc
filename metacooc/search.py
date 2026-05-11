#!/usr/bin/env python3
"""
search.py

Combined search functionality for metacooc.

Exposes:
  1. search_data(
         mode, data_dir, output_dir, search_string,
         ranks_for_search_inclusion=None, column_names=None,
         strict=False, tag="", inverse=False
     )
     - File-based interface: loads files from disk, writes the matching
       accessions to output_dir/search_results.txt, and returns the set of
       matching accessions.

  2. search_data_obj(
         search_mode, search_string, data_dir=None,
         ranks_for_search_inclusion=None, strict=False,
         column_names=None, inverse=False
     )
     - Object-based interface:
         * For mode "taxa_context", "focal_taxa", and "biome", data is an
           Ingredients object (or path resolvable by load_ingredients).
         * For mode "metadata", data is a metadata file.
       By default it returns the set of matching accessions.
       For focal_taxa mode, `return_details=True` returns a structured
       FocalSearchResolution object instead.

Modes
-----
  - taxa_context:
      Cohort / context-defining taxon search using boolean query grammar:
        '|' = OR across groups
        '+' = AND within a group

  - focal_taxa:
      Focal taxon resolution for cooccurrence workflows.
      Supports comma-separated focal taxa queries.
      Endpoint focal queries are exact-only.
      Non-terminal focal queries expand to self + descendants + aggregated self
      (if present).
      Optional 'LHS -> RHS' syntax is supported:
        - LHS defines the focal cohort
        - RHS defines retrieval-target focal taxa within that cohort workflow

  - metadata:
      Searches a metadata file (e.g. sra_metadata.tsv) for the search token,
      optionally restricting the search to specified columns.

  - biome:
      Returns samples whose biome matches one or more comma-separated biome names.
"""

import os
import shutil
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional, Set

from metacooc.pantry import load_ingredients, threshold_ingredients_presence
from metacooc.utils import (
    _RANK_PREFIXES,
    _PREFIX_TO_RANK,
    _parse_tokens,
    _token_rank,
    _terminal_rank_prefix,
    _deepest_rank_token,
)
from metacooc._data_config import *


METADATA_SEARCH_BACKEND_ENV = "METACOOC_METADATA_SEARCH_BACKEND"
METADATA_SEARCH_BACKENDS = {"auto", "external", "python"}
_EXTERNAL_METADATA_SEARCH_PROBE = {}


@dataclass(frozen=True)
class ParsedSearchString:
    lhs: str
    rhs: Optional[str] = None


@dataclass(frozen=True)
class FocalSearchResolution:
    """
    Structured result for focal_taxa search mode.

    Attributes
    ----------
    focal_rows_by_query_lhs
        Mapping from each LHS focal query string to resolved taxon row indices.
    focal_sample_union
        Union of samples containing any resolved LHS focal row.
    focal_rows_by_query_rhs
        Optional mapping from each RHS focal query string to resolved taxon row
        indices. This does not alter focal_sample_union; it is intended for
        downstream retrieval / reporting restriction.
    """
    focal_rows_by_query_lhs: Dict[str, Set[int]]
    focal_sample_union: Set[str]
    focal_rows_by_query_rhs: Optional[Dict[str, Set[int]]] = None


def _split_search_arrow(search_string: str) -> ParsedSearchString:
    """
    Parse optional retrieval syntax:

        LHS
        LHS -> RHS

    Rules
    -----
    - At most one top-level '->'
    - Empty LHS or RHS is invalid
    """
    if search_string is None:
        raise ValueError("search_string cannot be None")

    raw = search_string.strip()
    if not raw:
        raise ValueError("search_string cannot be empty")

    parts = raw.split("->")
    if len(parts) > 2:
        raise ValueError("Only one '->' is allowed in a search string")

    if len(parts) == 1:
        lhs = parts[0].strip()
        if not lhs:
            raise ValueError("search_string cannot be empty")
        return ParsedSearchString(lhs=lhs, rhs=None)

    lhs, rhs = parts[0].strip(), parts[1].strip()
    if not lhs:
        raise ValueError("Left-hand side of '->' cannot be empty")
    if not rhs:
        raise ValueError("Right-hand side of '->' cannot be empty")

    return ParsedSearchString(lhs=lhs, rhs=rhs)


def _parse_query_groups(q: str) -> List[List[str]]:
    """
    Current cohort grammar:
      '|' = OR across groups
      '+' = AND within a group

    Example
    -------
    A|B+C -> [["A"], ["B", "C"]]
    """
    groups: List[List[str]] = []
    for or_part in q.split("|"):
        or_part = or_part.strip()
        if not or_part:
            continue
        terms = [t.strip() for t in or_part.split("+") if t.strip()]
        if terms:
            groups.append(terms)
    return groups


def _parse_comma_queries(q: str) -> List[str]:
    """
    Comma-separated independent queries.
    Used for focal_taxa multi-query mode.
    """
    queries = [part.strip() for part in q.split(",") if part.strip()]
    if not queries:
        raise ValueError("No queries were provided.")
    return queries


def _search_taxon_rows(
    ingredients,
    search_string: str,
    ranks_for_search_inclusion: Optional[str] = None,
) -> Set[int]:
    """
    Row-resolution logic for cohort/context taxon searches.

    Semantics
    ---------
    Uses the deepest ranked token from `search_string`, then expands via the
    cached rank lookup for that token. This is the existing taxa_context logic.
    """
    if not search_string or not search_string.strip():
        return set()

    ingredients._ensure_taxa_lookups()

    rank, token = _deepest_rank_token(search_string)
    if rank is None or token is None:
        return set()

    if rank == "root":
        num_taxa = ingredients.presence_matrix.shape[0]
        candidates: Set[int] = set(range(num_taxa))
    else:
        candidates = set(ingredients._rank_lookups.get(rank, {}).get(token, ()))

    if not candidates:
        return set()

    if ranks_for_search_inclusion:
        rp = ranks_for_search_inclusion.strip().lower()
        if rp not in _RANK_PREFIXES:
            raise ValueError(
                f"Unknown rank '{ranks_for_search_inclusion}'. "
                f"Expected one of: {', '.join(_RANK_PREFIXES.keys())}"
            )

        rank_lookup = ingredients._rank_lookups.get(rp, {})
        if not rank_lookup:
            return set()

        allowed_rows = set().union(*rank_lookup.values())
        candidates &= allowed_rows

        if not candidates:
            return set()

    return candidates


def _parse_focal_query(q: str) -> List[str]:
    if "->" in q:
        raise ValueError(
            "Internal parsing error: focal query parser received unsplit '->' syntax."
        )
    if "|" in q or "+" in q:
        raise ValueError(
            "focal_taxa mode does not support '|' or '+'. "
            "Use commas to specify multiple focal taxa queries."
        )
    return _parse_comma_queries(q)


def _terminal_component(taxon: str) -> str:
    return taxon.split("; ")[-1]


def _exact_terminal_taxon_rows(ingredients, query: str) -> Set[int]:
    query = query.strip()
    matches = {
        i
        for i, taxon in enumerate(ingredients.taxa)
        if taxon == query or _terminal_component(taxon) == query
    }

    return matches


def _is_endpoint_focal_query(query: str) -> bool:
    """
    Endpoint focal queries are resolved by exact matching only.

    Endpoint forms:
      - any query ending in ' AGGREGATED'
      - species-level queries
    """
    q = query.strip()
    if not q:
        return False

    if q.endswith(" AGGREGATED"):
        return True

    rank, token = _deepest_rank_token(q)
    if rank is None or token is None:
        return False

    return rank == "species"

def _resolve_single_focal_taxa_query(ingredients, query: str) -> Set[int]:
    """
    Resolve one focal query.

    Endpoint queries:
      - exact endpoint row(s) only

    Non-terminal ranked queries:
      - descendant/context expansion from _search_taxon_rows()
      - exact row itself, if present
      - aggregated endpoint row, if present
    """
    query = query.strip()
    if not query:
        return set()

    ingredients._ensure_taxa_lookups()

    if _is_endpoint_focal_query(query):
        rows = _exact_terminal_taxon_rows(ingredients, query)

        if not rows:
            raise ValueError(f"No exact endpoint taxon matched focal query: {query!r}")
        return rows

    rank, token = _deepest_rank_token(query)
    if rank is None or token is None:
        raise ValueError(f"Could not parse ranked focal query: {query!r}")

    token = token.strip()
    rows: Set[int] = set()

    rows |= _search_taxon_rows(ingredients, token)
    rows |= _exact_terminal_taxon_rows(ingredients, token)
    rows |= _exact_terminal_taxon_rows(ingredients, f"{token} AGGREGATED")

    # temporary debug
    # print(f"DEBUG expanded focal query={query!r}")
    # print(f"DEBUG expanded token={token!r}")
    # print(f"DEBUG expanded n_matches={len(rows)}")

    if not rows:
        raise ValueError(f"No taxa matched focal query after expansion: {query!r}")

    return rows


def resolve_focal_taxa_queries(
    ingredients,
    search_string: str,
) -> Dict[str, Set[int]]:
    """
    Resolve comma-separated focal taxa queries independently.
    """
    out: Dict[str, Set[int]] = {}

    for query in _parse_focal_query(search_string):
        rows = _resolve_single_focal_taxa_query(ingredients, query)
        if not rows:
            raise ValueError(f"No taxa matched focal query: {query!r}")
        out[query] = rows

    return out


def resolve_rhs_taxa_context_rows(
    ingredients,
    rhs_string: str,
    ranks_for_search_inclusion: Optional[str] = None,
) -> Set[int]:
    """
    Resolve RHS retrieval taxa for context-style grammar.

    RHS grammar here follows the existing cohort grammar:
      '|' = OR
      '+' = AND

    Practical interpretation for row retrieval:
      - each term resolves to a row set
      - '+' intersects row sets
      - '|' unions row sets
    """
    groups = _parse_query_groups(rhs_string)
    total_rows: Set[int] = set()

    for terms in groups:
        row_hits = _search_taxon_rows(
            ingredients,
            terms[0],
            ranks_for_search_inclusion=ranks_for_search_inclusion,
        )
        for term in terms[1:]:
            row_hits &= _search_taxon_rows(
                ingredients,
                term,
                ranks_for_search_inclusion=ranks_for_search_inclusion,
            )
        total_rows |= row_hits

    return total_rows


def resolve_rhs_focal_rows(
    ingredients,
    rhs_string: str,
) -> Dict[str, Set[int]]:
    """
    Resolve RHS retrieval taxa in focal syntax:
      comma-separated independent focal queries
    """
    out: Dict[str, Set[int]] = {}
    for query in _parse_comma_queries(rhs_string):
        out[query] = _resolve_single_focal_taxa_query(ingredients, query)
    return out


def resolve_focal_taxa_search(
    ingredients,
    search_string: str,
) -> FocalSearchResolution:
    """
    Structured resolver for focal_taxa mode.

    LHS defines:
      - focal taxa to seed the focal cohort
      - focal sample union

    RHS, if present, defines:
      - an optional retrieval subset for downstream cooccurrence reporting
    """
    parsed = _split_search_arrow(search_string)

    focal_rows_by_query_lhs = resolve_focal_taxa_queries(ingredients, parsed.lhs)

    all_focal_rows: Set[int] = set()
    for rows in focal_rows_by_query_lhs.values():
        all_focal_rows |= rows

    if all_focal_rows:
        sub = ingredients.presence_matrix[sorted(all_focal_rows), :]
        _, cols = sub.nonzero()
        focal_sample_union = {ingredients.samples[c] for c in cols}
    else:
        focal_sample_union = set()

    focal_rows_by_query_rhs = None
    if parsed.rhs is not None:
        focal_rows_by_query_rhs = resolve_rhs_focal_rows(ingredients, parsed.rhs)

    return FocalSearchResolution(
        focal_rows_by_query_lhs=focal_rows_by_query_lhs,
        focal_sample_union=focal_sample_union,
        focal_rows_by_query_rhs=focal_rows_by_query_rhs,
    )


def search_by_focal_taxa(
    ingredients,
    search_string: str,
) -> Set[str]:
    """
    Resolve focal taxa queries and return the union of samples containing any
    resolved focal taxon row from the LHS query.
    """
    resolved = resolve_focal_taxa_search(ingredients, search_string)
    return resolved.focal_sample_union


def search_by_taxon(
    ingredients,
    search_string: str,
    ranks_for_search_inclusion: Optional[str] = None,
) -> Set[str]:
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
        list: The 1-based indices of the columns for AWK.
    """
    with open(metadata_file, "r") as f:
        headers = f.readline().strip().split(delimiter)

    indices = []
    for column_name in column_names:
        if column_name not in headers:
            print(f"Column '{column_name}' not found in metadata file.")
            raise ValueError(f"Column '{column_name}' not found in metadata file.")
        indices.append(headers.index(column_name) + 1)

    return indices


def _external_metadata_search_mode(column_names=None) -> str:
    return "columns" if column_names else "full"


def _external_metadata_search_tools(mode):
    return ("awk",) if mode == "columns" else ("grep",)


def _metadata_search_env():
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    return env


def _external_metadata_search_probe(mode):
    if mode == "columns":
        cmd = [
            "awk",
            "-F",
            "\t",
            "-v",
            "needle=soil",
            "NR==1{next} index(tolower($2), needle)>0 {print $1}",
        ]
    else:
        cmd = ["grep", "-i", "-F", "-e", "soil"]

    probe = subprocess.run(
        cmd,
        input="acc\tfield\nS001\tSOIL\n",
        capture_output=True,
        text=True,
        env=_metadata_search_env(),
    )
    return (
        probe.returncode == 0
        and "S001" in {line.split("\t", 1)[0] for line in probe.stdout.splitlines()}
    )


def _metadata_search_backend(column_names=None) -> str:
    global _EXTERNAL_METADATA_SEARCH_PROBE

    backend = os.environ.get(METADATA_SEARCH_BACKEND_ENV, "auto").strip().lower()
    if backend not in METADATA_SEARCH_BACKENDS:
        raise ValueError(
            f"{METADATA_SEARCH_BACKEND_ENV} must be one of "
            f"{', '.join(sorted(METADATA_SEARCH_BACKENDS))}; got {backend!r}"
        )

    if backend == "python":
        return "python"

    mode = _external_metadata_search_mode(column_names)
    required = _external_metadata_search_tools(mode)
    missing = [tool for tool in required if shutil.which(tool) is None]
    if missing:
        if backend == "external":
            raise RuntimeError(
                "Metadata search backend 'external' was requested, but required "
                f"tool(s) are missing: {', '.join(missing)}. Install the tools or "
                f"set {METADATA_SEARCH_BACKEND_ENV}=python."
            )
        return "python"

    if mode not in _EXTERNAL_METADATA_SEARCH_PROBE:
        _EXTERNAL_METADATA_SEARCH_PROBE[mode] = _external_metadata_search_probe(mode)

    if _EXTERNAL_METADATA_SEARCH_PROBE[mode]:
        return "external"

    if backend == "external":
        raise RuntimeError(
            "Metadata search backend 'external' was requested, but the external "
            "metadata search command probe failed. Install a POSIX-compatible "
            "shell environment with the required tools or "
            f"set {METADATA_SEARCH_BACKEND_ENV}=python."
        )

    return "python"


def _field_matches(fields, indices, needle_lower, inverse):
    matched = any(
        idx < len(fields) and needle_lower in fields[idx].casefold()
        for idx in indices
    )
    return not matched if inverse else matched


def _python_metadata_search(
    search_string,
    metadata_file,
    column_names=None,
    delimiter="\t",
    inverse=False,
):
    needle_lower = search_string.strip().casefold()
    out = set()

    if not column_names:
        with open(metadata_file, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.rstrip("\n")
                matched = needle_lower in line.casefold()
                if inverse:
                    matched = not matched
                if matched:
                    first = line.split(delimiter, 1)[0]
                    if first:
                        out.add(first)
        return out

    zero_based_indices = [
        idx - 1 for idx in get_column_indices(metadata_file, column_names, delimiter)
    ]
    with open(metadata_file, "r", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        for line in handle:
            fields = line.rstrip("\n").split(delimiter)
            if fields and _field_matches(
                fields, zero_based_indices, needle_lower, inverse
            ):
                out.add(fields[0])
    return out


def _external_metadata_search(
    search_string,
    metadata_file,
    column_names=None,
    delimiter="\t",
    inverse=False,
):
    needle = search_string.strip()

    if not column_names:
        cmd = ["grep", "-i", "-F"]
        if inverse:
            cmd.append("-v")
        cmd.extend(["-e", needle, metadata_file])

        result = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=_metadata_search_env(),
        )
        out = set()
        for line in result.stdout:
            first = line.rstrip("\n").split(delimiter, 1)[0]
            if first:
                out.add(first)
        stderr = result.stderr.read()
        returncode = result.wait()
        if returncode not in (0, 1):
            raise RuntimeError(
                f"Metadata search command failed with exit code {returncode}: "
                f"{stderr.strip()}"
            )
        return out

    col_idxs = get_column_indices(metadata_file, column_names, delimiter)
    needle_lower = needle.lower()
    if inverse:
        conds = [f"index(tolower(${i}), needle)==0" for i in col_idxs]
        cond = " && ".join(conds)
    else:
        conds = [f"index(tolower(${i}), needle)>0" for i in col_idxs]
        cond = " || ".join(conds)

    cmd = [
        "awk",
        "-F",
        delimiter,
        "-v",
        f"needle={needle_lower}",
        f"NR==1{{next}} {cond} {{print $1}}",
        metadata_file,
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=_metadata_search_env(),
    )
    if result.returncode not in (0, 1):
        raise RuntimeError(
            f"Metadata search command failed with exit code {result.returncode}: "
            f"{result.stderr.strip()}"
        )
    return set(line for line in result.stdout.splitlines() if line)


def grep_metadata(search_string, metadata_file, column_names=None, delimiter="\t", inverse=False):
    if not os.path.exists(metadata_file):
        raise FileNotFoundError(f"Metadata file '{metadata_file}' not found.")

    backend = _metadata_search_backend(column_names)
    if backend == "external":
        return _external_metadata_search(
            search_string,
            metadata_file,
            column_names=column_names,
            delimiter=delimiter,
            inverse=inverse,
        )
    return _python_metadata_search(
        search_string,
        metadata_file,
        column_names=column_names,
        delimiter=delimiter,
        inverse=inverse,
    )


def search_in_metadata(metadata, search_string, strict=False, column_names=None, inverse=False):
    """
    Search for a token in a metadata file with optional column restriction
    and inverse search capability.
    """
    if strict:
        column_names = [
            "acc",
            "organism",
            "env_biome_sam",
            "env_feature_sam",
            "env_material_sam",
            "biosamplemodel_sam",
        ]
    return grep_metadata(search_string, metadata, column_names, inverse=inverse)


def search_by_biome(ingredients, biome_names):
    """
    Return sample IDs whose biome (level_1 or level_2) matches any of `biome_names`.
    """
    if not hasattr(ingredients, "biomes_order"):
        ingredients._allocate_biomes()

    requested = (
        [b.strip() for b in biome_names.split(",")]
        if isinstance(biome_names, str)
        else list(biome_names)
    )

    available = set(ingredients.biomes_order.get("level_1", [])) | set(
        ingredients.biomes_order.get("level_2", [])
    )
    bad = [b for b in requested if b not in available]
    if bad:
        raise ValueError(f"Unknown biome(s): {bad}. Available: {sorted(available)}")

    out = set()
    for sample, (b1, b2) in ingredients.sample_to_biome.items():
        if (b1 in requested) or (b2 in requested):
            out.add(sample)
    return out


def _parse_query(q: str) -> List[List[str]]:
    return _parse_query_groups(q)


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
    metadata_file=None,
    return_details: bool = False,
    min_coverage=None,
    min_coverage_by_rank=None,
    min_relative_abundance=None,
    min_relative_abundance_by_rank=None,
):
    """
    Object-based search interface.

    Returns
    -------
    Default:
        set of matching accession IDs

    For focal_taxa with return_details=True:
        FocalSearchResolution
    """
    if search_mode == "taxon":
        raise ValueError(
            "search_mode='taxon' has been removed. "
            "Use 'taxa_context' for cohort/context taxon searches or "
            "'focal_taxa' for focal cooccurrence mode."
        )

    parsed = _split_search_arrow(search_string)

    if search_mode == "metadata":
        if parsed.rhs is not None:
            raise ValueError("'->' syntax is not supported in metadata mode")

        if metadata_file is None:
            defaulted_data_version = data_version is None
            data_version = data_version or LATEST_VERSION
            filenames, _ = get_file_info(data_version)
            if not data_dir:
                raise ValueError("data_dir must be provided if searching metadata without metadata_file")
            metadata_file = os.path.join(data_dir, filenames["sra_metadata"])
            if not os.path.exists(metadata_file):
                raise DataVersionError(
                    missing_data_version_message(
                        data_dir=data_dir,
                        data_version=data_version,
                        missing_path=metadata_file,
                        file_kind="metadata",
                        defaulted=defaulted_data_version,
                    )
                )
        elif not os.path.exists(metadata_file):
            raise FileNotFoundError(f"Missing '{metadata_file}'")

        if return_details:
            raise ValueError("return_details=True is only supported for focal_taxa mode")

        return search_in_metadata(
            metadata_file, parsed.lhs, strict, column_names, inverse
        )

    loader, search_fn = {
        "taxa_context": (load_ingredients, search_by_taxon),
        "focal_taxa": (load_ingredients, search_by_focal_taxa),
        "biome": (load_ingredients, search_by_biome),
    }.get(search_mode, (None, None))

    if loader is None:
        raise ValueError(
            "search_mode must be one of: "
            "'focal_taxa', 'taxa_context', 'metadata', 'biome'"
        )

    if search_mode == "focal_taxa" and inverse:
        raise ValueError("inverse searches are not supported in focal_taxa mode")

    if search_mode == "biome" and parsed.rhs is not None:
        raise ValueError("'->' syntax is not supported in biome mode")

    if search_mode == "taxa_context" and parsed.rhs is not None:
        raise ValueError("'->' syntax is not supported in taxa_context mode")

    if return_details and search_mode != "focal_taxa":
        raise ValueError("return_details=True is only supported for focal_taxa mode")

    ingredients = loader(
        data_dir,
        aggregated=aggregated,
        custom_ingredients=custom_ingredients,
        data_version=data_version,
    )
    ingredients = threshold_ingredients_presence(
        ingredients,
        min_coverage=min_coverage,
        min_coverage_by_rank=min_coverage_by_rank,
        min_relative_abundance=min_relative_abundance,
        min_relative_abundance_by_rank=min_relative_abundance_by_rank,
    )

    if search_mode == "focal_taxa":
        resolved = resolve_focal_taxa_search(ingredients, search_string)
        if return_details:
            return resolved
        return resolved.focal_sample_union

    groups = _parse_query_groups(parsed.lhs)

    if search_mode == "biome":
        for terms in groups:
            if len(terms) > 1:
                raise ValueError(
                    "AND queries (using '+') are not supported in biome mode; "
                    f"cannot process group: {terms!r}"
                )

    total_hits: Set[str] = set()

    for terms in groups:
        if search_mode == "taxa_context":
            hits = search_fn(
                ingredients,
                terms[0],
                ranks_for_search_inclusion,
            )
            for term in terms[1:]:
                hits &= search_fn(
                    ingredients,
                    term,
                    ranks_for_search_inclusion,
                )
        elif search_mode == "biome":
            hits = search_fn(ingredients, terms[0])
            for term in terms[1:]:
                hits &= search_fn(ingredients, term)
        else:
            raise ValueError(f"Unexpected search_mode during query resolution: {search_mode!r}")

        total_hits |= hits

    if inverse:
        return set(ingredients.samples) - total_hits
    return total_hits


def search_data(
    mode,
    data_dir,
    output_dir,
    search_string,
    ranks_for_search_inclusion=None,
    column_names=None,
    strict=False,
    tag="",
    inverse=False,
    custom_ingredients=None,
    data_version=None,
    list_column_names=False,
    list_biomes=False,
    aggregated=False,
    metadata_file=None,
    min_coverage=None,
    min_coverage_by_rank=None,
    min_relative_abundance=None,
    min_relative_abundance_by_rank=None,
):
    """
    File-based search wrapper for metacooc.

    Supports four modes:
      * 'taxa_context':
            loads an Ingredients object and resolves a cohort/sample set using
            taxon query grammar ('|' and '+').
      * 'focal_taxa':
            resolves focal taxa queries (comma-separated), returning the union
            of samples containing the resolved focal taxa from the LHS.
            This is primarily for cooccurrence workflows.
      * 'metadata':
            searches the metadata table.
      * 'biome':
            returns samples whose biome matches one or more requested biome names.

    The results (one accession per line) are written to:
        {output_dir}/{tag}search_results.txt
    """
    if list_biomes:
        ingredients = load_ingredients(
            data_dir=data_dir,
            aggregated=aggregated,
            custom_ingredients=custom_ingredients,
            data_version=data_version,
        )
        print(format_available_biomes(ingredients))
        return

    if list_column_names:
        if metadata_file is None:
            defaulted_data_version = data_version is None
            data_version = data_version or LATEST_VERSION
            filenames, _ = get_file_info(data_version)
            if not data_dir:
                raise ValueError("data_dir must be provided if listing metadata columns without metadata_file")
            metadata_file = os.path.join(data_dir, filenames["sra_metadata"])
            if not os.path.exists(metadata_file):
                raise DataVersionError(
                    missing_data_version_message(
                        data_dir=data_dir,
                        data_version=data_version,
                        missing_path=metadata_file,
                        file_kind="metadata",
                        defaulted=defaulted_data_version,
                    )
                )
        elif not os.path.exists(metadata_file):
            raise FileNotFoundError(f"Missing '{metadata_file}'")
        with open(metadata_file, "r") as f:
            headers = f.readline().strip().split("\t")
            print(headers)
        return

    matching_accessions = search_data_obj(
        search_mode=mode,
        search_string=search_string,
        data_dir=data_dir,
        ranks_for_search_inclusion=ranks_for_search_inclusion,
        strict=strict,
        column_names=column_names,
        inverse=inverse,
        custom_ingredients=custom_ingredients,
        data_version=data_version,
        aggregated=aggregated,
        metadata_file=metadata_file,
        min_coverage=min_coverage,
        min_coverage_by_rank=min_coverage_by_rank,
        min_relative_abundance=min_relative_abundance,
        min_relative_abundance_by_rank=min_relative_abundance_by_rank,
    )

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, f"{tag}search_results.txt")
    with open(output_file, "w") as f:
        for acc in sorted(matching_accessions):
            f.write(f"{acc}\n")

    print(
        f"Mode={mode!r}: Found {len(matching_accessions)} matching accessions. "
        f"Results saved to {output_file}"
    )


def format_available_biomes(ingredients) -> str:
    """
    Format available biome query terms for CLI output.
    """
    available = ingredients.available_biomes()
    if not any(available.values()):
        return "No biome annotations are available for this Ingredients object."

    lines = []
    for level in ("level_1", "level_2"):
        lines.append(f"{level}:")
        rows = available.get(level, [])
        if rows:
            for row in rows:
                sample_word = "sample" if int(row["n_samples"]) == 1 else "samples"
                lines.append(f"  {row['biome']} ({row['n_samples']} {sample_word})")
        else:
            lines.append("  none")
        if level != "level_2":
            lines.append("")
    return "\n".join(lines)
