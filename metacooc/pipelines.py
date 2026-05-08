#!/usr/bin/env python3
"""
pipelines.py

Pipeline functions for metacooc.

This module provides shared setup and command-specific execution for:
  - cooccurrence
  - association
  - structure
  - biome distribution

Design notes
------------
- taxa_context / metadata / biome are cohort-defining search modes
- focal_taxa is a cooccurrence-specific focal search mode
- focal_taxa may use:
      LHS
      LHS -> RHS
  where:
      * LHS defines the focal cohort and focal taxa
      * RHS optionally defines retrieval-target taxa for downstream reporting
"""

import os
import numpy as np
import pandas as pd

from metacooc.search import resolve_focal_taxa_queries, search_data_obj
from metacooc.filter import filter_data_obj
from metacooc.pantry import load_ingredients
from metacooc.analysis import (
    association_obj,
    cooccurrence_obj,
    select_taxa_universe,
    export_cooccurrence_outputs,
)
from metacooc.plot import plot_analysis_obj
from metacooc.clustering import determine_taxa_context
from metacooc.structure import structure_obj


def _rows_to_taxa_in_universe(
    ingredients,
    rows_by_query,
    taxa_universe,
):
    """
    Convert query -> row-index mapping into query -> taxa-name mapping,
    restricted to taxa remaining in taxa_universe.

    Output taxa are deduplicated while preserving order.
    """
    if not rows_by_query:
        return None

    taxa_arr = np.asarray(ingredients.taxa, dtype=object)
    taxa_universe_set = set(taxa_universe)

    out = {}
    for query, row_idx in rows_by_query.items():
        row_idx_arr = np.asarray(sorted(row_idx), dtype=np.int64)
        if row_idx_arr.size == 0:
            continue

        resolved_taxa = []
        seen = set()
        for taxon in taxa_arr[row_idx_arr].tolist():
            if taxon in taxa_universe_set and taxon not in seen:
                seen.add(taxon)
                resolved_taxa.append(taxon)

        if resolved_taxa:
            out[query] = resolved_taxa

    return out if out else None


def _taxa_indices_from_query(ingredients, taxa_query):
    rows_by_query = resolve_focal_taxa_queries(ingredients, taxa_query)
    rows = set()
    for query_rows in rows_by_query.values():
        rows.update(query_rows)
    return sorted(rows)


def run_shared_pipeline_setup(args):
    """
    Shared setup for cooccurrence, association, and structure pipelines.

    Returns
    -------
    tuple
        (
            null_ingredients,
            filtered_ingredients,
            taxa_universe,
            output_dir,
            focal_query_to_taxa,
            rhs_query_to_taxa,
        )

    Notes
    -----
    - For non-focal modes, focal_query_to_taxa and rhs_query_to_taxa are None.
    - For focal_taxa mode:
        * LHS determines the accession cohort
        * LHS is re-resolved after filtering to determine focal taxa remaining
        * RHS, if present, is resolved as taxon names using the same Ingredients
          object that produced the RHS row indices, and then restricted to the
          final analysis universe.
    """
    ingredients = load_ingredients(
        args.data_dir,
        args.aggregated,
        args.custom_ingredients,
        args.data_version,
    )

    focal_query_to_taxa = None
    rhs_query_to_taxa = None

    if args.search_mode == "focal_taxa":
        focal_resolution = search_data_obj(
            search_mode="focal_taxa",
            data_dir=args.data_dir,
            search_string=args.search_string,
            ranks_for_search_inclusion=args.ranks_for_search_inclusion,
            strict=args.strict,
            column_names=args.column_names,
            inverse=args.inverse,
            custom_ingredients=ingredients,
            data_version=args.data_version,
            aggregated=args.aggregated,
            metadata_file=getattr(args, "metadata_file", None),
            return_details=True,
        )
        matching_accessions = focal_resolution.focal_sample_union
    else:
        matching_accessions = search_data_obj(
            search_mode=args.search_mode,
            data_dir=args.data_dir,
            search_string=args.search_string,
            ranks_for_search_inclusion=args.ranks_for_search_inclusion,
            strict=args.strict,
            column_names=args.column_names,
            inverse=args.inverse,
            custom_ingredients=ingredients,
            data_version=args.data_version,
            aggregated=args.aggregated,
            metadata_file=getattr(args, "metadata_file", None),
        )

    if not matching_accessions:
        print("Pipeline: No matching accessions found. Exiting pipeline.")
        return None, None, None, None, None, None

    print(f"Pipeline: Found {len(matching_accessions)} matching accessions after initial filtering.")

    # LHS cohort filtering
    int_ingredients, is_successful = filter_data_obj(
        ingredients,
        accession_set=matching_accessions,
        min_taxa_count=args.min_taxa_count,
        min_sample_count=args.min_sample_count,
        filter_rank=args.filter_rank,
        taxa_count_rank=args.taxa_count_rank,
    )
    if not is_successful:
        return None, None, None, None, None, None

    sub_samples = int_ingredients.samples

    # Build null/background side first
    if args.null_scope is None:
        null_ingredients, is_successful = filter_data_obj(
            ingredients,
            accession_set=None,
            min_taxa_count=0 if args.remove_null_threshold else args.min_taxa_count,
            min_sample_count=0 if args.remove_null_threshold else args.min_sample_count,
            filter_rank=args.filter_rank,
            taxa_count_rank=args.taxa_count_rank,
        )
        if not is_successful:
            return None, None, None, None, None, None

    elif args.null_scope == "taxa":
        null_ingredients, is_successful = filter_data_obj(
            ingredients,
            accession_set=None,
            min_taxa_count=0 if args.remove_null_threshold else args.min_taxa_count,
            min_sample_count=0 if args.remove_null_threshold else args.min_sample_count,
            filter_rank=args.filter_rank,
            taxa_count_rank=args.taxa_count_rank,
        )
        if not is_successful:
            return None, None, None, None, None, None

        null_ingredients = determine_taxa_context(
            null_ingredients,
            focal_taxa=args.null_taxa_query,
            degree=args.taxa_degree,
            min_shared_samples_between_taxa=args.min_shared_samples_between_taxa,
        )

    elif args.null_scope in {"biome", "metadata"}:
        search_string = args.null_biome_query if args.null_scope == "biome" else args.null_metadata_query

        null_matching_accessions = search_data_obj(
            search_mode=args.null_scope,
            search_string=search_string,
            data_dir=args.data_dir,
            custom_ingredients=ingredients,
            data_version=args.data_version,
            aggregated=args.aggregated,
            metadata_file=getattr(args, "metadata_file", None),
        )

        null_ingredients, is_successful = filter_data_obj(
            ingredients,
            accession_set=null_matching_accessions,
            min_taxa_count=0 if args.remove_null_threshold else args.min_taxa_count,
            min_sample_count=0 if args.remove_null_threshold else args.min_sample_count,
            filter_rank=args.filter_rank,
            taxa_count_rank=args.taxa_count_rank,
        )
        if not is_successful:
            return None, None, None, None, None, None

    elif args.null_scope in {"biome_taxa", "metadata_taxa"}:
        search_mode = "biome" if args.null_scope == "biome_taxa" else "metadata"
        search_string = args.null_biome_query if search_mode == "biome" else args.null_metadata_query

        null_matching_accessions = search_data_obj(
            search_mode=search_mode,
            search_string=search_string,
            data_dir=args.data_dir,
            custom_ingredients=ingredients,
            data_version=args.data_version,
            aggregated=args.aggregated,
            metadata_file=getattr(args, "metadata_file", None),
        )

        null_ingredients, is_successful = filter_data_obj(
            ingredients,
            accession_set=null_matching_accessions,
            min_taxa_count=0 if args.remove_null_threshold else args.min_taxa_count,
            min_sample_count=0 if args.remove_null_threshold else args.min_sample_count,
            filter_rank=args.filter_rank,
            taxa_count_rank=args.taxa_count_rank,
        )
        if not is_successful:
            return None, None, None, None, None, None

        null_ingredients = determine_taxa_context(
            null_ingredients,
            focal_taxa=args.null_taxa_query,
            degree=args.taxa_degree,
            min_shared_samples_between_taxa=args.min_shared_samples_between_taxa,
        )

    else:
        raise ValueError(f"Unknown null_scope: {args.null_scope!r}")

    # Embed focal samples inside null/background space
    filtered_ingredients, is_successful = filter_data_obj(
        null_ingredients,
        accession_set=sub_samples,
    )
    if not is_successful:
        return None, None, None, None, None, None

    # Final analysis universe comes from the actual background
    taxa_universe = select_taxa_universe(null_ingredients, rank=args.filter_rank)

    if args.search_mode == "focal_taxa":
        # Re-resolve LHS on the filtered focal cohort
        filtered_focal_resolution = search_data_obj(
            search_mode="focal_taxa",
            data_dir=args.data_dir,
            search_string=args.search_string,
            ranks_for_search_inclusion=args.ranks_for_search_inclusion,
            strict=args.strict,
            column_names=args.column_names,
            inverse=False,
            custom_ingredients=filtered_ingredients,
            data_version=args.data_version,
            aggregated=args.aggregated,
            metadata_file=getattr(args, "metadata_file", None),
            return_details=True,
        )

        focal_query_to_taxa = _rows_to_taxa_in_universe(
            ingredients=filtered_ingredients,
            rows_by_query=filtered_focal_resolution.focal_rows_by_query_lhs,
            taxa_universe=taxa_universe,
        )

        if not focal_query_to_taxa:
            raise ValueError("No focal taxa remain after sample/count/rank filtering.")

        # IMPORTANT:
        # RHS row indices came from `focal_resolution`, which was resolved on `ingredients`.
        # Therefore convert them to taxon names using `ingredients`, not null_ingredients.
        if focal_resolution.focal_rows_by_query_rhs is not None:
            raw_rhs_query_to_taxa = _rows_to_taxa_in_universe(
                ingredients=ingredients,
                rows_by_query=focal_resolution.focal_rows_by_query_rhs,
                taxa_universe=taxa_universe,
            )

            if raw_rhs_query_to_taxa is None:
                print(
                    "RHS retrieval taxa were requested, but none could be resolved in the final "
                    "analysis universe after filtering. "
                    "No RHS-restricted co-occurrence output will be produced under the current settings."
                )
                rhs_query_to_taxa = None
            else:
                rhs_taxa_flat = []
                seen = set()
                for taxa_list in raw_rhs_query_to_taxa.values():
                    for taxon in taxa_list:
                        if taxon not in seen:
                            seen.add(taxon)
                            rhs_taxa_flat.append(taxon)

                rhs_query_to_taxa = {
                    focal_query: list(rhs_taxa_flat)
                    for focal_query in focal_query_to_taxa
                }

    if (
        getattr(args, "command", None) == "association"
        and set(filtered_ingredients.samples) == set(null_ingredients.samples)
    ):
        raise ValueError(
            "Term and null cohorts are identical — association requires a broader null.\n"
            "Fix by widening null_biome/null_scope or narrowing search_string."
        )

    os.makedirs(args.output_dir, exist_ok=True)
    return (
        null_ingredients,
        filtered_ingredients,
        taxa_universe,
        args.output_dir,
        focal_query_to_taxa,
        rhs_query_to_taxa,
    )


def run_structure(args):
    """
    Execute the structure analysis pipeline.
    """
    null_ing, filt_ing, taxa_universe, out_dir, _, _ = run_shared_pipeline_setup(args)

    if null_ing is None:
        return

    print(
        "Pipeline: Structure analysis being performed on Filtered Ingredients Presence Matrix with "
        f"{filt_ing.presence_matrix.shape[0]} taxa & {filt_ing.presence_matrix.shape[1]} samples"
    )

    structure_df = structure_obj(
        filt_ing,
        null_model=args.null_model,
        nm_n_reps=args.nm_n_reps,
        nm_seed=args.nm_seed,
        nm_n_workers=getattr(args, "nm_n_workers", None),
        nm_mp_start=getattr(args, "nm_mp_start", None),
        nm_sort_indices=getattr(args, "nm_sort_indices", False),
        nm_burn_in_steps=getattr(args, "nm_burn_in_steps", None),
        nm_steps_per_rep=getattr(args, "nm_steps_per_rep", None),
        nm_progress_every=getattr(args, "nm_progress_every", 25),
    )

    null_scope_prefix = "global" if args.null_scope is None else str(args.null_scope)
    output_path = os.path.join(out_dir, f"{args.tag}{null_scope_prefix}_structure.tsv")
    structure_df.to_csv(output_path, sep="\t", index=False)
    print(f"Pipeline: {null_scope_prefix} structure analysis saved to {output_path}")


def run_association(args):
    """
    Execute the association analysis pipeline.
    """
    null_ing, filt_ing, taxa_universe, out_dir, _, _ = run_shared_pipeline_setup(args)

    if null_ing is None:
        return

    print(
        "Pipeline: Association analysis being performed with Null Ingredients Presence Matrix with "
        f"{null_ing.presence_matrix.shape[0]} taxa & {null_ing.presence_matrix.shape[1]} samples"
    )

    single_df = association_obj(
        null_ing,
        filt_ing,
        threshold=args.threshold,
        null_model=args.null_model,
        nm_n_reps=args.nm_n_reps,
        compute_fisher=args.compute_fisher,
        nm_seed=args.nm_seed,
        nm_n_workers=getattr(args, "nm_n_workers", None),
        nm_mp_start=getattr(args, "nm_mp_start", None),
        nm_sort_indices=getattr(args, "nm_sort_indices", False),
        nm_burn_in_steps=getattr(args, "nm_burn_in_steps", None),
        nm_steps_per_rep=getattr(args, "nm_steps_per_rep", None),
        nm_progress_every=getattr(args, "nm_progress_every", 25),
    )

    null_scope_prefix = "global" if args.null_scope is None else str(args.null_scope)
    output_path = os.path.join(out_dir, f"{args.tag}{null_scope_prefix}_association.tsv")
    single_df.to_csv(output_path, sep="\t", index=False)
    print(f"Pipeline: {null_scope_prefix} association analysis saved to {output_path}")

    output_plot_file = os.path.join(out_dir, f"{args.tag}{null_scope_prefix}_plot.png")
    plot_analysis_obj(single_df, out_file=output_plot_file)
    print(f"Pipeline: Plotting {output_plot_file} complete.")


def run_cooccurrence(args):
    """
    Execute the co-occurrence analysis pipeline.
    """
    (
        null_ing,
        filt_ing,
        taxa_universe,
        out_dir,
        focal_query_to_taxa,
        rhs_query_to_taxa,
    ) = run_shared_pipeline_setup(args)

    if null_ing is None:
        return

    print(
        "Pipeline: Cooccurrence analysis being performed with Null Ingredients Presence Matrix with "
        f"{null_ing.presence_matrix.shape[0]} taxa & {null_ing.presence_matrix.shape[1]} samples"
    )

    # print("DEBUG focal_query_to_taxa =", focal_query_to_taxa)
    # print("DEBUG rhs_query_to_taxa =", rhs_query_to_taxa)

    edge_arrays, nodes_df = cooccurrence_obj(
        null_ing,
        taxa_universe,
        large=args.large,
        max_pairs=args.max_pairs,
        threshold=args.threshold,
        null_model=args.null_model,
        nm_n_reps=args.nm_n_reps,
        nm_seed=args.nm_seed,
        nm_n_workers=getattr(args, "nm_n_workers", None),
        nm_mp_start=getattr(args, "nm_mp_start", None),
        nm_sort_indices=getattr(args, "nm_sort_indices", False),
        nm_burn_in_steps=getattr(args, "nm_burn_in_steps", None),
        nm_steps_per_rep=getattr(args, "nm_steps_per_rep", None),
        nm_progress_every=getattr(args, "nm_progress_every", 25),
        focal_query_to_taxa=focal_query_to_taxa,
        rhs_query_to_taxa=rhs_query_to_taxa,
    )

    null_scope_prefix = "global" if args.null_scope is None else str(args.null_scope)

    export_cooccurrence_outputs(
        edge_arrays=edge_arrays,
        nodes_df=nodes_df,
        taxa_universe=taxa_universe,
        output_dir=out_dir,
        edges_base=f"{args.tag}{null_scope_prefix}_edges",
        nodes_base=f"{args.tag}{null_scope_prefix}_nodes",
        null_model=args.null_model,
        summary_n=100_000,
    )


def run_biome_distribution(args):
    """
    Determine and export the biome distribution of taxa across all annotated metagenomes.
    """
    os.makedirs(args.output_dir, exist_ok=True)

    ingredients = load_ingredients(args.data_dir, args.aggregated, args.custom_ingredients, args.data_version)
    ingredients, is_successful = filter_data_obj(
        ingredients,
        min_taxa_count=getattr(args, "min_taxa_count", None),
        min_sample_count=getattr(args, "min_sample_count", None),
        filter_rank=getattr(args, "filter_rank", None),
        taxa_count_rank=getattr(args, "taxa_count_rank", "species"),
    )
    if not is_successful:
        return

    biome_level = getattr(args, "biome_level", "level_1")
    biomes, presence, coverage, n_dropped = ingredients.biome_distribution(level=biome_level)
    biome_by_taxa_df = pd.DataFrame(data=presence.todense(), columns=ingredients.taxa, index=biomes)

    taxa_query = getattr(args, "taxa_query", None)
    if taxa_query:
        indices = _taxa_indices_from_query(ingredients, taxa_query)
        biome_by_query_df = biome_by_taxa_df.iloc[:, indices].T
        output_path = os.path.join(args.output_dir, f"{args.tag}taxa_biome_distribution.tsv")
        biome_by_query_df.to_csv(output_path, sep="\t")

    elif args.return_all_taxa:
        output_path = os.path.join(args.output_dir, f"{args.tag}taxa_biome_distribution.tsv")
        biome_by_taxa_df.to_csv(output_path, sep="\t")

    elif args.aggregated:
        if not [i for i in biome_by_taxa_df.columns if "AGGREGATED" in i]:
            print("WARNING: Ingredients did not contain aggregated taxa. Only species will be output")
        indices = [i for i, v in enumerate(biome_by_taxa_df.columns) if "s__" in v or "AGGREGATED" in v]
        biome_by_agg_df = biome_by_taxa_df.iloc[:, indices].T
        output_path = os.path.join(args.output_dir, f"{args.tag}taxa_biome_distribution.tsv")
        biome_by_agg_df.to_csv(output_path, sep="\t")

    else:
        indices = [i for i, v in enumerate(biome_by_taxa_df.columns) if "s__" in v]
        biome_by_species_df = biome_by_taxa_df.iloc[:, indices].T
        output_path = os.path.join(args.output_dir, f"{args.tag}taxa_biome_distribution_species.tsv")
        biome_by_species_df.to_csv(output_path, sep="\t")
