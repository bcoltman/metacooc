#!/usr/bin/env python3
"""
pipelines.py

This module defines pipeline functions for metacooc.

The primary function is:
    run_cooccurrence(args)

This function runs the full in-memory co-occurrence pipeline, which includes:
  1. Loading an Ingredients object (aggregated or raw) from args.data_dir.
  2. Performing a search (taxon or metadata mode) to obtain matching accessions.
  3. Filtering the Ingredients object to create:
         - A null object filtered solely by count thresholds.
         - A second filtered object based on both count thresholds and the matching accessions.
  4. Calculating ratios by comparing the filtered object to the null.
         Optionally, if args.ratio_threshold is provided, threshold filtering is applied.
  5. Saving the ratios (as a TSV file) and plotting them.

If file paths for ingredients or metadata indices are not explicitly provided, default files in args.data_dir are used.
"""

import os
import numpy as np
import pandas as pd

from metacooc.search import search_data_obj, resolve_focal_taxa_queries
from metacooc.filter import filter_data_obj
from metacooc.pantry import load_ingredients
from metacooc.analysis import (
    association_obj,
    cooccurrence_obj,
    select_taxa_universe,
    export_cooccurrence_outputs,
)
# from metacooc.analysis import association_obj, cooccurrence_obj, select_taxa_universe
from metacooc.plot import plot_analysis_obj
from metacooc.clustering import determine_taxa_context
from metacooc.structure import structure_obj

def _samples_from_taxon_rows(ingredients, row_idx: np.ndarray) -> set[str]:
    """
    Resolve a set of taxon row indices to the set of samples containing >=1 of them.
    """
    row_idx = np.asarray(row_idx, dtype=np.int64)
    if row_idx.size == 0:
        return set()
        
    hit = np.asarray(
        ingredients.presence_matrix[row_idx, :].getnnz(axis=0)
    ).ravel() > 0
    
    sample_arr = np.asarray(ingredients.samples, dtype=object)
    return set(sample_arr[hit].tolist())

def _resolve_focal_mode_initial_accessions(args, ingredients):
    """
    Resolve the initial accession cohort for focal_taxa mode directly from focal
    taxon row resolution, rather than via search_data_obj().
    
    Returns
    -------
    matching_accessions : set[str]
    focal_query_to_taxa_raw : dict[str, list[str]]
        Raw focal resolution against the unfiltered loaded Ingredients object.
        This will be re-resolved again after filtering to build the final
        focal_query_to_taxa used by cooccurrence.
    """
    focal_rows = resolve_focal_taxa_queries(ingredients, args.search_string)
    
    taxa_arr = np.asarray(ingredients.taxa, dtype=object)
    focal_query_to_taxa_raw = {}
    matching_accessions = set()
    
    for focal_query, row_idx in focal_rows.items():
        row_idx_arr = np.asarray(sorted(row_idx), dtype=np.int64)
        if row_idx_arr.size == 0:
            continue
            
        resolved_taxa = taxa_arr[row_idx_arr].tolist()
        focal_query_to_taxa_raw[focal_query] = resolved_taxa
        matching_accessions.update(_samples_from_taxon_rows(ingredients, row_idx_arr))
        
    return matching_accessions, focal_query_to_taxa_raw


def run_shared_pipeline_setup(args):
    """
    Shared setup for cooccurrence, association, and structure pipelines.
    
    Returns
    -------
    tuple
        (null_ingredients, filtered_ingredients, taxa_universe, output_dir, focal_query_to_taxa)
    """
    ingredients = load_ingredients(
        args.data_dir,
        args.aggregated,
        args.custom_ingredients,
        args.data_version,
    )
    
    focal_query_to_taxa = None
    
    if args.search_mode == "focal_taxa":
        matching_accessions, _ = _resolve_focal_mode_initial_accessions(args, ingredients)
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
        )
        
    if not matching_accessions:
        print("Pipeline: No matching accessions found. Exiting pipeline.")
        return None, None, None, None, None
        
    print(f"Pipeline: Found {len(matching_accessions)} matching accessions after initial filtering.")
    
    int_ingredients, is_successful = filter_data_obj(
        ingredients,
        accession_set=matching_accessions,
        min_taxa_count=args.min_taxa_count,
        min_sample_count=args.min_sample_count,
        filter_rank=args.filter_rank,
        taxa_count_rank=args.taxa_count_rank,
    )
    if not is_successful:
        return None, None, None, None, None
        
    sub_samples = int_ingredients.samples
    taxa_universe = select_taxa_universe(int_ingredients, rank=args.filter_rank)
    
    if args.search_mode == "focal_taxa":
        focal_rows = resolve_focal_taxa_queries(int_ingredients, args.search_string)
        taxa_arr = np.asarray(int_ingredients.taxa, dtype=object)
        taxa_universe_set = set(taxa_universe)
        
        focal_query_to_taxa = {}
        for focal_query, row_idx in focal_rows.items():
            row_idx_arr = np.asarray(sorted(row_idx), dtype=np.int64)
            if row_idx_arr.size == 0:
                continue
                
            resolved_taxa = [
                taxon
                for taxon in taxa_arr[row_idx_arr].tolist()
                if taxon in taxa_universe_set
            ]
            if resolved_taxa:
                focal_query_to_taxa[focal_query] = resolved_taxa
                
        if not focal_query_to_taxa:
            raise ValueError("No focal taxa remain after sample/count/rank filtering.")
            
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
            return None, None, None, None, None
            
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
            return None, None, None, None, None
            
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
            return None, None, None, None, None
            
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
            return None, None, None, None, None
            
        null_ingredients = determine_taxa_context(
            null_ingredients,
            focal_taxa=args.null_taxa_query,
            degree=args.taxa_degree,
            min_shared_samples_between_taxa=args.min_shared_samples_between_taxa,
        )
    else:
        raise ValueError(f"Unknown null_scope: {args.null_scope!r}")
        
    filtered_ingredients, is_successful = filter_data_obj(
        null_ingredients,
        accession_set=sub_samples,
    )
    if not is_successful:
        return None, None, None, None, None
        
    if (
        getattr(args, "command", None) == "association"
        and set(filtered_ingredients.samples) == set(null_ingredients.samples)
    ):
        raise ValueError(
            "Term and null cohorts are identical — association requires a broader null.\n"
            "Fix by widening null_biome/null_scope or narrowing search_string."
        )
        
    os.makedirs(args.output_dir, exist_ok=True)
    return null_ingredients, filtered_ingredients, taxa_universe, args.output_dir, focal_query_to_taxa



def run_structure(args):
    """
    Executes the association analysis pipeline using pre-processed ingredients data.
    
    The pipeline calculates either taxon or metadata enrichment metrics, depending on the specified mode.
    Results are saved as tab-separated files, and visualisations are generated.
    
    Args:
        args (argparse.Namespace): Command-line arguments or equivalent object containing:
            - tag (str): Prefix for output file names.
            - mode (str): Analysis mode, either "taxon" or "metadata".
            - filter_rank (str, optional): Taxonomic rank at which to apply filters.
            - large (bool): If True, optimises for large datasets.
            - max_pairs (int, optional): Maximum number of taxon pairs to consider.
            - threshold (float, optional): Threshold for filtering association metrics.
            
    Outputs:
        Saves the following files to the specified output directory:
            - If mode is "taxon":
                - {tag}taxon_enrichment.tsv: Tab-separated file of taxon enrichment results.
                - {tag}taxon_plot.png: Visualisation of the taxon enrichment analysis.
            - If mode is "metadata":
                - {tag}{mode}_enrichment.tsv: Tab-separated file of metadata enrichment results.
                - {tag}{mode}_plot.png: Visualisation of the metadata enrichment analysis.
    """
    
    
    null_ing, filt_ing, taxa_universe, out_dir, _ = run_shared_pipeline_setup(args)
    
    if null_ing is None:  # Early exit if setup failed
        return
    
    print(f"Pipeline: Structure analysis being performed on Filtered Ingredients Presence Matrix with "
          f"{filt_ing.presence_matrix.shape[0]} taxa & {filt_ing.presence_matrix.shape[1]} samples")
    
    structure_df = structure_obj(
            filt_ing,
            null_model=args.null_model,
            nm_n_reps=args.nm_n_reps,
            nm_random_state=args.nm_random_state
        )
    
    null_scope_prefix = "global" if args.null_scope is None else str(args.null_scope)
    output_path = os.path.join(out_dir, f"{args.tag}{null_scope_prefix}_structure.tsv")
    structure_df.to_csv(output_path, sep="\t", index=False)
    print(f"Pipeline: {null_scope_prefix} structure analysis saved to {output_path}")
    

def run_association(args):
    """
    Executes the association analysis pipeline using pre-processed ingredients data.
    """
    
    null_ing, filt_ing, taxa_universe, out_dir, _ = run_shared_pipeline_setup(args)
    
    if null_ing is None:  # Early exit if setup failed
        return
    
    print(f"Pipeline: Association analysis being performed with Null Ingredients Presence Matrix with "
          f"{null_ing.presence_matrix.shape[0]} taxa & {null_ing.presence_matrix.shape[1]} samples")
      
    single_df = association_obj(
            null_ing,
            filt_ing,
            threshold=args.threshold,
            null_model=args.null_model,
            nm_n_reps=args.nm_n_reps,
            compute_fisher=args.compute_fisher,
            nm_random_state=args.nm_random_state
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
    null_ing, filt_ing, taxa_universe, out_dir, focal_query_to_taxa = run_shared_pipeline_setup(args)
    
    if null_ing is None:
        return
        
    print(
        "Pipeline: Cooccurrence analysis being performed with Null Ingredients Presence Matrix with "
        f"{null_ing.presence_matrix.shape[0]} taxa & {null_ing.presence_matrix.shape[1]} samples"
    )
    
    edge_arrays, nodes_df = cooccurrence_obj(
        null_ing,
        taxa_universe,
        large=args.large,
        max_pairs=args.max_pairs,
        threshold=args.threshold,
        null_model=args.null_model,
        nm_n_reps=args.nm_n_reps,
        nm_random_state=args.nm_random_state,
        focal_query_to_taxa=focal_query_to_taxa,
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
    Determine and export the biome dsitribution of taxa across all annotated metagenomes
    """
    
    ingredients = load_ingredients(args.data_dir, args.aggregated, args.custom_ingredients, args.data_version)
    biomes, presence, coverage, n_dropped = ingredients.biome_distribution()
    biome_by_taxa_df = pd.DataFrame(data=presence.todense(), columns=ingredients.taxa, index=biomes)
    
    if args.return_all_taxa:
        output_path = os.path.join(args.output_dir, f"{args.tag}taxa_biome_distribution.tsv")
        # biome_by_taxa_df.T.to_csv(output_path, sep="\t")
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