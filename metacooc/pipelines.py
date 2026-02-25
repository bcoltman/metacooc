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
import pandas as pd

from metacooc.search import search_data_obj
from metacooc.filter import filter_data_obj
from metacooc.pantry import load_ingredients
from metacooc.analysis import association_obj, cooccurrence_obj, select_taxa_universe
from metacooc.plot import plot_analysis_obj
from metacooc.clustering import determine_taxa_context
from metacooc.structure import structure_obj

def run_shared_pipeline_setup(args):
    """
    Performs the shared setup steps for both co-occurrence and association pipelines.
    
    This function handles the loading of data, accession searching, and count-based filtering.
    It is designed to be called by both `run_cooccurrence` and `run_association` to avoid code duplication.
    
    Args:
        args (argparse.Namespace): Command-line arguments or equivalent object containing:
            - data_dir (str): Directory containing input data files.
            - aggregated (bool): If True, loads aggregated ingredients; otherwise, loads raw ingredients.
            - custom_ingredients (object, optional): Pre-loaded ingredients object.
            - data_version (str, optional): Version of the Sandpiper data used for data processing.
            - mode (str): Analysis mode, either "taxon" or "metadata".
            - search_string (str): Search term for filtering accessions.
            - ranks_for_search_inclusion (list, optional): Taxonomic ranks to restrict search.
            - strict (bool): If True, applies strict matching criteria.
            - column_names (list, optional): Column names for metadata mode.
            - inverse (bool): If True, inverts the search logic.
            - min_taxa_count (int): Minimum number of taxa required for inclusion.
            - min_sample_count (int): Minimum number of samples required for inclusion.
            - filter_rank (str, optional): Taxonomic rank at which to apply filters.
            - output_dir (str): Directory for saving output files.
            - remove_null_threshold (bool): If true, then remove min_taxa_count and min_sample_count threshold on null model
            - null_scope (str): Which scope to define null model from
            - null_biome_query (str): Biomes to include in the null. Default None is all
            - null_taxa_query (str): Taxa to include in local search. Default None includes all.
            - taxa_degree (int): Neighbourhood radius measured in taxa→sample expansions. Default 1 includes only those specified in null_taxa_query.
            - min_shared_samples_between_taxa (int):  Minimum number of samples in which two taxa must co-occur for a new taxon to be included during null scope expansion.
            
    Returns:
        tuple: A tuple containing the following elements:
            - null_ingredients (object): Ingredients object filtered by count thresholds.
            - filtered_ingredients (object): Ingredients object further filtered by matching accessions.
            - matching_accessions (set): Set of accessions matching the search criteria.
            - output_dir (str): Path to the directory for saving output files.
            If any step fails, returns (None, None, None, None).
    """
    # Step 1. Load the Ingredients object.
    ingredients = load_ingredients(args.data_dir, args.aggregated, args.custom_ingredients, args.data_version)
    
    # Step 2. Perform search.
    matching_accessions = search_data_obj(
        search_mode=args.search_mode,
        data_dir=args.data_dir,
        search_string=args.search_string,
        ranks_for_search_inclusion=args.ranks_for_search_inclusion,
        strict=args.strict,
        column_names=args.column_names,
        inverse=args.inverse,
        custom_ingredients=ingredients
    )
    
    if not matching_accessions:
        print("Pipeline: No matching accessions found. Exiting pipeline.")
        return None, None, None, None
        
    print(f"Pipeline: Found {len(matching_accessions)} matching accessions after initial filtering.")
    
    int_ingredients, is_successful = filter_data_obj(
        ingredients,
        accession_set=matching_accessions,
        min_taxa_count=args.min_taxa_count,
        min_sample_count=args.min_sample_count,
        filter_rank=args.filter_rank,
        taxa_count_rank=args.taxa_count_rank
    )
    if not is_successful:
        return None, None, None, None
    
    sub_samples = int_ingredients.samples
    
    taxa_universe = select_taxa_universe(int_ingredients, rank=args.filter_rank)
    
    
    if args.null_scope is None:
        
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=None, 
                                              min_taxa_count=0 if args.remove_null_threshold else args.min_taxa_count, 
                                              min_sample_count=0 if args.remove_null_threshold else args.min_sample_count, 
                                              filter_rank=args.filter_rank,
                                              taxa_count_rank=args.taxa_count_rank)
                                              
        if not is_successful:
            return None, None, None, None
        
    elif args.null_scope == "taxa":
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=None, 
                                              min_taxa_count=0 if args.remove_null_threshold else args.min_taxa_count, 
                                              min_sample_count=0 if args.remove_null_threshold else args.min_sample_count, 
                                              filter_rank=args.filter_rank,
                                              taxa_count_rank=args.taxa_count_rank)
        
        if not is_successful:
            return None, None, None, None
                                              
        null_ingredients = determine_taxa_context(null_ingredients,
                                       focal_taxa=args.null_taxa_query,
                                       degree=args.taxa_degree,
                                       min_shared_samples_between_taxa=args.min_shared_samples_between_taxa)
        
    elif args.null_scope == "biome" or args.null_scope == "metadata":
        search_string = args.null_biome_query if args.null_scope == "biome" else args.null_metadata_query
        
        null_matching_accessions = search_data_obj(search_mode=args.null_scope,
                                              search_string=search_string,
                                              custom_ingredients=ingredients)
        
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=null_matching_accessions, 
                                              min_taxa_count=0 if args.remove_null_threshold else args.min_taxa_count, 
                                              min_sample_count=0 if args.remove_null_threshold else args.min_sample_count, 
                                              filter_rank=args.filter_rank,
                                              taxa_count_rank=args.taxa_count_rank)
        
        if not is_successful:
            return None, None, None, None
    
    elif args.null_scope == "biome_taxa" or args.null_scope == "metadata_taxa" :
        
        search_mode = "biome" if args.null_scope == "biome_taxa" else "metadata"
        search_string = args.null_biome_query if search_mode == "biome" else args.null_metadata_query
        
        null_matching_accessions = search_data_obj(search_mode=search_mode,
                                              search_string=search_string,
                                              custom_ingredients=ingredients)
        
        null_ingredients, is_successful = filter_data_obj(ingredients, 
                                              accession_set=null_matching_accessions, 
                                              min_taxa_count=0 if args.remove_null_threshold else args.min_taxa_count, 
                                              min_sample_count=0 if args.remove_null_threshold else args.min_sample_count, 
                                              filter_rank=args.filter_rank,
                                              taxa_count_rank=args.taxa_count_rank)
        
        if not is_successful:
            return None, None, None, None
            
        null_ingredients = determine_taxa_context(null_ingredients,
                                       focal_taxa=args.null_taxa_query,
                                       degree=args.taxa_degree,
                                       min_shared_samples_between_taxa=args.min_shared_samples_between_taxa)
                                       
    
    filtered_ingredients, is_successful = filter_data_obj(null_ingredients, 
                                                          accession_set=sub_samples)
    
    if set(filtered_ingredients.samples) == set(null_ingredients.samples):
        raise ValueError(
            "Term and null cohorts are identical — association requires a broader null.\n"
            "Fix by widening null_biome/null_scope or narrowing search_string."
        )
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    return null_ingredients, filtered_ingredients, taxa_universe, args.output_dir



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
    
    
    null_ing, filt_ing, taxa_universe, out_dir = run_shared_pipeline_setup(args)
    
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
    
    null_ing, filt_ing, taxa_universe, out_dir = run_shared_pipeline_setup(args)
    
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
    Executes the co-occurrence analysis pipeline using pre-processed ingredients data.
    
    The pipeline calculates taxon co-occurrence metrics and saves the results as tab-separated files.
    Visualisations of the analysis are generated and saved as image files.
    
    Args:
        args (argparse.Namespace): Command-line arguments or equivalent object containing:
            - tag (str): Prefix for output file names.
            - mode (str): Analysis mode, either "taxon" or "metadata".
            - filter_rank (str, optional): Taxonomic rank at which to apply filters.
            - large (bool): If True, optimises for large datasets.
            - max_pairs (int, optional): Maximum number of taxon pairs to consider.
            - threshold (float, optional): Threshold for filtering co-occurrence ratios.
            - SIM9 (bool): If True, shuffles Null using SIM9 to determine Null distribution of J
            
    Outputs:
        Saves the following files to the specified output directory:
            - {tag}{mode}_edges.tsv: Tab-separated file of taxon co-occurrence edges.
            - {tag}{mode}_nodes.tsv: Tab-separated file of taxon co-occurrence nodes.
            - {tag}{mode}_plot.png: Visualisation of the co-occurrence analysis.
    """
    
    # null_ing, taxa_universe, matching_accs, out_dir = run_shared_pipeline_setup(args)
    # null_ing, filt_ing, matching_accs, out_dir = run_shared_pipeline_setup(args)
    null_ing, filt_ing, taxa_universe, out_dir = run_shared_pipeline_setup(args)
    
    if null_ing is None:  # Early exit if setup failed
        return
    
    print(f"Pipeline: Cooccurrence analysis being performed with Null Ingredients Presence Matrix with "
          f"{null_ing.presence_matrix.shape[0]} taxa & {null_ing.presence_matrix.shape[1]} samples")
      
    edges_df, nodes_df = cooccurrence_obj(
        null_ing,
        taxa_universe,
        large=args.large,
        max_pairs=args.max_pairs,
        threshold=args.threshold,
        null_model=args.null_model,
        nm_n_reps=args.nm_n_reps,
        nm_random_state=args.nm_random_state
    )
    # filter_rank=args.filter_rank,
    if edges_df is not None:
        null_scope_prefix = "global" if args.null_scope is None else str(args.null_scope)
        output_path = os.path.join(out_dir, f"{args.tag}{null_scope_prefix}_edges.tsv")
        edges_df.to_csv(output_path, sep="\t", index=False)
        print(f"Pipeline: Taxon edges analysis saved to {output_path}")
        
        output_path = os.path.join(out_dir, f"{args.tag}{null_scope_prefix}_nodes.tsv")
        nodes_df.to_csv(output_path, sep="\t", index=False)
        print(f"Pipeline: Taxon nodes analysis saved to {output_path}")


def run_biome_distribution(args):
    """
    Determine and export the biome dsitribution of taxa across all annotated metagenomes
    """
    
    ingredients = load_ingredients(args.data_dir, args.aggregated, args.custom_ingredients, args.data_version)
    biomes, presence, coverage, n_dropped = ingredients.biome_distribution()
    biome_by_taxa_df = pd.DataFrame(data=presence.todense(), columns=ingredients.taxa, index=biomes)
    
    if args.return_all_taxa:
        output_path = os.path.join(args.output_dir, f"{args.tag}taxa_biome_distribution.tsv")
        biome_by_taxa_df.T.to_csv(output_path, sep="\t")
    
    elif args.aggregated:
        if not [i for i in biome_by_taxa_df.columns if "AGGREGATED" in i]:
            print("WARNING: Ingredients did not contain aggregated taxa. Only species will be output")
        indices = [i for i, v in enumerate(biome_by_taxa_df.columns) if "s__" in v or "AGGREGATED" in v]
        biome_by_agg_df = biome_by_taxa_df.iloc[:, indices].T
        output_path = os.path.join(args.output_dir, f"{args.tag}taxa_biome_distribution{"_aggregated" if args.aggregated else ""}.tsv")
        biome_by_agg_df.T.to_csv(output_path, sep="\t")
    
    else:
        indices = [i for i, v in enumerate(biome_by_taxa_df.columns) if "s__" in v]
        biome_by_species_df = biome_by_taxa_df.iloc[:, indices].T
        output_path = os.path.join(args.output_dir, f"{args.tag}taxa_biome_distribution_species.tsv")
        biome_by_species_df.to_csv(output_path, sep="\t")