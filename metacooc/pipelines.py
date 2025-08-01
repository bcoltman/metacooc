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
         - A reference object filtered solely by count thresholds.
         - A second filtered object based on both count thresholds and the matching accessions.
  4. Calculating ratios by comparing the filtered object to the reference.
         Optionally, if args.ratio_threshold is provided, threshold filtering is applied.
  5. Saving the ratios (as a TSV file) and plotting them.

If file paths for ingredients or metadata indices are not explicitly provided, default files in args.data_dir are used.
"""

import os
import pandas as pd

from metacooc.search import search_data_obj
from metacooc.filter import filter_data_obj
from metacooc.ratios import calculate_ratios_obj
from metacooc.plot import plot_ratios_obj
from metacooc.pantry import load_ingredients

def run_cooccurrence(args):
    """
    Run the full co-occurrence pipeline in memory.
    
    Expected attributes in args:
      - data_dir, output_dir, aggregated (boolean)
      - mode: "taxon" or "metadata"
      - search_string, ranks_for_search_inclusion (for taxon mode)
      - index_file (optional, for metadata mode), strict
      - min_taxa_count, min_sample_count
      - ratio_threshold (optional)
    
    The pipeline:
      1. Loads the Ingredients object from args.data_dir. If args.aggregated is True, loads
         "ingredients_aggregated.pkl"; otherwise, loads "ingredients_raw.pkl".
      2. Runs a search using search_data_obj():
             - For taxon mode, the Ingredients object is searched (optionally restricting by ranks_for_search_inclusion).
             - For metadata mode, a default metadata index is loaded from args.data_dir (broad or broad_strict)
               unless args.index_file is provided.
             - For biome mode, filters based on the pre-defined biomes specified by manual allocation of accessions 
               to broad metadata categories
      3. Prints the number of matching accessions.
      4. Applies count-based filtering via filter_data_obj() to create two objects:
             - reference_ingredients: filtered only by min_taxa_count and min_sample_count.
             - filtered_ingredients: additionally filtered by the matching accessions.
      5. Calculates ratios using calculate_ratios_obj(). If args.ratio_threshold is provided,
         threshold filtering is applied.
      6. Saves the ratios DataFrame as "ratios.tsv" in args.output_dir and plots the results.
    """
    # Step 1. Load the Ingredients object.
    
    ingredients = load_ingredients(args.data_dir, args.aggregated, args.custom_ingredients, args.sandpiper_version)
    
    # Step 2. Perform search. Pass the preloaded ingredients as a custom_ingredients
    matching_accessions = search_data_obj(mode=args.mode, 
                                          data_dir=args.data_dir, 
                                          search_string=args.search_string, 
                                          ranks_for_search_inclusion=args.ranks_for_search_inclusion, 
                                          strict=args.strict, 
                                          column_names=args.column_names,
                                          inverse=args.inverse,
                                          custom_ingredients=ingredients)
                                          
    if not matching_accessions:
        print("Pipeline: No matching accessions found. Exiting pipeline.")
        return
    
    print(f"Pipeline: Found {len(matching_accessions)} matching accessions.")
    
    # Step 3. Filter the Ingredients object.
    # Create a 'reference' object filtered solely by count thresholds.
    reference_ingredients = filter_data_obj(
        ingredients,
        accession_set=None,
        min_taxa_count=args.min_taxa_count,
        min_sample_count=args.min_sample_count,
        filter_rank=args.filter_rank
    )
    # Create a 'filtered' object further filtered by matching accessions.
    filtered_ingredients = filter_data_obj(
        ingredients,
        accession_set=matching_accessions,
        min_taxa_count=args.min_taxa_count,
        min_sample_count=args.min_sample_count,
        filter_rank=args.filter_rank
    )

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    # Step 4. Calculate ratios.
    if args.ratio_threshold is None:
        ratios_df = calculate_ratios_obj(filtered_ingredients, reference_ingredients)
        output_path = os.path.join(args.output_dir, f"ratios{args.tag}.tsv")
        ratios_df.to_csv(output_path, sep="\t", index=False)
        print(f"Pipeline: Ratios saved to {output_path}")
    else:
        print(f"Pipeline: Applying threshold filtering: ratio >= {args.ratio_threshold}.")
        ratios_df, filtered_ratios_df = calculate_ratios_obj(filtered_ingredients, reference_ingredients, args.ratio_threshold)
        
        output_path = os.path.join(args.output_dir, f"ratios{args.tag}.tsv")
        ratios_df.to_csv(output_path, sep="\t", index=False)
        print(f"Pipeline: Ratios saved to {output_path}")
        
        filtered_path = os.path.join(args.output_dir, f"filtered_ratios{args.tag}.tsv")
        filtered_ratios_df.to_csv(filtered_path, sep="\t", index=False)
        print(f"Pipeline: Filtered ratios saved to {filtered_path}")
    
    # Step 5. Plot ratios.
    output_plot_file = os.path.join(args.output_dir, f"ratios_plot{args.tag}.png")
    plot_ratios_obj(ratios_df, output_plot_file=output_plot_file, ratio_threshold=args.ratio_threshold)
    print("Pipeline: Plotting complete.")

def run_biome_distribution(args):
    """
    Determine and export the biome dsitribution of taxa across all annotated metagenomes
    """
    
    ingredients = load_ingredients(args.data_dir, args.aggregated, args.custom_ingredients, args.sandpiper_version)
    biomes, presence, coverage, n_dropped = ingredients.biome_distribution()
    biome_by_taxa_df = pd.DataFrame(data=presence.todense(), columns=ingredients.taxa, index=biomes)
    
    if args.return_all_taxa:
        output_path = os.path.join(args.output_dir, f"taxa_biome_distribution{args.tag}.tsv")
        biome_by_taxa_df.T.to_csv(output_path, sep="\t")
    
    elif args.aggregated:
        if not [i for i in biome_by_taxa_df.columns if "agg" in i]:
            print("WARNING: Ingredients did not contain aggregated taxa. Only species will be output")
        indices = [i for i, v in enumerate(biome_by_taxa_df.columns) if "s__" in v or "AGG" in v]
        biome_by_agg_df = biome_by_taxa_df.iloc[:, indices].T
        output_path = os.path.join(args.output_dir, f"taxa_biome_distribution{"_aggregated" if args.aggregated else ""}{args.tag}.tsv")
        biome_by_agg_df.T.to_csv(output_path, sep="\t")
    
    else:
        indices = [i for i, v in enumerate(biome_by_taxa_df.columns) if "s__" in v]
        biome_by_species_df = biome_by_taxa_df.iloc[:, indices].T
        output_path = os.path.join(args.output_dir, f"taxa_biome_distribution{"_aggregated" if args.aggregated else ""}{args.tag}_species.tsv")
        biome_by_species_df.to_csv(output_path, sep="\t")