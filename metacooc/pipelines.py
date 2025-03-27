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
import pickle

from metacooc.search import search_data_obj
from metacooc.filter import filter_data_obj
from metacooc.ratios import calculate_ratios_obj
from metacooc.plot import plot_ratios_obj

def run_cooccurrence(args):
    """
    Run the full co-occurrence pipeline in memory.

    Expected attributes in args:
      - data_dir, output_dir, aggregated (boolean)
      - mode: "taxon" or "metadata"
      - search_string, rank (for taxon mode)
      - index_file (optional, for metadata mode), strict
      - min_taxa_count, min_sample_count
      - ratio_threshold (optional)

    The pipeline:
      1. Loads the Ingredients object from args.data_dir. If args.aggregated is True, loads
         "ingredients_aggregated.pkl"; otherwise, loads "ingredients_raw.pkl".
      2. Runs a search using search_data_obj():
             - For taxon mode, the Ingredients object is searched (optionally restricting by rank).
             - For metadata mode, a default metadata index is loaded from args.data_dir (broad or broad_strict)
               unless args.index_file is provided.
      3. Prints the number of matching accessions.
      4. Applies count-based filtering via filter_data_obj() to create two objects:
             - reference_ingredients: filtered only by min_taxa_count and min_sample_count.
             - filtered_ingredients: additionally filtered by the matching accessions.
      5. Calculates ratios using calculate_ratios_obj(). If args.ratio_threshold is provided,
         threshold filtering is applied.
      6. Saves the ratios DataFrame as "ratios.tsv" in args.output_dir and plots the results.
    """
    # Step 1. Load the Ingredients object.
    
    if args.aggregated:
        ingredients_file = os.path.join(args.data_dir, "ingredients_aggregated_genus.pkl")
    else:
        ingredients_file = os.path.join(args.data_dir, "ingredients_raw.pkl")
    if not os.path.exists(ingredients_file):
        raise FileNotFoundError(f"Ingredients file '{ingredients_file}' not found.")
    with open(ingredients_file, "rb") as f:
        ingredients = pickle.load(f)
    
    # Step 2. Perform search.
    matching_accessions = search_data_obj(args.mode, 
                                          args.data_dir, 
                                          args.search_string, 
                                          args.rank, 
                                          args.strict, 
                                          args.column_names,
                                          args.inverse)
                                          
    print(f"Pipeline: Found {len(matching_accessions)} matching accessions.")

    # Step 3. Filter the Ingredients object.
    # Create a 'reference' object filtered solely by count thresholds.
    reference_ingredients = filter_data_obj(
        ingredients,
        accession_set=None,
        min_taxa_count=args.min_taxa_count,
        min_sample_count=args.min_sample_count
    )
    # Create a 'filtered' object further filtered by matching accessions.
    filtered_ingredients = filter_data_obj(
        ingredients,
        accession_set=matching_accessions,
        min_taxa_count=args.min_taxa_count,
        min_sample_count=args.min_sample_count
    )

    # Step 4. Calculate ratios.
    if args.ratio_threshold is None:
        ratios_df = calculate_ratios_obj(filtered_ingredients, reference_ingredients)
        output_path = os.path.join(args.output_dir, f"ratios{args.tag}.tsv")
        ratios_df.to_csv(output_path, sep="\t", index=False)
        print(f"Pipeline: Ratios saved to {output_path}")
    else:
        print(f"Pipeline: Applying threshold filtering: ratio >= {args.ratio_threshold}.")
        # Assume calculate_ratios_obj() returns a tuple when threshold filtering is applied.
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