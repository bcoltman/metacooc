#!/usr/bin/env python3
"""
ratio.py

Calculate co-occurrence ratios for taxa by comparing two Ingredients objects:
a filtered Ingredients object and a reference Ingredients object.

Both Ingredients objects are assumed to have a cached attribute `total_counts`
computed at initialization. For each taxon common to both objects, the ratio is computed as:
    ratio = (filtered_count) / (reference_count)
with proper handling for taxa with zero reference counts.

This module also provides threshold filtering:
    filter_ratios_obj(ratios_df, ratio_threshold, counts_threshold)
which filters the ratios DataFrame to retain only taxa with a ratio >= ratio_threshold.

Usage (file-based):
    metacooc ratio --filtered path/to/ingredients_filtered.pkl \
       --reference path/to/ingredients_reference.pkl --output_dir /path/to/out \
       [--ratio_threshold 0.5 ]
"""

import os
import argparse
import pickle
import numpy as np
import pandas as pd

from metacooc.pantry import load_ingredients 


def calculate_ratios_obj(filtered, reference, ratio_threshold=None):
    """
    In-memory ratio calculation.

    For each taxon common to both filtered and reference, compute:
        ratio = filtered_count / reference_count
    using the cached total_counts attribute.

    Returns a DataFrame with columns: taxon, filtered_counts, reference_counts, and ratio.
    """
    filt_taxa = np.array(filtered.taxa)
    ref_taxa = np.array(reference.taxa)
    common_taxa, filt_idx, ref_idx = np.intersect1d(filt_taxa, ref_taxa, return_indices=True)
    if common_taxa.size == 0:
        raise ValueError("No common taxa found between filtered and reference datasets.")
    filtered_counts = filtered.total_counts[filt_idx]
    reference_counts = reference.total_counts[ref_idx]
    ratios = np.divide(
        filtered_counts,
        reference_counts,
        out=np.zeros_like(filtered_counts, dtype=float),
        where=reference_counts != 0
    )
    ratios_df = pd.DataFrame({
        'taxon': common_taxa,
        'filtered_counts': filtered_counts,
        'reference_counts': reference_counts,
        'ratio': ratios
    })
    
    ratios_df = ratios_df.sort_values(by='ratio', ascending=False)
    
    if ratio_threshold is not None:
        
        filtered_df = ratios_df[ratios_df['ratio'] >= ratio_threshold]
        return ratios_df, filtered_df
    else:
        return ratios_df

def calculate_ratios(reference_ingredients,
                     filtered_ingredients, 
                     output_dir, 
                     ratio_threshold=None, 
                     tag=None):
    """
    File-based ratio calculation.

    Loads the Ingredients objects (either from specified file paths or from defaults in data_dir),
    computes ratios, optionally applies threshold filtering, and saves the ratios DataFrame as TSV.

    Parameters:
        output_dir (str): Directory where ratios.tsv will be saved.
        data_dir (str): Directory to locate default Ingredients files if file paths are not provided.
        filtered_file (str): Path to the filtered Ingredients pickle.
        reference_file (str): Path to the reference Ingredients pickle.
        ratio_threshold (float, optional): Minimum ratio value to keep.
        counts_threshold (int, optional): Minimum reference count to keep.

    Returns:
        pd.DataFrame: The (optionally filtered) ratios DataFrame.
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
    reference = load_ingredients(data_dir=None, custom_ingredients=reference_ingredients)
    filtered = load_ingredients(data_dir=None, custom_ingredients=filtered_ingredients)
    # , filtered = load_ingredients(data_dir, filtered_file, reference_file)
    if ratio_threshold is None:
        ratios_df = calculate_ratios_obj(filtered, reference)
        
        output_path = os.path.join(output_dir, f"{tag}ratios.tsv")
        ratios_df.to_csv(output_path, sep="\t", index=False)
        print(f"Ratios saved to {output_path}")
        
    else:
        print(f"Applied threshold filtering: ratio >= {ratio_threshold}.")
        ratios_df, filtered_df = calculate_ratios_obj(filtered, reference, ratio_threshold)
        
        output_path = os.path.join(output_dir, f"{tag}ratios.tsv")
        ratios_df.to_csv(output_path, sep="\t", index=False)
        print(f"Ratios saved to {output_path}")
        
        output_path = os.path.join(output_dir, f"{tag}filtered_ratios.tsv")
        filtered_df.to_csv(output_path, sep="\t", index=False)
        print(f"Filtered ratios saved to {output_path}")