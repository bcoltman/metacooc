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

def load_ingredients(data_dir=None, filtered_file=None, reference_file=None):
    """
    Load the filtered and reference Ingredients objects.

    Parameters:
        data_dir (str): Directory to search for default files if file paths are not provided.
        filtered_file (str): Path to the filtered Ingredients pickle.
        reference_file (str): Path to the reference Ingredients pickle.

    Returns:
        tuple: (reference, filtered) Ingredients objects.
    """
    if filtered_file:
        if not os.path.exists(filtered_file):
            raise FileNotFoundError(f"Filtered Ingredients file '{filtered_file}' not found.")
    elif data_dir:
        filtered_file = os.path.join(data_dir, "ingredients_all_filtered.pkl")
    else:
        raise ValueError("Either 'data_dir' or 'filtered_file' must be provided.")
        
    if reference_file:
        if not os.path.exists(reference_file):
            raise FileNotFoundError(f"Reference Ingredients file '{reference_file}' not found.")
    elif data_dir:
        reference_file = os.path.join(data_dir, "ingredients_counts_filtered.pkl")
    else:
        raise ValueError("Either 'data_dir' or 'reference_file' must be provided.")
    
    with open(filtered_file, "rb") as f:
        filtered = pickle.load(f)
    with open(reference_file, "rb") as f:
        reference = pickle.load(f)
    return reference, filtered

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
    
    if ratio_threshold is not None:
        filtered_df = ratios_df[ratios_df['ratio'] >= ratio_threshold]
        return ratios_df, filtered_df
    else:
        return ratios_df

def calculate_ratios(output_dir, data_dir=None, filtered_file=None, reference_file=None,
                     ratio_threshold=None, tag=None):
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
    reference, filtered = load_ingredients(data_dir, filtered_file, reference_file)
    if ratio_threshold is None:
        ratios_df = calculate_ratios_obj(filtered, reference)
        
        output_path = os.path.join(output_dir, f"ratios{tag}.tsv")
        ratios_df.to_csv(output_path, sep="\t", index=False)
        print(f"Ratios saved to {output_path}")
        
    else:
        print(f"Applied threshold filtering: ratio >= {ratio_threshold}.")
        ratios_df, filtered_df = calculate_ratios_obj(filtered, reference, ratio_threshold)
        
        output_path = os.path.join(output_dir, f"ratios{tag}.tsv")
        ratios_df.to_csv(output_path, sep="\t", index=False)
        print(f"Ratios saved to {output_path}")
        
        output_path = os.path.join(output_dir, f"filtered_ratios{tag}.tsv")
        ratios_df.to_csv(output_path, sep="\t", index=False)
        print(f"Filtered ratios saved to {output_path}")
