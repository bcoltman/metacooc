#!/usr/bin/env python3
"""
plot.py

Plot co-occurrence ratios for taxa.

This module provides two interfaces:

1. plot_ratios_obj(ratios_df, output_plot_file=None, threshold=None)
   - Accepts an in-memory pandas DataFrame (ratios_df) and plots the ratios.
   - If output_plot_file is provided, the figure is saved; otherwise, it is displayed.

2. plot_ratios(ratios_file, output_dir, threshold=None)
   - File-based interface: loads the ratios TSV file (generated by calculate_ratios),
     calls plot_ratios_obj(), and saves the plot in output_dir (or displays it).

Usage (CLI):
    metacooc plot --ratios path/to/ratios.tsv --output_dir /path/to/out [--threshold 0.5]
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def plot_ratios_obj(ratios_df, output_plot_file, ratio_threshold=None):
    """
    Create a plot of taxon ratios from a ratios DataFrame.

    Parameters:
        ratios_df (pd.DataFrame): DataFrame with columns 'taxon', 'filtered_counts', 
                                  'reference_counts', and 'ratio'.
        output_plot_file (str, optional): If provided, save the figure to this file.
        threshold (float, optional): If provided, a vertical line will be drawn where the ratio
                                     falls below this threshold.
    """
    ratios_df = ratios_df[ratios_df['ratio'] != 0]
    
    # Sort the DataFrame by ratio in descending order.
    ratios_df_sorted = ratios_df.sort_values(by='ratio', ascending=False).reset_index(drop=True)

    # Create the plot.
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(ratios_df_sorted['ratio'].values, marker='o', markersize=5)
    ax.set_xlabel('Taxa (ordered by ratio)')
    ax.set_ylabel('Ratio')
    ax.set_title('Taxon Ratios')
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True)

    # Mark the threshold if provided.
    if ratio_threshold is not None:
        # Identify the last index where ratio is >= threshold.
        threshold_indices = ratios_df_sorted[ratios_df_sorted['ratio'] >= ratio_threshold].index
        if len(threshold_indices) > 0:
            cutoff_index = threshold_indices[-1]
            ax.axvline(x=cutoff_index, color='red', linestyle='--', label=f'Threshold = {ratio_threshold}')
            ax.legend()
        else:
            print(f"No taxa found with ratio >= {threshold}")

    fig.tight_layout()
    fig.savefig(output_plot_file)
    plt.close(fig)
    print(f"Plot saved to {output_plot_file}")

def plot_ratios(ratios_file, output_dir, ratio_threshold=None, tag=None):
    """
    File-based plotting function.

    Loads a ratios TSV file from disk, calls plot_ratios_obj() to create the plot,
    and saves the figure in output_dir as "ratios_plot.png".

    Parameters:
        ratios_file (str): Path to the ratios TSV file.
        output_dir (str): Directory where the plot image will be saved.
        threshold (float, optional): Threshold for marking the plot.
    """
    if not os.path.exists(ratios_file):
        raise FileNotFoundError(f"Ratios file '{ratios_file}' not found.")
    ratios_df = pd.read_csv(ratios_file, sep="\t")
    output_plot_file = os.path.join(output_dir, f"ratios_plot{tag}.png")
    plot_ratios_obj(ratios_df, output_plot_file=output_plot_file, ratio_threshold=ratio_threshold)