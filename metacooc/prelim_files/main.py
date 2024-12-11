#!/usr/bin/env python3

# main.py
import argparse
import os
import pandas as pd
from loader import download_table
from searcher import load_table1, process_tables, download_table
from searcher.__version__ import __version__

def main():
    parser = argparse.ArgumentParser(description="Search tables and process results.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--word1', type=str, help="Word to search in Table 1")
    group.add_argument('--word2', type=str, help="Word to search in Table 2")
    parser.add_argument('--sandpiper', type=str, default='data/table1.tsv', help="Path to Table 1 file")
    parser.add_argument('--sra_metadata', type=str, required=True, help="Path to Table 2 file")
    parser.add_argument('--output_dir', type=str, required=True, help="Output directory")
    parser.add_argument('--cutoff', type=int, default=10, help="Cutoff for minimal value in 'count' column")
    parser.add_argument('--chunksize', type=int, default=10000, help="Number of rows per chunk for processing Table 2")
    parser.add_argument('--download_table2', type=str, help="URL to download Table 2 if not provided locally")
    parser.add_argument('--version', action='store_true', help="Show version information and exit")
    
    args = parser.parse_args()
    
    if args.version:
        print(f"Table Searcher version {__version__}")
        return
    
    if args.sandpiper_table:
        if not os.path.exists(args.sandpiper_table):
            download_table(https://zenodo.org/api/records/11516218/files/sandpiper0.3.0.condensed.csv.gz/content, args.sandpiper)
        else:
            print(f"Sandpiper table already exists at {args.sandpiper}. Skipping download.")
    
    if args.sra_metadata:
        if not os.path.exists(args.sra_metadata):
            download_table(args.download_table2, sra_metadata)
        else:
            print(f"Table 2 already exists at {args.sra_metadata}. Skipping download.")
    
    if not os.path.exists(args.table2):
        print(f"Table 2 not found at {args.table2}. Please provide a valid path or URL to download.")
        return
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Load the tables
    sandpiper_df = pd.read_csv(args.sandpiper, sep="\t")
    
    sra_metadata_df = pd.read_csv(args.sra_metadata, sep="\t")
    
    # Determine the search word and table
    if args.word1:
        search_word = args.word1
        search_in_table1 = True
    else:
        search_word = args.word2
        search_in_table1 = False
    
    # Process the tables
    complete_merged, filtered_merged = process_tables(table1, args.table2, search_word, args.cutoff, search_in_table1, args.chunksize)
    
    # Write the results to files
    complete_merged_path = os.path.join(args.output_dir, "complete_merged.csv")
    filtered_merged_path = os.path.join(args.output_dir, "filtered_merged.csv")
    
    complete_merged.to_csv(complete_merged_path, index=False)
    filtered_merged.to_csv(filtered_merged_path, index=False)
    
    print(f"Complete merged table written to: {complete_merged_path}")
    print(f"Filtered merged table written to: {filtered_merged_path}")

if __name__ == "__main__":
    main()
