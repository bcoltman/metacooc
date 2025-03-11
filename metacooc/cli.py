#!/usr/bin/env python3
import argparse
import os
import pickle
from pathlib import Path

# make data directory in case it is needed
with Path("metacooc") as data_path:
    DEFAULT_DATA_DIR = str(data_path)+"/data"

if not os.path.exists(DEFAULT_DATA_DIR):
    os.makedirs(DEFAULT_DATA_DIR)

def parse_cli():
    # ------------------ Parent Parsers ------------------
    # Common arguments for all commands.
    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument("--data_dir", default=DEFAULT_DATA_DIR,
        help="Directory containing data files (default: %(default)s)")
    common_parser.add_argument("--output_dir", required=True,
        help="Directory where output files will be saved.")
    common_parser.add_argument("--tag", default="",
        help="Optional tag to append to output filenames for distinction.")


    # Arguments for formatting (used by 'format' and optionally in pipeline defaults).
    format_parent = argparse.ArgumentParser(add_help=False)
    format_parent.add_argument("--tax_profile",
        help="Taxonomic profile TSV file.")
    format_parent.add_argument("--metadata_file",
        help="Metadata TSV file.")
    format_parent.add_argument("--aggregation_pattern", choices=["Root", "d__", "p__", "c__", "o__", "f__", "g__", "s__"],
        help="Taxonomic level up to which to aggregate (e.g., 'g__').")

    # Arguments for metadata index (used by 'format' and 'search').
    metadata_common = argparse.ArgumentParser(add_help=False)
    metadata_common.add_argument("--index_type", choices=["broad", "exact"], default="broad",
        help="Index type for metadata search (default: 'broad').")
    metadata_common.add_argument("--strict", action="store_true",
        help="For metadata searches in broad mode, restrict index to a reduced set of columns.")
    metadata_common.add_argument("--column", help="Column name for exact index (if required).")
    metadata_common.add_argument("--aggregated", action="store_true",
        help="Flag to indicate that the aggregated Ingredients object should be used.")


    # Arguments for search.
    search_parent = argparse.ArgumentParser(add_help=False)
    search_parent.add_argument("--mode", "-m", choices=["taxon", "metadata"], required=True,
        help="Search mode: 'taxon' or 'metadata'.")
    search_parent.add_argument("--search_string", required=True,
        help="Search string to query.")
    search_parent.add_argument("--rank", help="Taxonomic rank filter (only for taxon mode).")
    search_parent.add_argument("--index_file",
        help="(Optional) Path to a pre-built metadata index (pickle).")

    # Arguments for filter.
    filter_parent = argparse.ArgumentParser(add_help=False)
    filter_parent.add_argument("--accessions_file",
        help="File containing accession numbers (one per line) to filter samples.")
    filter_parent.add_argument("--min_taxa_count", type=int,
        help="Minimum number of taxa a sample must have to be included.")
    filter_parent.add_argument("--min_sample_count", type=int,
        help="Minimum number of samples in which a taxon must be present.")

    # Arguments for ratio.
    ratio_parent = argparse.ArgumentParser(add_help=False)
    ratio_parent.add_argument("--filtered_file",
        help="(Optional) Path to the filtered Ingredients pickle file.")
    ratio_parent.add_argument("--reference_file",
        help="(Optional) Path to the reference Ingredients pickle file.")
    ratio_parent.add_argument("--ratio_threshold", type=float, default=0.5,
        help="Minimum ratio value to keep.")

    # ------------------ Main Parser and Subcommands ------------------
    parser = argparse.ArgumentParser(
        description="Co-occurrence data of microorganisms based on metagenome detection")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # 'download' subcommand (format and download run separately).
    download_cmd = subparsers.add_parser("download", parents=[common_parser],
        help="Download initial files (raw and aggregated Ingredients, metadata indices).")
    download_cmd.add_argument("--force", action="store_true",
        help="Force re-download even if files exist.")

    # 'format' subcommand (file-based formatting).
    format_cmd = subparsers.add_parser("format", parents=[common_parser, format_parent, metadata_common],
        help="Format data to generate Ingredients objects and metadata indices (file-based).")

    # 'search' subcommand (file-based).
    search_cmd = subparsers.add_parser("search", parents=[common_parser, metadata_common, search_parent],
        help="Perform a file-based search (using default files from data_dir).")

    # 'filter' subcommand (file-based).
    filter_cmd = subparsers.add_parser("filter", parents=[common_parser, metadata_common, filter_parent],
        help="Filter data by accession numbers (file-based).")

    # 'ratio' subcommand (file-based).
    ratio_cmd = subparsers.add_parser("ratio", parents=[common_parser, ratio_parent],
        help="Calculate co-occurrence ratios (file-based).")

    # 'plot' subcommand (file-based).
    plot_cmd = subparsers.add_parser("plot", parents=[common_parser, ratio_parent],
        help="Plot co-occurrence ratios (file-based).")
    # plot_cmd.add_argument("--ratios", required=True,
        # help="Path to the ratios file to plot.")

    # 'cooccurrence' subcommand: the full pipeline (in-memory).
    cooc_parser = subparsers.add_parser("cooccurrence", parents=[
        common_parser, metadata_common, search_parent, filter_parent, ratio_parent
    ], help="Run the full co-occurrence workflow (in-memory pipeline: search, filter, ratio, plot).")
    # In pipeline mode, the ingredients and indices are loaded from defaults if not explicitly provided.
    # (For example, if --filtered_file is not provided, the pipeline will use data_dir/ingredients_all_filtered.pkl,
    # and if --reference_file is not provided, it will use data_dir/ingredients_counts_filtered.pkl.)

    args = parser.parse_args()
    
    args.tag = f"_{args.tag}" if args.tag else ""
    

    # ------------------ Dispatch ------------------
    if args.command == "download":
        from metacooc.download import download_data
        download_data(data_dir=args.data_dir, 
                      force=args.force,
                      tag=args.tag)
    elif args.command == "format":
        from metacooc.format import format_data
        format_data(tax_profile=args.tax_profile, 
                    output_dir=args.output_dir, 
                    metadata_file=args.metadata_file, 
                    index_type=args.index_type, 
                    strict=args.strict, 
                    column=args.column, 
                    aggregated=args.aggregated, 
                    aggregation_pattern=args.aggregation_pattern,
                    tag=args.tag)
    elif args.command == "search":
        from metacooc.search import search_data
        search_data(mode=args.mode, 
                    data_dir=args.data_dir, 
                    output_dir=args.output_dir, 
                    search_string=args.search_string,
                    rank=args.rank, 
                    index_type=args.index_type, 
                    column=args.column,
                    strict=args.strict, 
                    index_file=args.index_file,
                    tag=args.tag)
    elif args.command == "filter":
        from metacooc.filter import filter_data
        filter_data(accessions_file=args.accessions_file, 
                    data_dir=args.data_dir, 
                    output_dir=args.output_dir,
                    aggregated=args.aggregated, 
                    min_taxa_count=args.min_taxa_count, 
                    min_sample_count=args.min_sample_count,
                    tag=args.tag)
    elif args.command == "ratio":
        from metacooc.ratios import calculate_ratios
        calculate_ratios(output_dir=args.output_dir, 
                         data_dir=args.data_dir,
                         filtered_file=args.filtered_file, 
                         reference_file=args.reference_file, 
                         ratio_threshold=args.ratio_threshold,
                         tag=args.tag)
    elif args.command == "plot":
        from metacooc.plot import plot_ratios
        plot_ratios(ratios=args.ratios, 
                    output_dir=args.output_dir, 
                    ratio_threshold=args.ratio_threshold,
                    tag=args.tag)
    elif args.command == "cooccurrence":
        from metacooc.pipelines import run_cooccurrence
        run_cooccurrence(args)

if __name__ == "__main__":
    parse_cli()
