#!/usr/bin/env python3
import argparse
import os
import pickle
from pathlib import Path

# Get the directory where this script is located
BASE_DIR = Path(__file__).resolve().parent

# Build the data directory path relative to the script directory
DEFAULT_DATA_DIR = BASE_DIR / "data"

# Create the directory if it doesn't exist
DEFAULT_DATA_DIR.mkdir(parents=True, exist_ok=True)


def parse_cli():
    
    # ----------------------
    # Common Parent Parsers
    # ----------------------

    # For commands that require data directory.
    data_parser = argparse.ArgumentParser(add_help=False)
    data_parser.add_argument("--data_dir", default=DEFAULT_DATA_DIR,
                             help="Directory containing data files (default: %(default)s)")
                             
    # For commands that require an output directory.
    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("--output_dir", required=True,
                               help="Directory where output files will be saved.")
                               
    # For commands that use an optional tag.
    tag_parser = argparse.ArgumentParser(add_help=False)
    tag_parser.add_argument("--tag", default="",
                            help="Optional tag to append to output filenames for distinction.")
                            
    # For search-related options.
    search_parser = argparse.ArgumentParser(add_help=False)
    search_parser.add_argument("--mode", choices=["taxon", "metadata"], required=True,
                               help="Search mode: 'taxon' or 'metadata'.")
    search_parser.add_argument("--search_string", required=True,
                               help="Search string to query.")
    search_parser.add_argument("--rank", help="Taxonomic rank filter (only for taxon mode).")
    search_parser.add_argument("--column_names", help="Column name for exact index (if required).")
    search_parser.add_argument("--strict", action="store_true",
                               help="Restrict search to a reduced set of columns.")
    search_parser.add_argument("--inverse", action="store_true",
                               help="Return inverse of search string.")
                               
    # For filter-related options.
    filter_parser = argparse.ArgumentParser(add_help=False)
    filter_parser.add_argument("--min_taxa_count", type=int,
                               help="Minimum number of taxa a sample must have to be included.")
    filter_parser.add_argument("--min_sample_count", type=int,
                               help="Minimum number of samples in which a taxon must be present.")
    filter_parser.add_argument("--aggregated", action="store_true",
                            help="Use the aggregated Ingredients object.")
                            
    # For ratio-related options.
    ratio_parser = argparse.ArgumentParser(add_help=False)
    ratio_parser.add_argument("--ratio_threshold", type=float, default=0.5,
                              help="Minimum ratio value to keep.")
                              
    # ----------------------
    # Subcommand Handlers
    # ----------------------
    def download_command(args):
        from metacooc.download import download_data
        download_data(data_dir=args.data_dir, force=args.force)

    def format_command(args):
        from metacooc.format import format_data
        format_data(tax_profile=args.tax_profile,
                    output_dir=args.output_dir,
                    aggregated=args.aggregated,
                    aggregation_pattern=args.aggregated_pattern)

    def search_command(args):
        from metacooc.search import search_data
        search_data(mode=args.mode,
                    data_dir=args.data_dir,
                    output_dir=args.output_dir,
                    search_string=args.search_string,
                    rank=args.rank,
                    column_names=args.column_names,
                    strict=args.strict,
                    inverse=args.inverse,
                    tag=args.tag)

    def filter_command(args):
        from metacooc.filter import filter_data
        filter_data(accessions=args.accessions,
                    data_dir=args.data_dir,
                    output_dir=args.output_dir,
                    aggregated=args.aggregated,
                    min_taxa_count=args.min_taxa_count,
                    min_sample_count=args.min_sample_count,
                    tag=args.tag)

    def ratio_command(args):
        from metacooc.ratios import calculate_ratios
        calculate_ratios(output_dir=args.output_dir,
                         data_dir=args.data_dir,
                         filtered_file=args.filtered_file,
                         reference_file=args.reference_file,
                         ratio_threshold=args.ratio_threshold,
                         tag=args.tag)

    def plot_command(args):
        from metacooc.plot import plot_ratios
        plot_ratios(ratios_file=args.ratios_file,
                    output_dir=args.output_dir,
                    ratio_threshold=args.ratio_threshold,
                    tag=args.tag)

    def cooccurrence_command(args):
        from metacooc.pipelines import run_cooccurrence
        run_cooccurrence(args)
        
    parser = argparse.ArgumentParser(
        description="Co-occurrence data of microorganisms based on metagenome detection")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Download
    download_sub = subparsers.add_parser("download",
                                          parents=[data_parser],
                                          help="Download initial files.")
    download_sub.add_argument("--force", action="store_true",
                              help="Force re-download even if files exist.")
    download_sub.set_defaults(func=download_command)

    # Format
    format_sub = subparsers.add_parser("format",
                                       parents=[output_parser],
                                       help="Format data to generate Ingredients objects.")
    format_sub.add_argument("--tax_profile", help="Taxonomic profile TSV file.")
    format_sub.add_argument("--aggregated", action="store_true",
                            help="Use the aggregated Ingredients object.")
    format_sub.add_argument("--aggregated_pattern",
                            help="Taxonomic level up to which to aggregate (e.g., 'g__').")
    format_sub.set_defaults(func=format_command)

    # Search
    search_sub = subparsers.add_parser("search",
                                       parents=[data_parser, output_parser, tag_parser, search_parser],
                                       help="Perform a file-based search.")
    search_sub.set_defaults(func=search_command)

    # Filter
    filter_sub = subparsers.add_parser("filter",
                                       parents=[data_parser, output_parser, tag_parser, filter_parser],
                                       help="Filter data by accession numbers.")
    filter_sub.add_argument("--accessions", help="File containing accession numbers to filter samples.")
    filter_sub.set_defaults(func=filter_command)

    # Ratio
    ratio_sub = subparsers.add_parser("ratio",
                                      parents=[data_parser, output_parser, tag_parser, ratio_parser],
                                      help="Calculate co-occurrence ratios.")
    ratio_sub.add_argument("--filtered_file", help="Path to the filtered Ingredients pickle file.")
    ratio_sub.add_argument("--reference_file", help="Path to the reference Ingredients pickle file.")
    ratio_sub.set_defaults(func=ratio_command)

    # Plot
    plot_sub = subparsers.add_parser("plot",
                                     parents=[output_parser, tag_parser, ratio_parser],
                                     help="Plot co-occurrence ratios.")
    plot_sub.add_argument("--ratios_file", required=True, help="Path to the ratios file to plot.")
    plot_sub.set_defaults(func=plot_command)

    # Cooccurrence (full pipeline; note: no output_dir or tag)
    cooc_sub = subparsers.add_parser("cooccurrence",
                                     parents=[data_parser, output_parser, tag_parser, search_parser, filter_parser, ratio_parser],
                                     help="Run the full co-occurrence workflow (in-memory).")
    cooc_sub.set_defaults(func=cooccurrence_command)

    args = parser.parse_args()
    # Normalize tag if provided.
    if hasattr(args, "tag"):
        args.tag = f"_{args.tag}" if args.tag else ""
    args.func(args)

if __name__ == "__main__":
    parse_cli()