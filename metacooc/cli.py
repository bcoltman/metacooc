#!/usr/bin/env python3

import argparse
# from metacooc.format import format_data
# from metacooc.download import download_data
from metacooc.search import search_data
# from metacooc.filter import filter_data
# from metacooc.ratios import calculate_ratios
# from metacooc.plot import plot_ratios

def parse_cli():
    parser = argparse.ArgumentParser(description="Co-occurrence data of microorganisms based on metagenome detection")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand for 'cooccurrence'
    cooc_parser = subparsers.add_parser("cooccurrence", help="run the full co-occurrence workflow")

    # Subcommand for 'format'
    format_parser = subparsers.add_parser("format", help="format data for co-occurrence analysis")

    # Subcommand for 'download'
    download_parser = subparsers.add_parser("download", help="download required data for co-occurrence analysis")

    # Subcommand for 'search'
    search_parser = subparsers.add_parser("search", help="search metadata and taxonomy files for the occurrence of terms")

    # Subcommand for 'filter'
    filter_parser = subparsers.add_parser("filter", help="filter data by search term or list of accession numbers")

    # Subcommand for 'ratios'
    ratios_parser = subparsers.add_parser("ratios", help="calculate co-occurrence ratios for all taxa in the filtered set")

    # Subcommand for 'plot'
    plot_parser = subparsers.add_parser("plot", help="plot co-occurrence ratios")

    # Add arguments for 'cooccurrence'
    cooc_parser.add_argument("blast_tab", help="BLAST/DIAMOND tabular output file")
    cooc_parser.add_argument("metadata", help="Metadata file")
    cooc_parser.add_argument("taxonomy", help="Taxonomy file")
    cooc_parser.add_argument("output_dir", help="Output directory")
    cooc_parser.add_argument("--data_dir", help="Directory to download and set up data")

    # Add arguments for 'format'
    format_parser.add_argument("--data_dir", help="Directory to format data")

    # Add arguments for 'download'
    download_parser.add_argument("--data_dir", help="Directory to download data")

    # Add arguments for 'search'
    search_parser.add_argument("--mode", "-m", choices=["taxon", "metadata"], required=True, help="Search mode")
    search_parser.add_argument("--strict", action="store_true", help="Enable strict search mode")
    search_parser.add_argument("--rank", help="Taxonomic rank to search")
    search_parser.add_argument("output_dir", help="Output directory for search results")
    search_parser.add_argument("search_string", help="Search string")

    # Add arguments for 'filter'
    filter_parser.add_argument("accessions_file", help="File containing accession numbers")
    filter_parser.add_argument("--data_dir", help="Directory containing data to filter")
    filter_parser.add_argument("--output_dir", help="Output directory for filtered data")

    # Add arguments for 'ratios'
    ratios_parser.add_argument("--filtered", help="Filtered file to calculate ratios")
    ratios_parser.add_argument("--data_dir", help="Directory containing data to calculate ratios")
    ratios_parser.add_argument("--output_dir", help="Output directory for ratios")

    # Add arguments for 'plot'
    plot_parser.add_argument("--ratios", help="Ratios file to plot")
    plot_parser.add_argument("--output_dir", help="Output directory for plot")

    args = parser.parse_args()

    if args.command == "cooccurrence":
        format_data(args.data_dir)
        download_data(args.data_dir)
        search_data(args.mode, args.strict, args.rank, args.output_dir, args.search_string)
        filter_data(args.accessions_file, args.data_dir, args.output_dir)
        calculate_ratios(args.filtered, args.data_dir, args.output_dir)
        plot_ratios(args.ratios, args.output_dir)
    elif args.command == "format":
        format_data(args.data_dir)
    elif args.command == "download":
        download_data(args.data_dir)
    elif args.command == "search":
        search_data(args.mode, args.strict, args.rank, args.output_dir, args.search_string)
    elif args.command == "filter":
        filter_data(args.accessions_file, args.data_dir, args.output_dir)
    elif args.command == "ratios":
        calculate_ratios(args.filtered, args.data_dir, args.output_dir)
    elif args.command == "plot":
        plot_ratios(args.ratios, args.output_dir)

if __name__ == "__main__":
    parse_cli()
