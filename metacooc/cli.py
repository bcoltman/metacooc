#!/usr/bin/env python3


import argparse
# from metacooc.filter import xyz

def parse_cli():
    parser = argparse.ArgumentParser(description="Co-occurence data of microorganims based on metagenome detection")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Subcommand for 'setup'
    setup_parser = subparsers.add_parser("setup", help="download and set up required data co-occurence analysis")
# example line to add arguments
#   pasr_parser.add_argument("blast_tab", help="BLAST/DIAMOND tabular output file")

    # Subcommand for 'cooccurrence'
    cooc_parser = subparsers.add_parser("cooccurrence", help="run the full co-occurence workflow")

    # Subcommand for 'filter'
    filter_parser = subparsers.add_parser("filter", help="filter data by search term or list of accession numbers")

    # Subcommand for 'ratio'
    ratio_parser = subparsers.add_parser("ratio", help="calculate cooccurrence ratios for all taxa in the filtered set")

    # Subcommand for 'search'
    search_parser = subparsers.add_parser("search", help="search metadata and taxonomy files for the occurrence of terms")

    # Subcommand for 'list'
    list_parser = subparsers.add_parser("list", help="get a list of accession numbers based on filtering criteria")


    args = parser.parse_args()

    if args.command == "setup":
        # do stuff
    elif args.command == "cooccurrence":
        # do stuff
    elif args.command == "filter":
        # do stuff
    elif args.command == "ratio":
        # do stuff
    elif args.command == "search":
        # do stuff
    elif args.command == "list":
        # do stuff





