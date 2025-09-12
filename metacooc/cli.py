#!/usr/bin/env python3
import argparse
from pathlib import Path

# Get the directory where this script is located
BASE_DIR = Path(__file__).resolve().parent
# Build the data directory path relative to the script directory
DEFAULT_DATA_DIR = BASE_DIR / "data"
# Create the data directory if it doesn't exist
DEFAULT_DATA_DIR.mkdir(parents=True, exist_ok=True)


def positive_int(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not an integer")
    
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"{value} is not a positive integer")
    return ivalue

def validate_ratio_threshold(value):
    try:
        value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not a valid float.")
    
    if value < 0.0 or value > 1.0:
        raise argparse.ArgumentTypeError("Ratio threshold must be between 0 and 1.")
    return value

def parse_cli():
    parser = argparse.ArgumentParser(
        description="Co-occurrence data of microorganisms based on metagenome detection"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    # ----------------------------
    # DOWNLOAD SUBCOMMAND
    # ----------------------------
    download_sub = subparsers.add_parser("download", help="Download initial files.")
    
    # Optional arguments group
    opt = download_sub.add_argument_group("optional arguments")
    
    opt.add_argument(
        "--data_dir",
        default=DEFAULT_DATA_DIR,
        help="Directory containing data files (default: %(default)s)",
    )
    opt.add_argument(
        "--list-versions",
        action="store_true",
        help="List available version",
    )
    opt.add_argument(
        "--sandpiper_version",
        default=None,
        help="Specify which data version to load (default: latest). Versions available for download can be listed with 'metacooc download --list_versions'"
    )
    opt.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if files exist.",
    )
    
    def download_command(args):
        from metacooc.download import download_data
        
        download_data(data_dir=args.data_dir, force=args.force, list_versions=args.list_versions, sandpiper_version=args.sandpiper_version)
    
    download_sub.set_defaults(func=download_command)
    
    # ----------------------------
    # FORMAT SUBCOMMAND
    # ----------------------------
    format_sub = subparsers.add_parser("format", help="Format data to generate Ingredients objects.")
    
    # Required arguments group
    req = format_sub.add_argument_group("required arguments")
    req.add_argument(
        "--output_dir",
        required=True,
        help="Directory where output files will be saved.",
    )
    req.add_argument(
        "--tax_profile",
        required=True,
        help="Taxonomic profile TSV file.",
    )
    # Optional arguments group
    opt = format_sub.add_argument_group("optional arguments")
    
    opt.add_argument(
        "--sample_to_biome_file",
        help="A CSV file linking SRA accessions to biome classifications.",
    )
    opt.add_argument(
        "--tag",
        help="Optional tag to prepend to output filenames for distinction.",
    )
    opt.add_argument(
        "--aggregated",
        action="store_true",
        help="Use the aggregated Ingredients object.",
    )
    
    def format_command(args):
        from metacooc.format import format_data
        
        format_data(
            tax_profile=args.tax_profile,
            output_dir=args.output_dir,
            sample_to_biome_file=args.sample_to_biome_file,
            aggregated=args.aggregated,
            tag=args.tag
        )
    
    format_sub.set_defaults(func=format_command)
    
    # ----------------------------
    # SEARCH SUBCOMMAND
    # ----------------------------
    search_sub = subparsers.add_parser("search", help="Perform a file-based search.")
    
    # Required arguments group
    req = search_sub.add_argument_group("required arguments (unless --list_column_names is used)")
    req.add_argument(
        "--mode",
        choices=["taxon", "metadata", "biome"],
        help="Search mode: 'taxon', 'metadata' or 'biome'.",
    )
    req.add_argument(
        "--search_string",
        type=str,
        help=(
            "Search string to query as a single token. "
            "Use '|' to separate OR-terms and '+' to separate AND-terms. "
            "Examples: 'foo|bar', 'foo+baz', 'foo|bar+baz|qux'. "
            "Search string must be wrapped in single or double quotes if it "
            "contains a special character or space e.g. 's__Escherichia coli'."
        ),
    )
    req.add_argument(
        "--output_dir",
        help="Directory where output files will be saved.",
    )
    
    opt = search_sub.add_argument_group("optional arguments")
    opt.add_argument(
    "--ranks_for_search_inclusion",
    choices=["domain", "phylum", "class", "order", "family", "genus", "species"],
    help=(
        "Taxa identified at a rank higher than this rank are excluded "
        "in taxon search mode."
    ),
    )
    opt.add_argument(
        "--data_dir",
        default=DEFAULT_DATA_DIR,
        help="Directory containing data files (default: %(default)s)",
    )
    opt.add_argument(
        "--tag",
        default="",
        help="Optional tag to prepend to output filenames for distinction.",
    )
    opt.add_argument(
        "--column_names",
        nargs='+',
        help="meta mode only: Restrict metadata search to within specified columns. Multiple entries can be specified as e.g. --column_names 'hello' 'world'",
    )
    opt.add_argument(
        "--custom_ingredients",
        help="Path to an Ingredients file to use instead of default. This Ingredients will be used, regardless of whether data_dir is specified. "
    )
    opt.add_argument(
        "--sandpiper_version",
        default=None,
        help="Specify which data version to load (default: latest). Versions available for download can be listed with 'metacooc download --list_versions'"
    )
    opt.add_argument(
        "--strict",
        action="store_true",
        help="meta mode only: Restrict metadata search to within a pre-defined and reduced set of columns i.e. ['acc', 'organism', 'env_biome_sam', 'env_feature_sam','env_material_sam', 'biosamplemodel_sam']",
    )
    opt.add_argument(
        "--inverse",
        action="store_true",
        help="meta mode only: Return inverse of search e.g. if using search_term 'soil' in metadata mode, then it will return all entries without soil.",
    )
    
    # Standalone argument for listing column names
    search_sub.add_argument(
        "--list_column_names",
        action="store_true",
        help="meta mode only: WARNING: Will produce lots of output! List available column names from NCBI metadata",
    )
    
    def search_command(args):
        from metacooc.search import search_data
        
        # If --list_column_names is specified, skip checking for required arguments
        if not args.list_column_names:
            if not args.mode or not args.search_string or not args.output_dir:
                search_sub.error("The following arguments are required unless --list_column_names is used: --mode, --search_string, --output_dir")
                
        if args.aggregated:
            args.tag = f"{args.tag}_aggregated_" if args.tag else "aggregated_"
        else:
            args.tag = f"{args.tag}_" if args.tag else ""
        
        search_data(
            mode=args.mode,
            data_dir=args.data_dir,
            output_dir=args.output_dir,
            search_string=args.search_string,
            ranks_for_search_inclusion=args.ranks_for_search_inclusion,
            column_names=args.column_names,
            strict=args.strict,
            inverse=args.inverse,
            tag=args.tag,
            custom_ingredients=args.custom_ingredients,
            sandpiper_version=args.sandpiper_version,
            list_column_names=args.list_column_names,
        )
        
    search_sub.set_defaults(func=search_command)
    
    # ----------------------------
    # FILTER SUBCOMMAND
    # ----------------------------
    filter_sub = subparsers.add_parser("filter", help="Filter data by accession numbers or other criteria.")

    # Required arguments group
    req = filter_sub.add_argument_group("required arguments")
    req.add_argument(
        "--output_dir",
        required=True,
        help="Directory where output files will be saved.",
    )

    # Optional arguments group
    opt = filter_sub.add_argument_group("optional arguments")
    opt.add_argument(
        "--data_dir",
        default=DEFAULT_DATA_DIR,
        help="Directory containing data files (default: %(default)s)",
    )
    opt.add_argument(
        "--tag",
        default="",
        help="Optional tag to prepend to output filenames for distinction.",
    )
    opt.add_argument(
        "--min_taxa_count",
        type=positive_int,
        help="Minimum number of taxa a sample must have to be included.",
    )
    opt.add_argument(
        "--min_sample_count",
        type=positive_int,
        help="Minimum number of samples in which a taxon must be present.",
    )
    opt.add_argument(
        "--accessions_file",
        help="File containing accession numbers to filter by.",
    )
    opt.add_argument(
        "--filter_rank",
        choices=["domain", "phylum", "class", "order", "family", "genus", "species"],
        help=(
            "Taxa identified at a rank higher than this rank are filtered out "
            "of results."
        ),
    )
    opt.add_argument(
        "--custom_ingredients",
        help="Ingredients file to use instead of default",
    )
    opt.add_argument(
        "--sandpiper_version",
        default=None,
        help="Specify which data version to load (default: latest). Versions available for download can be listed with 'metacooc download --list_versions'",
    )
    opt.add_argument(
        "--aggregated",
        action="store_true",
        help="Use the aggregated Ingredients object.",
    )
    
    def filter_command(args):
        from metacooc.filter import filter_data
        
        # Validate that at least one filtering option is provided
        if not any([args.min_taxa_count, args.min_sample_count, args.accessions_file, args.filter_rank]):
            filter_sub.error("At least one of the following arguments is required: --min_taxa_count, --min_sample_count, --accessions_file, or --filter_rank")
        
        if args.aggregated:
            args.tag = f"{args.tag}_aggregated_" if args.tag else "aggregated_"
        else:
            args.tag = f"{args.tag}_" if args.tag else ""
        
        filter_data(
            accessions_file=args.accessions_file,
            data_dir=args.data_dir,
            output_dir=args.output_dir,
            aggregated=args.aggregated,
            min_taxa_count=args.min_taxa_count,
            min_sample_count=args.min_sample_count,
            filter_rank=args.filter_rank,
            tag=args.tag,
            custom_ingredients=args.custom_ingredients,
            sandpiper_version=args.sandpiper_version,
        )
    
    filter_sub.set_defaults(func=filter_command)
    
    # ----------------------------
    # RATIO SUBCOMMAND
    # ----------------------------
    ratio_sub = subparsers.add_parser("ratio", help="Calculate co-occurrence ratios.")
    
    req = ratio_sub.add_argument_group("required arguments")
    req.add_argument(
        "--output_dir",
        required=True,
        help="Directory where output files will be saved.",
    )
    req.add_argument(
        "--filtered_file",
        required=True,
        help="Path to the filtered Ingredients pickle file.",
    )
    req.add_argument(
        "--reference_file",
        required=True,
        help="Path to the reference Ingredients pickle file.",
    )
    
    opt = ratio_sub.add_argument_group("optional arguments")
    opt.add_argument(
        "--tag",
        default="",
        help="Optional tag to prepend to output filenames for distinction.",
    )
    opt.add_argument(
        "--ratio_threshold",
        type=validate_ratio_threshold,
        default=0.5,
        help="Minimum ratio value to keep (default: %(default)s).",
    )
    
    def ratio_command(args):
        from metacooc.ratios import calculate_ratios
        
        args.tag = f"{args.tag}_" if args.tag else ""
        calculate_ratios(
            reference_ingredients=args.reference_file,
            filtered_ingredients=args.filtered_file,
            output_dir=args.output_dir,
            ratio_threshold=args.ratio_threshold,
            tag=args.tag,
        )
    
    ratio_sub.set_defaults(func=ratio_command)
    
    # ----------------------------
    # PLOT SUBCOMMAND
    # ----------------------------
    plot_sub = subparsers.add_parser("plot", help="Plot co-occurrence ratios.")
    
    req = plot_sub.add_argument_group("required arguments")
    req.add_argument(
        "--output_dir",
        required=True,
        help="Directory where output files will be saved.",
    )
    req.add_argument(
        "--ratios_file",
        required=True,
        help="Path to the ratios file to plot.",
    )
    
    opt = plot_sub.add_argument_group("optional arguments")
    opt.add_argument(
        "--tag",
        default="",
        help="Optional tag to prepend to output filenames for distinction.",
    )
    opt.add_argument(
        "--ratio_threshold",
        type=validate_ratio_threshold,
        default=0.5,
        help="Ratio threshold to include on plot (default: %(default)s).",
    )
    
    def plot_command(args):
        from metacooc.plot import plot_ratios
        
        args.tag = f"{args.tag}_" if args.tag else ""
        plot_ratios(
            ratios_file=args.ratios_file,
            output_dir=args.output_dir,
            ratio_threshold=args.ratio_threshold,
            tag=args.tag,
        )
    
    plot_sub.set_defaults(func=plot_command)
    
    # ----------------------------
    # COOCCURRENCE SUBCOMMAND
    # ----------------------------
    cooc_sub = subparsers.add_parser(
        "cooccurrence", help="Run the full co-occurrence workflow (in-memory)."
    )
    
    req = cooc_sub.add_argument_group("required arguments")
    req.add_argument(
        "--mode",
        choices=["taxon", "metadata", "biome"],
        required=True,
        help="Search mode: 'taxon', 'metadata' or 'biome'.",
    )
    req.add_argument(
        "--search_string",
        required=True,
        type=str,
        help=(
        "Search string to query as a single token. "
        "Use '|' to separate OR‑terms and '+' to separate AND‑terms. "
        "Examples: 'foo|bar', 'foo+baz', 'foo|bar+baz|qux'."
        )
    )
    req.add_argument(
        "--output_dir",
        required=True,
        help="Directory where output files will be saved.",
    )
    
    opt = cooc_sub.add_argument_group("optional arguments")
    opt.add_argument(
        "--data_dir",
        default=DEFAULT_DATA_DIR,
        help="Directory containing data files (default: %(default)s)",
    )
    opt.add_argument(
        "--tag",
        default="",
        help="Optional tag to prepend to output filenames for distinction.",
    )
    # Search-related optional
    opt.add_argument(
        "--ranks_for_search_inclusion",
        choices=["domain", "phylum", "class", "order", "family", "genus", "species"],
        help=(
            "Taxa identified at a rank higher than this rank are excluded "
            "in taxon search mode."
        ),
    )
    opt.add_argument(
        "--custom_ingredients",
        help="Initial Ingredients file to use instead of default"
    )
    opt.add_argument(
        "--column_names",
        nargs='+',
        help="meta mode only: Restrict metadata search to within specified columns. Multiple entries can be specified as e.g. --column_names 'hello' 'world'",
    )
    opt.add_argument(
        "--strict",
        action="store_true",
        help="Restrict search to a reduced set of columns.",
    )
    opt.add_argument(
        "--inverse",
        action="store_true",
        help="Return inverse of search string.",
    )
    # Filter-related optional
    opt.add_argument(
        "--min_taxa_count",
        type=positive_int,
        help="Minimum number of taxa a sample must have to be included.",
    )
    opt.add_argument(
        "--min_sample_count",
        type=positive_int,
        help="Minimum number of samples in which a taxon must be present.",
    )
    opt.add_argument(
        "--aggregated",
        action="store_true",
        help="Use the aggregated Ingredients object.",
    )
    opt.add_argument(
        "--filter_rank",
        choices=["domain", "phylum", "class", "order", "family", "genus", "species"],
        help=(
            "Taxa identified at a rank higher than this rank are filtered out "
            "of results."
        ),
    )
    # Ratio-related optional
    opt.add_argument(
        "--ratio_threshold",
        type=validate_ratio_threshold,
        default=0.5,
        help="Minimum ratio value to keep in filtered output and to plot on the output plot (default: %(default)s).",
    )
    opt.add_argument(
        "--sandpiper_version",
        default=None,
        help="Specify which data version to load (default: latest). Versions available for download can be listed with 'metacooc download --list_versions'"
    )
    
    def cooccurrence_command(args):
        from metacooc.pipelines import run_cooccurrence
        
        if args.aggregated:
            args.tag = f"{args.tag}_aggregated_" if args.tag else "aggregated_"
        else:
            args.tag = f"{args.tag}_" if args.tag else ""
        # Bundle everything into one args object for the pipeline
        run_cooccurrence(args)
    
    cooc_sub.set_defaults(func=cooccurrence_command)
    
    # ----------------------------
    # BIOME_DISTRIBUTION SUBCOMMAND
    # ----------------------------
    biome_sub = subparsers.add_parser(
        "biome_distribution", help="Return the biome distribution of Ingredients."
    )
    
    req = biome_sub.add_argument_group("required arguments")
    req.add_argument(
        "--output_dir",
        required=True,
        help="Directory where output files will be saved.",
    )
    
    opt = biome_sub.add_argument_group("optional arguments")
    opt.add_argument(
        "--data_dir",
        default=DEFAULT_DATA_DIR,
        help="Directory containing data files (default: %(default)s)",
    )
    opt.add_argument(
        "--tag",
        default="",
        help="Optional tag to prepend to output filenames for distinction.",
    )
    opt.add_argument(
        "--custom_ingredients",
        help="Initial Ingredients file to use instead of default"
    )
    opt.add_argument(
        "--sandpiper_version",
        default=None,
        help="Specify which data version to load (default: latest). Versions available for download can be listed with 'metacooc download --list_versions'"
    )
    opt.add_argument(
        "--aggregated",
        action="store_true",
        help="Use the aggregated Ingredients object.",
    )
    opt.add_argument(
        "--return_all_taxa",
        action="store_true",
        help="Specify whether to return distributions of all taxa (Not aggregated - original values) (default: True)"
    )
    def biome_distribution_command(args):
        from metacooc.pipelines import run_biome_distribution
        
        if args.aggregated:
            args.tag = f"{args.tag}_aggregated_" if args.tag else "aggregated_"
        else:
            args.tag = f"{args.tag}_" if args.tag else ""
        # Bundle everything into one args object for the pipeline
        run_biome_distribution(args)
    
    biome_sub.set_defaults(func=biome_distribution_command)
    
    # --------------
    # Parse & Dispatch
    # --------------
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    parse_cli()
