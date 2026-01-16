#!/usr/bin/env python3
import argparse
from pathlib import Path

# Constants
BASE_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = BASE_DIR / "data"
DEFAULT_DATA_DIR.mkdir(parents=True, exist_ok=True)
RANK_CHOICES = ["domain", "phylum", "class", "order", "family", "genus", "species"]
NULL_SCOPE_CHOICES = ["biome", "taxa", "metadata", "biome_taxa", "metadata_taxa"]
NULL_MODEL_CHOICES = ["FF","FE","EF","EE","EP"]
ANALYSIS_TYPE_CHOICES = ["cooccurence", "association", "structure"]
SEARCH_MODE_CHOICES = ["taxon", "metadata", "biome"]

# Helper functions
def positive_int(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not an integer")
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"{value} is not a positive integer")
    return ivalue

def validate_threshold(value):
    try:
        value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not a valid float.")
    if value < 0.0 or value > 1.0:
        raise argparse.ArgumentTypeError("Threshold must be between 0 and 1.")
    return value

def format_tag(tag, aggregated):
    if aggregated:
        return f"{tag}_aggregated_" if tag else "aggregated_"
    return f"{tag}_" if tag else ""

def add_subcommand(subparsers, name, help_text, func):
    sub = subparsers.add_parser(name, help=help_text)
    sub.set_defaults(func=func)
    func.__subparser__ = sub  # Attach the subparser to the function
    return sub


def check_required_args(args, required_args, subparser):
    missing = [arg for arg in required_args if getattr(args, arg) is None]
    if missing:
        subparser.error(f"The following arguments are required: {', '.join(missing)}")

# Argument group helpers (with group support)
def add_data_dir_and_version(parser, group=None):
    kwargs = {
        "default": DEFAULT_DATA_DIR,
        "help": "Directory containing data files (default: %(default)s)",
    }
    if group:
        group.add_argument("--data_dir", **kwargs)
    else:
        parser.add_argument("--data_dir", **kwargs)
        
    kwargs = {
        "default": None,
        "help": "Specify which data version to load (default: latest). Versions available for download can be listed with 'metacooc download --list_versions'",
    }
    if group:
        group.add_argument("--sandpiper_version", **kwargs)
    else:
        parser.add_argument("--sandpiper_version", **kwargs)

def add_tag_and_aggregated(parser, group=None):
    kwargs = {
        "default": "",
        "help": "Optional tag to prepend to output filenames for distinction.",
    }
    if group:
        group.add_argument("--tag", **kwargs)
    else:
        parser.add_argument("--tag", **kwargs)
        
    kwargs = {
        "action": "store_true",
        "help": "Use the aggregated Ingredients object.",
    }
    if group:
        group.add_argument("--aggregated", **kwargs)
    else:
        parser.add_argument("--aggregated", **kwargs)

def add_custom_ingredients(parser, group=None):
    kwargs = {
        "help": "Path to an Ingredients file to use instead of default.",
    }
    if group:
        group.add_argument("--custom_ingredients", **kwargs)
    else:
        parser.add_argument("--custom_ingredients", **kwargs)

def add_search_mode_and_string(parser, required=False, group=None):
    kwargs = {
        "choices": SEARCH_MODE_CHOICES,
        "required": required,
        "help": "Search mode: 'taxon', 'metadata' or 'biome'.",
    }
    if group:
        group.add_argument("--search_mode", **kwargs)
    else:
        parser.add_argument("--search_mode", **kwargs)
        
    kwargs = {
        "type": str,
        "required": required,
        "help": (
            "Search string to query as a single token. "
            "Use '|' to separate OR-terms and '+' to separate AND-terms. "
            "Examples: 'foo|bar', 'foo+baz', 'foo|bar+baz|qux'. "
            "Search string must be wrapped in single or double quotes if it "
            "contains a special character or space e.g. 's__Escherichia coli'."
        ),
    }
    if group:
        group.add_argument("--search_string", **kwargs)
    else:
        parser.add_argument("--search_string", **kwargs)

def add_search_args(parser, group=None):
    kwargs = {
        "choices": RANK_CHOICES,
        "help": "Taxa identified at a rank higher than this rank are excluded in taxon search mode (default: %(default)s).",
    }
    if group:
        group.add_argument("--ranks_for_search_inclusion", **kwargs)
    else:
        parser.add_argument("--ranks_for_search_inclusion", **kwargs)
        
    kwargs = {
        "nargs": '+',
        "help": "metadata mode only: Restrict metadata search to within specified columns. Multiple entries can be specified as e.g. --column_names 'hello' 'world'",
    }
    if group:
        group.add_argument("--column_names", **kwargs)
    else:
        parser.add_argument("--column_names", **kwargs)
        
    kwargs = {
        "action": "store_true",
        "help": "metadata mode only: Restrict metadata search to within a pre-defined and reduced set of columns i.e. ['acc', 'organism', 'env_biome_sam', 'env_feature_sam', 'env_material_sam', 'biosamplemodel_sam']",
    }
    if group:
        group.add_argument("--strict", **kwargs)
    else:
        parser.add_argument("--strict", **kwargs)
        
    kwargs = {
        "action": "store_true",
        "help": "Return inverse of search e.g. if using search_term 'soil' in 'biome' search_mode, then it will return all non-soil samples.",
    }
    if group:
        group.add_argument("--inverse", **kwargs)
    else:
        parser.add_argument("--inverse", **kwargs)

def add_null_scope_args(parser, group=None):
    kwargs = {
        "choices": NULL_SCOPE_CHOICES,
        "default": None,
        "help": (
            "Scope used for null model generation. "
            "By default, uses all samples, 'biome' restricts to --null_biome_query, "
            "'taxa' uses a neighbourhood around --null_taxa_query taxa and 'metadata' restricts to --null_metadata_query. "
            "'biome_taxa' first restricts to --null_biome_query and then a neighbourhood "
            "around -null_taxa_query taxa. 'metadata_taxa' first to --null_biome_query, then around --null_taxa_query"
        ),
    }
    if group:
        group.add_argument("--null_scope", **kwargs)
    else:
        parser.add_argument("--null_scope", **kwargs)
        
    kwargs = {
        "default": None,
        "help": "Biomes used to restrict samples for null model generation when null_scope includes 'biome'.",
    }
    if group:
        group.add_argument("--null_biome_query", **kwargs)
    else:
        parser.add_argument("--null_biome_query", **kwargs)
        
    kwargs = {
        "default": None,
        "help": "Focal taxon or taxa used to define the taxa neighbourhood when null_scope includes 'taxa' (default: %(default)s).",
    }
    if group:
        group.add_argument("--null_taxa_query", **kwargs)
    else:
        parser.add_argument("--null_taxa_query", **kwargs)
        
    kwargs = {
        "default": None,
        "help": "Metadata search term used to restrict samples for null_model generation when null_scope includes 'metadata'.",
    }
    if group:
        group.add_argument("--null_metadata_query", **kwargs)
    else:
        parser.add_argument("--null_metadata_query", **kwargs)
    
    kwargs = {
        "type": positive_int,
        "default": 1,
        "help": "Minimum number of samples in which two taxa must co-occur for a new taxon to "
                "be included during BFS expansion (default: %(default)s).",
             }
    
    if group:
        group.add_argument("--min_shared_samples_between_taxa", **kwargs)
    else:
        parser.add_argument("--min_shared_samples_between_taxa", **kwargs)
    

def validate_null_scope_args(args, subparser):
    if args.null_scope in ['biome', 'biome_taxa'] and args.null_biome_query is None:
        subparser.error("--null_biome_query is required when null_scope is 'biome' or 'biome_taxa'")
    if args.null_scope in ['metadata', 'metadata_taxa'] and args.null_metadata_query is None:
        subparser.error("--null_metadata_query is required when null_scope is 'metadata' or 'metadata_taxa'")





def add_filter_args(parser, group=None):
    kwargs = {
        "type": positive_int,
        "default": 1,
        "help": "Minimum number of taxa a sample must have to be included (default: %(default)s).",
    }
    if group:
        group.add_argument("--min_taxa_count", **kwargs)
    else:
        parser.add_argument("--min_taxa_count", **kwargs)
        
    kwargs = {
        "type": positive_int,
        "default": 1,
        "help": "Minimum number of samples in which a taxon must be present (default: %(default)s).",
    }
    if group:
        group.add_argument("--min_sample_count", **kwargs)
    else:
        parser.add_argument("--min_sample_count", **kwargs)
        
    kwargs = {
        "choices": RANK_CHOICES,
        "help": "Taxa identified at a rank higher than this rank are filtered out of results (default: all included).",
    }
    if group:
        group.add_argument("--filter_rank", **kwargs)
    else:
        parser.add_argument("--filter_rank", **kwargs)
        
    kwargs = {
        "choices": RANK_CHOICES,
        "default": "species",
        "help": (
            "Taxa identified at a rank higher than this rank are not used to "
            "determine the number of unique taxa in a sample (default: %(default)s)"
        ),
    }
    if group:
        group.add_argument("--taxa_count_rank", **kwargs)
    else:
        parser.add_argument("--taxa_count_rank", **kwargs)
        
    kwargs = {
        "action": "store_true",
        "help": (
            "Remove the min_taxa_count and min_sample_count threshold on the matrix defining the null background. "
            "In all modes, this is applied before the final filtered ingredients are subsetted. "
        ),
    }
    if group:
        group.add_argument("--remove_null_threshold", **kwargs)
    else:
        parser.add_argument("--remove_null_threshold", **kwargs)
        
    kwargs = {
        "type": positive_int,
        "default": 1,
        "help": "Sets the k-degree neighbourhood that is returned from the taxa search (default: %(default)s)",
    }
    if group:
        group.add_argument("--taxa_degree", **kwargs)
    else:
        parser.add_argument("--taxa_degree", **kwargs)

def add_null_model_args(parser, group=None):
    kwargs = {
        "choices": NULL_MODEL_CHOICES,
        "default": "FE",
        "help": "Null model to use (default: %(default)s).",
    }
    if group:
        group.add_argument("--null_model", **kwargs)
    else:
        parser.add_argument("--null_model", **kwargs)
        
    kwargs = {
        "type": positive_int,
        "default": 10000,
        "help": "Number of shuffled Null models to generate (default: %(default)s).",
    }
    if group:
        group.add_argument("--nm_n_reps", **kwargs)
    else:
        parser.add_argument("--nm_n_reps", **kwargs)
    
    kwargs = {
        "type": positive_int,
        "default": 42,
        "help": "Random state for sampling (default: %(default)s).",
    }
    if group:
        group.add_argument("--nm_random_state", **kwargs)
    else:
        parser.add_argument("--nm_random_state", **kwargs)
    
def add_fisher_args(parser, group=None):
    kwargs = {
        "action": "store_true",
        "help": "Perform Fishers exact test",
    }
    if group:
        group.add_argument("--compute_fisher", **kwargs)
    else:
        parser.add_argument("--compute_fisher", **kwargs)

def add_large_and_max_pairs_args(parser, group=None):
    kwargs = {
        "action": "store_true",
        "help": "Regardless of RAM usage, or file size - calculate all coccurences",
    }
    if group:
        group.add_argument("--large", **kwargs)
    else:
        parser.add_argument("--large", **kwargs)
        
    kwargs = {
        "type": positive_int,
        "default": 100000,
        "help": "If the number of taxon pairs exceeds this value, then cooccurence will not be determined unless --large is used.",
    }
    if group:
        group.add_argument("--max_pairs", **kwargs)
    else:
        parser.add_argument("--max_pairs", **kwargs)

def add_threshold_arg(parser, group=None):
    kwargs = {
        "type": validate_threshold,
        "default": 0,
        "help": "Minimum ** value required to output entry (default: %(default)s).",
    }
    if group:
        group.add_argument("--threshold", **kwargs)
    else:
        parser.add_argument("--threshold", **kwargs)

def add_output_dir(parser, required=True, group=None):
    kwargs = {
        "required": required,
        "help": "Directory where output files will be saved.",
    }
    if group:
        group.add_argument("--output_dir", **kwargs)
    else:
        parser.add_argument("--output_dir", **kwargs)

def add_list_column_names(parser, group=None):
    kwargs = {
        "action": "store_true",
        "help": "metadata mode only: WARNING: Will produce lots of output! List available column names from NCBI metadata",
    }
    if group:
        group.add_argument("--list_column_names", **kwargs)
    else:
        parser.add_argument("--list_column_names", **kwargs)

def add_accessions_file(parser, group=None):
    kwargs = {
        "help": "File containing accession numbers to filter by.",
    }
    if group:
        group.add_argument("--accessions_file", **kwargs)
    else:
        parser.add_argument("--accessions_file", **kwargs)

def add_force(parser, group=None):
    kwargs = {
        "action": "store_true",
        "help": "Force re-download even if files exist.",
    }
    if group:
        group.add_argument("--force", **kwargs)
    else:
        parser.add_argument("--force", **kwargs)

def add_list_versions(parser, group=None):
    kwargs = {
        "action": "store_true",
        "help": "List available versions",
    }
    if group:
        group.add_argument("--list-versions", **kwargs)
    else:
        parser.add_argument("--list-versions", **kwargs)

def add_tax_profile(parser, group=None):
    kwargs = {
        "required": True,
        "help": "Taxonomic profile TSV file.",
    }
    if group:
        group.add_argument("--tax_profile", **kwargs)
    else:
        parser.add_argument("--tax_profile", **kwargs)

def add_sample_to_biome_file(parser, group=None):
    kwargs = {
        "help": "A CSV file linking SRA accessions to biome classifications.",
    }
    if group:
        group.add_argument("--sample_to_biome_file", **kwargs)
    else:
        parser.add_argument("--sample_to_biome_file", **kwargs)

def add_filtered_file(parser, group=None):
    kwargs = {
        "required": True,
        "help": "Path to the filtered Ingredients pickle file.",
    }
    if group:
        group.add_argument("--filtered_file", **kwargs)
    else:
        parser.add_argument("--filtered_file", **kwargs)

def add_null_file(parser, group=None):
    kwargs = {
        "required": True,
        "help": "Path to the null Ingredients pickle file.",
    }
    if group:
        group.add_argument("--null_file", **kwargs)
    else:
        parser.add_argument("--null_file", **kwargs)

def add_ratios_file(parser, group=None):
    kwargs = {
        "required": True,
        "help": "Path to the ratios file to plot.",
    }
    if group:
        group.add_argument("--ratios_file", **kwargs)
    else:
        parser.add_argument("--ratios_file", **kwargs)

def add_analysis_type(parser, group=None):
    kwargs = {
        "choices": ANALYSIS_TYPE_CHOICES,
        "required": True,
        "help": "Analysis mode: 'cooccurence', 'association' or 'structure'.",
    }
    if group:
        group.add_argument("--analysis_type", **kwargs)
    else:
        parser.add_argument("--analysis_type", **kwargs)

def add_return_all_taxa(parser, group=None):
    kwargs = {
        "action": "store_true",
        "help": "Specify whether to return distributions of all taxa (Not aggregated - original values) (default: True)",
    }
    if group:
        group.add_argument("--return_all_taxa", **kwargs)
    else:
        parser.add_argument("--return_all_taxa", **kwargs)

# Subcommand functions
def download_command(args):
    from metacooc.download import download_data
    download_data(data_dir=args.data_dir, force=args.force, list_versions=args.list_versions, sandpiper_version=args.sandpiper_version)

def format_command(args):
    from metacooc.format import format_data
    format_data(
        tax_profile=args.tax_profile,
        output_dir=args.output_dir,
        sample_to_biome_file=args.sample_to_biome_file,
        aggregated=args.aggregated,
        tag=args.tag
    )

def search_command(args, subparser):
    from metacooc.search import search_data
    if not args.list_column_names:
        missing = []
        if not args.search_mode:
            missing.append("--search_mode")
        if not args.search_string:
            missing.append("--search_string")
        if not args.output_dir:
            missing.append("--output_dir")
        if missing:
            subparser.error(f"The following arguments are required unless --list_column_names is used: {', '.join(missing)}")
    args.tag = format_tag(args.tag, args.aggregated)
    search_data(
        mode=args.search_mode,
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

def filter_command(args, subparser):
    from metacooc.filter import filter_data
    if not any([args.min_taxa_count, args.min_sample_count, args.accessions_file, args.filter_rank]):
        subparser.error("At least one of the following arguments is required: --min_taxa_count, --min_sample_count, --accessions_file, or --filter_rank")
    args.tag = format_tag(args.tag, args.aggregated)
    filter_data(
        accessions_file=args.accessions_file,
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        aggregated=args.aggregated,
        min_taxa_count=args.min_taxa_count,
        min_sample_count=args.min_sample_count,
        filter_rank=args.filter_rank,
        taxa_count_rank=args.taxa_count_rank,
        tag=args.tag,
        custom_ingredients=args.custom_ingredients,
        sandpiper_version=args.sandpiper_version,
        min_shared_samples_between_taxa=args.min_shared_samples_between_taxa
    )

def analysis_command(args):
    from metacooc.analysis import association, cooccurrence
    args.tag = format_tag(args.tag, False)
    if args.analysis_type == "cooccurence":
        cooccurrence(
            null_ingredients=args.null_file,
            filtered_ingredients=args.filtered_file,
            output_dir=args.output_dir,
            tag=args.tag,
            filter_rank=args.filter_rank,
            large=args.large,
            max_pairs=args.max_pairs,
            threshold=args.threshold
        )
    else:
        association(
            null_ingredients=args.null_file,
            filtered_ingredients=args.filtered_file,
            output_dir=args.output_dir,
            tag=args.tag,
            threshold=args.threshold,
            null_model=args.null_model,
            nm_n_reps=args.nm_n_reps,
            nm_random_state=args.nm_random_state,
            compute_fisher=args.compute_fisher)

def plot_command(args):
    from metacooc.plot import plot_ratios
    args.tag = format_tag(args.tag, False)
    plot_ratios(
        ratios_file=args.ratios_file,
        output_dir=args.output_dir,
        ratio_threshold=args.threshold,
        tag=args.tag,
    )

def cooccurrence_command(args):
    from metacooc.pipelines import run_cooccurrence
    args.tag = format_tag(args.tag, args.aggregated)
    run_cooccurrence(args)

def association_command(args):
    from metacooc.pipelines import run_association
    args.tag = format_tag(args.tag, args.aggregated)
    run_association(args)

def structure_command(args):
    from metacooc.pipelines import run_structure
    args.tag = format_tag(args.tag, args.aggregated)
    run_structure(args)

def biome_distribution_command(args):
    from metacooc.pipelines import run_biome_distribution
    args.tag = format_tag(args.tag, args.aggregated)
    run_biome_distribution(args)

def parse_cli():
    parser = argparse.ArgumentParser(
        description="Co-occurrence data of microorganisms based on metagenome detection"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    
    # Download subcommand
    download_sub = add_subcommand(
        subparsers,
        "download",
        "Download initial files.",
        download_command,
    )
    opt = download_sub.add_argument_group("optional arguments")
    add_data_dir_and_version(download_sub, group=opt)
    add_list_versions(download_sub, group=opt)
    add_force(download_sub, group=opt)
    
    # Format subcommand
    format_sub = add_subcommand(
        subparsers,
        "format",
        "Format data to generate Ingredients objects.",
        format_command,
    )
    req = format_sub.add_argument_group("required arguments")
    opt = format_sub.add_argument_group("optional arguments")
    add_output_dir(format_sub, group=req)
    add_tax_profile(format_sub, group=req)
    add_tag_and_aggregated(format_sub, group=opt)
    add_sample_to_biome_file(format_sub, group=opt)
    
    # Search subcommand
    search_sub = add_subcommand(
        subparsers,
        "search",
        "Perform a file-based search.",
        lambda args: search_command(args, search_sub),
    )
    req = search_sub.add_argument_group("required arguments (unless --list_column_names is used)")
    opt = search_sub.add_argument_group("optional arguments")
    add_search_mode_and_string(search_sub, group=req)
    add_output_dir(search_sub, required=False, group=req)
    add_data_dir_and_version(search_sub, group=opt)
    add_tag_and_aggregated(search_sub, group=opt)
    add_custom_ingredients(search_sub, group=opt)
    add_search_args(search_sub, group=opt)
    add_list_column_names(search_sub, group=opt)
    
    # Filter subcommand
    filter_sub = add_subcommand(
        subparsers,
        "filter",
        "Filter data by accession numbers or other criteria.",
        lambda args: filter_command(args, filter_sub),
    )
    req = filter_sub.add_argument_group("required arguments")
    opt = filter_sub.add_argument_group("optional arguments")
    add_output_dir(filter_sub, group=req)
    add_data_dir_and_version(filter_sub, group=opt)
    add_tag_and_aggregated(filter_sub, group=opt)
    add_custom_ingredients(filter_sub, group=opt)
    add_filter_args(filter_sub, group=opt)
    add_null_scope_args(filter_sub, group=opt)
    add_accessions_file(filter_sub, group=opt)
    
    # Analysis subcommand
    analysis_sub = add_subcommand(
        subparsers,
        "analysis",
        "Perform co-occurrence analysis.",
        analysis_command,
    )
    req = analysis_sub.add_argument_group("required arguments")
    opt = analysis_sub.add_argument_group("optional arguments")
    add_search_mode_and_string(analysis_sub, required=True, group=req)
    add_analysis_type(analysis_sub, group=req)
    add_output_dir(analysis_sub, group=req)
    add_filtered_file(analysis_sub, group=req)
    add_null_file(analysis_sub, group=req)
    add_threshold_arg(analysis_sub, group=opt)
    add_large_and_max_pairs_args(analysis_sub, group=opt)
    analysis_sub.add_argument(
        "--filter_rank",
        choices=RANK_CHOICES,
        help="Taxa identified at a rank higher than this rank are filtered out of results.",
    )
    
    # Plot subcommand
    plot_sub = add_subcommand(
        subparsers,
        "plot",
        "Plot co-occurrence ratios.",
        plot_command,
    )
    req = plot_sub.add_argument_group("required arguments")
    opt = plot_sub.add_argument_group("optional arguments")
    add_output_dir(plot_sub, group=req)
    add_ratios_file(plot_sub, group=req)
    add_threshold_arg(plot_sub, group=opt)
    add_tag_and_aggregated(plot_sub, group=opt)
    
    # Cooccurrence subcommand
    cooc_sub = add_subcommand(
        subparsers,
        "cooccurrence",
        "Run the full co-occurrence workflow (in-memory).",
        cooccurrence_command,
    )
    req = cooc_sub.add_argument_group("required arguments")
    opt = cooc_sub.add_argument_group("optional arguments")
    add_search_mode_and_string(cooc_sub, required=True, group=req)
    add_output_dir(cooc_sub, group=req)
    add_data_dir_and_version(cooc_sub, group=opt)
    add_tag_and_aggregated(cooc_sub, group=opt)
    add_custom_ingredients(cooc_sub, group=opt)
    add_search_args(cooc_sub, group=opt)
    add_null_scope_args(cooc_sub, group=opt)
    add_filter_args(cooc_sub, group=opt)
    add_null_model_args(cooc_sub, group=opt)
    add_fisher_args(cooc_sub, group=opt)
    add_large_and_max_pairs_args(cooc_sub, group=opt)
    add_threshold_arg(cooc_sub, group=opt)
    
    # Association subcommand
    assoc_sub = add_subcommand(
        subparsers,
        "association",
        "Run the full association workflow (in-memory).",
        association_command,
    )
    req = assoc_sub.add_argument_group("required arguments")
    opt = assoc_sub.add_argument_group("optional arguments")
    add_search_mode_and_string(assoc_sub, required=True, group=req)
    add_output_dir(assoc_sub, group=req)
    add_data_dir_and_version(assoc_sub, group=opt)
    add_tag_and_aggregated(assoc_sub, group=opt)
    add_custom_ingredients(assoc_sub, group=opt)
    add_search_args(assoc_sub, group=opt)
    add_null_scope_args(assoc_sub, group=opt)
    add_filter_args(assoc_sub, group=opt)
    add_null_model_args(assoc_sub, group=opt)
    add_fisher_args(assoc_sub, group=opt)
    add_threshold_arg(assoc_sub, group=opt)
    
    # Association subcommand
    structure_sub = add_subcommand(
        subparsers,
        "structure",
        "Run the structure analysis workflow (in-memory).",
        structure_command,
    )
    req = structure_sub.add_argument_group("required arguments")
    opt = structure_sub.add_argument_group("optional arguments")
    add_search_mode_and_string(structure_sub, required=True, group=req)
    add_output_dir(structure_sub, group=req)
    add_data_dir_and_version(structure_sub, group=opt)
    add_tag_and_aggregated(structure_sub, group=opt)
    add_custom_ingredients(structure_sub, group=opt)
    add_search_args(structure_sub, group=opt)
    add_null_scope_args(structure_sub, group=opt)
    add_filter_args(structure_sub, group=opt)
    add_null_model_args(structure_sub, group=opt)
    
    # Biome distribution subcommand
    biome_sub = add_subcommand(
        subparsers,
        "biome_distribution",
        "Return the biome distribution of Ingredients.",
        biome_distribution_command,
    )
    req = biome_sub.add_argument_group("required arguments")
    opt = biome_sub.add_argument_group("optional arguments")
    add_output_dir(biome_sub, group=req)
    add_data_dir_and_version(biome_sub, group=opt)
    add_tag_and_aggregated(biome_sub, group=opt)
    add_custom_ingredients(biome_sub, group=opt)
    add_return_all_taxa(biome_sub, group=opt)
    
    args = parser.parse_args()
    if hasattr(args, 'null_scope'):
        validate_null_scope_args(args, args.func.__subparser__)
    args.func(args)


if __name__ == "__main__":
    parse_cli()