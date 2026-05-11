#!/usr/bin/env python3
import argparse

from metacooc._paths import DATA_DIR_ENV_VAR, default_data_dir

RANK_CHOICES = ["domain", "phylum", "class", "order", "family", "genus", "species"]
NULL_SCOPE_CHOICES = ["biome", "taxa", "metadata", "biome_taxa", "metadata_taxa"]
NULL_MODEL_CHOICES = ["FF", "FE", "EF", "EE"]
ANALYSIS_TYPE_CHOICES = ["cooccurrence", "association", "structure"]

COOCCURRENCE_SEARCH_MODE_CHOICES = ["focal_taxa", "taxa_context", "metadata", "biome"]
COHORT_SEARCH_MODE_CHOICES = ["taxa_context", "metadata", "biome"]
SEARCH_SUBCOMMAND_MODE_CHOICES = ["focal_taxa", "taxa_context", "metadata", "biome"]


class StackedHelpFormatter(argparse.HelpFormatter):
    """Place each help description below the option invocation."""

    def _format_action(self, action):
        help_position = min(self._action_max_length + 2, self._max_help_position)
        help_width = max(self._width - help_position, 11)
        action_header = self._format_action_invocation(action)
        parts = ["%*s%s\n" % (self._current_indent, "", action_header)]

        if action.help and action.help.strip():
            help_text = self._expand_help(action)
            if help_text:
                for line in self._split_lines(help_text, help_width):
                    parts.append("%*s%s\n" % (help_position, "", line))

        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))

        return self._join_parts(parts)


# Helper functions
def positive_int(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not an integer")
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"{value} is not a positive integer")
    return ivalue


def nonnegative_int(value):
    try:
        ivalue = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not an integer")
    if ivalue < 0:
        raise argparse.ArgumentTypeError(f"{value} is not a non-negative integer")
    return ivalue


def validate_mp_start(value):
    import multiprocessing as mp

    methods = mp.get_all_start_methods()
    if value not in methods:
        raise argparse.ArgumentTypeError(
            f"{value!r} is not a valid multiprocessing start method. "
            f"Expected one of: {', '.join(methods)}"
        )
    return value


def validate_threshold(value):
    try:
        value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not a valid float.")
    if value < 0.0 or value > 1.0:
        raise argparse.ArgumentTypeError("Threshold must be between 0 and 1.")
    return value


def validate_nonnegative_float(value):
    try:
        value = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value} is not a valid float.")
    if value < 0.0:
        raise argparse.ArgumentTypeError("Value must be non-negative.")
    return value


def _rank_threshold_pair(value, value_parser, flag_name):
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            f"{flag_name} values must have the form rank=value"
        )
    rank, raw = value.split("=", 1)
    rank = rank.strip().lower()
    if rank not in RANK_CHOICES:
        raise argparse.ArgumentTypeError(
            f"Unknown rank '{rank}'. Expected one of: {', '.join(RANK_CHOICES)}"
        )
    return rank, value_parser(raw)


def validate_coverage_rank_threshold(value):
    return _rank_threshold_pair(
        value,
        validate_nonnegative_float,
        "--min_coverage_by_rank",
    )


def validate_relative_abundance_rank_threshold(value):
    return _rank_threshold_pair(
        value,
        validate_threshold,
        "--min_relative_abundance_by_rank",
    )


def format_tag(tag, aggregated):
    if aggregated:
        return f"{tag}_aggregated_" if tag else "aggregated_"
    return f"{tag}_" if tag else ""


def add_subcommand(subparsers, name, help_text, func):
    sub = subparsers.add_parser(
        name,
        help=help_text,
        add_help=False,
        formatter_class=StackedHelpFormatter,
    )
    sub.set_defaults(func=func)
    func.__subparser__ = sub
    return sub


def check_required_args(args, required_args, subparser):
    missing = [arg for arg in required_args if getattr(args, arg) is None]
    if missing:
        subparser.error(f"The following arguments are required: {', '.join(missing)}")


# Argument group helpers
def add_help(parser, group=None):
    (group or parser).add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="show this help message and exit",
    )


def required_group(parser, description=None):
    return parser.add_argument_group("REQUIRED ARGUMENTS", description)


def argument_group(parser, title, description):
    return parser.add_argument_group(f"{title} ARGUMENTS", description)


def add_data_dir(parser, group=None):
    default = str(default_data_dir())
    kwargs = {
        "default": default,
        "help": (
            "Directory containing data files (default: %(default)s; "
            f"overridden by ${DATA_DIR_ENV_VAR})"
        ),
    }
    (group or parser).add_argument("--data_dir", **kwargs)


def add_data_version(parser, group=None, mode: str = "load"):
    """
    Add --data_version CLI option.

    mode:
        "load"   -> select which version to load (default: latest)
        "format" -> label/version to stamp into generated files (default: none)
    """
    if mode == "load":
        help_text = (
            "Specify which data version to load (default: latest). "
            "Versions available for download can be listed with "
            "'metacooc download --list_data_versions'."
        )
    elif mode == "format":
        help_text = (
            "Optional version label to embed in the generated Ingredients files. "
            "No default is applied."
        )
    else:
        raise ValueError("mode must be 'load' or 'format'")

    kwargs = {
        "default": None,
        "help": help_text,
    }
    (group or parser).add_argument("--data_version", **kwargs)


def add_tag_and_aggregated(parser, group=None):
    add_tag(parser, group=group)
    add_aggregated(parser, group=group)


def add_tag(parser, group=None):
    (group or parser).add_argument(
        "--tag",
        default="",
        help="Optional tag to prepend to output filenames for distinction.",
    )


def add_aggregated(parser, group=None):
    (group or parser).add_argument(
        "--aggregated",
        action="store_true",
        help="Use the aggregated Ingredients object.",
    )


def add_custom_ingredients(parser, group=None):
    (group or parser).add_argument(
        "--custom_ingredients",
        help="Path to an Ingredients directory to use instead of default.",
    )


def add_archive_ingredients(parser, group=None):
    (group or parser).add_argument(
        "--archive_ingredients",
        action="store_true",
        help="Also write .tar.gz archives for generated Ingredients directories.",
    )


def add_metadata_file(parser, group=None):
    (group or parser).add_argument(
        "--metadata_file",
        help=(
            "Explicit metadata TSV file to use for metadata search instead of "
            "resolving one from --data_dir/--data_version."
        ),
    )


def add_search_mode_and_string(parser, required=False, group=None, choices=None):
    choices = choices or COHORT_SEARCH_MODE_CHOICES

    (group or parser).add_argument(
        "--search_mode",
        choices=choices,
        required=required,
        help=f"Search mode. Allowed values: {', '.join(choices)}.",
    )

    (group or parser).add_argument(
        "--search_string",
        type=str,
        required=required,
        help=(
            "Search string. "
            "In taxa_context mode, use '|' to separate OR groups and '+' to separate AND terms. "
            "In focal_taxa mode, use commas to separate independent focal taxa queries. "
            "focal_taxa does not accept '|' or '+'. "
            "focal_taxa also supports optional 'LHS -> RHS' syntax, where LHS defines the focal cohort "
            "and RHS defines downstream retrieval-target taxa for cooccurrence reporting. "
            "Quote search strings containing spaces or special characters, e.g. "
            "'s__Escherichia coli', 'g__VFJL01 AGGREGATED', or "
            "'g__VFJL01 -> s__Ignicoccus hospitalis'."
        ),
    )


def add_search_args(parser, group=None):
    (group or parser).add_argument(
        "--ranks_for_search_inclusion",
        choices=RANK_CHOICES,
        help=(
            "Taxa identified at a rank higher than this rank are excluded in taxa_context "
            "search mode (default: no additional exclusion)."
        ),
    )

    (group or parser).add_argument(
        "--column_names",
        nargs="+",
        help=(
            "metadata mode only: Restrict metadata search to specified columns. "
            "Example: --column_names organism env_biome_sam"
        ),
    )

    (group or parser).add_argument(
        "--strict",
        action="store_true",
        help=(
            "metadata mode only: Restrict metadata search to a predefined reduced "
            "set of columns: ['acc', 'organism', 'env_biome_sam', 'env_feature_sam', "
            "'env_material_sam', 'biosamplemodel_sam']"
        ),
    )

    (group or parser).add_argument(
        "--inverse",
        action="store_true",
        help=(
            "Return inverse of search. Not supported in focal_taxa mode."
        ),
    )


def add_null_scope_args(parser, group=None):
    (group or parser).add_argument(
        "--null_scope",
        choices=NULL_SCOPE_CHOICES,
        default=None,
        help=(
            "Scope used for null model generation. "
            "By default, uses all samples. "
            "'biome' restricts to --null_biome_query, "
            "'taxa' uses a neighbourhood around --null_taxa_query, "
            "'metadata' restricts to --null_metadata_query, "
            "'biome_taxa' first restricts to --null_biome_query and then to a taxa neighbourhood, "
            "'metadata_taxa' first restricts to --null_metadata_query and then to a taxa neighbourhood."
        ),
    )

    (group or parser).add_argument(
        "--null_biome_query",
        default=None,
        help="Biomes used to restrict samples for null model generation when null_scope includes 'biome'.",
    )

    (group or parser).add_argument(
        "--null_taxa_query",
        default=None,
        help="Focal taxon or taxa used to define the taxa neighbourhood when null_scope includes 'taxa'.",
    )

    (group or parser).add_argument(
        "--null_metadata_query",
        default=None,
        help="Metadata search term used to restrict samples for null model generation when null_scope includes 'metadata'.",
    )

    (group or parser).add_argument(
        "--min_shared_samples_between_taxa",
        type=positive_int,
        default=1,
        help=(
            "Minimum number of samples in which two taxa must co-occur for a new taxon "
            "to be included during BFS expansion (default: %(default)s)."
        ),
    )


def validate_null_scope_args(args, subparser):
    if args.null_scope in ["biome", "biome_taxa"] and args.null_biome_query is None:
        subparser.error("--null_biome_query is required when null_scope is 'biome' or 'biome_taxa'")
    if args.null_scope in ["metadata", "metadata_taxa"] and args.null_metadata_query is None:
        subparser.error("--null_metadata_query is required when null_scope is 'metadata' or 'metadata_taxa'")


def add_count_filter_args(
    parser,
    group=None,
    min_taxa_count_default=1,
    min_sample_count_default=1,
):
    (group or parser).add_argument(
        "--min_taxa_count",
        type=positive_int,
        default=min_taxa_count_default,
        help="Minimum number of taxa a sample must have to be included (default: %(default)s).",
    )

    (group or parser).add_argument(
        "--min_sample_count",
        type=positive_int,
        default=min_sample_count_default,
        help="Minimum number of samples in which a taxon must be present (default: %(default)s).",
    )

    (group or parser).add_argument(
        "--filter_rank",
        choices=RANK_CHOICES,
        help="Taxa identified at a rank higher than this rank are filtered out of results (default: all included).",
    )

    (group or parser).add_argument(
        "--taxa_count_rank",
        choices=RANK_CHOICES,
        default="species",
        help=(
            "Taxa identified at a rank higher than this rank are not used to determine "
            "the number of unique taxa in a sample (default: %(default)s)"
        ),
    )


def add_filter_args(parser, group=None):
    add_count_filter_args(parser, group=group)

    (group or parser).add_argument(
        "--remove_null_threshold",
        action="store_true",
        help=(
            "Remove the min_taxa_count and min_sample_count threshold on the matrix "
            "defining the null background."
        ),
    )

    (group or parser).add_argument(
        "--taxa_degree",
        type=positive_int,
        default=1,
        help="Sets the k-degree neighbourhood returned from taxa neighbourhood expansion (default: %(default)s).",
    )


def add_presence_threshold_args(parser, group=None):
    target = group or parser
    target.add_argument(
        "--min_coverage",
        type=validate_nonnegative_float,
        default=None,
        help=(
            "Minimum raw coverage required for a taxon-sample cell to count as "
            "present. Applied before binary workflows."
        ),
    )
    target.add_argument(
        "--min_coverage_by_rank",
        nargs="+",
        type=validate_coverage_rank_threshold,
        default=None,
        metavar="rank=value",
        help=(
            "Rank-specific raw coverage thresholds, for example "
            "--min_coverage_by_rank species=0.02 genus=0.01."
        ),
    )
    target.add_argument(
        "--min_relative_abundance",
        type=validate_threshold,
        default=None,
        help=(
            "Minimum relative abundance fraction required for a taxon-sample "
            "cell to count as present, e.g. 0.02 for 2%%."
        ),
    )
    target.add_argument(
        "--min_relative_abundance_by_rank",
        nargs="+",
        type=validate_relative_abundance_rank_threshold,
        default=None,
        metavar="rank=value",
        help=(
            "Rank-specific relative abundance thresholds as fractions, for example "
            "--min_relative_abundance_by_rank species=0.02 genus=0.01."
        ),
    )


def add_null_model_args(parser, group=None):
    (group or parser).add_argument(
        "--null_model",
        choices=NULL_MODEL_CHOICES,
        default="FE",
        help="Null model to use (default: %(default)s).",
    )

    (group or parser).add_argument(
        "--nm_n_reps",
        type=positive_int,
        default=10000,
        help="Number of shuffled null models to generate (default: %(default)s).",
    )

    (group or parser).add_argument(
        "--nm_seed",
        type=nonnegative_int,
        default=None,
        help="Seed for null sampling. If omitted, a seed is generated and reported.",
    )

    (group or parser).add_argument(
        "--nm_n_workers",
        type=positive_int,
        default=None,
        help="Number of null simulation worker processes (default: one worker).",
    )

    (group or parser).add_argument(
        "--nm_mp_start",
        type=validate_mp_start,
        default=None,
        help="Multiprocessing start method for null simulations.",
    )

    (group or parser).add_argument(
        "--nm_burn_in_steps",
        type=nonnegative_int,
        default=None,
        help="FF null model burn-in steps per chain.",
    )

    (group or parser).add_argument(
        "--nm_steps_per_rep",
        type=nonnegative_int,
        default=None,
        help="FF null model Curveball steps between retained replicates.",
    )

    (group or parser).add_argument(
        "--nm_progress_every",
        type=positive_int,
        default=25,
        help="Null simulation progress update interval in replicates (default: %(default)s).",
    )


def add_fisher_args(parser, group=None):
    (group or parser).add_argument(
        "--compute_fisher",
        action="store_true",
        help="Perform Fisher's exact test.",
    )


def add_large_and_max_pairs_args(parser, group=None):
    (group or parser).add_argument(
        "--large",
        action="store_true",
        help="Calculate all cooccurrences regardless of RAM usage or output size.",
    )

    (group or parser).add_argument(
        "--max_pairs",
        type=positive_int,
        default=100000,
        help="If the number of taxon pairs exceeds this value, cooccurrence will not be determined unless --large is used.",
    )


def add_min_conditional_probability_arg(parser, group=None, default=0.0):
    (group or parser).add_argument(
        "--min_conditional_probability",
        type=validate_threshold,
        default=default,
        help=(
            "Minimum conditional probability required to output an association "
            "entry or co-occurrence edge "
            "(default: %(default)s)."
        ),
    )


def add_q_threshold_arg(parser, group=None):
    (group or parser).add_argument(
        "--q_threshold",
        type=validate_threshold,
        default=0.0,
        help="Maximum q-value to highlight in plots (default: %(default)s).",
    )


def add_output_dir(parser, required=True, group=None):
    (group or parser).add_argument(
        "--output_dir",
        required=required,
        help="Directory where output files will be saved.",
    )


def add_list_column_names(parser, group=None):
    (group or parser).add_argument(
        "--list_column_names",
        action="store_true",
        help="metadata mode only: WARNING: May produce lots of output. List available column names from NCBI metadata.",
    )


def add_list_biomes(parser, group=None):
    (group or parser).add_argument(
        "--list_biomes",
        action="store_true",
        help="List available biome query terms from the selected Ingredients object.",
    )


def add_accessions_file(parser, group=None):
    (group or parser).add_argument(
        "--accessions_file",
        help="File containing accession numbers to filter by.",
    )


def add_force(parser, group=None):
    (group or parser).add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if files exist.",
    )


def add_list_data_versions(parser, group=None):
    (group or parser).add_argument(
        "--list_data_versions",
        action="store_true",
        help="List available versions.",
    )


def add_tax_profile(parser, group=None):
    (group or parser).add_argument(
        "--tax_profile",
        required=True,
        help="Taxonomic profile TSV file.",
    )


def add_sample_to_biome_file(parser, group=None):
    (group or parser).add_argument(
        "--sample_to_biome_file",
        help="A CSV file linking SRA accessions to biome classifications.",
    )


def add_filtered_file(parser, group=None):
    (group or parser).add_argument(
        "--filtered_file",
        required=True,
        help="Path to the filtered Ingredients directory.",
    )


def add_null_file(parser, group=None, required=True):
    (group or parser).add_argument(
        "--null_file",
        required=required,
        help=(
            "Null Ingredients directory required for cooccurrence and association "
            "analyses; not used for structure."
        ),
    )


def add_analysis_file(parser, group=None):
    (group or parser).add_argument(
        "--analysis_file",
        required=True,
        help="Path to the analysis file to plot.",
    )


def add_analysis_type(parser, group=None):
    (group or parser).add_argument(
        "--analysis_type",
        choices=ANALYSIS_TYPE_CHOICES,
        required=True,
        help="Analysis mode: 'cooccurrence', 'association' or 'structure'.",
    )


def add_return_all_taxa(parser, group=None):
    (group or parser).add_argument(
        "--return_all_taxa",
        action="store_true",
        help="Return distributions of all taxa (not aggregated/original values).",
    )


def add_biome_distribution_args(parser, group=None):
    (group or parser).add_argument(
        "--taxa_query",
        help=(
            "Restrict output rows to taxa matching comma-separated focal-taxa-style "
            "queries, for example 'g__Nitrosomonas,g__Nitrosospira'."
        ),
    )
    (group or parser).add_argument(
        "--biome_level",
        choices=["level_1", "level_2"],
        default="level_1",
        help="Biome level to report (default: %(default)s).",
    )


# Subcommand functions
def download_command(args):
    from metacooc.download import download_data
    download_data(
        data_dir=args.data_dir,
        force=args.force,
        list_data_versions=args.list_data_versions,
        data_version=args.data_version,
    )


def format_command(args):
    from metacooc.format import format_data
    format_data(
        tax_profile=args.tax_profile,
        output_dir=args.output_dir,
        sample_to_biome_file=args.sample_to_biome_file,
        aggregated=args.aggregated,
        tag=args.tag,
        data_version=args.data_version,
        archive_ingredients=args.archive_ingredients,
    )


def search_command(args, subparser):
    from metacooc.search import search_data

    if not args.list_column_names and not args.list_biomes:
        missing = []
        if not args.search_mode:
            missing.append("--search_mode")
        if not args.search_string:
            missing.append("--search_string")
        if not args.output_dir:
            missing.append("--output_dir")
        if missing:
            subparser.error(
                "The following arguments are required unless --list_column_names "
                f"or --list_biomes is used: {', '.join(missing)}"
            )

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
        data_version=args.data_version,
        list_column_names=args.list_column_names,
        list_biomes=args.list_biomes,
        aggregated=args.aggregated,
        metadata_file=args.metadata_file,
        min_coverage=args.min_coverage,
        min_coverage_by_rank=args.min_coverage_by_rank,
        min_relative_abundance=args.min_relative_abundance,
        min_relative_abundance_by_rank=args.min_relative_abundance_by_rank,
    )


def filter_command(args, subparser):
    from metacooc.filter import filter_data

    if (
        args.min_taxa_count == 1
        and args.min_sample_count == 1
        and args.accessions_file is None
        and args.filter_rank is None
        and args.min_coverage is None
        and not args.min_coverage_by_rank
        and args.min_relative_abundance is None
        and not args.min_relative_abundance_by_rank
    ):
        subparser.error(
            "At least one of the following must be provided: "
            "--min_taxa_count (not 1), --min_sample_count (not 1), "
            "--accessions_file, --filter_rank, --min_coverage, "
            "--min_coverage_by_rank, --min_relative_abundance, or "
            "--min_relative_abundance_by_rank"
        )

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
        null_scope=args.null_scope,
        null_taxa_query=args.null_taxa_query,
        null_biome_query=args.null_biome_query,
        null_metadata_query=args.null_metadata_query,
        remove_null_threshold=args.remove_null_threshold,
        taxa_degree=args.taxa_degree,
        min_shared_samples_between_taxa=args.min_shared_samples_between_taxa,
        custom_ingredients=args.custom_ingredients,
        data_version=args.data_version,
        metadata_file=args.metadata_file,
        min_coverage=args.min_coverage,
        min_coverage_by_rank=args.min_coverage_by_rank,
        min_relative_abundance=args.min_relative_abundance,
        min_relative_abundance_by_rank=args.min_relative_abundance_by_rank,
    )


def analysis_command(args):
    args.tag = format_tag(args.tag, False)
    subparser = analysis_command.__subparser__

    if args.analysis_type in {"cooccurrence", "association"} and not args.null_file:
        subparser.error("--null_file is required when --analysis_type is cooccurrence or association")

    if args.analysis_type == "cooccurrence":
        from metacooc.analysis import cooccurrence
        cooccurrence(
            null_ingredients=args.null_file,
            filtered_ingredients=args.filtered_file,
            output_dir=args.output_dir,
            tag=args.tag,
            filter_rank=args.filter_rank,
            large=args.large,
            max_pairs=args.max_pairs,
            min_conditional_probability=(
                0.0
                if args.min_conditional_probability is None
                else args.min_conditional_probability
            ),
            null_model=args.null_model,
            nm_n_reps=args.nm_n_reps,
            nm_seed=args.nm_seed,
            nm_n_workers=args.nm_n_workers,
            nm_mp_start=args.nm_mp_start,
            nm_burn_in_steps=args.nm_burn_in_steps,
            nm_steps_per_rep=args.nm_steps_per_rep,
            nm_progress_every=args.nm_progress_every,
        )

    elif args.analysis_type == "structure":
        if args.min_conditional_probability is not None:
            subparser.error("--min_conditional_probability is only valid for --analysis_type cooccurrence and association")
        from metacooc.structure import structure
        structure(
            ingredients=args.filtered_file,
            output_dir=args.output_dir,
            tag=args.tag,
            null_model=args.null_model,
            nm_n_reps=args.nm_n_reps,
            compute_null=True,
            nm_seed=args.nm_seed,
            nm_n_workers=args.nm_n_workers,
            nm_mp_start=args.nm_mp_start,
            nm_burn_in_steps=args.nm_burn_in_steps,
            nm_steps_per_rep=args.nm_steps_per_rep,
            nm_progress_every=args.nm_progress_every,
        )

    elif args.analysis_type == "association":
        from metacooc.analysis import association
        association(
            null_ingredients=args.null_file,
            filtered_ingredients=args.filtered_file,
            output_dir=args.output_dir,
            tag=args.tag,
            min_conditional_probability=(
                0.0
                if args.min_conditional_probability is None
                else args.min_conditional_probability
            ),
            null_model=args.null_model,
            nm_n_reps=args.nm_n_reps,
            nm_seed=args.nm_seed,
            compute_fisher=args.compute_fisher,
            nm_n_workers=args.nm_n_workers,
            nm_mp_start=args.nm_mp_start,
            nm_burn_in_steps=args.nm_burn_in_steps,
            nm_steps_per_rep=args.nm_steps_per_rep,
            nm_progress_every=args.nm_progress_every,
        )


def plot_command(args):
    from metacooc.plot import plot_analysis
    args.tag = format_tag(args.tag, False)
    plot_analysis(
        df_file=args.analysis_file,
        output_dir=args.output_dir,
        tag=args.tag,
        q_thresh=args.q_threshold,
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


def build_parser():
    parser = argparse.ArgumentParser(
        description="Co-occurrence data of microorganisms based on metagenome detection",
        formatter_class=StackedHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Download
    download_sub = add_subcommand(
        subparsers,
        "download",
        "Download initial files.",
        download_command,
    )
    general = argument_group(download_sub, "GENERAL", None)
    data = argument_group(
        download_sub,
        "DATA SOURCE",
        "Choose where MetaCoOc data are stored and which data version to use.",
    )
    download_opts = argument_group(
        download_sub,
        "DOWNLOAD",
        "List available versions or force a fresh download.",
    )
    add_help(download_sub, group=general)
    add_data_dir(download_sub, group=data)
    add_data_version(download_sub, group=data)
    add_list_data_versions(download_sub, group=download_opts)
    add_force(download_sub, group=download_opts)

    # Format
    format_sub = add_subcommand(
        subparsers,
        "format",
        "Format data to generate Ingredients objects.",
        format_command,
    )
    req = required_group(format_sub, "Inputs that must be provided to format Ingredients.")
    general = argument_group(format_sub, "GENERAL", None)
    output = argument_group(
        format_sub,
        "OUTPUT",
        "Set output labels, generated Ingredients variants, and archive behavior.",
    )
    metadata = argument_group(
        format_sub,
        "METADATA",
        "Provide biome labels or metadata inputs to include in formatted Ingredients.",
    )
    add_output_dir(format_sub, group=req)
    add_tax_profile(format_sub, group=req)
    add_help(format_sub, group=general)
    add_tag_and_aggregated(format_sub, group=output)
    add_data_version(format_sub, group=output, mode="format")
    add_archive_ingredients(format_sub, group=output)
    add_sample_to_biome_file(format_sub, group=metadata)

    # Search
    search_sub = add_subcommand(
        subparsers,
        "search",
        "Perform a file-based search.",
        lambda args: search_command(args, search_sub),
    )
    req = required_group(
        search_sub,
        "Required unless --list_column_names or --list_biomes is used.",
    )
    general = argument_group(search_sub, "GENERAL", None)
    data = argument_group(
        search_sub,
        "DATA SOURCE",
        "Choose the Ingredients and metadata sources used for searching.",
    )
    output = argument_group(search_sub, "OUTPUT", "Set output filename labels.")
    search_opts = argument_group(
        search_sub,
        "SEARCH",
        "Refine how taxa, metadata, and biome searches are interpreted.",
    )
    presence = argument_group(
        search_sub,
        "PRESENCE THRESHOLD",
        "Convert abundance or coverage values to presence before searching.",
    )
    listing = argument_group(
        search_sub,
        "LISTING",
        "List metadata columns or available biome terms instead of running a search.",
    )
    add_search_mode_and_string(search_sub, group=req, choices=SEARCH_SUBCOMMAND_MODE_CHOICES)
    add_output_dir(search_sub, required=False, group=req)
    add_help(search_sub, group=general)
    add_data_dir(search_sub, group=data)
    add_data_version(search_sub, group=data)
    add_aggregated(search_sub, group=data)
    add_custom_ingredients(search_sub, group=data)
    add_metadata_file(search_sub, group=data)
    add_tag(search_sub, group=output)
    add_search_args(search_sub, group=search_opts)
    add_presence_threshold_args(search_sub, group=presence)
    add_list_column_names(search_sub, group=listing)
    add_list_biomes(search_sub, group=listing)

    # Filter
    filter_sub = add_subcommand(
        subparsers,
        "filter",
        "Filter data by accession numbers or other criteria.",
        lambda args: filter_command(args, filter_sub),
    )
    req = required_group(filter_sub, "Output location for filtered Ingredients.")
    general = argument_group(filter_sub, "GENERAL", None)
    data = argument_group(
        filter_sub,
        "DATA SOURCE",
        "Choose the Ingredients and metadata sources to filter.",
    )
    output = argument_group(filter_sub, "OUTPUT", "Set output filename labels.")
    selection = argument_group(
        filter_sub,
        "SAMPLE SELECTION",
        "Provide explicit accession inputs for sample filtering.",
    )
    matrix = argument_group(
        filter_sub,
        "MATRIX FILTER",
        "Filter samples and taxa by counts, rank, and taxa neighbourhood settings.",
    )
    presence = argument_group(
        filter_sub,
        "PRESENCE THRESHOLD",
        "Convert abundance or coverage values to presence before filtering.",
    )
    null_scope = argument_group(
        filter_sub,
        "NULL SCOPE",
        "Restrict the background sample set used for null model generation.",
    )
    add_output_dir(filter_sub, group=req)
    add_help(filter_sub, group=general)
    add_data_dir(filter_sub, group=data)
    add_data_version(filter_sub, group=data)
    add_aggregated(filter_sub, group=data)
    add_custom_ingredients(filter_sub, group=data)
    add_metadata_file(filter_sub, group=data)
    add_tag(filter_sub, group=output)
    add_accessions_file(filter_sub, group=selection)
    add_filter_args(filter_sub, group=matrix)
    add_presence_threshold_args(filter_sub, group=presence)
    add_null_scope_args(filter_sub, group=null_scope)

    # Analysis
    analysis_sub = add_subcommand(
        subparsers,
        "analysis",
        "Perform co-occurrence, association or structure analysis.",
        analysis_command,
    )
    req = required_group(analysis_sub, "Inputs that define the analysis to run.")
    general = argument_group(analysis_sub, "GENERAL", None)
    mode_inputs = argument_group(
        analysis_sub,
        "ANALYSIS TYPE SPECIFIC",
        "Inputs required only by selected --analysis_type values.",
    )
    output = argument_group(analysis_sub, "OUTPUT", "Set output filename labels and Ingredients variant labels.")
    analysis_opts = argument_group(
        analysis_sub,
        "ANALYSIS",
        "Control statistical tests, reporting thresholds, and co-occurrence limits.",
    )
    null_model = argument_group(
        analysis_sub,
        "NULL MODEL",
        "Configure null model selection, replicates, seeds, and multiprocessing.",
    )
    add_analysis_type(analysis_sub, group=req)
    add_output_dir(analysis_sub, group=req)
    add_filtered_file(analysis_sub, group=req)
    add_null_file(analysis_sub, group=mode_inputs, required=False)
    add_help(analysis_sub, group=general)
    add_tag_and_aggregated(analysis_sub, group=output)
    add_min_conditional_probability_arg(analysis_sub, group=analysis_opts, default=None)
    add_large_and_max_pairs_args(analysis_sub, group=analysis_opts)
    analysis_opts.add_argument(
        "--filter_rank",
        choices=RANK_CHOICES,
        help="Taxa identified at a rank higher than this rank are filtered out of results.",
    )
    add_fisher_args(analysis_sub, group=analysis_opts)
    add_null_model_args(analysis_sub, group=null_model)

    # Plot
    plot_sub = add_subcommand(
        subparsers,
        "plot",
        "Plot analysis.",
        plot_command,
    )
    req = required_group(plot_sub, "Inputs and output location for plotting.")
    general = argument_group(plot_sub, "GENERAL", None)
    output = argument_group(plot_sub, "OUTPUT", "Set output filename labels and Ingredients variant labels.")
    plot_opts = argument_group(plot_sub, "PLOT", "Control plot highlighting thresholds.")
    add_output_dir(plot_sub, group=req)
    add_analysis_file(plot_sub, group=req)
    add_help(plot_sub, group=general)
    add_tag_and_aggregated(plot_sub, group=output)
    add_q_threshold_arg(plot_sub, group=plot_opts)

    # Cooccurrence
    cooc_sub = add_subcommand(
        subparsers,
        "cooccurrence",
        "Run the full co-occurrence workflow (in-memory).",
        cooccurrence_command,
    )
    req = required_group(cooc_sub, "Query and output settings required for the full workflow.")
    general = argument_group(cooc_sub, "GENERAL", None)
    data = argument_group(
        cooc_sub,
        "DATA SOURCE",
        "Choose the Ingredients and metadata sources used by the workflow.",
    )
    output = argument_group(cooc_sub, "OUTPUT", "Set output filename labels.")
    search_opts = argument_group(
        cooc_sub,
        "SEARCH",
        "Refine how taxa, metadata, biome, and focal-taxa searches are interpreted.",
    )
    matrix = argument_group(
        cooc_sub,
        "MATRIX FILTER",
        "Filter samples and taxa by counts, rank, and taxa neighbourhood settings.",
    )
    presence = argument_group(
        cooc_sub,
        "PRESENCE THRESHOLD",
        "Convert abundance or coverage values to presence before binary analysis.",
    )
    null_scope = argument_group(
        cooc_sub,
        "NULL SCOPE",
        "Restrict the background sample set used for null model generation.",
    )
    null_model = argument_group(
        cooc_sub,
        "NULL MODEL",
        "Configure null model selection, replicates, seeds, and multiprocessing.",
    )
    analysis_opts = argument_group(
        cooc_sub,
        "ANALYSIS",
        "Control statistical tests, reporting thresholds, and co-occurrence limits.",
    )
    add_search_mode_and_string(cooc_sub, required=True, group=req, choices=COOCCURRENCE_SEARCH_MODE_CHOICES)
    add_output_dir(cooc_sub, group=req)
    add_help(cooc_sub, group=general)
    add_data_dir(cooc_sub, group=data)
    add_data_version(cooc_sub, group=data)
    add_aggregated(cooc_sub, group=data)
    add_custom_ingredients(cooc_sub, group=data)
    add_metadata_file(cooc_sub, group=data)
    add_tag(cooc_sub, group=output)
    add_search_args(cooc_sub, group=search_opts)
    add_filter_args(cooc_sub, group=matrix)
    add_presence_threshold_args(cooc_sub, group=presence)
    add_null_scope_args(cooc_sub, group=null_scope)
    add_null_model_args(cooc_sub, group=null_model)
    add_fisher_args(cooc_sub, group=analysis_opts)
    add_large_and_max_pairs_args(cooc_sub, group=analysis_opts)
    add_min_conditional_probability_arg(cooc_sub, group=analysis_opts)

    # Association
    assoc_sub = add_subcommand(
        subparsers,
        "association",
        "Run the full association workflow (in-memory).",
        association_command,
    )
    req = required_group(assoc_sub, "Query and output settings required for the full workflow.")
    general = argument_group(assoc_sub, "GENERAL", None)
    data = argument_group(
        assoc_sub,
        "DATA SOURCE",
        "Choose the Ingredients and metadata sources used by the workflow.",
    )
    output = argument_group(assoc_sub, "OUTPUT", "Set output filename labels.")
    search_opts = argument_group(
        assoc_sub,
        "SEARCH",
        "Refine how taxa, metadata, and biome searches are interpreted.",
    )
    matrix = argument_group(
        assoc_sub,
        "MATRIX FILTER",
        "Filter samples and taxa by counts, rank, and taxa neighbourhood settings.",
    )
    presence = argument_group(
        assoc_sub,
        "PRESENCE THRESHOLD",
        "Convert abundance or coverage values to presence before binary analysis.",
    )
    null_scope = argument_group(
        assoc_sub,
        "NULL SCOPE",
        "Restrict the background sample set used for null model generation.",
    )
    null_model = argument_group(
        assoc_sub,
        "NULL MODEL",
        "Configure null model selection, replicates, seeds, and multiprocessing.",
    )
    analysis_opts = argument_group(
        assoc_sub,
        "ANALYSIS",
        "Control statistical tests and reporting thresholds.",
    )
    add_search_mode_and_string(assoc_sub, required=True, group=req, choices=COHORT_SEARCH_MODE_CHOICES)
    add_output_dir(assoc_sub, group=req)
    add_help(assoc_sub, group=general)
    add_data_dir(assoc_sub, group=data)
    add_data_version(assoc_sub, group=data)
    add_aggregated(assoc_sub, group=data)
    add_custom_ingredients(assoc_sub, group=data)
    add_metadata_file(assoc_sub, group=data)
    add_tag(assoc_sub, group=output)
    add_search_args(assoc_sub, group=search_opts)
    add_filter_args(assoc_sub, group=matrix)
    add_presence_threshold_args(assoc_sub, group=presence)
    add_null_scope_args(assoc_sub, group=null_scope)
    add_null_model_args(assoc_sub, group=null_model)
    add_fisher_args(assoc_sub, group=analysis_opts)
    add_min_conditional_probability_arg(assoc_sub, group=analysis_opts)

    # Structure
    structure_sub = add_subcommand(
        subparsers,
        "structure",
        "Run the structure analysis workflow (in-memory).",
        structure_command,
    )
    req = required_group(structure_sub, "Query and output settings required for the full workflow.")
    general = argument_group(structure_sub, "GENERAL", None)
    data = argument_group(
        structure_sub,
        "DATA SOURCE",
        "Choose the Ingredients and metadata sources used by the workflow.",
    )
    output = argument_group(structure_sub, "OUTPUT", "Set output filename labels.")
    search_opts = argument_group(
        structure_sub,
        "SEARCH",
        "Refine how taxa, metadata, and biome searches are interpreted.",
    )
    matrix = argument_group(
        structure_sub,
        "MATRIX FILTER",
        "Filter samples and taxa by counts, rank, and taxa neighbourhood settings.",
    )
    presence = argument_group(
        structure_sub,
        "PRESENCE THRESHOLD",
        "Convert abundance or coverage values to presence before structure metrics.",
    )
    null_scope = argument_group(
        structure_sub,
        "NULL SCOPE",
        "Restrict the background sample set used for null model generation.",
    )
    null_model = argument_group(
        structure_sub,
        "NULL MODEL",
        "Configure null model selection, replicates, seeds, and multiprocessing.",
    )
    add_search_mode_and_string(structure_sub, required=True, group=req, choices=COHORT_SEARCH_MODE_CHOICES)
    add_output_dir(structure_sub, group=req)
    add_help(structure_sub, group=general)
    add_data_dir(structure_sub, group=data)
    add_data_version(structure_sub, group=data)
    add_aggregated(structure_sub, group=data)
    add_custom_ingredients(structure_sub, group=data)
    add_metadata_file(structure_sub, group=data)
    add_tag(structure_sub, group=output)
    add_search_args(structure_sub, group=search_opts)
    add_filter_args(structure_sub, group=matrix)
    add_presence_threshold_args(structure_sub, group=presence)
    add_null_scope_args(structure_sub, group=null_scope)
    add_null_model_args(structure_sub, group=null_model)

    # Biome distribution
    biome_sub = add_subcommand(
        subparsers,
        "biome_distribution",
        "Return the biome distribution of Ingredients.",
        biome_distribution_command,
    )
    req = required_group(biome_sub, "Output location for biome distribution tables.")
    general = argument_group(biome_sub, "GENERAL", None)
    data = argument_group(
        biome_sub,
        "DATA SOURCE",
        "Choose the Ingredients source used for biome summaries.",
    )
    output = argument_group(biome_sub, "OUTPUT", "Set output filename labels.")
    matrix = argument_group(
        biome_sub,
        "MATRIX FILTER",
        "Filter samples and taxa by counts and rank before reporting.",
    )
    presence = argument_group(
        biome_sub,
        "PRESENCE THRESHOLD",
        "Convert abundance or coverage values to presence before reporting.",
    )
    biome_opts = argument_group(
        biome_sub,
        "BIOME DISTRIBUTION",
        "Control taxa selection and biome level reporting.",
    )
    add_output_dir(biome_sub, group=req)
    add_help(biome_sub, group=general)
    add_data_dir(biome_sub, group=data)
    add_data_version(biome_sub, group=data)
    add_aggregated(biome_sub, group=data)
    add_custom_ingredients(biome_sub, group=data)
    add_tag(biome_sub, group=output)
    add_count_filter_args(
        biome_sub,
        group=matrix,
        min_taxa_count_default=None,
        min_sample_count_default=None,
    )
    add_presence_threshold_args(biome_sub, group=presence)
    add_return_all_taxa(biome_sub, group=biome_opts)
    add_biome_distribution_args(biome_sub, group=biome_opts)

    return parser


def parse_cli():
    from metacooc._data_config import DataVersionError

    parser = build_parser()
    args = parser.parse_args()

    if hasattr(args, "null_scope"):
        validate_null_scope_args(args, args.func.__subparser__)

    try:
        args.func(args)
    except DataVersionError as e:
        parser.error(str(e))


if __name__ == "__main__":
    parse_cli()
