"""Conditional argument validation for the MetaCoOc CLI."""

import os


def validate_null_scope_args(args, subparser):
    scope = args.null_scope
    requirements = {
        "biome": {"null_biome_query"},
        "taxa": {"null_taxa_query"},
        "metadata": {"null_metadata_query"},
        "biome_taxa": {"null_biome_query", "null_taxa_query"},
        "metadata_taxa": {"null_metadata_query", "null_taxa_query"},
    }
    required = requirements.get(scope, set())
    query_flags = {
        "null_biome_query": "--null_biome_query",
        "null_taxa_query": "--null_taxa_query",
        "null_metadata_query": "--null_metadata_query",
    }
    missing = [query_flags[name] for name in required if getattr(args, name) is None]
    if missing:
        subparser.error(
            f"{', '.join(missing)} required when --null_scope is {scope!r}"
        )
    unused = [
        flag
        for name, flag in query_flags.items()
        if getattr(args, name) is not None and name not in required
    ]
    if unused:
        subparser.error(
            f"{', '.join(unused)} cannot be used with --null_scope {scope!r}"
        )


def _subparser_error(args, message):
    args.func.__subparser__.error(message)


def _apply_default_data_release(args):
    from metacooc._data_config import get_default_data_release, load_registry

    registry = load_registry()
    args.data_release = get_default_data_release(registry)
    print(f"Using default data release: {args.data_release}")


def _validate_search_mode_options(args):
    mode = args.search_mode
    if mode != "metadata" and (args.column_names or args.strict):
        _subparser_error(
            args,
            "--column_names and --strict are only valid with --search_mode metadata",
        )
    if mode != "taxa_context" and args.ranks_for_search_inclusion is not None:
        _subparser_error(
            args,
            "--ranks_for_search_inclusion is only valid with "
            "--search_mode taxa_context",
        )
    if mode == "focal_taxa" and args.inverse:
        _subparser_error(args, "--inverse is not valid with --search_mode focal_taxa")


def _validate_search_cli(args):
    if args.list_column_names and args.list_biomes:
        _subparser_error(
            args,
            "--list_column_names and --list_biomes are mutually exclusive",
        )

    listing = args.list_column_names or args.list_biomes
    if listing:
        ignored = []
        for name, flag in (
            ("search_mode", "--search_mode"),
            ("search_string", "--search_string"),
            ("output_dir", "--output_dir"),
            ("ranks_for_search_inclusion", "--ranks_for_search_inclusion"),
            ("column_names", "--column_names"),
        ):
            if getattr(args, name) is not None:
                ignored.append(flag)
        for name, flag in (
            ("strict", "--strict"),
            ("inverse", "--inverse"),
        ):
            if getattr(args, name):
                ignored.append(flag)
        if args.tag:
            ignored.append("--tag")
        if any(
            getattr(args, name) is not None
            for name in (
                "min_coverage",
                "min_coverage_by_rank",
                "min_relative_abundance",
                "min_relative_abundance_by_rank",
            )
        ):
            ignored.append("presence-threshold options")
        if ignored:
            _subparser_error(
                args,
                f"{', '.join(ignored)} cannot be used with listing options",
            )

        if args.list_column_names:
            if args.custom_ingredients:
                _subparser_error(
                    args,
                    "--custom_ingredients cannot be used with --list_column_names",
                )
            if args.aggregated:
                _subparser_error(
                    args,
                    "--aggregated cannot be used with --list_column_names",
                )
            if args.metadata_file and args.data_release:
                _subparser_error(
                    args,
                    "--metadata_file and --data-release are alternative metadata "
                    "sources for --list_column_names",
                )
            if args.metadata_file is None and args.data_release is None:
                _apply_default_data_release(args)
        else:
            if args.metadata_file:
                _subparser_error(
                    args,
                    "--metadata_file cannot be used with --list_biomes",
                )
            if args.custom_ingredients is None and args.data_release is None:
                _apply_default_data_release(args)
        return

    missing = [
        flag
        for name, flag in (
            ("search_mode", "--search_mode"),
            ("search_string", "--search_string"),
            ("output_dir", "--output_dir"),
        )
        if getattr(args, name) is None
    ]
    if missing:
        _subparser_error(
            args,
            "The following arguments are required unless a listing option is "
            f"used: {', '.join(missing)}",
        )

    _validate_search_mode_options(args)
    if args.search_mode == "metadata":
        if args.custom_ingredients:
            _subparser_error(
                args,
                "--custom_ingredients is not used by metadata-only searches",
            )
        if args.aggregated:
            _subparser_error(
                args,
                "--aggregated is not used by metadata-only searches",
            )
        if args.metadata_file and args.data_release:
            _subparser_error(
                args,
                "--metadata_file and --data-release are alternative metadata sources",
            )
        if args.metadata_file is None and args.data_release is None:
            _apply_default_data_release(args)
    else:
        if args.metadata_file:
            _subparser_error(
                args,
                "--metadata_file is only valid with --search_mode metadata",
            )
        if args.custom_ingredients is None and args.data_release is None:
            _apply_default_data_release(args)


def _validate_pipeline_cli(args):
    _validate_search_mode_options(args)
    uses_metadata = args.search_mode == "metadata" or args.null_scope in {
        "metadata",
        "metadata_taxa",
    }
    if args.metadata_file and not uses_metadata:
        _subparser_error(
            args,
            "--metadata_file requires a metadata search mode or metadata null scope",
        )
    if args.custom_ingredients and uses_metadata and not args.metadata_file:
        _subparser_error(
            args,
            "--metadata_file is required for metadata searches with "
            "--custom_ingredients",
        )
    if args.custom_ingredients is None and args.data_release is None:
        _apply_default_data_release(args)


def _validate_analysis_cli(args):
    if args.analysis_type in {"cooccurrence", "association"}:
        if not args.null_file:
            _subparser_error(
                args,
                "--null_file is required when --analysis_type is "
                f"{args.analysis_type}",
            )
    elif args.null_file:
        _subparser_error(
            args,
            "--null_file is not used when --analysis_type is structure",
        )

    if args.analysis_type != "cooccurrence":
        incompatible = []
        if args.large:
            incompatible.append("--large")
        if args.max_pairs is not None:
            incompatible.append("--max_pairs")
        if args.filter_rank is not None:
            incompatible.append("--filter_rank")
        if incompatible:
            _subparser_error(
                args,
                f"{', '.join(incompatible)} only apply to cooccurrence analysis",
            )
    if args.analysis_type == "structure":
        incompatible = []
        if args.min_conditional_probability is not None:
            incompatible.append("--min_conditional_probability")
        if args.compute_fisher:
            incompatible.append("--compute_fisher")
        if incompatible:
            _subparser_error(
                args,
                f"{', '.join(incompatible)} do not apply to structure analysis",
            )


def _validate_input_paths(args):
    directory_options = {
        "custom_ingredients": "--custom_ingredients",
        "filtered_file": "--filtered_file",
        "null_file": "--null_file",
    }
    file_options = {
        "tax_profile": "--tax_profile",
        "sample_to_biome_file": "--sample_to_biome_file",
        "metadata_file": "--metadata_file",
        "accessions_file": "--accessions_file",
        "analysis_file": "--analysis_file",
    }
    for name, flag in directory_options.items():
        value = getattr(args, name, None)
        if value is not None and not os.path.isdir(value):
            _subparser_error(args, f"{flag} must be an existing directory: {value}")
    for name, flag in file_options.items():
        value = getattr(args, name, None)
        if value is not None and not os.path.isfile(value):
            _subparser_error(args, f"{flag} must be an existing file: {value}")
    output_dir = getattr(args, "output_dir", None)
    if output_dir is not None and os.path.exists(output_dir) and not os.path.isdir(output_dir):
        _subparser_error(args, f"--output_dir is not a directory: {output_dir}")


def prepare_cli_args(args):
    """Validate conditional CLI contracts and apply the published-data default."""
    command = args.command
    if getattr(args, "custom_ingredients", None) and getattr(
        args, "data_release", None
    ):
        _subparser_error(
            args,
            "--custom_ingredients and --data-release are mutually exclusive; "
            "custom Ingredients use the identity stored in their manifest",
        )
    if hasattr(args, "null_scope"):
        validate_null_scope_args(args, args.func.__subparser__)

    if command == "download":
        if args.list_data_releases:
            incompatible = []
            if args.data_release:
                incompatible.append("--data-release")
            if args.include_metadata:
                incompatible.append("--include-metadata")
            if args.force:
                incompatible.append("--force")
            if incompatible:
                _subparser_error(
                    args,
                    f"{', '.join(incompatible)} cannot be used with "
                    "--list-data-releases",
                )
        elif args.data_release is None:
            _apply_default_data_release(args)
    elif command == "search":
        _validate_search_cli(args)
    elif command == "filter":
        uses_metadata = args.null_scope in {"metadata", "metadata_taxa"}
        if args.metadata_file and not uses_metadata:
            _subparser_error(
                args,
                "--metadata_file requires a metadata null scope",
            )
        if args.custom_ingredients and uses_metadata and not args.metadata_file:
            _subparser_error(
                args,
                "--metadata_file is required for metadata null scopes with "
                "--custom_ingredients",
            )
        if args.custom_ingredients is None and args.data_release is None:
            _apply_default_data_release(args)
    elif command in {"cooccurrence", "association", "structure"}:
        _validate_pipeline_cli(args)
    elif command == "biome_distribution":
        if args.taxa_query and args.return_all_taxa:
            _subparser_error(
                args,
                "--taxa_query and --return_all_taxa are mutually exclusive",
            )
        if args.custom_ingredients is None and args.data_release is None:
            _apply_default_data_release(args)
    elif command == "analysis":
        _validate_analysis_cli(args)

    _validate_input_paths(args)
