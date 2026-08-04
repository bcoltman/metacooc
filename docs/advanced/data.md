# Data releases and Ingredients

## Published releases

Published data identifiers have the form `R<database>_<variant>_rev<revision>`,
for example `R226_globdb_rev1` or `R226_gtdb_rev1`. A revision is an immutable
scientific snapshot and is separate from both the MetaCoOc package version and
the Ingredients `format_version`.

The CLI can select the registry's current GlobDB release when
`--data_release` is omitted. Pin an exact release when an analysis needs to be
repeated, compared, or reported. Python APIs require an explicit release or
custom path.

Use `metacooc download --list_data_releases` to inspect the registry. Metadata
is optional because it is large. `--include_metadata` downloads the shared
metadata bundle; alternatively, provide a local TSV with `--metadata_file`.

## Why raw and aggregated Ingredients exist

MetaCoOc's published profiles come from
[Sandpiper](https://sandpiper.qut.edu.au/About) and follow SingleM's
rank-specific, or "unfilled", coverage convention. Coverage assigned directly
to a species is not automatically included in the coverage reported for its
genus. This preserves the rank at which each lineage was identified, including
lineages that could only be assigned to a higher taxonomic rank. The
[SingleM glossary](https://wwood.github.io/singlem/Glossary) describes this
profile convention in more detail.

Raw Ingredients retain those original profile rows. Aggregated Ingredients
retain the raw rows and add marked `AGGREGATED` rows whose coverage is summed
from all descendants. For example, if a sample contains coverage assigned to
two named species of a genus and additional coverage assigned only to that
genus, the aggregated genus represents the whole detected clade. The raw rows
preserve the three original assignments separately.

Use raw Ingredients when the rank at which a lineage was originally assigned
is part of the question, or when analysing the original species-level rows.
Use `--aggregated` when the question concerns the occurrence or abundance of a
whole higher-rank clade. Raw and aggregated rows are related views of the same
profile and should not be treated as independent observations.

## Custom Ingredients

A custom Ingredients directory can be supplied with `--custom_ingredients`.
Custom profile TSV files use the columns `sample`, `taxonomy`, and `coverage`.
A sample-to-biome mapping uses `accession`, `level_1`, and `level_2`. Keep the
input files, data identity, and preparation command with the analysis record.
