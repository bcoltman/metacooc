# Search and filtering

## Cohort versus null background

MetaCoOc separates two populations:

- **Cohort samples** are matched by `--search_mode` and `--search_string`, then
  filtered by the count and presence settings.
- **Null/background samples** are used for comparison and are controlled by
  `--null_scope` and its related query options.

For association, the cohort must be a strict subset of the background so that
there are non-cohort samples for comparison. The background can be global, or
can be restricted to a biome, metadata query, or taxon neighbourhood.

`--search_mode` defines how `--search_string` is interpreted:

- `biome` uses exact biome names; comma or `|` separates alternatives.
- `metadata` searches a literal substring, optionally restricted with `--column_names` or `--strict`.
- `taxa_context` uses `|` for OR groups and `+` for terms that must occur together.
- `focal_taxa` accepts comma-separated focal queries. Its optional `LHS -> RHS` form uses the left side to define the cohort and the right side to restrict the reported taxa.

Quote queries containing spaces or shell characters. The cohort query and the null/background query are independent.

Use `--min_taxa_count` and `--min_sample_count` to remove very small samples or rare taxa. `--min_coverage` and `--min_relative_abundance` convert abundance or coverage to presence before binary analyses; rank-specific variants are available when one threshold is not suitable for every rank. `--min_conditional_probability` removes weakly directional results from the written analysis output.

Increasing `--min_conditional_probability` can reduce runtime and output size
for large analyses. It is an output filter, not a replacement for checking
whether the cohort and background are appropriate.

## Metadata search

Metadata searches stream over the SRA metadata table. MetaCoOc uses an external
`grep`/`awk` backend when it passes a small availability check and falls back to
a pure-Python streaming search otherwise. Set
`METACOOC_METADATA_SEARCH_BACKEND=auto`, `external`, or `python` to choose the
backend explicitly.
