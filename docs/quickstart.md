# Quick setup

## Install and download Ingredients

MetaCoOc requires Python 3.10 or newer. Install the package in the environment where you will run the analyses:

```bash
python -m pip install metacooc
```

List the published snapshots before downloading if you need to choose a release:

```bash
metacooc download --list_data_releases
```

The commands below pin the current GlobDB snapshot so that the run is reproducible. The registry-selected default may change in a later package release.

```bash
metacooc download \
  --data_release R226_globdb_rev1
```

The download contains raw and aggregated Ingredients. Metadata is an optional, larger download; add `--include_metadata` only when you plan to search metadata columns.

## First association: a biome cohort

This example does not require the metadata bundle:

```bash
metacooc association \
  --data_release R226_globdb_rev1 \
  --search_mode biome \
  --search_string soil \
  --output_dir results/soil-association
```

The output directory contains:

```text
global_association_summary.tsv  # compact table for inspection
global_association.tsv          # detailed metrics
global_association_plot.png     # automatic association plot
```

The exact prefix can include a null-scope or tag. Open the summary first, then use the detailed table when you need additional effect sizes or counts.

## A metadata cohort

Download metadata explicitly before using `--search_mode metadata`:

```bash
metacooc download \
  --data_release R226_globdb_rev1 \
  --include_metadata

metacooc association \
  --data_release R226_globdb_rev1 \
  --search_mode metadata \
  --search_string "marine" \
  --output_dir results/marine-association
```

## First co-occurrence run

Once the metadata bundle is available, a focal-taxon query can be used to
inspect a taxon's neighbourhood within a defined background:

```bash
metacooc cooccurrence \
  --data_release R226_globdb_rev1 \
  --search_mode focal_taxa \
  --search_string "s__Methanocatella smithii" \
  --null_scope metadata \
  --null_metadata_query "human gut metagenome" \
  --filter_rank species \
  --min_conditional_probability 0.10 \
  --output_dir results/smithii-cooccurrence
```

The normal output is `global_edges_summary.tsv`, `global_edges.tsv`,
`global_nodes.tsv`, and `global_cooccurrence_plot.png`. Very large edge tables
use Parquet files instead of the detailed TSV; see [Outputs and scaling](advanced/outputs.md).
