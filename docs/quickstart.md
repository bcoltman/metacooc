# Quick setup

## Install MetaCoOc

Choose either of the following installation routes.

### Bioconda

```bash
conda create --name metacooc \
  --channel conda-forge \
  --channel bioconda \
  metacooc
conda activate metacooc
```

### PyPI

```bash
python -m pip install metacooc
```

Confirm that the command is available:

```bash
metacooc --help
```

See [Installation](installation.md) for a source checkout and documentation
development setup.

## Select and download data

MetaCoOc CLI commands can use the registry's current GlobDB release when
`--data_release` is omitted. This is convenient for an initial exploration.
Pinning an exact release is preferable when a run must be repeated, compared,
or reported because the registry default can change.

List the published snapshots when you need to inspect the available choices:

```bash
metacooc download --list_data_releases
```

The examples below pin `R226_globdb_rev1` and use metadata-defined cohorts, so
the optional metadata bundle is downloaded with the taxonomic profiles:

```bash
metacooc download \
  --data_release R226_globdb_rev1 \
  --include_metadata
```

Metadata is the larger, optional part of the download. Analyses based only on
biome labels or taxa do not require it.

## How the workflow fits together

```{figure} _images/metacooc-workflow.svg
:alt: MetaCoOc workflow from published profiles through cohort selection, filtering, analysis, and reporting, with background and null-model choices feeding the analysis.
:width: 100%

The main MetaCoOc workflow and the choices that define its comparison.
```

The search defines the cohort and therefore the biological question. Presence
thresholds and count filters decide which observations enter the analysis. The
background is selected separately because it defines the relevant comparison
population. For an important result, empirical null constraints can then test
whether the pattern remains unusual after selected properties of the community
matrix are explained.

## First association: bee gut

Biome labels provide broad, consistently named environments. Metadata is more
useful when the cohort depends on a specific ecological or host description.
Here, `bee gut` selects samples whose archive metadata contains that phrase:

```bash
metacooc association \
  --data_release R226_globdb_rev1 \
  --search_mode metadata \
  --search_string "bee gut" \
  --filter_rank species \
  --output_dir results/bee-gut-association
```

The default FE analysis provides a fast first comparison between the cohort
and its background. The output directory contains:

```text
global_association_summary.tsv  # compact table for inspection
global_association.tsv          # detailed metrics
global_association_plot.png     # automatic association plot
```

In the summary, positive phi indicates that a species occurs with the cohort
more often than expected from their separate frequencies. High
`p_cohort_given_taxon` means that samples containing the species are often in
the bee-gut cohort; high `p_taxon_given_cohort` means that the species is found
across much of the cohort. Read these effect sizes alongside the adjusted
q-value and support counts. They describe an archive association, not a direct
host relationship or mechanism.

The [association examples](examples/association.md) provide the complete
bee-gut and human-gut workflows, plots, and interpretation limits.

## First co-occurrence: an archaeal neighbourhood

A focal-taxon search asks which taxa occur with one microorganism. The cohort
contains samples with the focal taxon, while `--null_scope metadata` restricts
the comparison background to samples described as human-gut metagenomes:

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

The normal outputs are `global_edges_summary.tsv`, `global_edges.tsv`,
`global_nodes.tsv`, and `global_cooccurrence_plot.png`. A strong positive edge
has positive phi, suitable support, adjusted significance, and a conditional
probability that is meaningful in the reported source-to-target direction.
Such an edge prioritises a relationship for follow-up; it does not establish
physical contact, metabolite transfer, or sample independence.

The [co-occurrence example](examples/cooccurrence.md) develops this into a
three-stage *M. smithii* and *Christensenella* workflow. Very large edge tables
use Parquet instead of a detailed TSV; see [Outputs and
scaling](advanced/outputs.md).
