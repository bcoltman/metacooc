# MetaCoOc

**MetaCoOc** is a Python package for large-scale co-occurrence and association analysis of microorganisms in shotgun metagenomes. It provides a set of command-line tools and in-memory pipelines to:

- subset metagenomes by **taxon**, **metadata terms**, or **biome labels**
- define an explicit **null/background cohort** (global or local)
- compute **taxon–taxon co-occurrence** networks
- compute **taxon–term association** (enrichment/specificity) statistics
- quantify **community structure** (C-score, mean Jaccard, NODF) with optional null distributions
- export tidy TSV outputs suitable for downstream analysis and visualisation

---

## Table of Contents

- [Installation](#installation)
- [Data and download](#data-and-download)
- [Quick start](#quick-start)
- [Core pipelines](#core-pipelines)
  - [Co-occurrence](#co-occurrence)
  - [Association](#association)
  - [Structure](#structure)
  - [Biome distribution](#biome-distribution)
  - [Output schemas and metadata sidecars](#output-schemas-and-metadata-sidecars)
- [Key concepts](#key-concepts)
  - [Cohort vs null/background](#cohort-vs-nullbackground)
  - [Output cutoff semantics](#output-cutoff-semantics)
  - [Null scopes](#null-scopes)
  - [Null models](#null-models)
- [Command reference](#command-reference)
- [Advanced: generating Ingredients and metadata](#advanced-generating-ingredients-and-metadata)
- [License](#license)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)

---

## Installation

You will soon be able to install **MetaCoOc** via **pip** or **conda-forge**. For now, the repository can be cloned and installed like so:

```bash
git clone https://github.com/bcoltman/metacooc.git
cd metacooc
pip install -e .
```


---

## Data and download

MetaCoOc uses prebuilt datasets hosted on **Zenodo**. Data releases identify the underlying Sandpiper/database release and are available as two variants:

* `gtdb` — SingleM generated taxonomic profiles based on **GTDB** genomes
* `globdb` — SingleM generated taxonomic profiles based on **GlobDB** genomes


A complete data-release identifier must be:

```
R<database_release>_<variant>_rev<revision>
```

Examples:

* `R226_gtdb_rev1`
* `R226_globdb_rev1`

When `--data-release` is omitted from the CLI, MetaCoOc uses the exact default selected by the release registry: the latest current GlobDB snapshot. The current default is `R226_globdb_rev1`. MetaCoOc prints the selected identifier so that it can be recorded or supplied explicitly in later reproducible runs. Use `metacooc download --list-data-releases` to list all releases and identify the current default.

Supplying an exact identifier pins a run to that immutable snapshot. Explicit custom Ingredients and metadata paths suppress the corresponding default; the Python APIs continue to require explicit release or path selection.

A revision is a complete publication snapshot. If any scientific input changes, both variants and the shared files are published under the next revision, even when some contents are unchanged. Previous revisions remain immutable and downloadable.

The data release is separate from both the MetaCoOc Python package version and the internal Ingredients storage `format_version`.
Older semver-style identifiers such as `2.0.0_gtdb` are not accepted as published-data selectors.

Newly formatted Ingredients directories store `"format_version": 1` in their manifest. This integer identifies the on-disk schema, and MetaCoOc rejects directories whose format version is missing or unsupported. Official filenames also include `format1`. The exact source release is stored separately as `data_release` in both the manifest and the in-memory `Ingredients` object.

The release registry is refreshed from the MetaCoOc GitHub repository at most once per day, with a packaged registry as an offline fallback. The registry records the exact CLI default, allowing compatible data-only releases and the default selection to change without requiring a new MetaCoOc package release.

---

### Downloading data

#### Ingredients download

If you run:

```bash
metacooc download
```

MetaCoOc will:

1. Resolve and report the registry-selected latest GlobDB snapshot
2. Download and SHA-256 verify the raw and aggregated Ingredients archives
3. Install them into the default user data directory for your operating system

You can override the location using `--data_dir`:

```bash
metacooc download --data_dir ./my_data
```

To download GTDB instead, or to pin any exact snapshot, supply it explicitly:

```bash
metacooc download --data-release R226_gtdb_rev1
```

You can also set `METACOOC_DATA_DIR` to choose a persistent default location
for all MetaCoOc commands.

The download command only needs to be run once per data release and data directory.

After the files are downloaded, all subsequent analyses reuse the local copies.

---

#### What gets downloaded

For the current default data release:

```
R226_globdb_rev1
```

The following files are retrieved from Zenodo and installed locally:

##### Variant-specific

* `ingredients_raw_R226_globdb_rev1_format1/`
* `ingredients_aggregated_R226_globdb_rev1_format1/`

These are prebuilt **Ingredients** directories containing sparse matrices, labels, manifests, biome annotations, and cached taxonomic lookups.

##### Ingredients: raw vs aggregated

* **`ingredients_raw`** — uses *unfilled coverage* (as reported by [SingleM](https://wwood.github.io/singlem/tools/summarise)).  
  Coverage assigned to a taxon does **not** include coverage from its descendant taxa.

* **`ingredients_aggregated`** — uses *filled coverage* (as reported by [SingleM](https://wwood.github.io/singlem/tools/summarise)).  
  Coverage is propagated up the taxonomy, so each taxon includes the coverage of all taxa beneath it.

In short:
`raw` = coverage at that exact rank only.
`aggregated` = total coverage across the full subtree.

##### Optional base-release metadata

* `sra_metadata_R226_rev1.tsv`

This potentially large file is not downloaded by default. Download it when you need metadata searches or metadata-based cohort construction:

```bash
metacooc download --include-metadata
```

It is used for:

* metadata searches
* cohort construction

Biome classification data are stored inside each Ingredients directory. The publication also contains `sample_to_biome_R226_rev1.tsv.gz` for reproducibility, but MetaCoOc does not download it automatically.

---

## Quick start

### 1) Download a dataset

```bash
metacooc download
```


### 2) Run a complete pipeline

Example: association of taxa with the embedded biome term “soil”, using a global null background.
(This assumes the above download command has been used and therefore uses the default data directory.)

```bash
metacooc association \
  --search_mode biome \
  --search_string soil \
  --output_dir results/soil_assoc
```

This writes:

* `global_association.tsv`
* `global_plot.png`

With the default `FE` association null, no metadata sidecar is normally written because the association Jaccard null is handled analytically rather than by storing empirical null-run metadata.

---

## Core pipelines

MetaCoOc is designed around three analysis pipelines plus a biome summary export:

* `cooccurrence` — taxon–taxon edges/nodes
* `association` — taxon enrichment/specificity for a cohort vs a null/background
* `structure` — matrix-level structure metrics
* `biome_distribution` — taxon counts across annotated biome labels

The three analysis pipelines (`association`, `cooccurrence`, and `structure`) follow the same high-level pattern:

1. load Ingredients (raw or aggregated)
2. build a **cohort** of samples via `--search_mode` + `--search_string`
3. apply count-based filtering (`min_taxa_count`, `min_sample_count`, ranks)
4. define a **null/background** population (`--null_scope`, `--null_model`, …)
5. compute statistics and write TSV outputs (and plots where applicable)

The examples below use the registry-selected default release. Add an exact `--data-release` to pin a run or select GTDB.

### Co-occurrence

**Purpose:** build a directed taxon–taxon network where edges represent conditional co-occurrence above `--min_conditional_probability`.

```bash
metacooc cooccurrence \
  --search_mode taxa_context \
  --search_string "g__Nitrospira" \
  --ranks_for_search_inclusion genus \
  --output_dir results/nitrospira_cooc \
  --filter_rank species \
  --min_taxa_count 5 \
  --min_sample_count 5
```

Outputs:

* `global_edges.tsv`
* `global_nodes.tsv`
* for large edge tables: `global_edges.parquet`, `global_edges_taxa.parquet`, and `global_edges_summary.tsv`
* if an empirical null was run: `global_edges_metadata.tsv` and, for large edge tables, `global_edges_summary_metadata.tsv`

Notes:

* The **taxa universe** is determined from the cohort (after filtering), optionally restricted by `--filter_rank`.
* The co-occurrence statistics are computed on the chosen null/background Ingredients matrix, restricted to that taxa universe.
* Large universes can exceed pair limits; use `--max_pairs` or override with `--large`.
* Edge rows are directed: `source_taxon -> target_taxon`. The primary ranking metric is `p_target_given_source`.

---

### Association

**Purpose:** test which taxa are associated with a cohort (term) relative to a null/background population.

```bash
metacooc association \
  --search_mode metadata \
  --search_string soil \
  --output_dir results/soil_assoc_globdb \
  --filter_rank species \
  --min_taxa_count 50 \
  --min_sample_count 20
```

This metadata-based example requires the optional metadata table; install it first with `metacooc download --include-metadata` for the same selected release.

Outputs:

* `global_association.tsv`
* `global_plot.png`
* if an empirical null was run: `global_association_metadata.tsv`

The primary ranking columns are `p_cohort_given_taxon` (specificity: among samples containing the taxon, how often is the cohort present?) and `p_taxon_given_cohort` (sensitivity: among cohort samples, how often is the taxon present?).

---

### Structure

**Purpose:** quantify community structure in the cohort matrix (presence/absence), optionally with null distributions.

```bash
metacooc structure \
  --search_mode biome \
  --search_string soil \
  --output_dir results/soil_structure \
  --filter_rank species \
  --min_taxa_count 50 \
  --min_sample_count 20 \
  --null_model FE \
  --nm_n_reps 1000
```

Output:

* `global_structure.tsv`
* if null computation is enabled: `global_structure_metadata.tsv`

Metrics include:

* `c_score`
* `mean_jaccard`
* `nodf`

If null computation is enabled (via `nm_n_reps > 0`), the TSV also includes null mean/sd, SES, and empirical p-values.

---

### Biome distribution

**Purpose:** export a taxa × biome presence table.

```bash
metacooc biome_distribution \
  --output_dir results/biomes
```

Outputs a TSV, with behaviour depending on flags:

* `--return_all_taxa` → export all taxa
* `--taxa_query` → export only taxa matching the comma-separated query terms
* `--aggregated` → export species and aggregated taxa
* otherwise → export species to `taxa_biome_distribution_species.tsv`

Biome distribution does not run a null model, so it does not write a metadata sidecar.

---

### Output schemas and metadata sidecars

Result columns are ordered so that the most useful fields are near the left: identifiers first, then primary ranking metrics, supporting counts, effect-size/statistical metrics, optional Fisher tests, and optional empirical null summaries. The old run-level null columns are no longer repeated on every result row. When null-run metadata exists, it is written once to a sidecar TSV with columns `key` and `value`.

Null metadata sidecars use the same basename as the result file:

* `global_association.tsv` → `global_association_metadata.tsv`
* `global_edges.tsv` → `global_edges_metadata.tsv`
* `global_structure.tsv` → `global_structure_metadata.tsv`
* `global_edges_summary.tsv` → `global_edges_summary_metadata.tsv`

Sidecar keys are:

```text
null_model
null_replicates_requested
null_replicates_completed
null_replicates_ok
null_replicates_error
null_seed
null_seed_source
```

#### Association TSV

`global_association.tsv` contains one row per reported taxon. Base columns, in order:

```text
taxon
p_cohort_given_taxon
p_taxon_given_cohort
log2_rr_cohort_taxon_vs_not_taxon
rr_cohort_taxon_vs_not_taxon
log2_rr_taxon_cohort_vs_not_cohort
rr_taxon_cohort_vs_not_cohort
delta_p_taxon_cohort_vs_not_cohort
lift_taxon_cohort
jaccard_taxon_cohort
phi_coefficient
rr_cohort_given_taxon_vs_without_taxon
rr_taxon_given_cohort_vs_without_cohort
ln_rr_cohort_given_taxon_vs_without_taxon
ln_rr_taxon_given_cohort_vs_without_cohort
chi2_statistic
chi2_p_value
chi2_q_value_bh
chi2_log_p_value
chi2_log_q_value_bh
taxon_in_cohort_count
taxon_in_background_not_cohort_count
cohort_without_taxon_count
neither_taxon_nor_cohort_count
cohort_sample_count
background_not_cohort_sample_count
background_sample_count
p_taxon_given_not_cohort
p_cohort_given_not_taxon
```

With `--compute_fisher`, these columns are appended:

```text
fisher_odds_ratio
fisher_p_value
fisher_log_p_value
```

With simulated null models for association Jaccard summaries, these row-level null statistics are appended:

```text
jaccard_null_mean
jaccard_null_sd
jaccard_null_ses
jaccard_null_p_empirical
```

#### Co-occurrence TSV

`global_edges.tsv` contains one row per reported directed edge. For focal workflows, `focal_query` and `focal_taxon` appear first. Base edge columns, in order:

```text
source_taxon
target_taxon
p_target_given_source
p_source_given_target
shared_sample_count
source_taxon_sample_count
target_taxon_sample_count
source_only_sample_count
target_only_sample_count
neither_source_nor_target_sample_count
background_sample_count
source_taxon_prevalence
target_taxon_prevalence
cooccurrence_prevalence
lift_taxon_pair
jaccard_taxon_pair
phi_coefficient
chi2_statistic
chi2_p_value
chi2_q_value_bh
chi2_log_p_value
chi2_log_q_value_bh
rr_target_given_source_vs_without_source
rr_source_given_target_vs_without_target
ln_rr_target_given_source_vs_without_source
ln_rr_source_given_target_vs_without_target
```

With `--compute_fisher`, Fisher columns are appended:

```text
fisher_odds_ratio
fisher_p_value
fisher_log_p_value
```

With simulated null models for pairwise Jaccard summaries, these columns are appended:

```text
jaccard_null_mean
jaccard_null_sd
jaccard_null_ses
jaccard_null_p_empirical
jaccard_null_q_value_bh
jaccard_null_log_q_value_bh
```

`global_nodes.tsv` contains:

```text
taxon
taxon_sample_count
out_degree_p_target_given_source_gt_<threshold>
```

For very large edge tables, the full output is written compactly as `global_edges.parquet` with integer taxon IDs:

```text
source_taxon_index
target_taxon_index
shared_sample_count
jaccard_taxon_pair
chi2_log_p_value
chi2_log_q_value_bh
```

When simulated nulls are computed, compact parquet also stores the row-level null columns held in the edge arrays, including `jaccard_null_mean`, `jaccard_null_sd`, `jaccard_null_ses`, `jaccard_null_p_empirical`, and `jaccard_null_log_q_value_bh`. Fisher columns and derived `jaccard_null_q_value_bh` are reconstructed in readable TSV outputs rather than stored in the compact parquet edge table. `global_edges_taxa.parquet` maps `taxon_id` to `taxon` and `total_count`; `global_edges_summary.tsv` contains the highest-priority readable subset using the normal co-occurrence TSV schema.

#### Structure TSV

`global_structure.tsv` contains one row for each structure metric:

```text
metric
observed_value
observed_error
```

With null computation enabled, these columns are appended:

```text
null_mean
null_sd
null_standardized_effect_size
null_p_empirical
```

#### Biome Distribution TSV

Biome distribution output is a taxon-by-biome table with no hidden index column:

```text
taxon
<biome_label_1>
<biome_label_2>
...
```

Each biome column contains the number of samples in that biome where the taxon is present.

---

## Key concepts

### Cohort vs null/background

MetaCoOc explicitly separates:

* **Cohort (term) samples**: defined by `--search_mode` + `--search_string`, then filtered by count thresholds
* **Null/background samples**: defined by `--null_scope` (global by default), optionally restricted to a biome/metadata subset and/or a taxa neighbourhood

Association requires the cohort to be a **strict subset** of the null (i.e. there must be non-term samples).

---

### Output cutoff semantics

MetaCoOc uses metric-specific cutoff names:

* **association**: `--min_conditional_probability` filters taxa by cohort conditional probability
  `p_cohort_given_taxon > min_conditional_probability`, i.e. `P(cohort | taxon)`
* **cooccurrence**: `--min_conditional_probability` is the minimum **conditional probability** for an edge
  include directed edges where `p_target_given_source > min_conditional_probability`, i.e. `P(target | source)`

For large datasets, increasing `--min_conditional_probability` can dramatically reduce runtime and output size.

---

### Null scopes

`--null_scope` controls how the null/background is defined:

* *(default / None)*: global null across all samples (subject to thresholds)
* `biome`: restrict null samples to `--null_biome_query`
* `metadata`: restrict null samples to `--null_metadata_query`
* `taxa`: restrict null taxa to a neighbourhood around `--null_taxa_query`
* `biome_taxa`: biome-restricted samples, then taxa neighbourhood restriction
* `metadata_taxa`: metadata-restricted samples, then taxa neighbourhood restriction

Taxa neighbourhood parameters:

* `--taxa_degree` (radius of taxa→sample expansions)
* `--min_shared_samples_between_taxa` (BFS expansion gate)

---

### Null models

Null models generate shuffled presence/absence matrices for empirical baselines.

Supported models:

* **FF**: fixed row and column totals (Curveball; preserves both taxon prevalence and sample richness)
* **FE**: fixed row totals; columns equiprobable (preserves taxon prevalence)
* **EF**: fixed column totals; rows equiprobable (preserves sample richness)
* **EE**: fixed fill only (preserves overall sparsity)

In `association` and `cooccurrence`:

* For `FE`, analytic Fisher/χ² tests are sufficient, as Jaccard-based enrichment under an FE null is analytically equivalent to Fisher’s exact test.
* non-`FE` models compute an empirical null **for Jaccard-based metrics** (mean/sd/SES/empirical p)

In `structure`:

* null models provide empirical baselines for `c_score`, `mean_jaccard`, and `nodf`.

---

## Command reference

| Command              | What it does                                                                    |
| -------------------- | ------------------------------------------------------------------------------- |
| `download`           | Download Zenodo-hosted Ingredients and metadata files                           |
| `format`             | Convert raw taxonomic profiles into Ingredients objects                         |
| `search`             | Search by taxon / metadata / biome; returns matching accessions                 |
| `filter`             | Apply count thresholds and/or subset by accession list                          |
| `cooccurrence`       | Full in-memory co-occurrence pipeline                                           |
| `association`        | Full in-memory association pipeline                                             |
| `structure`          | Full in-memory structure pipeline                                               |
| `analysis`           | Compute association/cooccurrence/structure from explicit input files (advanced) |
| `plot`               | Plot TSV outputs from pipelines/analysis                                        |
| `biome_distribution` | Export taxa occurrence by biome                                                 |

### Search string grammar

`--search_string` is parsed as a single token with:

* `|` = OR
* `+` = AND

Examples:

* `"soil|sediment"`
* `"soil+forest"`
* `"soil|forest+rhizosphere"`

If your string contains spaces or special characters, quote it:

* `"s__Escherichia coli"`

---

## Advanced: generating Ingredients and metadata

This section is for reproducing datasets or using MetaCoOc on custom profiles.

### Formatting Sandpiper (or other) profiles into Ingredients

```bash
metacooc format \
  --tax_profile path/to/profiles.tsv \
  --output_dir ./my_data \
  --sample_to_biome_file path/to/sample_to_biome.tsv \
  --aggregated \
  --tag custom
```

This assumes `profiles.tsv` is organised like concatenated SingleM profiles (watch out for headers) and that `sample_to_biome.tsv` is a two-level labelling of the respective accessions.

For an official publication, omit `--tag` and provide the exact revision. The Ingredients schema version is added to the output names automatically:

```bash
metacooc format \
  --tax_profile path/to/sandpiper.R226.gtdb.tsv \
  --output_dir ./release \
  --sample_to_biome_file path/to/sample_to_biome_R226_rev1.tsv \
  --aggregated \
  --data-release R226_gtdb_rev1 \
  --archive_ingredients
```


### Using custom Ingredients

Most commands accept `--custom_ingredients` to bypass the default downloaded Ingredients:

```bash
metacooc association \
  --search_mode metadata \
  --search_string soil \
  --custom_ingredients path/to/ingredients_raw_custom \
  --metadata_file path/to/sra_metadata.tsv \
  --output_dir results/custom_assoc
```

### Notes on metadata

MetaCoOc can optionally download the parsed SRA metadata table belonging to an exact snapshot (`sra_metadata_<base>_rev<revision>.tsv`) by passing `--include-metadata` to `metacooc download`. Biome mappings are stored in Ingredients directories. If you build your own datasets, provide metadata explicitly with `--metadata_file` when needed.

---

## License

GNU GPL v3 or later (GPLv3+). See `LICENSE`.

---

## Acknowledgements

The NCBI SRA metadata parsing workflow was adapted from `public_sequencing_metadata_corrections` by W. Wood.

---

## Contact

Benjamin Coltman — [benjamin.coltman@univie.ac.at](mailto:benjamin.coltman@univie.ac.at)
Daan Speth — [daan.speth@univie.ac.at](mailto:daan.speth@univie.ac.at)
