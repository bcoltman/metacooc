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
- [Data and variants](#data-and-variants)
- [Quick start](#quick-start)
- [Core pipelines](#core-pipelines)
  - [Co-occurrence](#co-occurrence)
  - [Association](#association)
  - [Structure](#structure)
  - [Biome distribution](#biome-distribution)
- [Key concepts](#key-concepts)
  - [Cohort vs null/background](#cohort-vs-nullbackground)
  - [Threshold semantics](#threshold-semantics)
  - [Null scopes](#null-scopes)
  - [Null models](#null-models)
- [Command reference](#command-reference)
- [Advanced: generating Ingredients and metadata](#advanced-generating-ingredients-and-metadata)
- [License](#license)
- [Acknowledgements](#acknowledgements)
- [Contact](#contact)

---

## Installation

MetaCoOc can be installed via **pip** or **conda-forge**.

### pip

```bash
pip install metacooc
````

### conda-forge

```bash
conda install -c conda-forge metacooc
```

---

## Data and download

MetaCoOc uses prebuilt datasets hosted on **Zenodo**. Currently, the data release is available as **two Sandpiper database variants**:

* `gtdb` — SingleM generated taxonomic profiles based on **GTDB** genomes
* `globdb` — SingleM generated taxonomic profiles based on **GlobDB** genomes


A complete version string must be:

```
<base_version>_<variant>
```

Examples:

* `1.1.0_gtdb`
* `1.1.0_globdb`

If no version is specified, MetaCoOc defaults to the latest available base release with the default variant (`gtdb`).

---

### Downloading data

#### Default behaviour

If you run:

```bash
metacooc download
```

MetaCoOc will:

1. Use the **latest available version** (e.g. `1.1.0_gtdb`)
2. Download all required files
3. Install them into the **package data directory** (inside the installed metacooc environment)

You can override the location using:

```bash
metacooc download --data_dir ./my_data
```

The download command only needs to be run once per data version and data directory.

After the files are downloaded, all subsequent analyses reuse the local copies.

---

#### What gets downloaded

For a version such as:

```
1.1.0_gtdb
```

The following files are retrieved from Zenodo and decompressed locally:

##### Variant-specific

* `ingredients_raw_1.1.0_gtdb.pkl`
* `ingredients_aggregated_1.1.0_gtdb.pkl`

These are prebuilt **Ingredients** objects containing presence/absence matrices and cached taxonomic lookups.

##### Ingredients: raw vs aggregated

* **`ingredients_raw`** — uses *unfilled coverage* (as reported by [SingleM](https://wwood.github.io/singlem/tools/summarise)).  
  Coverage assigned to a taxon does **not** include coverage from its descendant taxa.

* **`ingredients_aggregated`** — uses *filled coverage* (as reported by [SingleM](https://wwood.github.io/singlem/tools/summarise)).  
  Coverage is propagated up the taxonomy, so each taxon includes the coverage of all taxa beneath it.

In short:
`raw` = coverage at that exact rank only.
`aggregated` = total coverage across the full subtree.

##### Base-version shared files

* `sra_metadata_1.1.0.tsv`
* `sample_to_biome_1.1.0.tsv`

These are used for:

* metadata searches
* biome classification
* cohort construction

---

## Quick start

### 1) Download a dataset

```bash
metacooc download
```


### 2) Run a complete pipeline

Example: association of taxa with the metadata term “soil”, using a global null background.
(This assumes the above download command has been used and therefore uses the default data directory and latest data version.)

```bash
metacooc association \
  --search_mode metadata \
  --search_string soil \
  --output_dir results/soil_assoc
```

This writes:

* `global_association.tsv`
* `global_plot.png`

---

## Core pipelines

MetaCoOc is designed around three “normal use” pipelines:

* `cooccurrence` — taxon–taxon edges/nodes
* `association` — taxon enrichment/specificity for a cohort vs a null/background
* `structure` — matrix-level structure metrics

All pipelines follow the same high-level pattern:

1. load Ingredients (raw or aggregated)
2. build a **cohort** of samples via `--search_mode` + `--search_string`
3. apply count-based filtering (`min_taxa_count`, `min_sample_count`, ranks)
4. define a **null/background** population (`--null_scope`, `--null_model`, …)
5. compute statistics and write TSV outputs (and plots where applicable)

### Co-occurrence

**Purpose:** build a directed taxon–taxon network where edges represent conditional co-occurrence above a threshold.

```bash
metacooc cooccurrence \
  --search_mode taxon \
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

Notes:

* The **taxa universe** is determined from the cohort (after filtering), optionally restricted by `--filter_rank`.
* The co-occurrence statistics are computed on the chosen null/background Ingredients matrix, restricted to that taxa universe.
* Large universes can exceed pair limits; use `--max_pairs` or override with `--large`.

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

Outputs:

* `global_association.tsv`
* `global_plot.png`

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
* otherwise → export species

---

## Key concepts

### Cohort vs null/background

MetaCoOc explicitly separates:

* **Cohort (term) samples**: defined by `--search_mode` + `--search_string`, then filtered by count thresholds
* **Null/background samples**: defined by `--null_scope` (global by default), optionally restricted to a biome/metadata subset and/or a taxa neighbourhood

Association requires the cohort to be a **strict subset** of the null (i.e. there must be non-term samples).

---

### Threshold semantics

The flag name is the same across pipelines, but the meaning depends on the analysis:

* **association**: `--threshold` filters taxa by **specificity**
  `p(T | X) > threshold`
* **cooccurrence**: `--threshold` is the minimum **conditional probability** for an edge
  include edges where `P(B | A) > threshold`

For large datasets, increasing `--threshold` can dramatically reduce runtime and output size.

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

THis assumes 'profiles.tsv' is organised like concatenated singleM profiles (watch out for headers) and that sample_to_biome.tsv is a 2-level labelling of the respective accessions


### Using custom Ingredients

Most commands accept `--custom_ingredients` to bypass the default downloaded Ingredients:

```bash
metacooc association \
  --search_mode metadata \
  --search_string soil \
  --custom_ingredients path/to/ingredients_raw_custom.pkl \
  --output_dir results/custom_assoc
```

### Notes on metadata

MetaCoOc ships a parsed SRA metadata table (`sra_metadata_<base>.tsv`) and a `sample_to_biome_<base>.tsv` mapping for biome-level queries. If you build your own datasets, ensure these files are present (or adjust your workflow to use only taxon-based searches).

---

## License

GNU GPL v3 or later (GPLv3+). See `LICENCE.txt`.

---

## Acknowledgements

The NCBI SRA metadata parsing workflow was adapted from `public_sequencing_metadata_corrections` by W. Wood.

---

## Contact

Benjamin Coltman — [benjamin.coltman@univie.ac.at](mailto:benjamin.coltman@univie.ac.at)
Daan Speth — [daan.speth@univie.ac.at](mailto:daan.speth@univie.ac.at)