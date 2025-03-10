# MetaCoOc

**MetaCoOc** is a Python package for analyzing co‑occurrence of microorganisms in metagenomic datasets. It provides a suite of command‑line tools to download data, format raw profiles, search metadata and taxonomic profiles, filter data based on user‐defined criteria, calculate co‑occurrence ratios, and visualize results. The package is designed to work both as a step‑by‑step workflow (via individual commands) and as a full in‑memory pipeline.

## Features

- **Download:** Retrieve default data files (e.g., Ingredients objects and metadata indices) that serve as inputs for subsequent analyses.
- **Format:** Convert raw taxonomic profiles (e.g., from Sandpiper output or SRA data) and metadata into a standardized, indexed format. This step generates the Ingredients objects and pre‑built metadata indices for fast lookups.
- **Search:** Query either the taxonomic profiles or the metadata for specific terms (e.g., to identify samples containing a given taxon or metadata keyword).
- **Filter:** Subset the data based on a list of accession numbers and/or count thresholds (e.g., minimum number of taxa per sample).
- **Ratios:** Calculate co‑occurrence ratios by comparing a filtered dataset to a reference, with options for threshold‑based filtering.
- **Plot:** Visualize co‑occurrence ratios as plots.

Additionally, **MetaCoOc** offers a full pipeline command (`cooccurrence`) that integrates the search, filter, ratio, and plot steps in memory, reducing I/O overhead while saving intermediate files if desired.

## Installation

Install **MetaCoOc** via pip or conda. If you use pip, clone the repository and install in editable mode:

```bash
git clone https://github.com/yourusername/metacooc.git
cd metacooc
pip install -e .
```

Ensure that your package includes the default data directory (or create one) so that the package can use relative paths if a data directory is not specified.

## Usage

Before running the main pipeline, you must download or generate the required input files. By default, the package expects the necessary data files to reside in a data directory (relative to the installation or specified by the user).

### Download

Download the default data files (Ingredients objects and metadata indices):

```bash
metacooc download --data_dir /path/to/data [--force]
```

The `--force` flag forces re‑download of files even if they already exist.

### Format

Generate Ingredients objects and metadata indices from your raw data files:

```bash
metacooc format --tax_profile metacooc/data/sandpiper0.3.0.condensed_intersection.csv \
                --metadata_file metacooc/data/sra_metadata_parsed_intersection.tsv \
                --output_dir metacooc/data/ \
                --index_type broad \
                --strict \
                --aggregated \
                --aggregation_pattern "g__"
```

This command creates the Ingredients files (both raw and aggregated) as well as pre‑built metadata indices that will be used by subsequent steps.

### Full Pipeline (Cooccurrence)

Run the complete in‑memory co‑occurrence pipeline. This step integrates search, filter, ratio, and plot. For example, to search for samples containing the taxon “Nitrospira”:

```bash
metacooc cooccurrence --mode taxon \
                      --output_dir metacooc/data/ \
                      --search_string "Nitrospira" \
                      --data_dir metacooc/data/ \
                      --aggregated \
                      --min_taxa_count 5 \
                      --min_sample_count 5 \
					  --tag Nitrospira
```

The tag option appends the tag e.g. Nitrospira to all files. 

If you wish to apply a ratio threshold filtering (e.g., only include taxa with a ratio ≥ 0.5), add the flag:

```bash
metacooc cooccurrence --mode taxon \
                      --output_dir metacooc/data/ \
                      --search_string "Nitrospira" \
                      --data_dir metacooc/data/ \
                      --aggregated \
                      --min_taxa_count 5 \
                      --min_sample_count 5 \
                      --ratio_threshold 0.5
```

### Individual Steps

You can also run each step of the workflow separately:

#### Search (File-Based)

```bash
metacooc search --mode metadata \
                --data_dir metacooc/data/ \
                --output_dir metacooc/data/ \
                --search_string "soil" \
                --index_type broad \
                --strict
```

#### Filter

```bash
metacooc filter --accessions_file metacooc/data/search_results.txt \
                --data_dir metacooc/data/ \
                --output_dir metacooc/data/ \
                --aggregated \
                --min_taxa_count 5 \
                --min_sample_count 5
```

#### Ratios

```bash
metacooc ratio --output_dir metacooc/data/ \
               --data_dir metacooc/data/ 
```

#### Plot

```bash
metacooc plot --ratios metacooc/data/ratios.tsv \
              --output_dir metacooc/data/
```

## Data Directory Defaults

If you do not specify `--data_dir` when running a command, **metacooc** uses a default data directory relative to the installation (set via your package’s entry point). You can override this behavior by providing your own data directory using the `--data_dir` option.

## License

This project is licensed under the GNU General Public License (GPLv3+). See [LICENCE.txt](LICENCE.txt) for details.

## Author

Your Name  
[Your Contact Info or Website]

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.