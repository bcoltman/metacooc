# MetaCoOc

**MetaCoOc** is a Python package for co-occurrence analysis of microorganisms in shotgun metagenomes. It provides a suite of command‑line tools that allow you to download necessary files, format raw input data, perform metadata and taxonomic searches, filter datasets, calculate co‑occurrence ratios, and generate plots for data visualization. The package is designed to work both as a step‑by‑step workflow (with individual commands) and as a full in‑memory pipeline (using the `cooccurrence` subcommand).

## Overview

The typical workflow with **MetaCoOc** involves:

1. **Download**: Retrieve the default data files (e.g., raw and aggregated Ingredients pickles, metadata index files) that serve as inputs for subsequent steps.
2. **Format**: Convert raw taxonomic profiles (e.g. from Sandpiper output) and metadata (e.g., from NCBI SRA) into standardized formats. This step generates the Ingredients objects and pre‑built metadata indices for fast lookups.
3. **Search**: Query the data for specific terms using either taxon-based or metadata-based search strategies.
4. **Filter**: Subset the dataset based on accession lists and count thresholds (e.g., a minimum number of taxa per sample or a minimum number of samples per taxon).
5. **Ratios**: Calculate co‑occurrence ratios by comparing a filtered dataset to a reference.
6. **Plot**: Visualize the ratios with customizable plots.
7. **Pipeline**: Optionally, run the entire search–filter–ratio–plot workflow in one command using the `cooccurrence` subcommand.

## Installation

You can install **MetaCoOc** via pip (or conda) after cloning the repository:

```bash
git clone https://github.com/bcoltman/metacooc.git
cd metacooc
pip install -e .
```

Make sure that your data files are included in the package (or available separately), as the package sets a default data directory relative to its installation.

## Usage

### Download

Before any analysis, download the default data files. These include:
- `ingredients_raw.pkl`
- `ingredients_aggregated.pkl`
- `broad_strict_metadata_index.pkl` - a pre‑built metadata index based on select metadata columns \["organism", "env_biome_sam", "env_feature_sam", "env_material_sam", "biosamplemodel_sam"\]

Run:

```bash
metacooc download --data_dir /path/to/data [--force]
```

The `--force` flag forces re‑download even if the files already exist.

Custom files can be processed and formatted using the format option - please see below

### Full Pipeline (Cooccurrence)

Once the files are downloaded and formatted, you can run the full in‑memory pipeline. For example, to search for samples containing the taxon "Nitrospira" (in taxon mode):

```bash
metacooc cooccurrence --mode taxon \
                      --data_dir metacooc/data/ \
                      --output_dir metacooc/data/ \
                      --aggregated \
                      --search_string "Nitrospira" \
                      --rank genus \
                      --min_taxa_count 5 \
                      --min_sample_count 5
```

You can also include a ratio threshold (e.g., only include taxa with ratio ≥ 0.5) by adding `--ratio_threshold 0.5`.

### Running Individual Steps

You can also run each step separately. For example:

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
metacooc ratio --data_dir metacooc/data/ \
               --output_dir metacooc/data/ \
               --ratio_threshold 0.5
```

#### Plot

```bash
metacooc plot --ratios metacooc/data/ratios.tsv \
              --output_dir metacooc/data/ \
              --ratio_threshold 0.5
```

## Data Preparation and Custom Parsing

### Parsing NCBI SRA Metadata

A key component of **MetaCoOc** is the standardized parsing of SRA metadata. The parsing workflow was adapted from [public_sequencing_metadata_corrections](https://github.com/wwood/public_sequencing_metadata_corrections). The process involves:

1. **Exporting Data from BigQuery**:  
   Run a query on [Google BigQuery](https://console.cloud.google.com/bigquery) such as:
   ```sql
	EXPORT DATA
	  OPTIONS( uri='gs://my-bucket-name/sra_metadata/*.json',
		format='JSON') AS
	SELECT
	  *
	FROM
	  nih-sra-datastore.sra.metadata
	WHERE
	  (librarysource = 'METAGENOMIC' OR organism = 'metagenome')
	  AND platform = 'ILLUMINA'
	  AND consent = 'public'
	  AND (mbases > 1000 OR (libraryselection = 'RANDOM' AND mbases > 100))
	  AND mbases <= 200000
	  AND librarysource NOT IN ('VIRAL RNA', 'METATRANSCRIPTOMIC', 'TRANSCRIPTOMIC')
   ```
   Replace `my-bucket-name` with your bucket name.

2. **Downloading the Data**:  
   Use `gsutil` to download the JSON files locally:
   ```bash
   gsutil -m cp -r gs://my-bucket-name/sra_metadata/*.json .
   ```
   Once downloaded, you can delete the files from the bucket to avoid storage costs.

3. **Setting Up the Parsing Environment**:  
   Create and activate a conda environment:
   ```bash
   mamba env create -n public_sequencing_metadata_corrections -f env.yml
   conda activate public_sequencing_metadata_corrections
   ```

4. **Parsing the Metadata**:  
   Run the provided parsing script (adapted and distributed in the `extras` directory of **metacooc**) to generate a cleaned TSV file:
   ```bash
   utils/parse_merge_sra_metadata.py --json-dir path/to_json_dir \
									 --output-file metacooc/data/sra_metadata_parsed.tsv \
									 --threads 4
   ```

### Downloading Sandpiper Data

The sandpiper file is downloaded from Zenodo. For example, use the following URL:
```
https://zenodo.org/records/11516218/files/sandpiper0.3.0.condensed.csv.gz?download=1
```
After downloading, decompress and process it as needed.

### Format

Format your raw data into the standardized Ingredients objects and build metadata indices. For example, to process Sandpiper output and NCBI SRA metadata:

```bash
metacooc format --tax_profile metacooc/data/sandpiper0.3.0.condensed_intersection.csv \
                --metadata_file metacooc/data/sra_metadata_parsed.tsv \
                --output_dir metacooc/data/ \
                --index_type broad \
                --strict \
                --aggregated \
                --aggregation_pattern "g__"
```

This command processes your raw files, generating both the raw and aggregated Ingredients objects as well as metadata indices that will be used by subsequent commands.

## Default Data Directory

If you do not supply `--data_dir`, **MetaCoOc** uses a default directory relative to the package installation. This is set using `importlib.resources` in the CLI, so that your data files are always located relative to the installed package.

## License

**MetaCoOc** is distributed under the GNU General Public License v3 or later (GPLv3+). See [LICENCE.txt](LICENCE.txt) for details.

## Acknowledgements

The NCBI SRA metadata parsing script is adapted from [public_sequencing_metadata_corrections](https://github.com/wwood/public_sequencing_metadata_corrections) by [W. Wood](https://github.com/wwood).

## Contact

Benjamin Coltman - [benjamin.coltman@univie.ac.at](benjamin.coltman@univie.ac.at)
Daan Speth - [daan.speth@univie.ac.at](daan.speth@univie.ac.at)    

