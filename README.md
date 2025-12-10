# MetaCoOc

**MetaCoOc** is a Python package for co-occurrence analysis of microorganisms in shotgun metagenomes. It provides a suite of command‑line tools that allow you to download necessary files, format raw input data, perform metadata and taxonomic searches, filter datasets, calculate co‑occurrence ratios, and generate plots for data visualization. The package is designed to work both as a step‑by‑step workflow (with individual commands) and as a full in‑memory pipeline (using the `cooccurrence` subcommand).

## WARNING: This repository is in active development. Feel free to explore the package and please share any feedback or issues you encounter


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

## Usage

### Download

Before any analysis, download the default data files. MetaCoOc sets a default data directory relative to its installation. These include:
- `ingredients_raw_1.0.0.pkl`
- `ingredients_aggregated_1.0.0.pkl`
- `sra_metadata_1.0.0.tsv`

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
                      --output_dir /path/to/output_dir \
                      --aggregated \
                      --search_string "Nitrospira" \
                      --ranks_for_search_inclusion genus \
                      --min_taxa_count 5 \
                      --min_sample_count 5
```

You can also include a ratio threshold (e.g., only include taxa with ratio ≥ 0.5) by adding `--ratio_threshold 0.5`.

### Running Individual Steps

You can also run each step separately. For example:

#### Search (File-Based)

We filter all metagenomes based on the word soil being listed in their metadata
```bash
metacooc search --mode metadata \
                --output_dir /path/to/output_dir \
                --search_string "soil" \
                --strict
```
					 
#### Filter

We then filter our Ingredients based on the resulting accession list

```bash
metacooc filter --accessions_file /path/to/output_dir/search_results.txt \
                --output_dir /path/to/output_dir \
                --aggregated \
                --min_taxa_count 5 \
                --min_sample_count 5
```
					 
#### Ratios

The last function created two Ingredients. One is an intermediate file, filtered according to:
- minimum taxa counts - which is a threshold for sample inclusion based on the minimum unique taxa in sample
- minimum sample counts - which is a threshold for taxa inclusion based on the number of samples it is detected in 

The other is a further filtered version of this file, filtering out the samples that didn't contain "soil" in their metadata


We use this intermediate file as our reference and compare our filtered file to it
```bash
metacooc ratio --output_dir /path/to/output_dir  \
               --ratio_threshold 0.5 \
			   --reference_file /path/to/output_dir/ingredients_counts_filtered.pkl \
			   --filtered_file /path/to/output_dir/ingredients_all_filtered.pkl
```

#### Plot
Finally, plot it.

```bash
metacooc plot --ratios_file /path/to/output_dir/ratios.tsv \
              --output_dir /path/to/output_dir \
              --ratio_threshold 0.5
```

## Some alternative ways to use MetaCoOc

### Using previous versions of Ingredients

Users can use previous versions of the zenodo hosted Ingredients by specifying the version number with the --sandpiper_version flag e.g.

```
metacooc cooccurrence --mode taxon \
                      --output_dir /path/to/output_dir \
                      --search_string "Nitrospira" \
                      --sandpiper_version "0.3.0"
```

### Custom directory for Ingredients files

If you have downloaded and setup the Ingredients files in a custom directory, then specify this location with --data_dir. This can be specified for all subfunctions e.g.

By default, it will still search for the Ingredients files which are specified within _data_config.py, which is updated upon a new release of metacooc

```bash
metacooc search --mode metadata \
                --data_dir /path/to/data_dir \
				--output_dir /path/to/output_dir \
                --search_string "soil" \
                --strict
```

### Using own custom_ingredients files

Cooccurrence, search and filter can all take a custom_ingredients file. This can be specified using the --custom_ingredients flag
```bash
metacooc cooccurrence --mode taxon \
                      --output_dir /path/to/output_dir \
                      --aggregated \
                      --search_string "Nitrospira" \
                      --ranks_for_search_inclusion genus \
                      --min_taxa_count 5 \
                      --min_sample_count 5 \
                      --custom_ingredients /path/to/custom_ingredients.pkl

metacooc search --mode metadata \
                --output_dir /path/to/output_dir \
                --search_string "soil" \
                --strict \
                --custom_ingredients /path/to/custom_ingredients.pkl

metacooc filter --accessions_file /path/to/output_dir/search_results.txt \
                --output_dir /path/to/output_dir \
                --aggregated \
                --min_taxa_count 5 \
                --min_sample_count 5 \
                --custom_ingredients /path/to/custom_ingredients.pkl
```


## How we generated and prepared the data on [Zenodo](https://doi.org/10.5281/zenodo.15283587)

### Community profiles
#### Downloading Community profiles (Sandpiper Data)

Sandpiper hosts community profiles generated by annotating public metagenome datasets using SingleM.

The data is available to download on Zenodo and the below code highlights how the data was formatted to work with MetaCoOc. 

1. **Download and unzip data**
```bash
curl https://zenodo.org/records/15233363/files/sandpiper1.0.0.condensed.csv.gz?download=1 -o sandpiper1.0.0.condensed.csv.gz
gzip -d sandpiper1.0.0.condensed.csv.gz
```

2. **Format profiles**

Format your raw data into the standardized Ingredients objects. The **format** command processes your raw files, generating both the raw and aggregated Ingredients objects that will be used by **MetaCoOc** commands.For example, to process Sandpiper output (two column csv):

```bash

metacooc format --tax_profile metacooc/data/sandpiper1.0.0.condensed.csv \
					--output_dir metacooc/data/ \
					--aggregated \
					--sample_to_biome_file metacooc/data/sample_to_biome_1.0.2.tsv \
					--tag "1.0.2""
```

### Metadata associated with metagenomes

A key component of **MetaCoOc** is the standardized parsing of SRA metadata. First, we identify the accessions reported in Sandpiper and then use these accessions to retrieve their submission metadata from the NCBI SRA

```bash
tail -n+2 metacooc/data/sandpiper1.0.0.condensed.csv | cut -f1 - | uniq > metacooc/data/sandpiper1.0.0.uniqueaccessions.txt
```

##### A. Approach 1
This parsing workflow was adapted from [public_sequencing_metadata_corrections](https://github.com/wwood/public_sequencing_metadata_corrections). The process involves:

1. **Exporting Data from BigQuery**:
	```sql
	-- 1) Define the external table (one‐time setup; put this in your project/dataset)
	CREATE EXTERNAL TABLE `project.dataset.acc_list` (
	  acc STRING
	)
	OPTIONS (
	  format       = 'CSV',           
	  uris         = ['gs://my-bucket-name/sandpiper1.0.0.uniqueaccessions.txt'],
	  skip_leading_rows = 0                         -- no header row
	);

	-- 2) Use it in your EXPORT DATA
	EXPORT DATA
	  OPTIONS(
		uri    = 'gs://my-bucket-name/sra_metadata/*.json',
		format = 'JSON'
	  ) AS
	SELECT
	  m.*
	FROM
	  `nih-sra-datastore.sra.metadata` AS m
	  INNER JOIN `project.dataset.acc_list` AS a
		ON m.acc = a.acc
	WHERE
	  m.platform = 'ILLUMINA';
	```

2. **Downloading the Data**:  
   Use `gsutil` to download the JSON files locally:
   ```bash
   gsutil -m cp -r gs://my-bucket-name/sra_metadata/*.json .
   ```
   Once downloaded, you can delete the files from the bucket to avoid storage costs.

3. **Setting Up the Parsing Environment**:  
   Create and activate a conda environment:
   ```bash
   conda create -n metacooc
   conda activate metacooc
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   conda install python==3.13.2 pandas==2.2.3 scipy==1.15.2 matplotlib==3.10.1 iso8601==2.1.0 requests==2.32.3 
   conda install -c bioconda kingfisher==0.4.1
   ```

4. **Parsing the Metadata**:  
   Run the provided parsing script (adapted and distributed in the `extras` directory of **metacooc**) to generate a cleaned TSV file. 
   WARNING: Can take a while as we process and keep all atributes from the NCBI SRA. We chose to use all attributes due to the inconsistency of the metadata attached with samples uploaded to the SRA.
   We found that relevant information for filtering is often included in a wide range of non-standardised columns and therefore chose to retain all of it.
   ```bash
   utils/parse_merge_sra_metadata.py --json-dir path/to_json_dir \
									 --output-file metacooc/data/sra_metadata.tsv \
									 --threads 4
   ```

##### B. Approach 2.

a. **Query the SRA for the annotation using [kingfisher](https://wwood.github.io/kingfisher-download/)**
	This can take a while
	```bash
	tail -n+2 metacooc/data/sandpiper1.0.0.condensed.csv | cut -f1 - | uniq > metacooc/data/sandpiper1.0.0.uniqueaccessions.txt

	kingfisher annotate --run-identifiers-list metacooc/data/sandpiper1.0.0.uniqueaccessions.txt \
		-f csv \
		--all-columns \
		-o metacooc/data/sra_metadata_1.0.0.tsv
	```
	
	Can run multiple processes but anymore than 2 and you'll likely hit the API limits:
	```bash
	TMPDIR=$(mktemp -d)
	split -l 10000 metacooc/data/sandpiper1.0.0.uniqueaccessions.txt $TMPDIR/batch_
	
	export TMPDIR  # So it's available to the subshells
	
	find "$TMPDIR" -name 'batch_*' | parallel -j 2 '
		batch_file={}
		out_file="${batch_file}.csv"
		kingfisher annotate --run-identifiers-list "$batch_file" -f csv --all-columns -o "$out_file"'
	```

##### C. Approach 3 - original

Approach 3 retrieves metadata of all submissions passing its filters and was also adapted from [public_sequencing_metadata_corrections](https://github.com/wwood/public_sequencing_metadata_corrections). This process differs as it requires later filtering of the sra_metadata to remove accessions not included in sandpiper.

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
   conda create -n metacooc
   conda activate metacooc
   conda config --add channels conda-forge
   conda config --set channel_priority strict
   conda install python==3.13.2 pandas==2.2.3 scipy==1.15.2 matplotlib==3.10.1 iso8601==2.1.0 requests==2.32.3 
   conda install -c bioconda kingfisher==0.4.1
   ```

4. **Parsing the Metadata**:  
   Run the provided parsing script (adapted and distributed in the `extras` directory of **metacooc**) to generate a cleaned TSV file. 
   WARNING: Can take a while as we process and keep all atributes from the NCBI SRA. We chose to use all attributes due to the inconsistency of the metadata attached with samples uploaded to the SRA.
   We found that relevant information for filtering is often included in a wide range of non-standardised columns and therefore chose to retain all of it.
   ```bash
   utils/parse_merge_sra_metadata.py --json-dir path/to_json_dir \
									 --output-file metacooc/data/sra_metadata.tsv \
									 --threads 4
   ```

5. **Filtering the Metadata**
	Community profiles are not included in Sandpiper for all SRA entries. To reduce the size of the metacooc/data/sra_metadata_parsed.tsv file, uploaded to Zenodo, we excluded accessions without an entry in Sandpiper.
	```bash
	awk -F'\t' 'NR > 1 { print $1 }' metacooc/data/sra_metadata.tsv | uniq > metacooc/data/sra_metadata_parsed_uniqueaccessions.txt
	
	awk 'NR==FNR { a[$1]; next } $1 in a' \
		metacooc/data/sandpiper1.0.0.uniqueaccessions.txt \
		metacooc/data/sra_metadata_parsed_uniqueaccessions.txt \
		> metacooc/data/intersection_uniqueaccessions.txt
		
	awk 'NR==FNR { a[$1]; next } $1 in a' \
		metacooc/data/sandpiper1.0.0.uniqueaccessions.txt \
		metacooc/data/sra_metadata_parsed_uniqueaccessions.txt \
		> metacooc/data/intersection_uniqueaccessions.txt
		
	awk -F'\t' 'NR == FNR { lookup[$1]; next } FNR == 1 || ($1 in lookup)' \
		metacooc/data/intersection_uniqueaccessions.txt \
		metacooc/data/sra_metadata_parsed.tsv \
		> metacooc/data/sra_metadata_1.0.0.tsv
	```

## Default Data Directory

If you do not supply `--data_dir`, **MetaCoOc** uses a default directory relative to the package installation. This is set using `importlib.resources` in the CLI, so that your data files are always located relative to the installed package.

## License

**MetaCoOc** is distributed under the GNU General Public License v3 or later (GPLv3+). See [LICENCE.txt](LICENCE.txt) for details.

## Acknowledgements

The NCBI SRA metadata parsing script is adapted from [public_sequencing_metadata_corrections](https://github.com/wwood/public_sequencing_metadata_corrections) by [W. Wood](https://github.com/wwood).

## Contact

Benjamin Coltman - [benjamin.coltman@univie.ac.at](benjamin.coltman@univie.ac.at)
Daan Speth - [daan.speth@univie.ac.at](daan.speth@univie.ac.at)    

