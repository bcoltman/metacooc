# MetaCoOccur

A package to search and process tables.

## Installation

### Using pip

pip install table_searcher






# Step 1: Extract unique values and save them to a temporary file
awk -F'\t' 'NR > 1 { print $1 }' data/sandpiper0.3.0.condensed.csv | uniq > data/sandpiper0.3.0.uniqueaccessions.txt

awk -F'\t' 'NR > 1 { print $1 }' data/sra_metadata_parsed.tsv | uniq > data/sra_metadata_parsed_uniqueaccessions.txt

awk 'NR==FNR { a[$1]; next } $1 in a' \
	data/sandpiper0.3.0.uniqueaccessions.csv \
	data/sra_metadata_parsed_uniqueaccessions.txt \
	> data/intersection_uniqueaccessions.txt
	
	
# comm -12 data/sandpiper0.3.0.uniqueaccessions.txt data/sra_metadata_parsed_uniqueaccessions.txt > data/intersection_uniqueaccessions.txt

# Step 2: Filter file2.txt based on values in a specific column
awk -F'\t' 'NR == FNR { lookup[$1]; next } FNR == 1 || ($1 in lookup)' \
	data/intersection_uniqueaccessions.txt \
	data/sra_metadata_parsed.tsv \
	> data/sra_metadata_parsed_intersection.tsv
	
awk -F'\t' '{ for (i=1; i<=44; i++) printf "%s\t", $i; print "" }' \
	data/sra_metadata_parsed_intersection.tsv \
	> data/sra_metadata_parsed_intersection_reduced.tsv


# Step 2: Filter file2.txt based on values in a specific column
awk -F'\t' 'NR == FNR { lookup[$1]; next } FNR == 1 || ($1 in lookup)' \
	data/intersection_uniqueaccessions.txt \
	data/sandpiper0.3.0.condensed.csv \
	> data/sandpiper0.3.0.condensed_intersection.csv
	


module load conda
conda activate /lisc/scratch/dome/coltman/metacooccur/env/metacooccur/
cd /lisc/scratch/dome/coltman/metacooccur/
# jupyter notebook --no-browser --port=8080
jupyter notebook MetaCoOccur.ipynb --no-browser --port=8080





# create script for executing on SLURM to calculate all ratios
load sandpiper and metadata
look at general filters e.g. over 5



# currently doesn't prevent sample sform time series etc. Doesn't exclude time series etc
# could implement something using the union of samples to determine uniqueness of samples


# # Step 2: Filter file2.txt based on values in a specific column
# awk -F'\t' 'NR == FNR { lookup[$1]; next } FNR == 1 || ($1 in lookup)' \
	# data/sandpiper0.3.0.uniqueaccessions.csv \
	# data/sra_metadata_parsed.tsv \
	# > data/sra_metadata_parsed_sandpiper.tsv

```python
import copy
import pandas as pd
from multiprocessing import Pool, Manager

# Define a global variable to hold the RatioAnalyzer instance per pool worker
worker_ratio_analyzer = None

# Initializer function for each worker to set up its own copy of RatioAnalyzer
def init_worker(base_ratio_analyzer):
    global worker_ratio_analyzer
    worker_ratio_analyzer = copy.deepcopy(base_ratio_analyzer)  # Create a single copy per worker

# Define the processing function to use the worker-specific RatioAnalyzer
def process_species(species):
    global worker_ratio_analyzer
    # Filter samples by the current species
    worker_ratio_analyzer.sandpiper_data.filter_samples_by_taxa_string(taxa_list=[species], combine=True)
    
    # Compute taxon ratios and filter them based on thresholds
    worker_ratio_analyzer.compute_taxon_ratios(filter_type='combined_filter')
    filtered_ratios = worker_ratio_analyzer.filter_ratios(ratio_threshold=0.85, counts_threshold=5)
    
    # Add a column to identify the comparator species and return the result
    filtered_ratios.loc[:, "comparator"] = species
    return filtered_ratios

# Initialize and filter sandpiper_data as in your original code
sd = copy.deepcopy(sandpiper_data)
sd.add_taxa_levels_to_presence_matrix(pattern="g__")
sd.filter_taxa_by_sample_count(min_count_samples=10)

all_species = [taxon for taxon in sd.get_filtered_taxa() if "s__" in taxon]
print(f"{len(all_species)} species")

sd.filter_sample_by_taxa_count(min_count_taxa=10)
base_ratio_analyzer = RatioAnalyzer(sd, total_counts_filter='filtered_by_taxa_count')

# Set up the multiprocessing pool with the initializer function
num_workers = 12  # Number of workers you want in the pool
with Pool(processes=num_workers, initializer=init_worker, initargs=(base_ratio_analyzer,)) as pool:
    # Process each species using the pool
    all_ratios = pool.map(process_species, all_species)

# Concatenate all the results from each worker into a single DataFrame
all_ratios = pd.concat(all_ratios)
all_ratios
```
