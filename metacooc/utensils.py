#!/usr/bin/env python3

class Recipes:
    def __init__(self, metadata_tsv):
        self.metadata_df = pd.read_csv(metadata_tsv, delimiter='\t', dtype=str, engine="c")
        self.filtered_metadata_df = self.metadata_df.copy()
    
    def filter_by_accessions(self, accessions):
        self.filtered_metadata_df = self.metadata_df[self.metadata_df['acc'].isin(accessions)]
    
    def filter_by_search_terms(self, search_terms, columns_to_search=None):
        search_columns = columns_to_search if columns_to_search else self.metadata_df.columns.tolist()
        
        mask = pd.Series(False, index=self.metadata_df.index)
        for col in search_columns:
            for term in search_terms:
                mask |= self.metadata_df[col].str.contains(term, case=False, na=False)
        
        self.filtered_metadata_df = self.metadata_df[mask]
    
    def get_filtered_accessions(self):
        return self.filtered_metadata_df['acc'].unique().tolist()

class Ingredients:
    def __init__(self, sandpiper_tsv):
        # Load the data and create the sparse matrices
        self.samples, self.taxa = self._build_indices(sandpiper_tsv)
        self.presence_matrix, self.coverage_matrix = self._create_sparse_matrices(sandpiper_tsv)
        
        # Dictionary to store filter indices
        self.filters = {
            'unfiltered_samples': np.arange(len(self.samples)),
            'filtered_by_accessions': np.arange(len(self.samples)),
            'filtered_by_taxa_count': np.arange(len(self.samples)),
            'combined_filter': np.arange(len(self.samples)),  
            'filtered_taxa': np.arange(len(self.taxa))
            
        }
        
    def _build_indices(self, sandpiper_tsv):
        # Build sample and taxon indices
        chunk_iterator = pd.read_csv(sandpiper_tsv, delimiter='\t', dtype=str, usecols=['sample', 'taxonomy'], chunksize=100000, engine="c")
        
        sample_to_index = {}
        taxon_to_index = {}
        
        for chunk in chunk_iterator:
            for sample in chunk['sample'].unique():
                if sample not in sample_to_index:
                    sample_to_index[sample] = len(sample_to_index)
            for taxon in chunk['taxonomy'].unique():
                if taxon not in taxon_to_index:
                    taxon_to_index[taxon] = len(taxon_to_index)
        
        # Convert index mappings to lists
        samples = [None] * len(sample_to_index)
        for sample, idx in sample_to_index.items():
            samples[idx] = sample
        
        taxa = [None] * len(taxon_to_index)
        for taxon, idx in taxon_to_index.items():
            taxa[idx] = taxon
        
        # Store mappings for later use
        self.sample_to_index = sample_to_index
        self.taxon_to_index = taxon_to_index
        
        return samples, taxa
    
    def _create_sparse_matrices(self, sandpiper_tsv):
        # Initialize lists for sparse matrix construction
        data_presence = []
        data_coverage = []
        row_indices = []
        col_indices = []
        
        chunk_iterator = pd.read_csv(sandpiper_tsv, delimiter='\t', dtype=str, chunksize=100000, engine="c")
        
        for chunk in chunk_iterator:
            # Map samples and taxa to indices
            sample_indices = chunk['sample'].map(self.sample_to_index)
            taxon_indices = chunk['taxonomy'].map(self.taxon_to_index)
            
            # Convert coverage to numeric
            coverage_values = pd.to_numeric(chunk['coverage'], errors='coerce').fillna(0)
            
            # Append data
            data_presence.extend([1] * len(chunk))
            data_coverage.extend(coverage_values)
            row_indices.extend(sample_indices)
            col_indices.extend(taxon_indices)
        
        num_samples = len(self.samples)
        num_taxa = len(self.taxa)
        
        # Create presence/absence sparse matrix
        presence_matrix = csr_matrix((data_presence, (row_indices, col_indices)), shape=(num_samples, num_taxa), dtype=int)
        
        # Create coverage sparse matrix
        coverage_matrix = csr_matrix((data_coverage, (row_indices, col_indices)), shape=(num_samples, num_taxa), dtype=float)
        
        return presence_matrix, coverage_matrix
    
    def filter_samples_by_taxa_string(self, taxa_list, filter_type='unfiltered_samples', combine=False):
        """
        Filters samples based on the presence of taxa strings provided in `taxa_list`.
        
        Parameters:
        - taxa_list: List of taxa strings to search for in the taxa names.
        - filter_type: Specifies the filter type to apply before filtering by taxa string.
        """
        # Step 1: Identify indices of taxa that match any of the strings in taxa_list
        organism_indices = [
            idx for idx, taxon in enumerate(self.get_filtered_taxa())
            if any(org in taxon for org in taxa_list)
        ]
        
        # Step 2: Check presence of these taxa in each sample and filter samples accordingly
        presence_matrix = self.get_presence_matrix(filter_type=filter_type)
        organism_presence = presence_matrix[:, organism_indices].sum(axis=1).A.flatten() > 0
        filtered_sample_indices = [i for i, present in enumerate(organism_presence) if present]
    
        # Step 3: Update filters with the new sample indices
        self.filters['filtered_by_accessions'] = np.array(filtered_sample_indices)
        
        if combine:
            # Compute the combined filter as the intersection of filtered_by_accessions and filtered_by_taxa_count
            filtered_sample_set = set(self.filters['filtered_by_accessions'])
            valid_sample_set = set(self.filters['filtered_by_taxa_count'])
            self.filters['combined_filter'] = np.array(list(filtered_sample_set.intersection(valid_sample_set)))
        
    def filter_by_accessions(self, accessions, combine=False):
        # Map accessions to indices
        accession_to_index = self.sample_to_index
        self.filters['filtered_by_accessions'] = np.array([accession_to_index[acc] for acc in accessions if acc in accession_to_index])

        if combine:
            # Compute the combined filter as the intersection of filtered_by_accessions and filtered_by_taxa_count
            filtered_sample_set = set(self.filters['filtered_by_accessions'])
            valid_sample_set = set(self.filters['filtered_by_taxa_count'])
            self.filters['combined_filter'] = np.array(list(filtered_sample_set.intersection(valid_sample_set)))
        
    def filter_taxa_by_sample_count(self, filter_type='unfiltered_samples', min_count_samples=0):
        # Apply taxa threshold if min_sample_count > 0
        if min_count_samples > 0:
            presence_matrix = self.presence_matrix[self.filters[filter_type], :]
            taxa_counts = np.array((presence_matrix > 0).sum(axis=0)).flatten()
            valid_taxa_indices = np.where(taxa_counts >= min_count_samples)[0]
            self.filters['filtered_taxa'] = valid_taxa_indices
        else:
            self.filters['filtered_taxa'] = np.arange(len(self.taxa))
    
    def filter_taxa_level(self, level="s__", filter_type='unfiltered_samples', combine=False):

        # Determine which taxa to use based on the combine flag
        if combine: 
            filtered_taxa = self.get_filtered_taxa()
        else:
            filtered_taxa = self.taxa
        
        # Get indices of taxa containing the specified level
        valid_taxa_indices = np.array([
            idx for idx, taxon in enumerate(self.taxa)
            if level in taxon and (not combine or taxon in filtered_taxa)
        ])
        
        # Set 'filtered_taxa' filter to be the indices of valid taxa
        self.filters['filtered_taxa'] = valid_taxa_indices
        
    def filter_sample_by_taxa_count(self, min_count_taxa=0, combine=False):
        # Filter samples by taxa count
        presence_matrix = self.get_presence_matrix(filter_type='unfiltered_samples')
        sample_taxa_counts = np.array(presence_matrix.sum(axis=1)).flatten()
        valid_sample_indices = np.where(sample_taxa_counts >= min_count_taxa)[0]
        self.filters['filtered_by_taxa_count'] = valid_sample_indices
        
        if combine:
            # Compute the combined filter as the intersection of filtered_by_accessions and filtered_by_taxa_count
            filtered_sample_set = set(self.filters['filtered_by_accessions'])
            valid_sample_set = set(self.filters['filtered_by_taxa_count'])
            self.filters['combined_filter'] = np.array(list(filtered_sample_set.intersection(valid_sample_set)))
        
    def get_presence_matrix(self, filter_type='combined_filter'):
        sample_indices = self.filters[filter_type]
        return self.presence_matrix[sample_indices, :][:, self.filters['filtered_taxa']]
        
    def get_filtered_coverage_matrix(self, filter_type='combined_filter'):
        sample_indices = self.filters[filter_type]
        return self.coverage_matrix[sample_indices, :][:, self.filters['filtered_taxa']]
    
    def get_filtered_samples(self, filter_type='combined_filter'):
        indices = self.filters[filter_type]
        return [self.samples[i] for i in indices]
    
    def get_filtered_taxa(self):
        # Return the list of filtered taxa
        return [self.taxa[i] for i in self.filters['filtered_taxa']]
    
    def get_filtered_sandpiper_df(self, filter_type='combined_filter'):
        sample_indices = self.filters[filter_type]
        presence_matrix_filtered = self.get_presence_matrix(filter_type)
        coo = presence_matrix_filtered.tocoo()
        data = {
            'sample': [self.samples[sample_indices[i]] for i in coo.row],
            'taxonomy': [self.taxa[self.filters['filtered_taxa'][j]] for j in coo.col],
            'presence': coo.data
        }
        filtered_df = pd.DataFrame(data)
        return filtered_df
        
    def add_taxa_levels_to_presence_matrix(self, pattern="g__"):
        """
        Adds additional taxa for levels matching the specified pattern in the taxa strings.
        Updates the presence matrix to include higher taxonomic levels.
        
        Parameters:
        - pattern: Regex pattern to identify specific levels (e.g., "g__" for genus level).
        """
        
        # Define the columns for taxonomic levels
        level_names = ["Root", "d__", "p__", "c__", "o__", "f__", "g__", "s__"]
        
        # Split taxonomic strings into a DataFrame with each level as a separate column
        taxa_series = pd.Series(self.taxa)
        taxa_df = taxa_series.str.split("; ", expand=True)
        taxa_df.columns = level_names[:taxa_df.shape[1]]  # Limit columns to the existing levels
        
        # Identify the target level based on the pattern
        target_level = next((level for level in level_names if re.fullmatch(pattern, level)), None)
        if not target_level or target_level not in taxa_df.columns:
            print(f"No matching level found for pattern: {pattern}")
            return
        
        # Get the higher-level taxa for each taxon
        higher_level_taxa = taxa_df[target_level]
        
        # Identify unique higher-level taxa
        unique_levels = higher_level_taxa.dropna().unique()
        
        # Map higher-level taxa to indices
        level_to_index = {level: idx for idx, level in enumerate(unique_levels)}
        
        # Create the mapping matrix T (num_taxa x num_unique_levels)
        # For each taxon, set T[taxon_idx, higher_level_idx] = 1
        num_taxa = len(self.taxa)
        num_unique_levels = len(unique_levels)
        row_indices = []
        col_indices = []
        data = []
        
        for taxon_idx, level in enumerate(higher_level_taxa):
            if pd.isnull(level):
                continue  # Skip if level is NaN
            level_idx = level_to_index[level]
            row_indices.append(taxon_idx)
            col_indices.append(level_idx)
            data.append(1)
        
        T = csr_matrix((data, (row_indices, col_indices)), shape=(num_taxa, num_unique_levels), dtype=int)
        
        # Compute the new presence matrix L (num_samples x num_unique_levels)
        P = self.presence_matrix  # Original presence matrix (num_samples x num_taxa)
        L = P @ T  # Matrix multiplication
        
        # Binarize L (set all positive entries to 1)
        L.data = np.ones_like(L.data)
        
        # Extend self.taxa with new higher-level taxa
        self.taxa.extend(unique_levels)
        
        # Extend the presence matrix by concatenating L to P
        self.presence_matrix = hstack([P, L], format='csr')
        self.filters['filtered_taxa'] = np.arange(len(self.taxa))


class RatioAnalyzer:
    def __init__(self, sandpiper_data, total_counts_filter='unfiltered_samples'):
        self.sandpiper_data = sandpiper_data
        self.ratios_df = None
        self.total_counts = self.calculate_counts(filter_type=total_counts_filter) # Store total counts to avoid recomputation
        self.filtered_taxa = self.sandpiper_data.get_filtered_taxa()
        
    def calculate_counts(self, filter_type='unfiltered_samples'):
        # Calculate counts per taxon in the filtered dataset (filtered samples and filtered taxa)
        presence_matrix = self.sandpiper_data.get_presence_matrix(filter_type=filter_type)
        counts = np.array(presence_matrix.sum(axis=0)).flatten()
        return counts
        
    def compute_taxon_ratios(self, filter_type='filtered_by_accessions'):
        
        # Step 1: Calculate filtered counts
        filtered_counts = self.calculate_counts(filter_type=filter_type)

        # Step 2: Compute ratios
        with np.errstate(divide='ignore', invalid='ignore'):
            ratios = np.true_divide(filtered_counts, self.total_counts)
            ratios[~np.isfinite(ratios)] = 0  # Replace NaN and inf with 0
        
        # return # Step 3: Get the list of filtered taxa
        # filtered_taxa = self.sandpiper_data.get_filtered_taxa()
        
        # Step 4: Create DataFrame with the computed ratios
        self.ratios_df = pd.DataFrame({
            'taxon': self.filtered_taxa,
            'filtered_counts': filtered_counts,
            'all_counts': self.total_counts,
            'ratio': ratios
        })
        
        
    def plot_ratios(self, output_plot_file=None, threshold=None):
        # Sort the dataframe by ratio in descending order
        ratios_df_sorted = self.ratios_df.sort_values(by='ratio', ascending=False).reset_index(drop=True)
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(ratios_df_sorted['ratio'].values, marker='o', markersize=5)
        ax.set_xlabel('Taxa (ordered by ratio)')
        ax.set_ylabel('Ratio')
        ax.set_title('Taxon Ratios')
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True)
        
        # Mark the threshold if provided
        if threshold is not None:
            # Find the index where the ratio falls below the threshold
            threshold_indices = ratios_df_sorted[ratios_df_sorted['ratio'] >= threshold].index
            if len(threshold_indices) > 0:
                cutoff_index = threshold_indices[-1]
                ax.axvline(x=cutoff_index, color='red', linestyle='--', label=f'Threshold = {threshold}')
                ax.legend()
            else:
                print(f"No taxa found with ratio >= {threshold}")
        
        # Save or show the plot
        fig.tight_layout()
        if output_plot_file:
            fig.savefig(output_plot_file)
            plt.close(fig)
        else:
            plt.show()
        
    def filter_ratios(self, ratio_threshold, counts_threshold): #, taxa_to_exclude=None):
        # Filter the ratios based on the threshold
        filtered_ratios_df = self.ratios_df[self.ratios_df['ratio'] >= ratio_threshold]
        filtered_ratios_df = filtered_ratios_df[filtered_ratios_df['all_counts'] >= counts_threshold]
        return filtered_ratios_df
