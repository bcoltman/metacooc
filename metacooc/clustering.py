import numpy as np
from typing import Optional, Iterable, Tuple, List, Set
from metacooc.search import _search_taxon_columns

def cluster_taxa(
    ingredients,
    min_cooccurrence_count: int = 2,
    min_sample_partner_count: int = 0,
    sample_shared_taxa_threshold: int = 1,
    min_taxa_partner_count: int = 0
) -> List[List[str]]:
    """
    Partition taxa into clusters based on co-occurrence thresholds,
    removing low-overlap samples and excluding taxa absent in all samples.
    
    Workflow:
    1. Compute sample × sample co-occurrence counts (number of shared taxa).
    2. Optionally remove samples that:
       - Overlap with fewer than `min_sample_partner_count` other samples
         (sharing at least `sample_shared_taxa_threshold` taxa each).
       - Share fewer than `min_taxa_partner_count` taxa with any single sample.
    3. Remove taxa not present in any remaining samples.
    4. Build an undirected graph on the remaining taxa: an edge exists if
       two taxa co-occur in at least `min_cooccurrence_count` samples.
    5. Identify connected components as clusters.
    
    Parameters
    ----------
    min_cooccurrence_count : int, optional
        Minimum number of shared samples required to connect two taxa.
        Default is 2.
    min_sample_partner_count : int, optional
        Minimum number of other samples a sample must overlap with to be retained.
        Overlap here means sharing at least `sample_shared_taxa_threshold` taxa.
        Default is 0 (keep all samples).
    sample_shared_taxa_threshold : int, optional
        Minimum number of taxa two samples must share to count as overlapping.
        Default is 1 (original behavior).
    min_taxa_partner_count : int, optional
        Minimum number of taxa a sample must share with its best-matching partner
        sample to be retained. Default is 0.
        
    Returns
    -------
    List[List[str]]
        A list of taxon clusters (connected components) based on the specified
        co-occurrence threshold. Samples and taxa that fail filters are excluded
        prior to clustering.
    """
    # 1) Sample × sample co-occurrence (shared taxa counts)
    samp_coocc = ingredients._presence_matrix @ ingredients._presence_matrix.T
    samp_coocc.setdiag(0)
    
    # Count how many samples each sample overlaps with (shares ≥ threshold taxa)
    overlaps = (samp_coocc >= sample_shared_taxa_threshold).sum(axis=1).A1
    # Compute the maximum number of taxa shared with any single sample
    max_shared_taxa = samp_coocc.max(axis=1).toarray().flatten()
    
    # 2) Filter samples by partner count and taxa depth
    sample_mask = np.ones(len(ingredients.samples), dtype=bool)
    if min_sample_partner_count > 0:
        sample_mask &= overlaps >= min_sample_partner_count
    if min_taxa_partner_count > 0:
        sample_mask &= max_shared_taxa >= min_taxa_partner_count
    if not np.any(sample_mask):
        ingredients._clusters_cache = []
        ingredients._clusters_valid = True
        return []
    matrix = ingredients._presence_matrix[sample_mask, :]
    
    # 3) Filter taxa by presence
    taxa_counts = np.array((matrix > 0).sum(axis=0)).flatten()
    taxon_mask = taxa_counts > 0
    if not np.any(taxon_mask):
        ingredients._clusters_cache = []
        ingredients._clusters_valid = True
        return []
    matrix = matrix[:, taxon_mask]
    pruned_taxa = [t for t, keep in zip(ingredients.taxa, taxon_mask) if keep]
    
    # 4) Build taxa co-occurrence graph
    coocc = matrix.T @ matrix
    coocc.setdiag(0)
    adj = (coocc >= min_cooccurrence_count).astype(int)
    adj.eliminate_zeros()
    
    # 5) Extract connected components
    n_comp, labels = connected_components(csgraph=adj, directed=False)
    clusters = [[] for _ in range(n_comp)]
    for idx, lbl in enumerate(labels):
        clusters[lbl].append(pruned_taxa[idx])
    
    # Cache & mark valid
    ingredients._clusters_cache = clusters
    ingredients._clusters_valid = True
    return clusters
 

def _resolve_focal_taxa_indices(ingredients, focal_taxa):
    """
    Resolve focal_taxa (indices, full tax strings, or rank tokens) → np.ndarray[int].
    
    Uses the same deepest-token logic as search_by_taxon, via _search_taxon_columns.
    """
    # Normalise to list
    if isinstance(focal_taxa, (str, int, np.integer)):
        focal_taxa = [focal_taxa]
    elif isinstance(focal_taxa, np.ndarray):
        focal_taxa = focal_taxa.tolist()
        
    ingredients._ensure_taxa_lookups()
    n_taxa = len(ingredients.taxa)
    
    indices = set()
    not_found = []
    
    for t in focal_taxa:
        # --- numeric indices ---
        if isinstance(t, (int, np.integer)):
            i = int(t)
            if i < 0 or i >= n_taxa:
                raise IndexError(f"Taxon index {i} out of range [0, {n_taxa - 1}]")
            indices.add(i)
            continue
            
        if not isinstance(t, str):
            not_found.append(t)
            continue
            
        # 1) try exact taxa string match
        try:
            idx = ingredients.taxa.index(t)
            indices.add(idx)
            continue
        except ValueError:
            pass
            
        # 2) fall back to taxonomy-based search using deepest ranked token
        col_idxs = _search_taxon_columns(ingredients, t)
        if col_idxs:
            indices.update(col_idxs)
        else:
            not_found.append(t)
            
    if not_found:
        raise KeyError(f"Could not resolve taxa: {', '.join(map(str, not_found))}")
        
    return np.array(sorted(indices), dtype=int)


import numpy as np

def determine_taxa_context(
    ingredients,
    focal_taxa,
    degree: int = 1,
    min_shared_samples_between_taxa: int = 1,
):
    """
    Return a filtered Ingredients object containing only samples within a k-degree
    taxonomic neighbourhood around one or more focal taxa.
    
    The neighbourhood is defined by a breadth-first search (BFS) over the bipartite
    sample–taxon graph induced by the binary presence matrix:
    
        taxa → samples → taxa → samples → ...
        
    Each BFS iteration expands from the current set of taxa to samples in which they
    occur, then from those samples to additional taxa. Expansion proceeds for
    `degree` taxa→sample steps.
    
    A taxon is only added during the samples→taxa expansion if it co-occurs with at
    least one already-reached taxon in a minimum number of shared samples, enforcing
    a co-occurrence strength constraint on taxon–taxon connectivity.
    
    Parameters
    ----------
    ingredients : Ingredients
        Object containing a binary presence matrix (CSR), sample identifiers, and
        taxon identifiers.
        
    focal_taxa : int | str | list[int | str]
        One or more focal taxa, specified by column index, full taxon string, or
        taxonomy token.
        
    degree : int, optional (≥1)
        Neighbourhood radius measured in taxa→sample expansions.
        * degree = 1 → samples containing the focal taxa
        * degree = 2 → samples connected via co-occurring taxa
        * degree = k → iteratively expand outward
        
    min_shared_samples_between_taxa : int, optional (≥1)
        Minimum number of samples in which two taxa must co-occur for a new taxon to
        be included during BFS expansion.
        
    Returns
    -------
    Ingredients
        A filtered Ingredients subset containing only samples reachable within the
        specified neighbourhood.
    """
    
    if degree < 1:
        raise ValueError("degree must be ≥ 1")
    if min_shared_samples_between_taxa < 1:
        raise ValueError("min_shared_samples_between_taxa must be ≥ 1")
        
    P = ingredients.presence_matrix  # CSR binary
    n_samples, n_taxa = P.shape
    
    focal_idx = _resolve_focal_taxa_indices(ingredients, focal_taxa)
    if focal_idx.size == 0:
        return ingredients.filtered_samples(np.zeros(n_samples, dtype=bool))
        
    visited_samples = np.zeros(n_samples, dtype=bool)
    visited_taxa = np.zeros(n_taxa, dtype=bool)
    
    frontier_taxa = np.zeros(n_taxa, dtype=bool)
    frontier_taxa[focal_idx] = True
    visited_taxa[focal_idx] = True
    
    P_T = P.T
    
    for step in range(degree):
        if not frontier_taxa.any():
            break
            
        # taxa → samples
        sample_hits = (P @ frontier_taxa.astype(int)) > 0
        new_samples = np.asarray(sample_hits).ravel() & ~visited_samples
        visited_samples |= new_samples
        
        if (not new_samples.any()) or (step == degree - 1):
            break
            
        # samples → taxa (candidates seen in the new samples)
        candidate_taxa = (P_T @ new_samples.astype(int)) > 0
        candidate_taxa = np.asarray(candidate_taxa).ravel()
        
        # enforce minimum shared samples between "the two" taxa
        if min_shared_samples_between_taxa > 1:
            frontier_cols = np.where(frontier_taxa)[0]
            # shared sample counts between each frontier taxon and every taxon
            shared = (P[:, frontier_cols].T @ P)  # (n_frontier x n_taxa)
            # require co-occurrence with at least one frontier taxon in >= K samples
            shared_max = np.asarray(shared.max(axis=0)).ravel()
            candidate_taxa &= (shared_max >= min_shared_samples_between_taxa)
            
        new_taxa = candidate_taxa & ~visited_taxa
        visited_taxa |= new_taxa
        frontier_taxa = new_taxa
        
    print(
        f"{int(visited_taxa.sum())} taxa detected within {degree} degrees "
        f"of focal taxa across {int(visited_samples.sum())} samples."
    )
    
    return ingredients.filtered_samples(visited_samples)