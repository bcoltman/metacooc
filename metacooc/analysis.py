#!/usr/bin/env python3
# metacooc/analytics/analysis.py

from __future__ import annotations
from typing import Optional, Iterable, Tuple, List, Set, Callable
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import chi2 as chi2_dist, fisher_exact as _fisher_exact
from scipy.sparse import csr_matrix, csc_matrix

from metacooc.pantry import *
from metacooc.utils import _RANK_PREFIXES, stream_edges

from metacooc.null_models import (
    null_matrices,
    parallel_null_reduce_vector,
    stat_fn_association_jaccard,
    stat_fn_cooccurrence_jaccard,
    _best_mp_start
)


import multiprocessing as mp
import os

_SMOOTH = 0.5

def _chi2_phi_from_counts(a, b, c, d):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    c = np.asarray(c, float)
    d = np.asarray(d, float)

    n    = a + b + c + d
    row1 = a + b
    row2 = c + d
    col1 = a + c
    col2 = b + d

    with np.errstate(divide="ignore", invalid="ignore"):
        exp_a = row1 * col1 / n
        exp_b = row1 * col2 / n
        exp_c = row2 * col1 / n
        exp_d = row2 * col2 / n

        chi2 = ((a - exp_a) ** 2) / exp_a \
             + ((b - exp_b) ** 2) / exp_b \
             + ((c - exp_c) ** 2) / exp_c \
             + ((d - exp_d) ** 2) / exp_d

        denom = np.sqrt((a + b) * (a + c) * (b + d) * (c + d))
        phi   = (a * d - b * c) / denom
        
        invalid = (
        (a < 0) | (b < 0) | (c < 0) | (d < 0) |               # impossible tables
        np.isclose(denom, 0) |                                # φ denominator zero
        ~np.isfinite(chi2) | ~np.isfinite(phi) |              # any statistic undefined
        (a + b == 0) | (c + d == 0) |                         # taxon only absent or only present
        (a + c == 0) | (b + d == 0)                           # no contrast in T vs ¬T
        )

    return chi2, phi, invalid


def table_metrics(a: np.ndarray,
                  b: np.ndarray,
                  c: np.ndarray,
                  d: np.ndarray,
                  compute_fisher: bool = True) -> pd.DataFrame:
    """
    Vectorised χ², p, φ, directional risk ratios, and optional Fisher's exact test
    from 2×2 tables:
    
               B present   B absent
    A present      a           b
    A absent       c           d
    """
    
    chi2, phi, invalid = _chi2_phi_from_counts(a,b,c,d)
    log_p = chi2_dist.logsf(chi2, df=1)  # natural log, very negative but finite
    p = np.exp(log_p)                    # will underflow to 0 for extreme chi2    
        
    # directional RRs (_SMOOTH for stability)
    P_B_given_A    = (a + _SMOOTH) / (a + b + _SMOOTH)
    P_B_given_notA = (c + _SMOOTH) / (c + d + _SMOOTH)
    RR_A_to_B      = P_B_given_A / P_B_given_notA
    with np.errstate(divide="ignore", invalid="ignore"):
        logRR_A_to_B   = np.log(RR_A_to_B)
    
    P_A_given_B    = (a + _SMOOTH) / (a + c + _SMOOTH)
    P_A_given_notB = (b + _SMOOTH) / (b + d + _SMOOTH)
    RR_B_to_A      = P_A_given_B / P_A_given_notB
    with np.errstate(divide="ignore", invalid="ignore"):
        logRR_B_to_A   = np.log(RR_B_to_A)
    
    # --- Fisher's exact test (optional, loop but only on valid tables) ---
    n_tables = a.shape[0]
    fisher_odds = np.full(n_tables, np.nan, dtype=float)
    fisher_p = np.full(n_tables, np.nan, dtype=float)
    log_fisher_p = np.full(n_tables, np.nan, dtype=float)
    
    if compute_fisher:
        valid_idx = np.where(~invalid)[0]
        for i in valid_idx:
            # 2×2 table for this row
            odds, pval = _fisher_exact(
                [[a[i], b[i]],
                 [c[i], d[i]]],
                alternative="two-sided"
            )
            fisher_odds[i] = odds
            fisher_p[i] = pval
            
        with np.errstate(divide="ignore", invalid="ignore"):
            log_fisher_p = np.where(fisher_p > 0, np.log(fisher_p), -np.inf)
            
    return pd.DataFrame({
        "chi2": chi2,
        "p": p,
        "log_p": log_p,
        "phi": phi,
        "RR_A_to_B": RR_A_to_B,
        "RR_B_to_A": RR_B_to_A,
        "logRR_A_to_B": logRR_A_to_B,
        "logRR_B_to_A": logRR_B_to_A,
        "fisher_odds": fisher_odds,
        "fisher_p": fisher_p,
        "log_fisher_p": log_fisher_p,
        "invalid_table": invalid,
    })

def bh_qvalues_from_logp(log_p, m_total: Optional[int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Benjamini–Hochberg FDR q-values (monotone), computed from log p-values.

    Parameters
    ----------
    log_p : array-like
        Natural-log p-values (e.g. from scipy.stats.chi2.logsf),
        for the *tested* hypotheses.
    m_total : int, optional
        Total number of hypotheses you conceptually correct for.
        If None (default), uses len(log_p).
        If > len(log_p), behaves as if the missing (unseen) hypotheses had p=1.

    Returns
    -------
    q : np.ndarray
        BH-adjusted q-values on the original scale.
    log_q : np.ndarray
        Corresponding natural-log q-values.
    """
    log_p = np.asarray(log_p, dtype=float)
    m = log_p.size
    if m == 0:
        return np.array([], float), np.array([], float)
        
    n = int(m_total) if m_total is not None else m
    
    finite_mask = np.isfinite(log_p)
    if not finite_mask.any():
        # everything is NaN / inf: return all-NaN q's
        return np.full(m, np.nan), np.full(m, np.nan)
    
    order = np.argsort(log_p[finite_mask])
    ranks = np.arange(1, finite_mask.sum() + 1, dtype=float)
    log_q_finite = log_p[finite_mask][order] + np.log(n) - np.log(ranks)
    log_q_finite = np.minimum.accumulate(log_q_finite[::-1])[::-1]

    log_q_final = np.full(m, np.nan, dtype=float)
    idx = np.where(finite_mask)[0]
    log_q_final[idx[order]] = log_q_finite
    q_final = np.exp(log_q_final)
    # q_final[q_final > 1.0] = 1.0
    
    return q_final, log_q_final


def select_taxa_universe(
    ing: "Ingredients",
    rank: Optional[str] = None,
    terminal_only: bool = True,   # True -> terminal rank equals `rank`; False -> lineage contains `rank`
    require_present: bool = True  # keep only taxa with any presence in current sample set
) -> List[str]:
    """
    Build the taxa universe to analyse, leveraging Ingredients' caches and
    current presence.
    
    Base set
    --------
    The starting set of taxa is:
    
      - If ``require_present=True`` (default): all taxa with
        ``total_counts > 0`` in the current ``ing`` object (after any prior
        sample filtering).
      - If ``require_present=False``: all ``ing.taxa`` in their stored order.
      
    Rank-based restriction
    ----------------------
    If ``rank`` is provided (e.g. ``"species"``), the base set is further
    restricted using the cached rank information in ``Ingredients``:
    
      - ``terminal_only=True`` (default):
          Keep taxa whose *terminal* rank equals ``rank``. This uses the cached
          ``_terminal_rank_prefixes``.
      - ``terminal_only=False``:
          Keep taxa whose lineage string contains the rank prefix anywhere.
          
    If ``rank`` is ``None``, the base set is returned unchanged.
    
    Parameters
    ----------
    ing : Ingredients
        Ingredients object providing ``taxa``, ``total_counts`` and cached
        rank lookups.
    rank : str, optional
        Taxonomic rank name (case-insensitive), e.g. ``"species"``. Must be a
        key in ``_RANK_PREFIXES`` if provided.
    terminal_only : bool, default True
        Whether to require the terminal rank to match ``rank`` (True) or to
        accept taxa where the rank appears anywhere in the lineage (False).
    require_present : bool, default True
        If True, restrict to taxa with ``total_counts > 0`` in the current
        sample set. If False, include all taxa.
        
    Returns
    -------
    List[str]
        Ordered list of taxon names forming the analysis universe.
    """
    # Base availability
    if require_present:
        present_set = {t for t, cnt in zip(ing.taxa, ing.total_counts) if cnt > 0}
    else:
        present_set = set(ing.taxa)
        
    # Start from base (present-only or all)
    base = [t for t in ing.taxa if t in present_set]
    
    # 2) Rank-based selection
    if rank is None:
        return base
        
    rk = rank.lower().strip()
    if rk not in _RANK_PREFIXES:
        raise ValueError(
            f"Unknown rank '{rank}'. Expected one of: {', '.join(_RANK_PREFIXES.keys())}"
        )
    pref = _RANK_PREFIXES[rk]
    
    # Ensure caches exist
    ing._ensure_taxa_lookups()
    
    if terminal_only:
        # Build a taxon -> terminal prefix map once
        term_map = {t: tp for t, tp in zip(ing.taxa, ing._terminal_rank_prefixes)}
        return [t for t in base if term_map.get(t) == pref]
    else:
        # Looser: lineage contains this rank anywhere
        return [t for t in base if pref in t]



# ========== Workflow A: filtered vs null (single-taxon enrichment) ==========

def association_obj(
    null_ingredients: "Ingredients",
    filtered_ingredients: "Ingredients",
    threshold: float = 0.0,
    null_model: str = "FE",               # "FE" or "FF" 
    nm_n_reps: int = 1000,
    nm_random_state: int | None = None,
    compute_fisher: bool = False,
) -> pd.DataFrame:
    """
    Entry point for association / term-enrichment analysis.

    Parameters
    ----------
    null_ingredients : Ingredients
        Background cohort (depends on null_scope, biome/local, etc).
    filtered_ingredients : Ingredients
        Term cohort T: samples matching the search term (soil, Nitrospira_D+, ...).
        Must have samples ⊆ null_ingredients.samples.
    threshold : float
        Optional filter on p_T_given_X (taxon-centric specificity).
    null_model : {"FE", "FF"}
        "FE" -> analytic tests based on 2×2 tables (χ², Fisher, etc.) only.
        "FF" -> same per-taxon stats + (optionally) SIM9-based community metrics.
    community_structure : bool
        If True and null_model == "FF", run SIM9 on the null presence matrix
        and attach C-score, NODF, mean Jaccard (global + subset) in attrs.
    nm_n_reps : int
        Number of SIM9 null matrices if community_structure is True.
    nm_random_state : int or None
        Random seed for SIM9.

    Returns
    -------
    out : pd.DataFrame
        One row per taxon with enrichment metrics and p-values.
        If community_structure is True and null_model == "FF",
        out.attrs["sim9_matrix_metrics"] contains matrix-level metrics.
    """
    
    out = _association_core(
        null_ingredients=null_ingredients,
        filtered_ingredients=filtered_ingredients,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_random_state=nm_random_state,
    )

    if threshold is not None:
        out = out[out["p_T_given_X"] > threshold].copy()

    return out


def association(
    null_ingredients: "Ingredients",
    filtered_ingredients: "Ingredients",
    output_dir: str,
    tag: str | None = None,
    threshold: float = 0.0,
    null_model: str = "FE",
    nm_n_reps: int = 1000,
    nm_random_state: int | None = None,
    compute_fisher: bool = False,
) -> pd.DataFrame:
    """
    Entry point for association / term-enrichment analysis.

    Parameters
    ----------
    null_ingredients : Ingredients
        Background cohort (depends on null_scope, biome/local, etc).
    filtered_ingredients : Ingredients
        Term cohort T: samples matching the search term (soil, Nitrospira_D+, ...).
        Must have samples ⊆ null_ingredients.samples.
    threshold : float
        Optional filter on p_T_given_X (taxon-centric specificity).
    null_model : {"FE", "FF"}
        "FE" -> analytic tests based on 2×2 tables (χ², Fisher, etc.) only.
        "FF" -> same per-taxon stats + (optionally) SIM9-based community metrics.
    community_structure : bool
        If True and null_model == "FF", run SIM9 on the null presence matrix
        and attach C-score, NODF, mean Jaccard (global + subset) in attrs.
    nm_n_reps : int
        Number of SIM9 null matrices if community_structure is True.
    nm_random_state : int or None
        Random seed for SIM9.

    Returns
    -------
    out : pd.DataFrame
        One row per taxon with enrichment metrics and p-values.
        If community_structure is True and null_model == "FF",
        out.attrs["sim9_matrix_metrics"] contains matrix-level metrics.
    """
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    # Normalise / load Ingredients objects
    null = load_ingredients(data_dir=None, custom_ingredients=null_ingredients)
    filtered = load_ingredients(data_dir=None, custom_ingredients=filtered_ingredients)

    association_df = association_obj(
        null_ingredients=null,
        filtered_ingredients=filtered,
        threshold=threshold,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_random_state=nm_random_state,
        compute_fisher=compute_fisher
    )
    
    # Write node table if available
    if association_df is not None:
        output_path = os.path.join(output_dir, f"{tag}association.tsv")
        association_df.to_csv(output_path, sep="\t", index=False)
        print(f"Association analysis saved to {output_path}")


def presence_submatrix_by_taxa(ingredients: Ingredients, taxa_subset: List[str]) -> sp.csr_matrix:
    row_map = {t: i for i, t in enumerate(ingredients.taxa)}
    rows = [row_map[t] for t in taxa_subset]
    return ingredients.presence_matrix[rows,:].tocsr()

def symmetric_counts_from_matrix(ingredients: Ingredients,
                                 taxa_universe: List[str]) -> Tuple[np.ndarray, sp.csr_matrix]:
    X = presence_submatrix_by_taxa(ingredients, taxa_universe)
    totals = np.asarray(X.sum(axis=1)).ravel()
    co_mat = (X @ X.T).tocsr()
    return totals, co_mat

def conditional_probabilities(co_mat: sp.csr_matrix,
                              totals: np.ndarray) -> Tuple[sp.csr_matrix, sp.csr_matrix]:
    totals = totals.astype(float, copy=False)
    inv = np.zeros_like(totals, dtype=float)
    nz = totals > 0
    inv[nz] = 1.0 / totals[nz]
    D_inv = sp.diags(inv)
    return (D_inv @ co_mat).tocsr(), (co_mat @ D_inv).tocsr()


def build_edge_chunk(iAs: np.ndarray,
                     iBs: np.ndarray,
                     totals: np.ndarray,
                     co_csr: sp.csr_matrix,
                     P_BA_csr: sp.csr_matrix,
                     P_AB_csr: sp.csr_matrix,
                     taxa_universe: List[str],
                     N_total: int) -> pd.DataFrame:
    """
    Build a DataFrame of edge-level co-occurrence metrics for a chunk of (A,B) pairs.
    
    For each pair (A,B):
    
        A_total_count  = #(A present)
        B_total_count  = #(B present)
        A_B_intersection_count = #(A and B co-present)
        
        2×2 table (rows = A present/absent, cols = B present/absent):
        
                     B present      B absent
        A present        a              b
        A absent         c              d
        
        where:
            a = A_B_intersection_count
            b = A_total_count - a
            c = B_total_count - a
            d = N_total - (a + b + c)
            
    We compute:
        - chi2, p, phi, RR_A_to_B, RR_B_to_A, logRR_A_to_B, logRR_B_to_A (from table_metrics)
        - p_A, p_B, p_A_and_B, lift
        - P_B_given_A, P_A_given_B
        - jaccard = a / (A_total + B_total - a)
    """
    A_totals = totals[iAs]
    B_totals = totals[iBs]
    inter = co_csr[iAs, iBs].A1  # co-occurrence counts
    
    # 2×2 table for each (A,B)
    a = inter.astype(float, copy=False)
    b = (A_totals - inter).astype(float, copy=False)
    c = (B_totals - inter).astype(float, copy=False)
    d = (N_total - (a + b + c)).astype(float, copy=False)
    
    mets = table_metrics(a, b, c, d)
    
    Pba = P_BA_csr[iAs, iBs].A1
    Pab = P_AB_csr[iAs, iBs].A1
    
    # Marginal and joint prevalences + symmetric lift + Jaccard
    with np.errstate(divide="ignore", invalid="ignore"):
        p_A = A_totals / float(N_total)
        p_B = B_totals / float(N_total)
        p_A_and_B = inter / float(N_total)
        
        denom_lift = p_A * p_B
        lift = np.divide(
            p_A_and_B,
            denom_lift,
            out=np.zeros_like(p_A_and_B, dtype=float),
            where=denom_lift != 0
        )
        
        union = A_totals + B_totals - inter
        jaccard = np.divide(
            inter,
            union,
            out=np.zeros_like(inter, dtype=float),
            where=union > 0
        )
    
    taxa_arr = np.asarray(taxa_universe, dtype=object)
    df = pd.DataFrame({
        "A_taxon": taxa_arr[iAs],
        "B_taxon": taxa_arr[iBs],
        # raw counts
        "A_B_intersection_count": inter,
        "A_total_count": A_totals,
        "B_total_count": B_totals,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "N": N_total,
        # marginal / joint prevalences
        "p_A": p_A,
        "p_B": p_B,
        "p_A_and_B": p_A_and_B,
        "lift": lift,
        # directional conditionals
        "P_B_given_A": Pba,
        "P_A_given_B": Pab,
        # symmetric similarity
        "jaccard": jaccard,
    })
    
    df = df.join(mets)
    df = df.loc[~df["invalid_table"]].copy()
    df.drop(columns="invalid_table", inplace=True)
    return df


def _cooccur_core(
    ing: Ingredients,
    taxa_universe: List[str],
    threshold: float = 0.1,
    m_total: Optional[int] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, sp.csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    """
    Pairwise co-occurrence edges (within 'ing') and node summary.
    
    Returns
    -------
    edges_df : DataFrame or None
    nodes_df : DataFrame
    X_sub    : csr_matrix (len(taxa_universe) x samples)
    totals   : 1D array of taxon totals
    iA_all   : 1D array of A indices (into taxa_universe) for each edge
    iB_all   : 1D array of B indices (into taxa_universe) for each edge
    """
    # Build the restricted presence matrix once
    X_sub = presence_submatrix_by_taxa(ing, taxa_universe)  # sparse, binary
    totals = np.asarray(X_sub.sum(axis=1)).ravel()
    co_mat = (X_sub @ X_sub.T).tocsr()
    P_BA, P_AB = conditional_probabilities(co_mat, totals)
    
    N_total = len(ing.samples)
    co_csr   = co_mat.tocsr(copy=False)
    P_BA_csr = P_BA.tocsr(copy=False)
    P_AB_csr = P_AB.tocsr(copy=False)
    
    chunks: List[pd.DataFrame] = []
    iA_all_list: List[np.ndarray] = []
    iB_all_list: List[np.ndarray] = []
    
    for iA, iB in stream_edges(P_BA_csr, threshold):
        iA_all_list.append(iA)
        iB_all_list.append(iB)
        chunks.append(
            build_edge_chunk(
                iA, iB,
                totals=totals,
                co_csr=co_csr,
                P_BA_csr=P_BA_csr,
                P_AB_csr=P_AB_csr,
                taxa_universe=taxa_universe,
                N_total=N_total,
            )
        )
    
    if chunks:
        edges_df = pd.concat(chunks, ignore_index=True)
        
        q_bh, log_q_bh = bh_qvalues_from_logp(edges_df["log_p"].values, m_total=m_total)
        edges_df["q_bh"] = q_bh
        edges_df["log_q_bh"] = log_q_bh
        
        iA_all = np.concatenate(iA_all_list)
        iB_all = np.concatenate(iB_all_list)
    else:
        edges_df = None
        iA_all = np.array([], dtype=int)
        iB_all = np.array([], dtype=int)
    
    deg_fwd = np.asarray((P_BA_csr > threshold).sum(axis=1)).ravel()
    nodes_df = pd.DataFrame({
        "taxon": taxa_universe,
        "total_count": totals.astype(int, copy=False),
        f"degree_PBA_gt_{threshold}": deg_fwd
    })
    
    return edges_df, nodes_df, X_sub, totals, iA_all, iB_all


def should_run_cooccurrence(
    n_taxa: int,
    large: bool,
    max_pairs: int = 100_000,
) -> Tuple[bool, int]:
    """
    Decide whether to run co-occurrence given user intent + scale.
    Returns (run, estimated_pairs).
      - Else compute nC2; require (large==True) OR estimated_pairs <= max_pairs.
    """
    pairs = (n_taxa * (n_taxa - 1)) // 2
    if large:
        return True, pairs
    return (pairs <= max_pairs), pairs



def cooccurrence(
    null_ingredients,
    filtered_ingredients,
    output_dir: str,
    tag: str | None = None,
    filter_rank: Optional[str] = None,   # e.g. "species" ties to aggregated datasets
    large: bool = False,                 # --large to allow huge runs
    max_pairs: int = 100_000,            # soft cap unless --large
    threshold: float = 0.1,
    null_model: str = "FE",
    nm_n_reps: int = 10,
    nm_random_state: int = 42,
):
    """
    Run taxon–taxon co-occurrence analysis and write edge/node tables to disk.
    
    This is a convenience wrapper around :func:`taxon_cooccurrence_obj` that:
    
      1. Loads / normalises the provided Ingredients objects.
      2. Builds the taxa universe using :func:`select_taxa_universe` on the
         *filtered* Ingredients (optionally restricted by a taxonomic rank).
      3. Runs co-occurrence analysis for that taxa universe.
      4. Writes the resulting edge- and node-level tables to TSV files.
      
    Output files
    ------------
    Files are written into ``output_dir`` (created if needed) with optional
    ``tag`` prefix:
    
      - ``{tag}taxon_edges.tsv`` – edge-level co-occurrence metrics
        (one row per taxon–taxon pair above ``threshold``), if any edges exist.
      - ``{tag}taxon_nodes.tsv`` – node-level summary (one row per taxon),
        if the analysis is not skipped.
            
    Parameters
    ----------
    null_ingredients :
        Background cohort used as the null community for co-occurrence. This is
        passed to :func:`load_ingredients` via ``custom_ingredients``, so it can
        be either an ``Ingredients`` instance or whatever ``load_ingredients``
        knows how to handle.
    filtered_ingredients :
        Cohort used to define the taxon universe (after any filtering).
        Also passed through :func:`load_ingredients` in the same way.
    output_dir : str
        Directory where output TSV files will be written. Created if it does not
        exist.
    tag : str, optional
        Optional prefix for output filenames (e.g. ``"soil_"`` → ``"soil_taxon_edges.tsv"``).
        If ``None``, no prefix is used.
    filter_rank : str, optional
        If given, passed as ``rank`` to :func:`select_taxa_universe` to restrict
        the analysis to taxa whose terminal rank (or lineage) matches this rank.
    large : bool, default False
        If ``False``, co-occurrence analysis is skipped when the estimated number
        of taxon pairs exceeds ``max_pairs``. If ``True``, this safety check is
        overridden and all pairs are analysed.
    max_pairs : int, default 100_000
        Soft cap on the number of taxon pairs unless ``large=True``.
    threshold : float, default 0.1
        Minimum conditional probability (``P(B|A)``) required for an edge to be
        included in the edge table.
    null_model : {"FE", "FF"}, default "FE"
        Co-occurrence null model; passed through to :func:`taxon_cooccurrence_obj`.
    nm_n_reps : int, default 10
        Number of SIM9 null matrices for the "FF" null model.
        
    Returns
    -------
    edges_df : pandas.DataFrame or None
        Edge-level co-occurrence metrics. ``None`` if co-occurrence was skipped
        (too many pairs and ``large=False``) or if no edges exceeded the
        ``threshold``.
    nodes_df : pandas.DataFrame or None
        Node-level summary. ``None`` only if co-occurrence was skipped entirely.
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        
    # Normalise / load Ingredients objects
    null = load_ingredients(data_dir=None, custom_ingredients=null_ingredients)
    filtered = load_ingredients(data_dir=None, custom_ingredients=filtered_ingredients)
    
    # Define taxa universe on the *filtered* cohort
    taxa_universe = select_taxa_universe(filtered, rank=filter_rank)
    
    # Run co-occurrence
    edges_df, nodes_df = cooccurrence_obj(
        null,
        taxa_universe,
        large=large,
        max_pairs=max_pairs,
        threshold=threshold,
        null_model= null_model,
        nm_n_reps=nm_n_reps,
        nm_random_state=nm_random_state
    )
    
    
    # Write node table if available
    if nodes_df is not None:
        nodes_path = os.path.join(output_dir, f"{tag}taxon_nodes.tsv")
        nodes_df.to_csv(nodes_path, sep="\t", index=False)
        print(f"Taxon nodes analysis saved to {nodes_path}")
        
    # Write edge table only if we actually have edges
    if edges_df is not None:
        edges_path = os.path.join(output_dir, f"{tag}taxon_edges.tsv")
        edges_df.to_csv(edges_path, sep="\t", index=False)
        print(f"Taxon edges analysis saved to {edges_path}")
    else:
        print("No co-occurrence edges above threshold; no edge file written.")





def _association_core(
    null_ingredients: "Ingredients",
    filtered_ingredients: "Ingredients",
    null_model: str = "FE",
    nm_n_reps: int = 1000,
    nm_random_state: int | None = None,
    compute_fisher: bool = False,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_sort_indices: bool = False,
    nm_burn_in_steps: int | None = None,   # FF only; forwarded to null_matrices via helper
    nm_steps_per_rep: int | None = None,   # FF only; forwarded to null_matrices via helper
) -> pd.DataFrame:
    """
    Core taxon-term enrichment (no community-structure metrics).

    Observed: 2x2 association metrics + jaccard(taxon, term) = a / (b + N_T)

    Null (Jaccard only):
      - Uses null_matrices(X_full, model=..., n_reps=...) inside a parallel reducer.
      - The reducer parallelises across workers by giving each worker its own slice
        of replicates and independent RNG seed.
      - No null matrices are shipped back to the parent process; workers return
        only partial accumulators.

    Cross-platform:
      - Works on Linux/macOS/Windows; on spawn platforms, ensure your program entrypoint
        is guarded by `if __name__ == "__main__":`.
    """
    # --- sanity: sample set relationship (filtered ⊆ null) ---
    null_samples = set(null_ingredients.samples)
    filt_samples = set(filtered_ingredients.samples)
    if not filt_samples.issubset(null_samples):
        raise ValueError("filtered_ingredients.samples must be a subset of null_ingredients.samples")

    # --- align taxa ---
    null_taxa = np.array(null_ingredients.taxa, dtype=object)
    filt_taxa = np.array(filtered_ingredients.taxa, dtype=object)
    taxa, filt_idx, null_idx = np.intersect1d(filt_taxa, null_taxa, return_indices=True)
    if taxa.size == 0:
        raise ValueError("No taxa intersect the null and filtered ingredients.")

    # --- counts and cohort sizes ---
    N_T = float(len(filtered_ingredients.samples))
    N_null = float(len(null_ingredients.samples))
    if N_T <= 0 or N_null <= 0 or N_null < N_T:
        raise ValueError(f"Invalid cohort sizes: N_T={N_T}, N_null={N_null}")
    if N_null == N_T:
        raise ValueError("Null and term cohorts are identical (N_null == N_T): no non-term samples available.")

    # X present in T / non-T (counts over samples in each cohort)
    a_raw = filtered_ingredients.total_counts[filt_idx].astype(float, copy=False)  # X in term
    ref_counts = null_ingredients.total_counts[null_idx].astype(float, copy=False)
    b_raw = ref_counts - a_raw                                                     # X in non-term

    N_notT = N_null - N_T
    c_raw = N_T - a_raw
    d_raw = N_notT - b_raw

    with np.errstate(divide="ignore", invalid="ignore"):
        jaccard_obs_all = np.divide(
            a_raw,
            (b_raw + N_T),
            out=np.zeros_like(a_raw, dtype=float),
            where=(b_raw + N_T) > 0,
        )

    mets = table_metrics(a_raw, b_raw, c_raw, d_raw, compute_fisher=compute_fisher)

    # Smoothed counts for probability-based metrics
    a = a_raw + _SMOOTH
    b = b_raw + _SMOOTH
    c = c_raw + _SMOOTH
    d = d_raw + _SMOOTH

    with np.errstate(divide="ignore", invalid="ignore"):
        p_X_given_T = a / (a + c)
        p_X_given_notT = b / (b + d)
        RR_T = np.divide(p_X_given_T, p_X_given_notT, out=np.zeros_like(p_X_given_T), where=p_X_given_notT != 0)
        log2_RR_T = np.where(RR_T > 0, np.log2(RR_T), -np.inf)
        delta_p_T = p_X_given_T - p_X_given_notT

        p_T_given_X = a / (a + b)
        p_T_given_notX = c / (c + d)
        RR_X = np.divide(p_T_given_X, p_T_given_notX, out=np.zeros_like(p_T_given_X), where=p_T_given_notX != 0)
        log2_RR_X = np.where(RR_X > 0, np.log2(RR_X), -np.inf)

        P_X = (a + b) / N_null
        P_T = (a + c) / N_null
        P_XT = a / N_null
        lift = np.divide(P_XT, P_X * P_T, out=np.zeros_like(P_XT), where=(P_X * P_T) != 0)

    out = pd.DataFrame(
        {
            "taxon": taxa,
            "a": a_raw, "b": b_raw, "c": c_raw, "d": d_raw,
            "N_T": N_T, "N_notT": N_notT, "N_null": N_null,
            "p_X_given_T": p_X_given_T,
            "p_X_given_notT": p_X_given_notT,
            "RR_T": RR_T,
            "log2_RR_T": log2_RR_T,
            "delta_p_T": delta_p_T,
            "p_T_given_X": p_T_given_X,
            "p_T_given_notX": p_T_given_notX,
            "RR_X": RR_X,
            "log2_RR_X": log2_RR_X,
            "lift": lift,
            "jaccard": jaccard_obs_all,
        }
    )

    out = out.join(
        mets[
            [
                "chi2", "p", "log_p", "phi",
                "RR_A_to_B", "RR_B_to_A", "logRR_A_to_B", "logRR_B_to_A",
                "fisher_odds", "fisher_p", "log_fisher_p",
                "invalid_table",
            ]
        ]
    )

    q_bh, log_q_bh = bh_qvalues_from_logp(out["log_p"].values)
    out["q_bh"] = q_bh
    out["log_q_bh"] = log_q_bh

    valid_mask = ~out["invalid_table"].values
    out = out.loc[valid_mask].copy()
    out.drop(columns="invalid_table", inplace=True)

    null_idx_valid = null_idx[valid_mask]
    out.attrs["null_idx_valid"] = np.asarray(null_idx_valid, dtype=int)
    out.attrs["taxa_valid"] = out["taxon"].to_numpy(dtype=object, copy=False)

    # ---- null (parallel) ----
    n_reps = int(nm_n_reps) if nm_n_reps is not None else 0
    if n_reps <= 0 or out.shape[0] == 0:
        return out

    suffix = str(null_model).upper()
    
    # FE: association determined analyticlaly - no need for prbababilisticanalytic-only
    if suffix == "FE":
        print("FE: association determined analytically - no need for shuffling null and probabilistic approach")
        return out

    # Prepare full matrix once
    X_full = null_ingredients.presence_matrix.tocsr()
    X_full.eliminate_zeros()
    X_full.sum_duplicates()
    X_full.sort_indices()

    n_rows, n_cols = X_full.shape

    sample_index = {s: i for i, s in enumerate(null_ingredients.samples)}
    term_cols = np.array([sample_index[s] for s in filtered_ingredients.samples], dtype=np.int64)

    mask_cols = np.zeros(n_cols, dtype=bool)
    mask_cols[term_cols] = True
    nonterm_cols = np.where(~mask_cols)[0].astype(np.int64, copy=False)

    subset_idx = np.asarray(null_idx_valid, dtype=np.int64)
    # subset_idx.sort()

    obs_jacc = out["jaccard"].to_numpy(dtype=float, copy=False)

    mp_start = _best_mp_start() if nm_mp_start is None else str(nm_mp_start)

    j_res = parallel_null_reduce_vector(
        X=X_full,
        model=suffix,
        n_reps=n_reps,
        obs=obs_jacc,
        stat_fn=stat_fn_association_jaccard,
        random_state=nm_random_state,
        n_workers=nm_n_workers,
        mp_start=mp_start,
        term_cols=term_cols,
        nonterm_cols=nonterm_cols,
        subset_idx=subset_idx,
        n_rows=n_rows,
        N_T=N_T,
    )

    out[f"jaccard_null_mean_{suffix}"] = j_res["mean"]
    out[f"jaccard_null_sd_{suffix}"] = j_res["sd"]
    out[f"jaccard_ses_{suffix}"] = j_res["ses"]
    out[f"jaccard_p_{suffix}"] = j_res["p_emp"]
    
    out[f"n_ok_{suffix}"]      = int(j_res["n_ok"])
    out[f"n_err_{suffix}"]     = int(j_res["n_err"])
    out[f"n_done_{suffix}"]    = int(j_res["n_done"])
    out[f"n_requested_{suffix}"]    = int(j_res["n_target"])

    return out


def cooccurrence_obj(
    null_ingredients: "Ingredients",
    taxa_universe: List[str],
    large: bool = False,
    max_pairs: int = 100_000,
    threshold: float = 0.1,
    null_model: str = "FE",
    nm_n_reps: int = 10,
    nm_random_state: int | None = None,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_sort_indices: bool = False,
    nm_burn_in_steps: int | None = None,   # FF only
    nm_steps_per_rep: int | None = None,   # FF only
) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame]]:
    """
    Pairwise co-occurrence of taxa.
    
    Observed:
      - computed via _cooccur_core (χ², Fisher, φ, RRs, observed Jaccard, etc.)
      
    Null (Jaccard only):
      - parallel reduction over null replicates from null_matrices on FULL presence matrix
      - per-edge Jaccard under null
      - attaches mean/sd/SES/p with suffix: *_FF, *_FE, etc.
      - BH q-values from empirical p (with m_total=est_pairs)
    """
    run_co, est_pairs = should_run_cooccurrence(
        n_taxa=len(taxa_universe),
        large=large,
        max_pairs=max_pairs,
    )
    
    if not run_co:
        if not large and est_pairs > max_pairs:
            print(
                f"Estimated pairs ({est_pairs:,}) exceed max_pairs={max_pairs:,}. "
                "Skipping complete co-occurrence analysis.\n\n"
                "To perform complete co-occurrence analysis, consider the following options:\n"
                "  • Increase the --max_pairs parameter, or\n"
                "  • Use --large to override (note: memory usage may exceed 100 GB, and large output files will be generated).\n\n"
                "To reduce memory usage, you can:\n"
                "  • Increase the --threshold value [0, 1],\n"
                "  • Lower the --filter_rank (e.g. to 'species'),\n"
                "  • Adjust the input ingredients by:\n"
                "    - Raising the --min_taxa_count, or\n"
                "    - Raising the --min_sample_count."
            )
        return None, None
        
    edges_df, nodes_df, X_sub, totals, iA_all, iB_all = _cooccur_core(
        null_ingredients,
        taxa_universe,
        threshold=threshold,
        m_total=est_pairs,
    )
    
    if edges_df is None or len(edges_df) == 0:
        return edges_df, nodes_df
        
    n_reps = int(nm_n_reps) if nm_n_reps is not None else 0
    if n_reps <= 0:
        return edges_df, nodes_df
        
    suffix = str(null_model).upper()
    
    # FE: association determined analyticlaly - no need for prbababilisticanalytic-only
    if suffix == "FE":
        print("FE: association determined analytically - no need for shuffling null and probabilistic approach")
        return edges_df, nodes_df
        
    # FULL matrix (CSR normalised once)
    X_full = null_ingredients.presence_matrix.tocsr()
    X_full.eliminate_zeros()
    X_full.sum_duplicates()
    X_full.sort_indices()
    
    # map taxa_universe -> rows in X_full
    tax_map = {t: i for i, t in enumerate(null_ingredients.taxa)}
    subset_idx = np.array([tax_map[t] for t in taxa_universe], dtype=np.int64)
    
    obs_jacc = edges_df["jaccard"].to_numpy(dtype=float, copy=False)
    
    mp_start = _best_mp_start() if nm_mp_start is None else str(nm_mp_start)
    
    j_res = parallel_null_reduce_vector(
        X=X_full,
        model=suffix,
        n_reps=n_reps,
        obs=obs_jacc,
        stat_fn=stat_fn_cooccurrence_jaccard,
        random_state=nm_random_state,
        n_workers=nm_n_workers,
        mp_start=mp_start,
        subset_idx=subset_idx,
        iA=iA_all,
        iB=iB_all,
    )
    
    edges_df[f"jaccard_null_mean_{suffix}"] = j_res["mean"]
    edges_df[f"jaccard_null_sd_{suffix}"] = j_res["sd"]
    edges_df[f"jaccard_ses_{suffix}"] = j_res["ses"]
    edges_df[f"jaccard_p_{suffix}"] = j_res["p_emp"]
    
    log_p = np.log(edges_df[f"jaccard_p_{suffix}"].to_numpy(dtype=float, copy=False))
    q, log_q = bh_qvalues_from_logp(log_p, m_total=est_pairs)
    edges_df[f"jaccard_q_{suffix}"] = q
    edges_df[f"log_jaccard_q_{suffix}"] = log_q
    
    edges_df[f"n_ok_{suffix}"]      = int(j_res["n_ok"])
    edges_df[f"n_err_{suffix}"]     = int(j_res["n_err"])
    edges_df[f"n_done_{suffix}"]    = int(j_res["n_done"])
    edges_df[f"n_requested_{suffix}"]    = int(j_res["n_target"])
    
    return edges_df, nodes_df