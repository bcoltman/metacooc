import numpy as np
import scipy.sparse as sp
import pandas as pd

from metacooc.utils import stream_csr_upper_threshold
from metacooc.null_models import (
    # null_matrices,
    parallel_null_reduce_vector,
    stat_fn_structure_metrics,
    _best_mp_start
)

def _nodf_sum_from_overlap(
    O: sp.csr_matrix,
    deg: np.ndarray,
    chunk_rows: int = 50_000,
) -> tuple[float, int]:
    """
    Sum contributions over i<j where deg[i] > deg[j] > 0:
        (O[i,j]/deg[j]) * 100
    """
    deg = deg.astype(np.float64, copy=False)

    total = 0.0
    pairs = 0

    for i, j, ov in stream_csr_upper_threshold(O, threshold=0.0, chunk_rows=chunk_rows):
        ov = ov.astype(np.float64, copy=False)
        dj = deg[j]
        valid = (deg[i] > dj) & (dj > 0)
        if np.any(valid):
            total += float(np.sum((ov[valid] / dj[valid]) * 100.0))
            pairs += int(np.sum(valid))

    return total, pairs


def compute_nodf_streamed(
    X: sp.spmatrix,
    chunk_rows: int = 50_000,
) -> float:
    """
    Exact NODF with streamed *processing* (no COO conversion).
    Still requires forming overlap matrices (X@X.T and X.T@X).
    Returns 0..100 scale.
    """

    n_rows, n_cols = X.shape
    # if n_rows < 2 and n_cols < 2:
    if n_rows < 2 or n_cols < 2:
        return np.nan
        
    # row part
    row_deg = np.diff(X.indptr)
    Or = (X @ X.T).tocsr(copy=False)
    nodf_rows_sum, row_pairs = _nodf_sum_from_overlap(Or, row_deg, chunk_rows=chunk_rows)
    
    # col part
    Xc = X.tocsc(copy=False)
    col_deg = np.diff(Xc.indptr)
    Oc = (Xc.T @ Xc).tocsr(copy=False)
    nodf_cols_sum, col_pairs = _nodf_sum_from_overlap(Oc, col_deg, chunk_rows=chunk_rows)
    
    total_pairs = row_pairs + col_pairs
    if total_pairs == 0:
        return np.nan
        
    return float((nodf_rows_sum + nodf_cols_sum) / total_pairs)


def compute_c_score(X: sp.spmatrix, chunk_rows: int = 50_000) -> float:
    """
    C-score across taxa for X = taxa × samples.
    
    For taxa i,j:
        C_ij = (R_i - S_ij) * (R_j - S_ij)
    where R_i is row degree (#samples occupied), and S_ij is shared samples.
    
    Returns mean C-score across all i<j taxa pairs, or np.nan if undefined.
    """
    
    R = np.diff(X.indptr).astype(np.float64, copy=False)
    n_taxa = R.size
    if n_taxa < 2:
        return np.nan
        
    num_pairs = n_taxa * (n_taxa - 1) / 2.0
    sum_R = R.sum()
    sum_R2 = np.square(R).sum()
    total_RiRj = (sum_R * sum_R - sum_R2) / 2.0
    
    # Overlap matrix among taxa
    S = (X @ X.T).tocsr(copy=False)
    
    # We need sum_{i<j} S_ij*(R_i+R_j) and sum_{i<j} S_ij^2 over nonzero S_ij only
    Sij_Rsum = 0.0
    Sij_sqsum = 0.0
    
    for i, j, sij in stream_csr_upper_threshold(S, threshold=0.0, chunk_rows=chunk_rows):
        sij = sij.astype(np.float64, copy=False)
        Sij_Rsum  += float(np.sum(sij * (R[i] + R[j])))
        Sij_sqsum += float(np.sum(sij * sij))
        
    # Sum over all i<j:
    # C_ij = R_i R_j - S_ij(R_i+R_j) + S_ij^2
    total_C = total_RiRj - Sij_Rsum + Sij_sqsum
    
    return float(total_C / num_pairs)


def mean_jaccard_dot(X: sp.spmatrix, chunk_rows: int = 50_000) -> float:
    """
    Mean pairwise Jaccard across taxa for X = taxa × samples.

    Mean is over all pairs among non-empty taxa:
        mean_{i<j} |Ti∩Tj| / |Ti∪Tj|
    """
    
    deg_all = np.diff(X.indptr).astype(np.int64, copy=False)
    nonempty = deg_all > 0
    n = int(nonempty.sum())
    if n < 2:
        return np.nan

    # Restrict to non-empty taxa (keeps unions meaningful)
    X = X[nonempty, :].tocsr(copy=False)
    deg = np.diff(X.indptr).astype(np.float64, copy=False)

    total_pairs = n * (n - 1) / 2.0

    # intersections
    S = (X @ X.T).tocsr(copy=False)

    total = 0.0
    for i, j, inter in stream_csr_upper_threshold(S, threshold=0.0, chunk_rows=chunk_rows):
        inter = inter.astype(np.float64, copy=False)
        union = deg[i] + deg[j] - inter
        m = union > 0
        if np.any(m):
            total += float(np.sum(inter[m] / union[m]))

    return float(total / total_pairs)


def structure_obj(
    ingredients: "Ingredients",
    null_model: str = "FE",               # "FE" or "FF"
    nm_n_reps: int = 1000,
    compute_null: bool = True,
    nm_random_state: int | None = 42,
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
    return _structure_core(
        ingredients=ingredients,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_random_state=nm_random_state,
        compute_null=compute_null,
    )

def _structure_core(
    ingredients: "Ingredients",
    null_model: str = "FF",
    nm_n_reps: int = 0,
    nm_random_state: int | None = None,
    compute_null: bool = False,
    chunk_rows: int = 50_000,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_sort_indices: bool = False,
    nm_burn_in_steps: int | None = None,   # FF only
    nm_steps_per_rep: int | None = None,   # FF only
    nm_progress_every: int = 25,           # tqdm granularity inside workers
) -> pd.DataFrame:
    """
    Community-structure metrics on the observed presence/absence matrix, with optional null distributions.
    
    Observed metrics:
      - C-score
      - Mean Jaccard (dot-based)
      - NODF (streamed)
      
    Null handling:
      - Uses parallel_null_reduce_vector(...) over null_matrices(...).
      - For FF, each worker corresponds to one independent Markov chain (one seed),
        and processes reps_per_worker replicates from that chain.
      - tqdm progress updates are driven by worker-side progress events.
      
    Requirements
    ------------
    - stat_fn_structure_metrics must be defined at *module scope* (picklable) in metacooc.null_models
      and must compute a length-3 vector [c_score, mean_jaccard, nodf] using worker globals
      (e.g., chunk_rows) as needed.
      
    Returns
    -------
    pd.DataFrame with one row per metric. If compute_null, attaches null summary columns.
    """
    # --- observed matrix (CSR, normalised once) ---
    X_obs = ingredients.presence_matrix.tocsr()
    X_obs.eliminate_zeros()
    X_obs.sum_duplicates()
    X_obs.sort_indices()
    
    # ---- observed metrics ----
    rows: list[dict] = []
    
    # C-score
    cscore_obs = np.nan
    cscore_err = None
    try:
        cscore_obs = float(compute_c_score(X_obs, chunk_rows=int(chunk_rows)))
    except Exception as e:
        cscore_err = str(e)
        
    rows.append(
        {
            "metric": "c_score",
            "obs": cscore_obs if np.isfinite(cscore_obs) else np.nan,
            "obs_error": cscore_err,
        }
    )
    
    # Mean Jaccard
    mj_obs = np.nan
    mj_err = None
    try:
        mj_obs = float(mean_jaccard_dot(X_obs, chunk_rows=int(chunk_rows)))
    except Exception as e:
        mj_err = str(e)
        
    rows.append(
        {
            "metric": "mean_jaccard",
            "obs": mj_obs if np.isfinite(mj_obs) else np.nan,
            "obs_error": mj_err,
        }
    )
    
    # NODF
    nodf_obs = np.nan
    nodf_err = None
    try:
        nodf_obs = float(compute_nodf_streamed(X_obs, chunk_rows=int(chunk_rows)))
    except MemoryError as e:
        nodf_err = str(e)
    except Exception as e:
        nodf_err = str(e)
        
    rows.append(
        {
            "metric": "nodf",
            "obs": nodf_obs if np.isfinite(nodf_obs) else np.nan,
            "obs_error": nodf_err,
        }
    )
    
    # If no null requested, return observed-only
    if (not compute_null) or (nm_n_reps is None) or (int(nm_n_reps) <= 0):
        return pd.DataFrame(rows)
        
    # ---- null reduction (vectorised over the 3 metrics) ----
    suffix = str(null_model).upper()
    n_reps = int(nm_n_reps)
    
    obs_vec = np.array([cscore_obs, mj_obs, nodf_obs], dtype=float)
    if not np.all(np.isfinite(obs_vec)):
        return pd.DataFrame(rows)
        
    mp_start = _best_mp_start() if nm_mp_start is None else str(nm_mp_start)
    
    j_res = parallel_null_reduce_vector(
        X=X_obs,
        model=suffix,
        n_reps=n_reps,
        obs=obs_vec,
        stat_fn=stat_fn_structure_metrics,
        random_state=nm_random_state,
        n_workers=nm_n_workers,
        mp_start=mp_start,
        sort_indices=nm_sort_indices,
        burn_in_steps=nm_burn_in_steps,
        steps_per_rep=nm_steps_per_rep,
        progress_every=int(nm_progress_every),
        # spawn-safe init kwargs (simple types only)
        chunk_rows=int(chunk_rows),
        structure_do_nodf=True,
    )
    
    # ---- attach null summaries back to rows ----
    out_rows: list[dict] = []
    for i, r in enumerate(rows):
        payload = dict(r)
        payload[f"null_mean_{suffix}"] = float(j_res["mean"][i]) if np.isfinite(j_res["mean"][i]) else np.nan
        payload[f"null_sd_{suffix}"]   = float(j_res["sd"][i])   if np.isfinite(j_res["sd"][i])   else np.nan
        payload[f"ses_{suffix}"]       = float(j_res["ses"][i])  if np.isfinite(j_res["ses"][i])  else np.nan
        payload[f"p_emp_{suffix}"]     = float(j_res["p_emp"][i]) if np.isfinite(j_res["p_emp"][i]) else np.nan
        payload[f"n_ok_{suffix}"]      = int(j_res["n_ok"])
        payload[f"n_err_{suffix}"]     = int(j_res["n_err"])
        payload[f"n_done_{suffix}"]    = int(j_res["n_done"])
        payload[f"n_requested_{suffix}"]    = int(j_res["n_target"])
        out_rows.append(payload)

    return pd.DataFrame(out_rows)