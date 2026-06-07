#structure.py
import os
import warnings

import numpy as np
import scipy.sparse as sp
import pandas as pd

from metacooc.null_models import (
    _best_mp_start,
    parallel_null_reduce_vector,
    stat_fn_structure_metrics,
)

from metacooc.pantry import load_ingredients, presence_for_counts
from metacooc.output import null_metadata_from_reduction, write_result_metadata_sidecar


_LARGE_STRUCTURE_PAIR_WARNING = 100_000_000

STRUCTURE_BASE_COLUMNS = [
    "metric",
    "observed_value",
    "observed_error",
]

STRUCTURE_NULL_COLUMNS = [
    "null_mean",
    "null_sd",
    "null_standardized_effect_size",
    "null_p_empirical",
]


def _finalize_structure_output(df: pd.DataFrame, *, include_null: bool) -> pd.DataFrame:
    columns = list(STRUCTURE_BASE_COLUMNS)
    if include_null:
        columns.extend(STRUCTURE_NULL_COLUMNS)
    attrs = dict(df.attrs)
    out = df.loc[:, columns].copy()
    out.attrs.update(attrs)
    return out


def _warn_if_large_structure_matrix(X: sp.spmatrix, chunk_rows: int) -> None:
    n_rows, n_cols = X.shape
    row_pairs = n_rows * (n_rows - 1) // 2
    col_pairs = n_cols * (n_cols - 1) // 2
    if max(row_pairs, col_pairs) < _LARGE_STRUCTURE_PAIR_WARNING:
        return

    warnings.warn(
        "Structure metrics use bounded-memory blocking controlled by "
        f"chunk_rows={int(chunk_rows)}, but runtime may be substantial for "
        f"{n_rows} taxa x {n_cols} samples.",
        RuntimeWarning,
        stacklevel=2,
    )


def _check_chunk_rows(chunk_rows: int) -> int:
    chunk_rows = int(chunk_rows)
    if chunk_rows <= 0:
        raise ValueError("chunk_rows must be a positive integer")
    return chunk_rows


def _iter_upper_overlap_blocks(
    X: sp.csr_matrix,
    chunk_rows: int,
):
    """
    Yield strict upper-triangle entries from X @ X.T one row block at a time.

    The full overlap matrix is never materialized; each yielded block is backed
    only by the sparse product for X[start:end] @ X.T.
    """
    chunk_rows = _check_chunk_rows(chunk_rows)
    n_rows = X.shape[0]
    if n_rows < 2:
        return

    XT = X.T.tocsr(copy=False)
    index_dtype = np.result_type(X.indices.dtype, np.int64)

    for r0 in range(0, n_rows, chunk_rows):
        r1 = min(r0 + chunk_rows, n_rows)
        block = (X[r0:r1, :] @ XT).tocsr(copy=False)
        block.eliminate_zeros()

        indptr = block.indptr
        indices = block.indices
        data = block.data

        total_keep = 0
        for local_i in range(r1 - r0):
            global_i = r0 + local_i
            start = indptr[local_i]
            end = indptr[local_i + 1]
            if start == end:
                continue
            total_keep += int(np.count_nonzero(indices[start:end] > global_i))

        if total_keep == 0:
            continue

        rows = np.empty(total_keep, dtype=index_dtype)
        cols = np.empty(total_keep, dtype=index_dtype)
        vals = np.empty(total_keep, dtype=data.dtype)

        pos = 0
        for local_i in range(r1 - r0):
            global_i = r0 + local_i
            start = indptr[local_i]
            end = indptr[local_i + 1]
            if start == end:
                continue

            block_cols = indices[start:end]
            keep = block_cols > global_i
            k = int(np.count_nonzero(keep))
            if k == 0:
                continue

            rows[pos:pos + k] = global_i
            cols[pos:pos + k] = block_cols[keep]
            vals[pos:pos + k] = data[start:end][keep]
            pos += k

        yield rows, cols, vals


def _nodf_sum_from_blocks(
    X: sp.csr_matrix,
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

    for i, j, ov in _iter_upper_overlap_blocks(X, chunk_rows=chunk_rows):
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
    Exact NODF with bounded-memory blocked overlap reductions.
    Returns 0..100 scale.
    """

    X = presence_for_counts(X).tocsr(copy=True)
    X.eliminate_zeros()
    X.sum_duplicates()
    X.sort_indices()
    n_rows, n_cols = X.shape
    if n_rows < 2 or n_cols < 2:
        return np.nan
        
    # row part
    row_deg = np.diff(X.indptr)
    nodf_rows_sum, row_pairs = _nodf_sum_from_blocks(X, row_deg, chunk_rows=chunk_rows)
    
    # col part
    X_by_col = X.T.tocsr(copy=False)
    col_deg = np.diff(X_by_col.indptr)
    nodf_cols_sum, col_pairs = _nodf_sum_from_blocks(
        X_by_col,
        col_deg,
        chunk_rows=chunk_rows,
    )
    
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
    
    X = presence_for_counts(X).tocsr(copy=True)
    X.eliminate_zeros()
    X.sum_duplicates()
    X.sort_indices()
    R = np.diff(X.indptr).astype(np.float64, copy=False)
    n_taxa = R.size
    if n_taxa < 2:
        return np.nan
        
    num_pairs = n_taxa * (n_taxa - 1) / 2.0
    sum_R = R.sum()
    sum_R2 = np.square(R).sum()
    total_RiRj = (sum_R * sum_R - sum_R2) / 2.0
    
    # We need sum_{i<j} S_ij*(R_i+R_j) and sum_{i<j} S_ij^2.
    Sij_Rsum = 0.0
    Sij_sqsum = 0.0
    
    for i, j, sij in _iter_upper_overlap_blocks(X, chunk_rows=chunk_rows):
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
    
    X = presence_for_counts(X).tocsr(copy=True)
    X.eliminate_zeros()
    X.sum_duplicates()
    X.sort_indices()
    deg_all = np.diff(X.indptr).astype(np.int64, copy=False)
    nonempty = deg_all > 0
    n = int(nonempty.sum())
    if n < 2:
        return np.nan

    # Restrict to non-empty taxa (keeps unions meaningful)
    X = X[nonempty, :].tocsr(copy=False)
    deg = np.diff(X.indptr).astype(np.float64, copy=False)

    total_pairs = n * (n - 1) / 2.0

    total = 0.0
    for i, j, inter in _iter_upper_overlap_blocks(X, chunk_rows=chunk_rows):
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
    nm_seed: int | None = None,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_burn_in_steps: int | None = None,
    nm_steps_per_rep: int | None = None,
    nm_progress_every: int = 25,
) -> pd.DataFrame:
    """
    Entry point for community-structure analysis.

    Parameters
    ----------
    ingredients : Ingredients
        Cohort matrix used to compute C-score, mean Jaccard, and NODF.
    null_model : {"FE", "EF", "EE", "FF"}
        Null model passed to the structure null reduction when compute_null is True.
    compute_null : bool
        If True, compute null summaries for the observed structure metrics.
    nm_n_reps : int
        Number of null replicates if compute_null is True.
    nm_seed : int or None
        Random seed for null generation.

    Returns
    -------
    out : pd.DataFrame
        One row per structure metric. Run-level null metadata, if present, is
        stored in out.attrs["null_metadata"].
    """
    return _structure_core(
        ingredients=ingredients,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_seed=nm_seed,
        compute_null=compute_null,
        nm_n_workers=nm_n_workers,
        nm_mp_start=nm_mp_start,
        nm_burn_in_steps=nm_burn_in_steps,
        nm_steps_per_rep=nm_steps_per_rep,
        nm_progress_every=nm_progress_every,
    )

def _structure_core(
    ingredients: "Ingredients",
    null_model: str = "FF",
    nm_n_reps: int = 0,
    nm_seed: int | None = None,
    compute_null: bool = False,
    chunk_rows: int = 50_000,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
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
    X_obs = presence_for_counts(ingredients).tocsr()
    X_obs.eliminate_zeros()
    X_obs.sum_duplicates()
    X_obs.sort_indices()
    _warn_if_large_structure_matrix(X_obs, chunk_rows=int(chunk_rows))
    
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
            "observed_value": cscore_obs if np.isfinite(cscore_obs) else np.nan,
            "observed_error": cscore_err,
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
            "observed_value": mj_obs if np.isfinite(mj_obs) else np.nan,
            "observed_error": mj_err,
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
            "observed_value": nodf_obs if np.isfinite(nodf_obs) else np.nan,
            "observed_error": nodf_err,
        }
    )
    
    # If no null requested, return observed-only
    if (not compute_null) or (nm_n_reps is None) or (int(nm_n_reps) <= 0):
        return _finalize_structure_output(pd.DataFrame(rows), include_null=False)
        
    # ---- null reduction (vectorised over the 3 metrics) ----
    suffix = str(null_model).upper()
    n_reps = int(nm_n_reps)
    
    obs_vec = np.array([cscore_obs, mj_obs, nodf_obs], dtype=float)
    if not np.all(np.isfinite(obs_vec)):
        return _finalize_structure_output(pd.DataFrame(rows), include_null=False)
        
    mp_start = _best_mp_start() if nm_mp_start is None else str(nm_mp_start)
    
    j_res = parallel_null_reduce_vector(
        X=X_obs,
        model=suffix,
        n_reps=n_reps,
        obs=obs_vec,
        stat_fn=stat_fn_structure_metrics,
        seed=nm_seed,
        n_workers=nm_n_workers,
        mp_start=mp_start,
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
        payload["null_mean"] = float(j_res["mean"][i]) if np.isfinite(j_res["mean"][i]) else np.nan
        payload["null_sd"] = float(j_res["sd"][i]) if np.isfinite(j_res["sd"][i]) else np.nan
        payload["null_standardized_effect_size"] = (
            float(j_res["ses"][i]) if np.isfinite(j_res["ses"][i]) else np.nan
        )
        payload["null_p_empirical"] = (
            float(j_res["p_emp"][i]) if np.isfinite(j_res["p_emp"][i]) else np.nan
        )
        out_rows.append(payload)

    out = pd.DataFrame(out_rows)
    out.attrs["null_metadata"] = null_metadata_from_reduction(j_res)

    print(
        f"INFO: Null simulation seed {j_res['null_seed']} "
        f"({j_res['null_seed_source']}, model={j_res['null_model']})."
    )
        
    return _finalize_structure_output(out, include_null=True)


def structure(
    ingredients,
    output_dir: str,
    tag: str | None = None,
    null_model: str = "FE",
    nm_n_reps: int = 1000,
    compute_null: bool = True,
    nm_seed: int | None = None,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_burn_in_steps: int | None = None,
    nm_steps_per_rep: int | None = None,
    nm_progress_every: int = 25,
) -> pd.DataFrame:
    """
    
    This is the file-writing wrapper around :func:`structure_obj` that runs community-structure analysis and write results to disk.
    
    Parameters
    ----------
    ingredients :
        Input Ingredients object, or any value accepted by
        :func:`metacooc.pantry.load_ingredients` via ``custom_ingredients``.
    output_dir : str
        Directory where the output TSV file will be written.
    tag : str, optional
        Optional prefix for the output filename.
    null_model : str, default "FE"
        Null model passed through to :func:`structure_obj`.
    nm_n_reps : int, default 1000
        Number of null replicates.
    compute_null : bool, default True
        Whether to compute null distributions.
    nm_seed : int or None
        Random seed for null generation. If None, an entropy-backed seed is generated.
        
    Returns
    -------
    pd.DataFrame
        Structure metrics table.
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    ing = load_ingredients(data_dir=None, custom_ingredients=ingredients)
    
    structure_df = structure_obj(
        ingredients=ing,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        compute_null=compute_null,
        nm_seed=nm_seed,
        nm_n_workers=nm_n_workers,
        nm_mp_start=nm_mp_start,
        nm_burn_in_steps=nm_burn_in_steps,
        nm_steps_per_rep=nm_steps_per_rep,
        nm_progress_every=nm_progress_every,
    )
    
    output_path = os.path.join(output_dir, f"{tag}structure.tsv")
    structure_df.to_csv(output_path, sep="\t", index=False)
    print(f"Structure analysis saved to {output_path}")
    metadata_path = write_result_metadata_sidecar(
        output_path,
        structure_df.attrs.get("null_metadata"),
    )
    if metadata_path is not None:
        print(f"Structure metadata saved to {metadata_path}")
    
    return structure_df
