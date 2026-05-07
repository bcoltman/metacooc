#!/usr/bin/env python3
# metacooc/analytics/analysis.py

from __future__ import annotations
from typing import Optional, Iterable, Tuple, List, Set, Callable
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import chi2 as chi2_dist, fisher_exact as _fisher_exact
from scipy.sparse import csr_matrix, csc_matrix
from scipy import special
import multiprocessing as mp
import os
from dataclasses import dataclass, field


from metacooc.pantry import *
from metacooc.utils import _RANK_PREFIXES, stream_edge_values, _stream_csr_entries

from metacooc.null_models import (
    parallel_null_reduce_vector,
    stat_fn_association_jaccard,
    stat_fn_cooccurrence_jaccard,
    _best_mp_start
)


_SMOOTH = 0.5


@dataclass
class CooccurrenceArrays:
    cols: dict[str, np.ndarray]
    meta: dict[str, object] = field(default_factory=dict)

    @property
    def n_rows(self) -> int:
        if not self.cols:
            return 0
        first_key = next(iter(self.cols))
        return int(len(self.cols[first_key]))


def _chi2_phi_from_counts(a, b, c, d):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    c = np.asarray(c, dtype=np.float64)
    d = np.asarray(d, dtype=np.float64)

    row1 = a + b
    row2 = c + d
    col1 = a + c
    col2 = b + d
    n = row1 + row2

    cross = a * d
    cross -= b * c

    denom_sq = row1 * row2
    denom_sq *= col1
    denom_sq *= col2

    with np.errstate(divide="ignore", invalid="ignore"):
        phi = cross / np.sqrt(denom_sq)
        chi2 = n * (cross * cross) / denom_sq

    invalid = (
        (a < 0) | (b < 0) | (c < 0) | (d < 0) |
        (row1 == 0) | (row2 == 0) |
        (col1 == 0) | (col2 == 0) |
        ~np.isfinite(chi2) | ~np.isfinite(phi)
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

    chi2, phi, invalid = _chi2_phi_from_counts(a, b, c, d)
    log_p = chi2_dist.logsf(chi2, df=1)
    p = np.exp(log_p)

    P_B_given_A = (a + _SMOOTH) / (a + b + _SMOOTH)
    P_B_given_notA = (c + _SMOOTH) / (c + d + _SMOOTH)
    RR_A_to_B = P_B_given_A / P_B_given_notA
    with np.errstate(divide="ignore", invalid="ignore"):
        logRR_A_to_B = np.log(RR_A_to_B)

    P_A_given_B = (a + _SMOOTH) / (a + c + _SMOOTH)
    P_A_given_notB = (b + _SMOOTH) / (b + d + _SMOOTH)
    RR_B_to_A = P_A_given_B / P_A_given_notB
    with np.errstate(divide="ignore", invalid="ignore"):
        logRR_B_to_A = np.log(RR_B_to_A)

    n_tables = a.shape[0]
    fisher_odds = np.full(n_tables, np.nan, dtype=float)
    fisher_p = np.full(n_tables, np.nan, dtype=float)
    log_fisher_p = np.full(n_tables, np.nan, dtype=float)

    if compute_fisher:
        valid_idx = np.where(~invalid)[0]
        for i in valid_idx:
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


def select_taxa_universe(
    ing: "Ingredients",
    rank: Optional[str] = None,
    terminal_only: bool = True,
    require_present: bool = True
) -> List[str]:
    """
    Build the taxa universe to analyse, leveraging Ingredients' caches and
    current presence.
    """
    if require_present:
        present_set = {t for t, cnt in zip(ing.taxa, ing.total_counts) if cnt > 0}
    else:
        present_set = set(ing.taxa)

    base = [t for t in ing.taxa if t in present_set]

    if rank is None:
        return base

    rk = rank.lower().strip()
    if rk not in _RANK_PREFIXES:
        raise ValueError(
            f"Unknown rank '{rank}'. Expected one of: {', '.join(_RANK_PREFIXES.keys())}"
        )
    pref = _RANK_PREFIXES[rk]

    ing._ensure_taxa_lookups()

    if terminal_only:
        term_map = {t: tp for t, tp in zip(ing.taxa, ing._terminal_rank_prefixes)}
        return [t for t in base if term_map.get(t) == pref]
    else:
        return [t for t in base if pref in t]


# ========== Workflow A: filtered vs null (single-taxon enrichment) ==========

def association_obj(
    null_ingredients: "Ingredients",
    filtered_ingredients: "Ingredients",
    threshold: float = 0.0,
    null_model: str = "FE",
    nm_n_reps: int = 1000,
    nm_random_state: int | None = None,
    nm_memory_mb: float | None = None,
    compute_fisher: bool = False,
) -> pd.DataFrame:
    out = _association_core(
        null_ingredients=null_ingredients,
        filtered_ingredients=filtered_ingredients,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_random_state=nm_random_state,
        nm_memory_mb=nm_memory_mb,
    )

    if threshold is not None:
        out = out[out["p_T_given_X"] > threshold].copy()

    return out


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
    nm_burn_in_steps: int | None = None,
    nm_steps_per_rep: int | None = None,
    nm_memory_mb: float | None = None,
) -> pd.DataFrame:
    """
    Core taxon-term enrichment (no community-structure metrics).
    """
    null_samples = set(null_ingredients.samples)
    filt_samples = set(filtered_ingredients.samples)
    if not filt_samples.issubset(null_samples):
        raise ValueError("filtered_ingredients.samples must be a subset of null_ingredients.samples")

    null_taxa = np.array(null_ingredients.taxa, dtype=object)
    filt_taxa = np.array(filtered_ingredients.taxa, dtype=object)
    taxa, filt_idx, null_idx = np.intersect1d(filt_taxa, null_taxa, return_indices=True)
    if taxa.size == 0:
        raise ValueError("No taxa intersect the null and filtered ingredients.")

    N_T = float(len(filtered_ingredients.samples))
    N_null = float(len(null_ingredients.samples))
    if N_T <= 0 or N_null <= 0 or N_null < N_T:
        raise ValueError(f"Invalid cohort sizes: N_T={N_T}, N_null={N_null}")
    if N_null == N_T:
        raise ValueError("Null and term cohorts are identical (N_null == N_T): no non-term samples available.")

    a_raw = filtered_ingredients.total_counts[filt_idx].astype(float, copy=False)
    ref_counts = null_ingredients.total_counts[null_idx].astype(float, copy=False)
    b_raw = ref_counts - a_raw

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

    out["log_q_bh"] = bh_logq_from_logp(out["log_p"].to_numpy(dtype=np.float64, copy=False))
    out["q_bh"] = np.exp(out["log_q_bh"].to_numpy(dtype=np.float64, copy=False))

    valid_mask = ~out["invalid_table"].values
    out = out.loc[valid_mask].copy()
    out.drop(columns="invalid_table", inplace=True)

    null_idx_valid = null_idx[valid_mask]
    out.attrs["null_idx_valid"] = np.asarray(null_idx_valid, dtype=int)
    out.attrs["taxa_valid"] = out["taxon"].to_numpy(dtype=object, copy=False)

    n_reps = int(nm_n_reps) if nm_n_reps is not None else 0
    if n_reps <= 0 or out.shape[0] == 0:
        return out

    if null_model == "FE":
        print("FE: association determined analytically - no need for shuffling null and probabilistic approach")
        return out

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

    obs_jacc = out["jaccard"].to_numpy(dtype=float, copy=False)

    mp_start = _best_mp_start() if nm_mp_start is None else str(nm_mp_start)

    j_res = parallel_null_reduce_vector(
        X=X_full,
        model=null_model,
        n_reps=n_reps,
        obs=obs_jacc,
        stat_fn=stat_fn_association_jaccard,
        random_state=nm_random_state,
        n_workers=nm_n_workers,
        mp_start=mp_start,
        memory_mb=nm_memory_mb,
        term_cols=term_cols,
        nonterm_cols=nonterm_cols,
        subset_idx=subset_idx,
        n_rows=n_rows,
        N_T=N_T,
    )

    out[f"jaccard_null_mean_{null_model}"] = j_res["mean"]
    out[f"jaccard_null_sd_{null_model}"] = j_res["sd"]
    out[f"jaccard_ses_{null_model}"] = j_res["ses"]
    out[f"jaccard_p_{null_model}"] = j_res["p_emp"]

    out[f"n_ok_{null_model}"] = int(j_res["n_ok"])
    out[f"n_err_{null_model}"] = int(j_res["n_err"])
    out[f"n_done_{null_model}"] = int(j_res["n_done"])
    out[f"n_requested_{null_model}"] = int(j_res["n_target"])

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
    nm_memory_mb: float | None = None,
    compute_fisher: bool = False,
) -> pd.DataFrame:
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    null = load_ingredients(data_dir=None, custom_ingredients=null_ingredients)
    filtered = load_ingredients(data_dir=None, custom_ingredients=filtered_ingredients)

    association_df = association_obj(
        null_ingredients=null,
        filtered_ingredients=filtered,
        threshold=threshold,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_random_state=nm_random_state,
        nm_memory_mb=nm_memory_mb,
        compute_fisher=compute_fisher
    )

    if association_df is not None:
        output_path = os.path.join(output_dir, f"{tag}association.tsv")
        association_df.to_csv(output_path, sep="\t", index=False)
        print(f"Association analysis saved to {output_path}")


def presence_submatrix_by_taxa(ingredients: Ingredients, taxa_subset: List[str]) -> sp.csr_matrix:
    row_map = {t: i for i, t in enumerate(ingredients.taxa)}
    rows = [row_map[t] for t in taxa_subset]
    return ingredients.presence_matrix[rows, :].tocsr()


def _cooccur_core(
    ing: Ingredients,
    taxa_universe: List[str],
    threshold: float = 0.1,
    m_total: Optional[int] = None,
    compute_fisher: bool = False,
) -> Tuple[Optional[CooccurrenceArrays], pd.DataFrame]:
    """
    Pairwise co-occurrence edges (within 'ing') and node summary.
    """
    X_sub = presence_submatrix_by_taxa(ing, taxa_universe)
    totals = np.asarray(X_sub.sum(axis=1)).ravel().astype(np.int32, copy=False)
    co_mat = (X_sub @ X_sub.T).tocsr()
    N_total = int(len(ing.samples))

    n_taxa = len(totals)
    deg_fwd = np.zeros(n_taxa, dtype=np.int64)

    n_keep = 0
    for iA, iB, inter in stream_edge_values(co_mat, min_value=0):
        A_totals = totals[iA]

        if threshold > 0:
            keep = inter > (threshold * A_totals)
            if not np.any(keep):
                continue
            iA = iA[keep]
            iB = iB[keep]
            inter = inter[keep]
            A_totals = A_totals[keep]

        B_totals = totals[iB]
        valid = _valid_edge_mask_from_counts(inter, A_totals, B_totals, N_total)
        if not np.any(valid):
            continue

        iA_valid = iA[valid]
        n_keep += int(iA_valid.size)
        deg_fwd += np.bincount(iA_valid, minlength=n_taxa)

    nodes_df = pd.DataFrame({
        "taxon": taxa_universe,
        "total_count": totals.astype(int, copy=False),
        f"degree_PBA_gt_{threshold}": deg_fwd.astype(int, copy=False),
    })

    if n_keep == 0:
        return None, nodes_df

    iA_all = np.empty(n_keep, dtype=np.int32)
    iB_all = np.empty(n_keep, dtype=np.int32)
    inter_all = np.empty(n_keep, dtype=np.int32)
    log_p_all = np.empty(n_keep, dtype=np.float64)
    jaccard_all = np.empty(n_keep, dtype=np.float32)

    pos = 0
    for iA, iB, inter in stream_edge_values(co_mat, min_value=0):
        A_totals = totals[iA]

        if threshold > 0:
            keep = inter > (threshold * A_totals)
            if not np.any(keep):
                continue
            iA = iA[keep]
            iB = iB[keep]
            inter = inter[keep]
            A_totals = A_totals[keep]

        B_totals = totals[iB]
        valid = _valid_edge_mask_from_counts(inter, A_totals, B_totals, N_total)
        if not np.any(valid):
            continue

        iA = iA[valid].astype(np.int32, copy=False)
        iB = iB[valid].astype(np.int32, copy=False)
        inter = inter[valid].astype(np.int32, copy=False)
        A_totals = A_totals[valid]
        B_totals = B_totals[valid]

        log_p_chunk, jaccard_chunk = _compute_core_edge_metrics(
            inter=inter,
            A_totals=A_totals,
            B_totals=B_totals,
            N_total=N_total,
        )

        k = int(inter.size)
        sl = slice(pos, pos + k)

        iA_all[sl] = iA
        iB_all[sl] = iB
        inter_all[sl] = inter
        log_p_all[sl] = log_p_chunk
        jaccard_all[sl] = jaccard_chunk

        pos += k

    edge_arrays = CooccurrenceArrays(
        cols={
            "iA": iA_all,
            "iB": iB_all,
            "inter": inter_all,
            "jaccard": jaccard_all,
            "log_p": log_p_all,
        },
        meta={
            "totals": totals,
            "N_total": N_total,
            "compute_fisher": bool(compute_fisher),
        },
    )

    edge_arrays.cols["log_q_bh"] = bh_logq_from_logp(log_p_all, m_total=m_total)

    return edge_arrays, nodes_df


def _cooccur_core_focal(
    ing: Ingredients,
    taxa_universe: List[str],
    focal_local_idx: np.ndarray,
    threshold: float = 0.1,
    m_total: Optional[int] = None,
    compute_fisher: bool = False,
) -> Tuple[Optional[CooccurrenceArrays], pd.DataFrame]:
    """
    Focal-only co-occurrence core.

    Computes only focal x all intersections and emits edges where A_taxon is focal.
    Focal-focal pairs are emitted once with orientation normalised to iA < iB.
    """
    X_sub = presence_submatrix_by_taxa(ing, taxa_universe).tocsr()
    totals = np.asarray(X_sub.sum(axis=1)).ravel().astype(np.int32, copy=False)
    N_total = int(len(ing.samples))

    focal_local_idx = np.asarray(focal_local_idx, dtype=np.int64)
    focal_set = set(focal_local_idx.tolist())

    nodes_df = pd.DataFrame({
        "taxon": taxa_universe,
        "total_count": totals.astype(int, copy=False),
        f"degree_PBA_gt_{threshold}": np.zeros(len(taxa_universe), dtype=int),
    })

    if focal_local_idx.size == 0:
        return None, nodes_df

    X_focal = X_sub[focal_local_idx, :]
    co_focal = (X_focal @ X_sub.T).tocsr()

    n_taxa = len(taxa_universe)
    deg_fwd = np.zeros(n_taxa, dtype=np.int64)

    n_keep = 0
    for f_rows, j_cols, inter in _stream_csr_entries(co_focal, min_value=0):
        iA = focal_local_idx[f_rows]
        A_totals = totals[iA]

        nonself = iA != j_cols
        if not np.any(nonself):
            continue
        iA = iA[nonself]
        j_cols = j_cols[nonself]
        inter = inter[nonself]
        A_totals = A_totals[nonself]

        if threshold > 0:
            keep = inter > (threshold * A_totals)
            if not np.any(keep):
                continue
            iA = iA[keep]
            j_cols = j_cols[keep]
            inter = inter[keep]
            A_totals = A_totals[keep]

        B_totals = totals[j_cols]
        valid = _valid_edge_mask_from_counts(inter, A_totals, B_totals, N_total)
        if not np.any(valid):
            continue

        iA = iA[valid]
        j_cols = j_cols[valid]

        both_focal = np.fromiter((j in focal_set for j in j_cols), dtype=bool, count=j_cols.size)
        if np.any(both_focal):
            keep2 = ~both_focal | (iA < j_cols)
            if not np.any(keep2):
                continue
            iA = iA[keep2]
            j_cols = j_cols[keep2]

        n_keep += int(iA.size)
        deg_fwd += np.bincount(iA, minlength=n_taxa)

    nodes_df[f"degree_PBA_gt_{threshold}"] = deg_fwd.astype(int, copy=False)

    if n_keep == 0:
        return None, nodes_df

    iA_all = np.empty(n_keep, dtype=np.int32)
    iB_all = np.empty(n_keep, dtype=np.int32)
    inter_all = np.empty(n_keep, dtype=np.int32)
    log_p_all = np.empty(n_keep, dtype=np.float64)
    jaccard_all = np.empty(n_keep, dtype=np.float32)

    pos = 0
    for f_rows, j_cols, inter in _stream_csr_entries(co_focal, min_value=0):
        iA = focal_local_idx[f_rows]
        A_totals = totals[iA]

        nonself = iA != j_cols
        if not np.any(nonself):
            continue
        iA = iA[nonself]
        j_cols = j_cols[nonself]
        inter = inter[nonself]
        A_totals = A_totals[nonself]

        if threshold > 0:
            keep = inter > (threshold * A_totals)
            if not np.any(keep):
                continue
            iA = iA[keep]
            j_cols = j_cols[keep]
            inter = inter[keep]
            A_totals = A_totals[keep]

        B_totals = totals[j_cols]
        valid = _valid_edge_mask_from_counts(inter, A_totals, B_totals, N_total)
        if not np.any(valid):
            continue

        iA = iA[valid]
        iB = j_cols[valid]
        inter = inter[valid].astype(np.int32, copy=False)
        A_totals = A_totals[valid]
        B_totals = B_totals[valid]

        both_focal = np.fromiter((j in focal_set for j in iB), dtype=bool, count=iB.size)
        if np.any(both_focal):
            keep2 = ~both_focal | (iA < iB)
            if not np.any(keep2):
                continue
            iA = iA[keep2]
            iB = iB[keep2]
            inter = inter[keep2]
            A_totals = A_totals[keep2]
            B_totals = B_totals[keep2]

        log_p_chunk, jaccard_chunk = _compute_core_edge_metrics(
            inter=inter,
            A_totals=A_totals,
            B_totals=B_totals,
            N_total=N_total,
        )

        k = int(inter.size)
        sl = slice(pos, pos + k)

        iA_all[sl] = iA.astype(np.int32, copy=False)
        iB_all[sl] = iB.astype(np.int32, copy=False)
        inter_all[sl] = inter
        log_p_all[sl] = log_p_chunk
        jaccard_all[sl] = jaccard_chunk
        pos += k

    edge_arrays = CooccurrenceArrays(
        cols={
            "iA": iA_all,
            "iB": iB_all,
            "inter": inter_all,
            "jaccard": jaccard_all,
            "log_p": log_p_all,
        },
        meta={
            "totals": totals,
            "N_total": N_total,
            "compute_fisher": bool(compute_fisher),
        },
    )

    edge_arrays.cols["log_q_bh"] = bh_logq_from_logp(log_p_all, m_total=m_total)
    return edge_arrays, nodes_df


def estimate_focal_pairs(n_taxa: int, n_focal: int) -> int:
    """
    Maximum number of focal-anchored output pairs.
    """
    n_taxa = int(n_taxa)
    n_focal = int(n_focal)

    if n_focal <= 0 or n_taxa <= 1:
        return 0
    if n_focal > n_taxa:
        raise ValueError(f"n_focal ({n_focal}) cannot exceed n_taxa ({n_taxa})")

    return n_focal * (n_taxa - 1)


def should_run_cooccurrence(
    n_taxa: int,
    large: bool,
    max_pairs: int = 100_000,
    n_focal: int | None = None,
) -> Tuple[bool, int]:
    """
    Decide whether to run co-occurrence.
    """
    if n_focal is None:
        pairs = (n_taxa * (n_taxa - 1)) // 2
    else:
        pairs = estimate_focal_pairs(n_taxa=n_taxa, n_focal=n_focal)

    if large:
        return True, pairs
    return (pairs <= max_pairs), pairs

def _subset_focal_edge_arrays(
    edge_arrays: CooccurrenceArrays,
    taxa_universe: List[str],
    focal_query_to_taxa: dict[str, List[str]],
) -> CooccurrenceArrays:
    """
    Restrict compact edge arrays to focal-anchored edges and orient them so the
    focal taxon is always iA / A_taxon.

    Output columns added
    --------------------
    focal_query
        The focal query responsible for emitting the row.
    focal_taxon
        The concrete focal taxon anchoring the row.

    Notes
    -----
    - Focal–focal edges are duplicated once per focal taxon.
    - If the same focal taxon is included in multiple focal queries, its edges
      are duplicated once per focal query.
    """
    tax_to_idx = {t: i for i, t in enumerate(taxa_universe)}

    src_iA = np.asarray(edge_arrays.cols["iA"], dtype=np.int64)
    src_iB = np.asarray(edge_arrays.cols["iB"], dtype=np.int64)

    collected = {k: [] for k in edge_arrays.cols}
    focal_query_parts = []
    focal_taxon_parts = []

    for focal_query, focal_taxa in focal_query_to_taxa.items():
        for focal_taxon in focal_taxa:
            focal_idx = tax_to_idx.get(focal_taxon)
            if focal_idx is None:
                continue

            left_idx = np.flatnonzero(src_iA == focal_idx)
            if left_idx.size:
                for key, arr in edge_arrays.cols.items():
                    collected[key].append(np.asarray(arr)[left_idx])
                focal_query_parts.append(np.full(left_idx.size, focal_query, dtype=object))
                focal_taxon_parts.append(np.full(left_idx.size, focal_taxon, dtype=object))

            right_idx = np.flatnonzero(src_iB == focal_idx)
            if right_idx.size:
                for key, arr in edge_arrays.cols.items():
                    arr_np = np.asarray(arr)
                    if key == "iA":
                        collected[key].append(src_iB[right_idx].astype(arr_np.dtype, copy=False))
                    elif key == "iB":
                        collected[key].append(src_iA[right_idx].astype(arr_np.dtype, copy=False))
                    else:
                        collected[key].append(arr_np[right_idx])
                focal_query_parts.append(np.full(right_idx.size, focal_query, dtype=object))
                focal_taxon_parts.append(np.full(right_idx.size, focal_taxon, dtype=object))

    out_cols = {}
    for key, parts in collected.items():
        if parts:
            out_cols[key] = np.concatenate(parts)
        else:
            out_cols[key] = np.asarray(edge_arrays.cols[key])[:0].copy()

    out_cols["focal_query"] = (
        np.concatenate(focal_query_parts) if focal_query_parts else np.array([], dtype=object)
    )
    out_cols["focal_taxon"] = (
        np.concatenate(focal_taxon_parts) if focal_taxon_parts else np.array([], dtype=object)
    )

    return CooccurrenceArrays(cols=out_cols, meta=dict(edge_arrays.meta))


def _subset_edges_by_rhs_queries(
    edge_arrays: CooccurrenceArrays,
    taxa_universe: List[str],
    rhs_query_to_taxa: dict[str, List[str]],
) -> CooccurrenceArrays | None:
    """
    Restrict an already focal-expanded edge set using RHS query taxa.

    Expected input
    --------------
    `edge_arrays` must already have:
      - iA
      - iB
      - focal_query

    Semantics
    ---------
    Keep an edge iff:
      - its focal_query exists in rhs_query_to_taxa, and
      - its B_taxon name belongs to the RHS taxa mapped to that focal_query

    Notes
    -----
    - This assumes `_subset_focal_edge_arrays()` has already run, so focal_query
      is present and iA is already oriented as the focal side.
    - Returning None here means:
          RHS taxa resolved, but no focal->RHS co-occurrence edges survived.
      That is a valid biological/data outcome, not a parsing error.
    """
    if edge_arrays is None or edge_arrays.n_rows == 0:
        return edge_arrays

    if "focal_query" not in edge_arrays.cols:
        raise ValueError(
            "_subset_edges_by_rhs_queries requires edge_arrays.cols['focal_query']. "
            "Run _subset_focal_edge_arrays first."
        )

    taxa_arr = np.asarray(taxa_universe, dtype=object)
    focal_query_arr = np.asarray(edge_arrays.cols["focal_query"], dtype=object)
    iB_arr = np.asarray(edge_arrays.cols["iB"], dtype=np.int64)
    b_taxon_arr = taxa_arr[iB_arr]

    rhs_query_to_taxa_set = {
        q: set(taxa)
        for q, taxa in rhs_query_to_taxa.items()
        if taxa
    }

    keep = np.zeros(edge_arrays.n_rows, dtype=bool)

    unique_queries = np.unique(focal_query_arr)
    for q in unique_queries:
        rhs_taxa_set = rhs_query_to_taxa_set.get(q)
        if not rhs_taxa_set:
            continue

        q_mask = (focal_query_arr == q)
        if not np.any(q_mask):
            continue

        keep[q_mask] = np.isin(b_taxon_arr[q_mask], list(rhs_taxa_set))

    if not np.any(keep):
        return None

    out_cols = {
        key: np.asarray(arr)[keep]
        for key, arr in edge_arrays.cols.items()
    }

    return CooccurrenceArrays(
        cols=out_cols,
        meta=dict(edge_arrays.meta),
    )


def cooccurrence_obj(
    null_ingredients: "Ingredients",
    taxa_universe: List[str],
    focal_query_to_taxa: dict[str, List[str]] | None = None,
    rhs_query_to_taxa: dict[str, List[str]] | None = None,
    large: bool = False,
    max_pairs: int = 100_000,
    threshold: float = 0.1,
    null_model: str = "FE",
    nm_n_reps: int = 10,
    nm_random_state: int | None = None,
    nm_memory_mb: float | None = None,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_sort_indices: bool = False,
    nm_burn_in_steps: int | None = None,
    nm_steps_per_rep: int | None = None,
) -> Tuple[Optional[CooccurrenceArrays], Optional[pd.DataFrame]]:
    """
    Pairwise co-occurrence of taxa.

    If focal_query_to_taxa is provided, computation is restricted to focal-anchored
    pairs only, using the union of resolved focal taxa across queries.

    If rhs_query_to_taxa is also provided, the already focal-expanded edge set is
    further restricted so that B_taxon belongs to the RHS taxa for the matching
    focal_query. This does not trigger additional null generation.
    """
    focal_local_idx = None
    if focal_query_to_taxa:
        tax_to_idx = {t: i for i, t in enumerate(taxa_universe)}
        focal_idx_set = {
            tax_to_idx[t]
            for taxa in focal_query_to_taxa.values()
            for t in taxa
            if t in tax_to_idx
        }
        focal_local_idx = np.asarray(sorted(focal_idx_set), dtype=np.int64)
        n_focal = int(focal_local_idx.size)
    else:
        n_focal = None

    run_co, est_pairs = should_run_cooccurrence(
        n_taxa=len(taxa_universe),
        large=large,
        max_pairs=max_pairs,
        n_focal=n_focal,
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

    if focal_local_idx is None:
        edge_arrays, nodes_df = _cooccur_core(
            null_ingredients,
            taxa_universe,
            threshold=threshold,
            m_total=est_pairs,
            compute_fisher=False,
        )
    else:
        edge_arrays, nodes_df = _cooccur_core_focal(
            null_ingredients,
            taxa_universe,
            focal_local_idx=focal_local_idx,
            threshold=threshold,
            m_total=est_pairs,
            compute_fisher=False,
        )

    if edge_arrays is None or edge_arrays.n_rows == 0:
        return edge_arrays, nodes_df

    if focal_query_to_taxa:
        edge_arrays = _subset_focal_edge_arrays(
            edge_arrays=edge_arrays,
            taxa_universe=taxa_universe,
            focal_query_to_taxa=focal_query_to_taxa,
        )

        if edge_arrays is None or edge_arrays.n_rows == 0:
            return edge_arrays, nodes_df

        if rhs_query_to_taxa:
            edge_arrays = _subset_edges_by_rhs_queries(
                edge_arrays=edge_arrays,
                taxa_universe=taxa_universe,
                rhs_query_to_taxa=rhs_query_to_taxa,
            )

            if edge_arrays is None or edge_arrays.n_rows == 0:
                print(
                    "No focal->RHS co-occurrence edges survived after RHS restriction. "
                    "This means the RHS taxa were valid targets, but no edges to those "
                    "B-side taxa survived under the current focal cohort, taxa universe, "
                    "and threshold settings."
                )
                return None, nodes_df

        edge_arrays.cols["log_q_bh"] = bh_logq_from_logp(
            edge_arrays.cols["log_p"],
            m_total=edge_arrays.n_rows,
        )

    n_reps = int(nm_n_reps) if nm_n_reps is not None else 0
    if n_reps <= 0:
        return edge_arrays, nodes_df

    if null_model == "FE":
        print("FE: cooccurrence determined analytically - no need for shuffling null and probabilistic approach")
        return edge_arrays, nodes_df

    X_full = null_ingredients.presence_matrix.tocsr()
    X_full.eliminate_zeros()
    X_full.sum_duplicates()
    X_full.sort_indices()

    tax_map = {t: i for i, t in enumerate(null_ingredients.taxa)}
    subset_idx = np.array([tax_map[t] for t in taxa_universe], dtype=np.int64)

    obs_jacc = edge_arrays.cols["jaccard"].astype(np.float64, copy=False)
    iA_all = edge_arrays.cols["iA"].astype(np.int64, copy=False)
    iB_all = edge_arrays.cols["iB"].astype(np.int64, copy=False)

    mp_start = _best_mp_start() if nm_mp_start is None else str(nm_mp_start)

    j_res = parallel_null_reduce_vector(
        X=X_full,
        model=null_model,
        n_reps=n_reps,
        obs=obs_jacc,
        stat_fn=stat_fn_cooccurrence_jaccard,
        random_state=nm_random_state,
        n_workers=nm_n_workers,
        mp_start=mp_start,
        memory_mb=nm_memory_mb,
        subset_idx=subset_idx,
        iA=iA_all,
        iB=iB_all,
    )

    edge_arrays.cols[f"jaccard_null_mean_{null_model}"] = np.asarray(j_res["mean"], dtype=np.float32)
    edge_arrays.cols[f"jaccard_null_sd_{null_model}"] = np.asarray(j_res["sd"], dtype=np.float32)
    edge_arrays.cols[f"jaccard_ses_{null_model}"] = np.asarray(j_res["ses"], dtype=np.float32)
    edge_arrays.cols[f"jaccard_p_{null_model}"] = np.asarray(j_res["p_emp"], dtype=np.float64)

    edge_arrays.cols[f"log_jaccard_q_{null_model}"] = bh_logq_from_logp(
        np.log(edge_arrays.cols[f"jaccard_p_{null_model}"]),
        m_total=edge_arrays.n_rows if focal_query_to_taxa else est_pairs,
    )

    edge_arrays.meta[f"n_ok_{null_model}"] = int(j_res["n_ok"])
    edge_arrays.meta[f"n_err_{null_model}"] = int(j_res["n_err"])
    edge_arrays.meta[f"n_done_{null_model}"] = int(j_res["n_done"])
    edge_arrays.meta[f"n_requested_{null_model}"] = int(j_res["n_target"])

    return edge_arrays, nodes_df

def export_cooccurrence_outputs(
    *,
    edge_arrays: CooccurrenceArrays | None,
    nodes_df: pd.DataFrame | None,
    taxa_universe: List[str],
    output_dir: str,
    edges_base: str,
    nodes_base: str,
    null_model: str = "FE",
    summary_n: int = 100_000,
) -> None:
    """
    Export co-occurrence outputs with a user-friendly policy.
    """
    os.makedirs(output_dir, exist_ok=True)

    if nodes_df is not None:
        nodes_path = os.path.join(output_dir, f"{nodes_base}.tsv")
        nodes_df.to_csv(nodes_path, sep="\t", index=False)
        print(f"Taxon nodes analysis saved to {nodes_path}")

    if edge_arrays is None or edge_arrays.n_rows == 0:
        print("No co-occurrence edges above threshold; no edge file written.")
        return

    n_rows = edge_arrays.n_rows

    if n_rows <= summary_n:
        edges_tsv_path = os.path.join(output_dir, f"{edges_base}.tsv")
        summary_df = _build_summary_df(
            edge_arrays=edge_arrays,
            taxa_universe=taxa_universe,
            null_model=null_model,
            summary_n=summary_n,
        )
        summary_df.to_csv(
            edges_tsv_path,
            sep="\t",
            index=False,
            float_format="%.6g",
        )
        print(f"Taxon edges analysis saved to {edges_tsv_path}")
        return

    parquet_path = os.path.join(output_dir, f"{edges_base}")
    _write_full_edges_parquet_from_arrays(
        edge_arrays=edge_arrays,
        taxa_universe=taxa_universe,
        output_path=parquet_path,
    )
    print(
        f"Full compact taxon edges analysis saved to "
        f"{os.path.splitext(parquet_path)[0]}.parquet "
        f"and {os.path.splitext(parquet_path)[0]}_taxa.parquet"
    )

    summary_df = _build_summary_df(
        edge_arrays=edge_arrays,
        taxa_universe=taxa_universe,
        null_model=null_model,
        summary_n=summary_n,
    )

    summary_path = os.path.join(output_dir, f"{edges_base}_summary.tsv")
    summary_df.to_csv(
        summary_path,
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    print(
        f"Summary taxon edges analysis saved to {summary_path} "
        f"({len(summary_df):,} of {n_rows:,} rows)"
    )


def cooccurrence(
    null_ingredients,
    filtered_ingredients,
    output_dir: str,
    tag: str | None = None,
    filter_rank: Optional[str] = None,
    large: bool = False,
    max_pairs: int = 100_000,
    threshold: float = 0.1,
    null_model: str = "FE",
    nm_n_reps: int = 10,
    nm_random_state: int = 42,
    nm_memory_mb: float | None = None,
):
    """
    Run taxon–taxon co-occurrence analysis and write edge/node tables to disk.
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    null = load_ingredients(data_dir=None, custom_ingredients=null_ingredients)
    filtered = load_ingredients(data_dir=None, custom_ingredients=filtered_ingredients)

    taxa_universe = select_taxa_universe(filtered, rank=filter_rank)

    edge_arrays, nodes_df = cooccurrence_obj(
        null,
        taxa_universe,
        large=large,
        max_pairs=max_pairs,
        threshold=threshold,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_random_state=nm_random_state,
        nm_memory_mb=nm_memory_mb,
    )

    export_cooccurrence_outputs(
        edge_arrays=edge_arrays,
        nodes_df=nodes_df,
        taxa_universe=taxa_universe,
        output_dir=output_dir,
        edges_base=f"{tag}taxon_edges",
        nodes_base=f"{tag}taxon_nodes",
        null_model=null_model,
        summary_n=100_000,
    )



def table_metrics_arrays(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    d: np.ndarray,
    compute_fisher: bool = False,
) -> dict[str, np.ndarray]:
    """
    Vectorised χ², p, φ, directional risk ratios, and optional Fisher's exact test
    from 2×2 tables, returning NumPy arrays instead of a DataFrame.
    """
    chi2, phi, invalid = _chi2_phi_from_counts(a, b, c, d)
    log_p = chi2_dist.logsf(chi2, df=1)
    p = np.exp(log_p)

    P_B_given_A = (a + _SMOOTH) / (a + b + _SMOOTH)
    P_B_given_notA = (c + _SMOOTH) / (c + d + _SMOOTH)
    RR_A_to_B = P_B_given_A / P_B_given_notA
    with np.errstate(divide="ignore", invalid="ignore"):
        logRR_A_to_B = np.log(RR_A_to_B)

    P_A_given_B = (a + _SMOOTH) / (a + c + _SMOOTH)
    P_A_given_notB = (b + _SMOOTH) / (b + d + _SMOOTH)
    RR_B_to_A = P_A_given_B / P_A_given_notB
    with np.errstate(divide="ignore", invalid="ignore"):
        logRR_B_to_A = np.log(RR_B_to_A)

    n_tables = a.shape[0]
    fisher_odds = np.full(n_tables, np.nan, dtype=float)
    fisher_p = np.full(n_tables, np.nan, dtype=float)
    log_fisher_p = np.full(n_tables, np.nan, dtype=float)

    if compute_fisher:
        valid_idx = np.where(~invalid)[0]
        for i in valid_idx:
            odds, pval = _fisher_exact(
                [[a[i], b[i]],
                 [c[i], d[i]]],
                alternative="two-sided",
            )
            fisher_odds[i] = odds
            fisher_p[i] = pval

        with np.errstate(divide="ignore", invalid="ignore"):
            log_fisher_p = np.where(fisher_p > 0, np.log(fisher_p), -np.inf)

    return {
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
    }


def _edge_summary_indices(
    cols: dict[str, np.ndarray],
    null_model: str,
    summary_n: int,
) -> np.ndarray:
    log_null_q = f"log_jaccard_q_{null_model}"
    if log_null_q in cols:
        score = cols[log_null_q]
    else:
        score = cols["log_q_bh"]

    n = len(score)
    if n == 0:
        return np.array([], dtype=np.int64)

    if np.isnan(score).any():
        score_for_sort = np.array(score, copy=True)
        score_for_sort[np.isnan(score_for_sort)] = np.inf
    else:
        score_for_sort = score

    if n <= summary_n:
        return np.arange(n, dtype=np.int64)

    k = int(summary_n)
    idx = np.argpartition(score_for_sort, k - 1)[:k]
    idx = idx[np.argsort(score_for_sort[idx], kind="quicksort")]
    return idx


def _build_summary_df(
    edge_arrays: CooccurrenceArrays,
    taxa_universe: List[str],
    null_model: str,
    summary_n: int,
) -> pd.DataFrame:
    idx = _edge_summary_indices(
        edge_arrays.cols,
        null_model=null_model,
        summary_n=summary_n,
    )
    return _reconstruct_edge_frame(
        edge_arrays=edge_arrays,
        taxa_universe=taxa_universe,
        idx=idx,
    )


def _write_full_edges_parquet_from_arrays(
    edge_arrays: CooccurrenceArrays,
    taxa_universe: List[str],
    output_path: str,
) -> None:
    """
    Write full co-occurrence output in compact form.
    """
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as e:
        raise ImportError(
            "Parquet export requested for a large co-occurrence edge table, but "
            "no parquet engine is installed. Install 'pyarrow'."
        ) from e

    cols = edge_arrays.cols
    totals = np.asarray(edge_arrays.meta["totals"], dtype=np.int32)
    taxa_arr = np.asarray(taxa_universe, dtype=object)

    edges_path = f"{output_path}.parquet"
    taxa_path = f"{output_path}_taxa.parquet"

    edge_table_dict = {
        "iA": pa.array(np.asarray(cols["iA"], dtype=np.int32), type=pa.int32()),
        "iB": pa.array(np.asarray(cols["iB"], dtype=np.int32), type=pa.int32()),
        "A_B_intersection_count": pa.array(np.asarray(cols["inter"], dtype=np.int32), type=pa.int32()),
        "jaccard": pa.array(np.asarray(cols["jaccard"], dtype=np.float32), type=pa.float32()),
        "log_p": pa.array(np.asarray(cols["log_p"], dtype=np.float64), type=pa.float64()),
        "log_q_bh": pa.array(np.asarray(cols["log_q_bh"], dtype=np.float64), type=pa.float64()),
    }

    for key, arr in cols.items():
        if key in {"iA", "iB", "inter", "jaccard", "log_p", "log_q_bh"}:
            continue

        arr = np.asarray(arr)
        if arr.dtype == np.float32:
            edge_table_dict[key] = pa.array(arr, type=pa.float32())
        elif arr.dtype == np.float64:
            edge_table_dict[key] = pa.array(arr, type=pa.float64())
        elif arr.dtype == np.int32:
            edge_table_dict[key] = pa.array(arr, type=pa.int32())
        elif arr.dtype == np.int64:
            edge_table_dict[key] = pa.array(arr, type=pa.int64())
        else:
            edge_table_dict[key] = pa.array(arr)

    edge_table = pa.table(edge_table_dict)
    pq.write_table(edge_table, edges_path, compression="snappy")

    taxa_table = pa.table({
        "taxon_id": pa.array(np.arange(taxa_arr.shape[0], dtype=np.int32), type=pa.int32()),
        "taxon": pa.array(taxa_arr),
        "total_count": pa.array(totals, type=pa.int32()),
    })
    pq.write_table(taxa_table, taxa_path, compression="snappy")


def bh_logq_from_logp(
    log_p: np.ndarray,
    m_total: Optional[int] = None,
    *,
    block_size: int = 5_000_000,
) -> np.ndarray:
    """
    Exact Benjamini-Hochberg adjustment from natural-log p-values.
    """
    log_p = np.asarray(log_p, dtype=np.float64)
    m = log_p.size

    if m == 0:
        return np.array([], dtype=np.float64)

    n = int(m_total) if m_total is not None else m

    nan_mask = np.isnan(log_p)
    if nan_mask.any():
        work_src = log_p.copy()
        work_src[nan_mask] = np.inf
    else:
        work_src = log_p

    order = np.argsort(work_src, kind="quicksort")
    sorted_log_q = work_src[order].copy()
    sorted_log_q += np.log(float(n))

    for start in range(0, m, block_size):
        end = min(start + block_size, m)
        ranks = np.arange(start + 1, end + 1, dtype=np.float64)
        sorted_log_q[start:end] -= np.log(ranks)

    tail_min = np.inf
    for end in range(m, 0, -block_size):
        start = max(0, end - block_size)
        block = sorted_log_q[start:end]

        rev = block[::-1]
        np.minimum.accumulate(rev, out=rev)

        if np.isfinite(tail_min):
            np.minimum(block, tail_min, out=block)

        tail_min = block[0]

    log_q = np.empty_like(sorted_log_q)
    log_q[order] = sorted_log_q

    if nan_mask.any():
        log_q[nan_mask] = np.nan

    return log_q


def _valid_edge_mask_from_counts(
    inter: np.ndarray,
    A_totals: np.ndarray,
    B_totals: np.ndarray,
    N_total: int,
) -> np.ndarray:
    """
    Cheap validity mask from integer counts only.
    """
    d = N_total - A_totals - B_totals + inter

    return (
        (inter >= 0) &
        (A_totals > 0) & (A_totals < N_total) &
        (B_totals > 0) & (B_totals < N_total) &
        (d >= 0)
    )


def _compute_core_edge_metrics(
    inter: np.ndarray,
    A_totals: np.ndarray,
    B_totals: np.ndarray,
    N_total: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute only the core per-edge metrics we actually need to keep in memory.
    """
    a = inter.astype(np.float64, copy=False)
    row1 = A_totals.astype(np.float64, copy=False)
    col1 = B_totals.astype(np.float64, copy=False)
    n = float(N_total)

    cross = a * n - row1 * col1

    denom = row1 * (n - row1)
    denom *= col1
    denom *= (n - col1)

    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = n * (cross * cross) / denom
        log_p = np.log(special.gammaincc(0.5, 0.5 * chi2))

    union = row1 + col1 - a
    jaccard = np.divide(
        a,
        union,
        out=np.zeros(a.shape[0], dtype=np.float32),
        where=union > 0,
    )

    return log_p, jaccard


def _reconstruct_edge_frame(
    edge_arrays: CooccurrenceArrays,
    taxa_universe: List[str],
    idx: np.ndarray | slice,
) -> pd.DataFrame:
    cols = edge_arrays.cols
    totals = np.asarray(edge_arrays.meta["totals"], dtype=np.int32)
    N_total = int(edge_arrays.meta["N_total"])
    compute_fisher = bool(edge_arrays.meta.get("compute_fisher", False))

    iA = cols["iA"][idx]
    iB = cols["iB"][idx]
    inter = cols["inter"][idx]

    A_totals = totals[iA]
    B_totals = totals[iB]

    a = inter.astype(np.float64, copy=False)
    b = (A_totals - inter).astype(np.float64, copy=False)
    c = (B_totals - inter).astype(np.float64, copy=False)
    d = (N_total - A_totals - B_totals + inter).astype(np.float64, copy=False)

    mets = table_metrics_arrays(a, b, c, d, compute_fisher=compute_fisher)

    with np.errstate(divide="ignore", invalid="ignore"):
        p_A = A_totals / float(N_total)
        p_B = B_totals / float(N_total)
        p_A_and_B = inter / float(N_total)

        denom_lift = p_A * p_B
        lift = np.divide(
            p_A_and_B,
            denom_lift,
            out=np.zeros_like(p_A_and_B, dtype=float),
            where=denom_lift != 0,
        )

        P_B_given_A = np.divide(
            inter,
            A_totals,
            out=np.zeros_like(inter, dtype=float),
            where=A_totals > 0,
        )

        P_A_given_B = np.divide(
            inter,
            B_totals,
            out=np.zeros_like(inter, dtype=float),
            where=B_totals > 0,
        )

    taxa_arr = np.asarray(taxa_universe, dtype=object)

    out = {
        "A_taxon": taxa_arr[iA],
        "B_taxon": taxa_arr[iB],
        "A_B_intersection_count": inter,
        "A_total_count": A_totals,
        "B_total_count": B_totals,
        "a": a,
        "b": b,
        "c": c,
        "d": d,
        "N": np.full(inter.shape[0], N_total, dtype=np.int32),
        "p_A": p_A,
        "p_B": p_B,
        "p_A_and_B": p_A_and_B,
        "lift": lift,
        "P_B_given_A": P_B_given_A,
        "P_A_given_B": P_A_given_B,
        "jaccard": cols["jaccard"][idx],
        "chi2": mets["chi2"],
        "p": np.exp(cols["log_p"][idx]),
        "log_p": cols["log_p"][idx],
        "phi": mets["phi"],
        "RR_A_to_B": mets["RR_A_to_B"],
        "RR_B_to_A": mets["RR_B_to_A"],
        "logRR_A_to_B": mets["logRR_A_to_B"],
        "logRR_B_to_A": mets["logRR_B_to_A"],
        "fisher_odds": mets["fisher_odds"],
        "fisher_p": mets["fisher_p"],
        "log_fisher_p": mets["log_fisher_p"],
        "q_bh": np.exp(cols["log_q_bh"][idx]),
        "log_q_bh": cols["log_q_bh"][idx],
    }

    for key, arr in cols.items():
        if key in {"iA", "iB", "inter", "jaccard", "log_p", "log_q_bh"}:
            continue

        if key.startswith("log_jaccard_q_"):
            model = key.removeprefix("log_jaccard_q_")
            out[f"jaccard_q_{model}"] = np.exp(arr[idx])
            out[key] = arr[idx]
        else:
            out[key] = arr[idx]

    for key, value in edge_arrays.meta.items():
        if key in {"totals", "N_total", "compute_fisher"}:
            continue
        out[key] = np.full(len(iA), value)

    return pd.DataFrame(out)
