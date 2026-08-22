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
from metacooc.output import (
    null_metadata_from_reduction,
    with_compact_null_metadata,
)


_SMOOTH = 0.5
_COOCCURRENCE_NULL_EDGE_WARNING_THRESHOLD = 1_000_000
_COOCCURRENCE_PARQUET_BATCH_SIZE = 100_000

ASSOCIATION_BASE_COLUMNS = [
    "taxon",
    "p_cohort_given_taxon",
    "p_taxon_given_cohort",
    "log2_rr_cohort_taxon_vs_not_taxon",
    "rr_cohort_taxon_vs_not_taxon",
    "log2_rr_taxon_cohort_vs_not_cohort",
    "rr_taxon_cohort_vs_not_cohort",
    "delta_p_taxon_cohort_vs_not_cohort",
    "lift_taxon_cohort",
    "jaccard_taxon_cohort",
    "phi_coefficient",
    "rr_cohort_given_taxon_vs_without_taxon",
    "rr_taxon_given_cohort_vs_without_cohort",
    "ln_rr_cohort_given_taxon_vs_without_taxon",
    "ln_rr_taxon_given_cohort_vs_without_cohort",
    "chi2_statistic",
    "chi2_p_value",
    "chi2_q_value_bh",
    "chi2_log_p_value",
    "chi2_log_q_value_bh",
    "taxon_in_cohort_count",
    "taxon_in_background_not_cohort_count",
    "cohort_without_taxon_count",
    "neither_taxon_nor_cohort_count",
    "cohort_sample_count",
    "background_not_cohort_sample_count",
    "background_sample_count",
    "p_taxon_given_not_cohort",
    "p_cohort_given_not_taxon",
]

ASSOCIATION_FISHER_COLUMNS = [
    "fisher_odds_ratio",
    "fisher_p_value",
    "fisher_log_p_value",
]

ASSOCIATION_NULL_COLUMNS = [
    "jaccard_null_mean",
    "jaccard_null_sd",
    "jaccard_null_ses",
    "jaccard_null_p_empirical",
]

ASSOCIATION_SUMMARY_BASE_COLUMNS = [
    "taxon",
    "phi_coefficient",
    "p_cohort_given_taxon",
    "p_taxon_given_cohort",
    "taxon_in_cohort_count",
    "cohort_sample_count",
    "taxon_in_background_not_cohort_count",
    "background_not_cohort_sample_count",
    "chi2_q_value_bh",
]

ASSOCIATION_SUMMARY_NULL_COLUMNS = [
    "jaccard_null_ses",
    "jaccard_null_p_empirical",
]

COOCCURRENCE_BASE_COLUMNS = [
    "source_taxon",
    "target_taxon",
    "p_target_given_source",
    "p_source_given_target",
    "shared_sample_count",
    "source_taxon_sample_count",
    "target_taxon_sample_count",
    "source_only_sample_count",
    "target_only_sample_count",
    "neither_source_nor_target_sample_count",
    "background_sample_count",
    "source_taxon_prevalence",
    "target_taxon_prevalence",
    "cooccurrence_prevalence",
    "lift_taxon_pair",
    "jaccard_taxon_pair",
    "phi_coefficient",
    "chi2_statistic",
    "chi2_p_value",
    "chi2_q_value_bh",
    "chi2_log_p_value",
    "chi2_log_q_value_bh",
    "rr_target_given_source_vs_without_source",
    "rr_source_given_target_vs_without_target",
    "ln_rr_target_given_source_vs_without_source",
    "ln_rr_source_given_target_vs_without_target",
]

COOCCURRENCE_FISHER_COLUMNS = [
    "fisher_odds_ratio",
    "fisher_p_value",
    "fisher_log_p_value",
]

COOCCURRENCE_NULL_COLUMNS = [
    "jaccard_null_mean",
    "jaccard_null_sd",
    "jaccard_null_ses",
    "jaccard_null_p_empirical",
    "jaccard_null_q_value_bh",
    "jaccard_null_log_q_value_bh",
]

COOCCURRENCE_SUMMARY_BASE_COLUMNS = [
    "source_taxon",
    "target_taxon",
    "phi_coefficient",
    "p_target_given_source",
    "p_source_given_target",
    "shared_sample_count",
    "source_taxon_sample_count",
    "target_taxon_sample_count",
    "chi2_q_value_bh",
]

COOCCURRENCE_SUMMARY_NULL_COLUMNS = [
    "jaccard_null_ses",
    "jaccard_null_p_empirical",
    "jaccard_null_q_value_bh",
]


def _cooccurrence_degree_column(min_conditional_probability: float) -> str:
    return f"out_degree_p_target_given_source_gt_{min_conditional_probability}"


def _association_output_columns(
    *,
    compute_fisher: bool,
    include_null: bool,
) -> list[str]:
    columns = list(ASSOCIATION_BASE_COLUMNS)
    if compute_fisher:
        columns.extend(ASSOCIATION_FISHER_COLUMNS)
    if include_null:
        columns.extend(ASSOCIATION_NULL_COLUMNS)
    return columns


def _finalize_association_output(
    df: pd.DataFrame,
    *,
    compute_fisher: bool,
    include_null: bool,
) -> pd.DataFrame:
    columns = _association_output_columns(
        compute_fisher=compute_fisher,
        include_null=include_null,
    )
    attrs = dict(df.attrs)
    out = df.loc[:, columns].copy()
    out.attrs.update(attrs)
    return out


def _build_association_summary(df: pd.DataFrame) -> pd.DataFrame:
    """Build the interpretation-first association result table."""
    columns = list(ASSOCIATION_SUMMARY_BASE_COLUMNS)
    columns.extend(
        column for column in ASSOCIATION_SUMMARY_NULL_COLUMNS if column in df.columns
    )

    attrs = dict(df.attrs)
    out = df.loc[:, columns].copy()

    count_columns = [
        "taxon_in_cohort_count",
        "cohort_sample_count",
        "taxon_in_background_not_cohort_count",
        "background_not_cohort_sample_count",
    ]
    for column in count_columns:
        out[column] = out[column].astype(np.int64)

    out["_summary_abs_phi"] = out["phi_coefficient"].abs()
    out.sort_values(
        by=[
            "chi2_q_value_bh",
            "_summary_abs_phi",
            "taxon_in_cohort_count",
            "taxon",
        ],
        ascending=[True, False, False, True],
        na_position="last",
        kind="mergesort",
        inplace=True,
    )
    out.drop(columns="_summary_abs_phi", inplace=True)
    out.reset_index(drop=True, inplace=True)
    out.attrs.update(attrs)
    return out


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
    min_conditional_probability: float = 0.0,
    null_model: str = "FE",
    nm_n_reps: int = 1000,
    nm_seed: int | None = None,
    compute_fisher: bool = False,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_burn_in_steps: int | None = None,
    nm_steps_per_rep: int | None = None,
    nm_progress_every: int = 25,
) -> pd.DataFrame:
    """Calculate taxon-to-cohort association metrics in memory.

    ``filtered_ingredients`` supplies the cohort and ``null_ingredients``
    supplies the comparison background. The returned table contains detailed
    metrics; use :func:`export_association_outputs` when summary and detailed
    files are required.
    """
    out = _association_core(
        null_ingredients=null_ingredients,
        filtered_ingredients=filtered_ingredients,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_seed=nm_seed,
        compute_fisher=compute_fisher,
        nm_n_workers=nm_n_workers,
        nm_mp_start=nm_mp_start,
        nm_burn_in_steps=nm_burn_in_steps,
        nm_steps_per_rep=nm_steps_per_rep,
        nm_progress_every=nm_progress_every,
    )

    if min_conditional_probability is not None:
        attrs = dict(out.attrs)
        out = out[out["p_cohort_given_taxon"] > min_conditional_probability].copy()
        out.attrs.update(attrs)

    return out


def export_association_outputs(
    association_df: pd.DataFrame,
    output_path: str,
    *,
    null_model: str = "FE",
) -> None:
    """Write detailed and reduced association results for one analysis run."""
    metadata = association_df.attrs.get("null_metadata")
    fallback_null_model = "FE" if str(null_model).upper() == "FE" else None

    detailed_df = with_compact_null_metadata(
        association_df,
        metadata,
        fallback_null_model=fallback_null_model,
    )
    detailed_df.to_csv(output_path, sep="\t", index=False)
    print(f"Detailed association analysis saved to {output_path}")

    summary_path = f"{os.path.splitext(output_path)[0]}_summary.tsv"
    summary_df = with_compact_null_metadata(
        _build_association_summary(association_df),
        metadata,
        fallback_null_model=fallback_null_model,
    )
    summary_df.to_csv(
        summary_path,
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    print(f"Association summary saved to {summary_path}")


def _association_core(
    null_ingredients: "Ingredients",
    filtered_ingredients: "Ingredients",
    null_model: str = "FE",
    nm_n_reps: int = 1000,
    nm_seed: int | None = None,
    compute_fisher: bool = False,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_burn_in_steps: int | None = None,
    nm_steps_per_rep: int | None = None,
    nm_progress_every: int = 25,
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

    cohort_sample_count = float(len(filtered_ingredients.samples))
    background_sample_count = float(len(null_ingredients.samples))
    if cohort_sample_count <= 0 or background_sample_count <= 0 or background_sample_count < cohort_sample_count:
        raise ValueError(
            "Invalid cohort sizes: "
            f"cohort_sample_count={cohort_sample_count}, "
            f"background_sample_count={background_sample_count}"
        )
    if background_sample_count == cohort_sample_count:
        raise ValueError(
            "Null and term cohorts are identical: no non-term samples available."
        )

    taxon_in_cohort_count = filtered_ingredients.total_counts[filt_idx].astype(float, copy=False)
    taxon_in_background_total_count = null_ingredients.total_counts[null_idx].astype(float, copy=False)
    taxon_in_background_not_cohort_count = taxon_in_background_total_count - taxon_in_cohort_count

    background_not_cohort_sample_count = background_sample_count - cohort_sample_count
    cohort_without_taxon_count = cohort_sample_count - taxon_in_cohort_count
    neither_taxon_nor_cohort_count = (
        background_not_cohort_sample_count - taxon_in_background_not_cohort_count
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        jaccard_taxon_cohort = np.divide(
            taxon_in_cohort_count,
            taxon_in_background_not_cohort_count + cohort_sample_count,
            out=np.zeros_like(taxon_in_cohort_count, dtype=float),
            where=(taxon_in_background_not_cohort_count + cohort_sample_count) > 0,
        )

    mets = table_metrics(
        taxon_in_cohort_count,
        taxon_in_background_not_cohort_count,
        cohort_without_taxon_count,
        neither_taxon_nor_cohort_count,
        compute_fisher=compute_fisher,
    )

    taxon_in_cohort_smoothed = taxon_in_cohort_count + _SMOOTH
    taxon_in_background_not_cohort_smoothed = taxon_in_background_not_cohort_count + _SMOOTH
    cohort_without_taxon_smoothed = cohort_without_taxon_count + _SMOOTH
    neither_taxon_nor_cohort_smoothed = neither_taxon_nor_cohort_count + _SMOOTH

    with np.errstate(divide="ignore", invalid="ignore"):
        p_taxon_given_cohort = taxon_in_cohort_smoothed / (
            taxon_in_cohort_smoothed + cohort_without_taxon_smoothed
        )
        p_taxon_given_not_cohort = taxon_in_background_not_cohort_smoothed / (
            taxon_in_background_not_cohort_smoothed + neither_taxon_nor_cohort_smoothed
        )
        rr_taxon_cohort_vs_not_cohort = np.divide(
            p_taxon_given_cohort,
            p_taxon_given_not_cohort,
            out=np.zeros_like(p_taxon_given_cohort),
            where=p_taxon_given_not_cohort != 0,
        )
        log2_rr_taxon_cohort_vs_not_cohort = np.where(
            rr_taxon_cohort_vs_not_cohort > 0,
            np.log2(rr_taxon_cohort_vs_not_cohort),
            -np.inf,
        )
        delta_p_taxon_cohort_vs_not_cohort = (
            p_taxon_given_cohort - p_taxon_given_not_cohort
        )

        p_cohort_given_taxon = taxon_in_cohort_smoothed / (
            taxon_in_cohort_smoothed + taxon_in_background_not_cohort_smoothed
        )
        p_cohort_given_not_taxon = cohort_without_taxon_smoothed / (
            cohort_without_taxon_smoothed + neither_taxon_nor_cohort_smoothed
        )
        rr_cohort_taxon_vs_not_taxon = np.divide(
            p_cohort_given_taxon,
            p_cohort_given_not_taxon,
            out=np.zeros_like(p_cohort_given_taxon),
            where=p_cohort_given_not_taxon != 0,
        )
        log2_rr_cohort_taxon_vs_not_taxon = np.where(
            rr_cohort_taxon_vs_not_taxon > 0,
            np.log2(rr_cohort_taxon_vs_not_taxon),
            -np.inf,
        )

        p_taxon = (
            taxon_in_cohort_smoothed + taxon_in_background_not_cohort_smoothed
        ) / background_sample_count
        p_cohort = (
            taxon_in_cohort_smoothed + cohort_without_taxon_smoothed
        ) / background_sample_count
        p_taxon_and_cohort = taxon_in_cohort_smoothed / background_sample_count
        lift_taxon_cohort = np.divide(
            p_taxon_and_cohort,
            p_taxon * p_cohort,
            out=np.zeros_like(p_taxon_and_cohort),
            where=(p_taxon * p_cohort) != 0,
        )

    out = pd.DataFrame(
        {
            "taxon": taxa,
            "p_cohort_given_taxon": p_cohort_given_taxon,
            "p_taxon_given_cohort": p_taxon_given_cohort,
            "log2_rr_cohort_taxon_vs_not_taxon": log2_rr_cohort_taxon_vs_not_taxon,
            "rr_cohort_taxon_vs_not_taxon": rr_cohort_taxon_vs_not_taxon,
            "log2_rr_taxon_cohort_vs_not_cohort": log2_rr_taxon_cohort_vs_not_cohort,
            "rr_taxon_cohort_vs_not_cohort": rr_taxon_cohort_vs_not_cohort,
            "delta_p_taxon_cohort_vs_not_cohort": delta_p_taxon_cohort_vs_not_cohort,
            "lift_taxon_cohort": lift_taxon_cohort,
            "jaccard_taxon_cohort": jaccard_taxon_cohort,
            "taxon_in_cohort_count": taxon_in_cohort_count,
            "taxon_in_background_not_cohort_count": taxon_in_background_not_cohort_count,
            "cohort_without_taxon_count": cohort_without_taxon_count,
            "neither_taxon_nor_cohort_count": neither_taxon_nor_cohort_count,
            "cohort_sample_count": cohort_sample_count,
            "background_not_cohort_sample_count": background_not_cohort_sample_count,
            "background_sample_count": background_sample_count,
            "p_taxon_given_not_cohort": p_taxon_given_not_cohort,
            "p_cohort_given_not_taxon": p_cohort_given_not_taxon,
        }
    )

    metrics = pd.DataFrame(
        {
            "chi2_statistic": mets["chi2"],
            "chi2_p_value": mets["p"],
            "chi2_log_p_value": mets["log_p"],
            "phi_coefficient": mets["phi"],
            "rr_cohort_given_taxon_vs_without_taxon": mets["RR_A_to_B"],
            "rr_taxon_given_cohort_vs_without_cohort": mets["RR_B_to_A"],
            "ln_rr_cohort_given_taxon_vs_without_taxon": mets["logRR_A_to_B"],
            "ln_rr_taxon_given_cohort_vs_without_cohort": mets["logRR_B_to_A"],
            "invalid_table": mets["invalid_table"],
        }
    )
    if compute_fisher:
        metrics = metrics.join(
            mets[["fisher_odds", "fisher_p", "log_fisher_p"]].rename(
                columns={
                    "fisher_odds": "fisher_odds_ratio",
                    "fisher_p": "fisher_p_value",
                    "log_fisher_p": "fisher_log_p_value",
                }
            )
        )

    out = out.join(metrics)

    out["chi2_log_q_value_bh"] = bh_logq_from_logp(
        out["chi2_log_p_value"].to_numpy(dtype=np.float64, copy=False)
    )
    out["chi2_q_value_bh"] = np.exp(out["chi2_log_q_value_bh"].to_numpy(dtype=np.float64, copy=False))

    valid_mask = ~out["invalid_table"].values
    out = out.loc[valid_mask].copy()
    out.drop(columns="invalid_table", inplace=True)

    background_taxon_indices_valid = null_idx[valid_mask]
    out.attrs["background_taxon_indices_valid"] = np.asarray(background_taxon_indices_valid, dtype=int)
    out.attrs["taxa_valid"] = out["taxon"].to_numpy(dtype=object, copy=False)

    n_reps = int(nm_n_reps) if nm_n_reps is not None else 0
    if n_reps <= 0 or out.shape[0] == 0:
        return _finalize_association_output(
            out,
            compute_fisher=compute_fisher,
            include_null=False,
        )

    if null_model == "FE":
        print("FE: association determined analytically - no need for shuffling null and probabilistic approach")
        return _finalize_association_output(
            out,
            compute_fisher=compute_fisher,
            include_null=False,
        )

    X_full = presence_for_counts(null_ingredients).tocsr()
    X_full.eliminate_zeros()
    X_full.sum_duplicates()
    X_full.sort_indices()

    n_rows, n_cols = X_full.shape

    sample_index = {s: i for i, s in enumerate(null_ingredients.samples)}
    term_cols = np.array([sample_index[s] for s in filtered_ingredients.samples], dtype=np.int64)

    mask_cols = np.zeros(n_cols, dtype=bool)
    mask_cols[term_cols] = True
    nonterm_cols = np.where(~mask_cols)[0].astype(np.int64, copy=False)

    subset_idx = np.asarray(background_taxon_indices_valid, dtype=np.int64)

    observed_jaccard = out["jaccard_taxon_cohort"].to_numpy(dtype=float, copy=False)

    mp_start = _best_mp_start() if nm_mp_start is None else str(nm_mp_start)

    null_reduction = parallel_null_reduce_vector(
        X=X_full,
        model=null_model,
        n_reps=n_reps,
        obs=observed_jaccard,
        stat_fn=stat_fn_association_jaccard,
        seed=nm_seed,
        n_workers=nm_n_workers,
        mp_start=mp_start,
        burn_in_steps=nm_burn_in_steps,
        steps_per_rep=nm_steps_per_rep,
        progress_every=int(nm_progress_every),
        term_cols=term_cols,
        nonterm_cols=nonterm_cols,
        subset_idx=subset_idx,
        n_rows=n_rows,
        N_T=cohort_sample_count,
    )

    out["jaccard_null_mean"] = null_reduction["mean"]
    out["jaccard_null_sd"] = null_reduction["sd"]
    out["jaccard_null_ses"] = null_reduction["ses"]
    out["jaccard_null_p_empirical"] = null_reduction["p_emp"]
    out.attrs["null_metadata"] = null_metadata_from_reduction(null_reduction)
    print(
        f"INFO: Null simulation seed {null_reduction['null_seed']} "
        f"({null_reduction['null_seed_source']}, model={null_reduction['null_model']})."
    )

    return _finalize_association_output(
        out,
        compute_fisher=compute_fisher,
        include_null=True,
    )


def association(
    null_ingredients: "Ingredients",
    filtered_ingredients: "Ingredients",
    output_dir: str,
    tag: str | None = None,
    min_conditional_probability: float = 0.0,
    null_model: str = "FE",
    nm_n_reps: int = 1000,
    nm_seed: int | None = None,
    compute_fisher: bool = False,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_burn_in_steps: int | None = None,
    nm_steps_per_rep: int | None = None,
    nm_progress_every: int = 25,
) -> None:
    """Run association analysis and write detailed and summary TSV files.

    ``null_ingredients`` and ``filtered_ingredients`` may be Ingredients
    objects or values accepted by :func:`metacooc.pantry.load_ingredients`.
    The detailed and summary result files are written below ``output_dir``.
    The file-writing wrapper currently returns ``None``.
    """
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    null = load_ingredients(data_dir=None, custom_ingredients=null_ingredients)
    filtered = load_ingredients(data_dir=None, custom_ingredients=filtered_ingredients)

    association_df = association_obj(
        null_ingredients=null,
        filtered_ingredients=filtered,
        min_conditional_probability=min_conditional_probability,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_seed=nm_seed,
        compute_fisher=compute_fisher,
        nm_n_workers=nm_n_workers,
        nm_mp_start=nm_mp_start,
        nm_burn_in_steps=nm_burn_in_steps,
        nm_steps_per_rep=nm_steps_per_rep,
        nm_progress_every=nm_progress_every,
    )

    if association_df is not None:
        output_path = os.path.join(output_dir, f"{tag}association.tsv")
        export_association_outputs(
            association_df,
            output_path,
            null_model=null_model,
        )


def presence_submatrix_by_taxa(ingredients: Ingredients, taxa_subset: List[str]) -> sp.csr_matrix:
    row_map = {t: i for i, t in enumerate(ingredients.taxa)}
    rows = [row_map[t] for t in taxa_subset]
    return presence_for_counts(ingredients)[rows, :].tocsr()


def _cooccur_core(
    ing: Ingredients,
    taxa_universe: List[str],
    min_conditional_probability: float = 0.1,
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

        if min_conditional_probability > 0:
            keep = inter > (min_conditional_probability * A_totals)
            if not np.any(keep):
                continue
            iA = iA[keep]
            iB = iB[keep]
            inter = inter[keep]
            A_totals = A_totals[keep]

        B_totals = totals[iB]
        valid = _valid_edge_mask_from_counts(
            shared_sample_count=inter,
            source_taxon_sample_count=A_totals,
            target_taxon_sample_count=B_totals,
            N_total=N_total,
        )
        if not np.any(valid):
            continue

        iA_valid = iA[valid]
        n_keep += int(iA_valid.size)
        deg_fwd += np.bincount(iA_valid, minlength=n_taxa)

    nodes_df = pd.DataFrame({
        "taxon": taxa_universe,
        "taxon_sample_count": totals.astype(int, copy=False),
        _cooccurrence_degree_column(min_conditional_probability): deg_fwd.astype(int, copy=False),
    })

    if n_keep == 0:
        return None, nodes_df

    source_taxon_index_all = np.empty(n_keep, dtype=np.int32)
    target_taxon_index_all = np.empty(n_keep, dtype=np.int32)
    shared_sample_count_all = np.empty(n_keep, dtype=np.int32)
    chi2_log_p_value_all = np.empty(n_keep, dtype=np.float64)
    jaccard_taxon_pair_all = np.empty(n_keep, dtype=np.float32)

    pos = 0
    for iA, iB, inter in stream_edge_values(co_mat, min_value=0):
        A_totals = totals[iA]

        if min_conditional_probability > 0:
            keep = inter > (min_conditional_probability * A_totals)
            if not np.any(keep):
                continue
            iA = iA[keep]
            iB = iB[keep]
            inter = inter[keep]
            A_totals = A_totals[keep]

        B_totals = totals[iB]
        valid = _valid_edge_mask_from_counts(
            shared_sample_count=inter,
            source_taxon_sample_count=A_totals,
            target_taxon_sample_count=B_totals,
            N_total=N_total,
        )
        if not np.any(valid):
            continue

        iA = iA[valid].astype(np.int32, copy=False)
        iB = iB[valid].astype(np.int32, copy=False)
        inter = inter[valid].astype(np.int32, copy=False)
        A_totals = A_totals[valid]
        B_totals = B_totals[valid]

        chi2_log_p_value_chunk, jaccard_taxon_pair_chunk = _compute_core_edge_metrics(
            shared_sample_count=inter,
            source_taxon_sample_count=A_totals,
            target_taxon_sample_count=B_totals,
            N_total=N_total,
        )

        k = int(inter.size)
        sl = slice(pos, pos + k)

        source_taxon_index_all[sl] = iA
        target_taxon_index_all[sl] = iB
        shared_sample_count_all[sl] = inter
        chi2_log_p_value_all[sl] = chi2_log_p_value_chunk
        jaccard_taxon_pair_all[sl] = jaccard_taxon_pair_chunk

        pos += k

    edge_arrays = CooccurrenceArrays(
        cols={
            "source_taxon_index": source_taxon_index_all,
            "target_taxon_index": target_taxon_index_all,
            "shared_sample_count": shared_sample_count_all,
            "jaccard_taxon_pair": jaccard_taxon_pair_all,
            "chi2_log_p_value": chi2_log_p_value_all,
        },
        meta={
            "totals": totals,
            "N_total": N_total,
            "compute_fisher": bool(compute_fisher),
        },
    )

    edge_arrays.cols["chi2_log_q_value_bh"] = bh_logq_from_logp(
        chi2_log_p_value_all,
        m_total=m_total,
    )

    return edge_arrays, nodes_df


def _cooccur_core_focal(
    ing: Ingredients,
    taxa_universe: List[str],
    focal_local_idx: np.ndarray,
    min_conditional_probability: float = 0.1,
    m_total: Optional[int] = None,
    compute_fisher: bool = False,
) -> Tuple[Optional[CooccurrenceArrays], pd.DataFrame]:
    """
    Focal-only co-occurrence core.

    Computes only focal x all intersections and emits edges where source_taxon is focal.
    Focal-focal pairs are emitted once with orientation normalised to iA < iB.
    """
    X_sub = presence_submatrix_by_taxa(ing, taxa_universe).tocsr()
    totals = np.asarray(X_sub.sum(axis=1)).ravel().astype(np.int32, copy=False)
    N_total = int(len(ing.samples))

    focal_local_idx = np.asarray(focal_local_idx, dtype=np.int64)
    focal_set = set(focal_local_idx.tolist())

    nodes_df = pd.DataFrame({
        "taxon": taxa_universe,
        "taxon_sample_count": totals.astype(int, copy=False),
        _cooccurrence_degree_column(min_conditional_probability): np.zeros(len(taxa_universe), dtype=int),
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

        if min_conditional_probability > 0:
            keep = inter > (min_conditional_probability * A_totals)
            if not np.any(keep):
                continue
            iA = iA[keep]
            j_cols = j_cols[keep]
            inter = inter[keep]
            A_totals = A_totals[keep]

        B_totals = totals[j_cols]
        valid = _valid_edge_mask_from_counts(
            shared_sample_count=inter,
            source_taxon_sample_count=A_totals,
            target_taxon_sample_count=B_totals,
            N_total=N_total,
        )
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

    nodes_df[_cooccurrence_degree_column(min_conditional_probability)] = deg_fwd.astype(int, copy=False)

    if n_keep == 0:
        return None, nodes_df

    source_taxon_index_all = np.empty(n_keep, dtype=np.int32)
    target_taxon_index_all = np.empty(n_keep, dtype=np.int32)
    shared_sample_count_all = np.empty(n_keep, dtype=np.int32)
    chi2_log_p_value_all = np.empty(n_keep, dtype=np.float64)
    jaccard_taxon_pair_all = np.empty(n_keep, dtype=np.float32)

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

        if min_conditional_probability > 0:
            keep = inter > (min_conditional_probability * A_totals)
            if not np.any(keep):
                continue
            iA = iA[keep]
            j_cols = j_cols[keep]
            inter = inter[keep]
            A_totals = A_totals[keep]

        B_totals = totals[j_cols]
        valid = _valid_edge_mask_from_counts(
            shared_sample_count=inter,
            source_taxon_sample_count=A_totals,
            target_taxon_sample_count=B_totals,
            N_total=N_total,
        )
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

        chi2_log_p_value_chunk, jaccard_taxon_pair_chunk = _compute_core_edge_metrics(
            shared_sample_count=inter,
            source_taxon_sample_count=A_totals,
            target_taxon_sample_count=B_totals,
            N_total=N_total,
        )

        k = int(inter.size)
        sl = slice(pos, pos + k)

        source_taxon_index_all[sl] = iA.astype(np.int32, copy=False)
        target_taxon_index_all[sl] = iB.astype(np.int32, copy=False)
        shared_sample_count_all[sl] = inter
        chi2_log_p_value_all[sl] = chi2_log_p_value_chunk
        jaccard_taxon_pair_all[sl] = jaccard_taxon_pair_chunk
        pos += k

    edge_arrays = CooccurrenceArrays(
        cols={
            "source_taxon_index": source_taxon_index_all,
            "target_taxon_index": target_taxon_index_all,
            "shared_sample_count": shared_sample_count_all,
            "jaccard_taxon_pair": jaccard_taxon_pair_all,
            "chi2_log_p_value": chi2_log_p_value_all,
        },
        meta={
            "totals": totals,
            "N_total": N_total,
            "compute_fisher": bool(compute_fisher),
        },
    )

    edge_arrays.cols["chi2_log_q_value_bh"] = bh_logq_from_logp(
        chi2_log_p_value_all,
        m_total=m_total,
    )
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
    focal taxon is always source_taxon.

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

    source_indices = np.asarray(edge_arrays.cols["source_taxon_index"], dtype=np.int64)
    target_indices = np.asarray(edge_arrays.cols["target_taxon_index"], dtype=np.int64)

    collected = {k: [] for k in edge_arrays.cols}
    focal_query_parts = []
    focal_taxon_parts = []

    for focal_query, focal_taxa in focal_query_to_taxa.items():
        for focal_taxon in focal_taxa:
            focal_idx = tax_to_idx.get(focal_taxon)
            if focal_idx is None:
                continue

            left_idx = np.flatnonzero(source_indices == focal_idx)
            if left_idx.size:
                for key, arr in edge_arrays.cols.items():
                    collected[key].append(np.asarray(arr)[left_idx])
                focal_query_parts.append(np.full(left_idx.size, focal_query, dtype=object))
                focal_taxon_parts.append(np.full(left_idx.size, focal_taxon, dtype=object))

            right_idx = np.flatnonzero(target_indices == focal_idx)
            if right_idx.size:
                for key, arr in edge_arrays.cols.items():
                    arr_np = np.asarray(arr)
                    if key == "source_taxon_index":
                        collected[key].append(target_indices[right_idx].astype(arr_np.dtype, copy=False))
                    elif key == "target_taxon_index":
                        collected[key].append(source_indices[right_idx].astype(arr_np.dtype, copy=False))
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
      - source_taxon_index
      - target_taxon_index
      - focal_query

    Semantics
    ---------
    Keep an edge iff:
      - its focal_query exists in rhs_query_to_taxa, and
      - its target_taxon name belongs to the RHS taxa mapped to that focal_query

    Notes
    -----
    - This assumes `_subset_focal_edge_arrays()` has already run, so
      focal_query is present and source_taxon_index is already oriented as the
      focal side.
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
    target_index_arr = np.asarray(edge_arrays.cols["target_taxon_index"], dtype=np.int64)
    target_taxon_arr = taxa_arr[target_index_arr]

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

        keep[q_mask] = np.isin(target_taxon_arr[q_mask], list(rhs_taxa_set))

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
    min_conditional_probability: float = 0.1,
    null_model: str = "FE",
    nm_n_reps: int = 10,
    nm_seed: int | None = None,
    compute_fisher: bool = False,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_burn_in_steps: int | None = None,
    nm_steps_per_rep: int | None = None,
    nm_progress_every: int = 25,
) -> Tuple[Optional[CooccurrenceArrays], Optional[pd.DataFrame]]:
    """
    Pairwise co-occurrence of taxa.

    ``null_ingredients`` defines the sample/background matrix and
    ``taxa_universe`` defines the taxa considered. Focal and RHS mappings are
    optional restrictions for focal-taxa workflows.

    If focal_query_to_taxa is provided, computation is restricted to focal-anchored
    pairs only, using the union of resolved focal taxa across queries.

    If rhs_query_to_taxa is also provided, the already focal-expanded edge set is
    further restricted so that target_taxon belongs to the RHS taxa for the matching
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
                "  • Increase the --min_conditional_probability value [0, 1],\n"
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
            min_conditional_probability=min_conditional_probability,
            m_total=est_pairs,
            compute_fisher=compute_fisher,
        )
    else:
        edge_arrays, nodes_df = _cooccur_core_focal(
            null_ingredients,
            taxa_universe,
            focal_local_idx=focal_local_idx,
            min_conditional_probability=min_conditional_probability,
            m_total=est_pairs,
            compute_fisher=compute_fisher,
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
                    "target-side taxa survived under the current focal cohort, taxa universe, "
                    "and min_conditional_probability settings."
                )
                return None, nodes_df

        edge_arrays.cols["chi2_log_q_value_bh"] = bh_logq_from_logp(
            edge_arrays.cols["chi2_log_p_value"],
            m_total=edge_arrays.n_rows,
        )

    n_reps = int(nm_n_reps) if nm_n_reps is not None else 0
    if n_reps <= 0:
        return edge_arrays, nodes_df

    if null_model == "FE":
        print("FE: cooccurrence determined analytically - no need for shuffling null and probabilistic approach")
        return edge_arrays, nodes_df

    X_full = presence_for_counts(null_ingredients).tocsr()
    X_full.eliminate_zeros()
    X_full.sum_duplicates()
    X_full.sort_indices()

    tax_map = {t: i for i, t in enumerate(null_ingredients.taxa)}
    subset_idx = np.array([tax_map[t] for t in taxa_universe], dtype=np.int64)

    observed_jaccard = edge_arrays.cols["jaccard_taxon_pair"].astype(np.float64, copy=False)
    source_taxon_indices = edge_arrays.cols["source_taxon_index"].astype(np.int64, copy=False)
    target_taxon_indices = edge_arrays.cols["target_taxon_index"].astype(np.int64, copy=False)

    mp_start = _best_mp_start() if nm_mp_start is None else str(nm_mp_start)

    if edge_arrays.n_rows >= _COOCCURRENCE_NULL_EDGE_WARNING_THRESHOLD:
        print(
            "WARNING: Cooccurrence null simulation will evaluate "
            f"{edge_arrays.n_rows:,} edge pairs per replicate across {n_reps:,} replicates. "
            f"Estimated candidate pairs before filtering: {est_pairs:,}. "
            "The statistic will use bounded dot-product batching to limit RAM, "
            "but runtime and output size may be substantial."
        )

    j_res = parallel_null_reduce_vector(
        X=X_full,
        model=null_model,
        n_reps=n_reps,
        obs=observed_jaccard,
        stat_fn=stat_fn_cooccurrence_jaccard,
        seed=nm_seed,
        n_workers=nm_n_workers,
        mp_start=mp_start,
        burn_in_steps=nm_burn_in_steps,
        steps_per_rep=nm_steps_per_rep,
        progress_every=int(nm_progress_every),
        subset_idx=subset_idx,
        iA=source_taxon_indices,
        iB=target_taxon_indices,
    )

    edge_arrays.cols["jaccard_null_mean"] = np.asarray(j_res["mean"], dtype=np.float32)
    edge_arrays.cols["jaccard_null_sd"] = np.asarray(j_res["sd"], dtype=np.float32)
    edge_arrays.cols["jaccard_null_ses"] = np.asarray(j_res["ses"], dtype=np.float32)
    edge_arrays.cols["jaccard_null_p_empirical"] = np.asarray(j_res["p_emp"], dtype=np.float64)

    edge_arrays.cols["jaccard_null_log_q_value_bh"] = bh_logq_from_logp(
        np.log(edge_arrays.cols["jaccard_null_p_empirical"]),
        m_total=edge_arrays.n_rows if focal_query_to_taxa else est_pairs,
    )

    null_metadata = null_metadata_from_reduction(j_res)
    edge_arrays.meta["null_metadata"] = null_metadata
    edge_arrays.meta.update(null_metadata)
    print(
        f"INFO: Null simulation seed {j_res['null_seed']} "
        f"({j_res['null_seed_source']}, model={j_res['null_model']})."
    )

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
) -> str | None:
    """
    Export co-occurrence outputs and return the detailed edge result path.
    """
    os.makedirs(output_dir, exist_ok=True)

    if nodes_df is not None:
        nodes_path = os.path.join(output_dir, f"{nodes_base}.tsv")
        nodes_df.to_csv(nodes_path, sep="\t", index=False)
        print(f"Taxon nodes analysis saved to {nodes_path}")

    if edge_arrays is None or edge_arrays.n_rows == 0:
        print("No co-occurrence edges above min_conditional_probability; no edge file written.")
        return None

    n_rows = edge_arrays.n_rows
    metadata = edge_arrays.meta.get("null_metadata")
    fallback_null_model = "FE" if str(null_model).upper() == "FE" else None

    if n_rows <= summary_n:
        edges_tsv_path = os.path.join(output_dir, f"{edges_base}.tsv")
        detailed_df = _reconstruct_edge_frame(
            edge_arrays=edge_arrays,
            taxa_universe=taxa_universe,
            idx=slice(None),
        )
        summary_df = _build_cooccurrence_summary(
            edge_arrays=edge_arrays,
            taxa_universe=taxa_universe,
            summary_n=summary_n,
            detailed_df=detailed_df,
        )
        detailed_df = with_compact_null_metadata(
            detailed_df,
            metadata,
            fallback_null_model=fallback_null_model,
        )
        detailed_df.to_csv(
            edges_tsv_path,
            sep="\t",
            index=False,
            float_format="%.6g",
        )
        print(f"Detailed taxon edges analysis saved to {edges_tsv_path}")

        summary_df = with_compact_null_metadata(
            summary_df,
            metadata,
            fallback_null_model=fallback_null_model,
        )
        summary_path = os.path.join(output_dir, f"{edges_base}_summary.tsv")
        summary_df.to_csv(
            summary_path,
            sep="\t",
            index=False,
            float_format="%.6g",
        )
        print(f"Taxon edges summary saved to {summary_path}")
        return edges_tsv_path

    parquet_path = os.path.join(output_dir, f"{edges_base}")
    _write_full_edges_parquet_from_arrays(
        edge_arrays=edge_arrays,
        taxa_universe=taxa_universe,
        output_path=parquet_path,
        null_model=null_model,
    )
    print(
        f"Full detailed taxon edges analysis saved to "
        f"{os.path.splitext(parquet_path)[0]}.parquet "
        f"and {os.path.splitext(parquet_path)[0]}_taxa.parquet"
    )
    summary_df = _build_cooccurrence_summary(
        edge_arrays=edge_arrays,
        taxa_universe=taxa_universe,
        summary_n=summary_n,
    )
    summary_df = with_compact_null_metadata(
        summary_df,
        metadata,
        fallback_null_model=fallback_null_model,
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
    return f"{parquet_path}.parquet"


def cooccurrence(
    null_ingredients,
    filtered_ingredients,
    output_dir: str,
    tag: str | None = None,
    filter_rank: Optional[str] = None,
    large: bool = False,
    max_pairs: int = 100_000,
    min_conditional_probability: float = 0.1,
    null_model: str = "FE",
    nm_n_reps: int = 10,
    nm_seed: int | None = None,
    compute_fisher: bool = False,
    *,
    nm_n_workers: int | None = None,
    nm_mp_start: str | None = None,
    nm_burn_in_steps: int | None = None,
    nm_steps_per_rep: int | None = None,
    nm_progress_every: int = 25,
):
    """
    Run taxon–taxon co-occurrence analysis and write edge/node tables to disk.

    Inputs may be Ingredients objects or values accepted by
    :func:`metacooc.pantry.load_ingredients`. The output directory receives a
    reduced edge summary, detailed edges, and node results; very large edge
    tables use Parquet output. The file-writing wrapper returns ``None``.
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
        min_conditional_probability=min_conditional_probability,
        null_model=null_model,
        nm_n_reps=nm_n_reps,
        nm_seed=nm_seed,
        compute_fisher=compute_fisher,
        nm_n_workers=nm_n_workers,
        nm_mp_start=nm_mp_start,
        nm_burn_in_steps=nm_burn_in_steps,
        nm_steps_per_rep=nm_steps_per_rep,
        nm_progress_every=nm_progress_every,
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
    edge_arrays: CooccurrenceArrays,
    taxa_universe: List[str],
    summary_n: int,
) -> np.ndarray:
    cols = edge_arrays.cols
    score_column = _cooccurrence_summary_log_q_column(edge_arrays)
    score = np.asarray(cols[score_column], dtype=np.float64)
    n = len(score)
    if n == 0:
        return np.array([], dtype=np.int64)

    if np.isnan(score).any():
        score_for_sort = np.array(score, copy=True)
        score_for_sort[np.isnan(score_for_sort)] = np.inf
    else:
        score_for_sort = score

    k = min(n, int(summary_n))
    if k <= 0:
        return np.array([], dtype=np.int64)

    if n <= k:
        return np.arange(n, dtype=np.int64)

    cutoff = np.partition(score_for_sort, k - 1)[k - 1]
    better = np.flatnonzero(score_for_sort < cutoff)
    tied = np.flatnonzero(score_for_sort == cutoff)
    remaining = k - better.size

    if tied.size > remaining:
        phi = _edge_phi(edge_arrays, tied)
        phi = np.nan_to_num(np.abs(phi), nan=-np.inf)
        support = np.asarray(cols["shared_sample_count"])[tied]
        source_idx = np.asarray(cols["source_taxon_index"], dtype=np.int64)[tied]
        target_idx = np.asarray(cols["target_taxon_index"], dtype=np.int64)[tied]

        taxa_order = np.argsort(np.asarray(taxa_universe, dtype=object), kind="stable")
        taxa_rank = np.empty(taxa_order.size, dtype=np.int64)
        taxa_rank[taxa_order] = np.arange(taxa_order.size, dtype=np.int64)

        tied_order = np.lexsort(
            (
                tied,
                taxa_rank[target_idx],
                taxa_rank[source_idx],
                -support,
                -phi,
            )
        )
        tied = tied[tied_order[:remaining]]

    return np.concatenate((better, tied[:remaining])).astype(np.int64, copy=False)


def _cooccurrence_summary_log_q_column(
    edge_arrays: CooccurrenceArrays,
) -> str:
    if "jaccard_null_log_q_value_bh" in edge_arrays.cols:
        return "jaccard_null_log_q_value_bh"
    return "chi2_log_q_value_bh"


def _edge_phi(
    edge_arrays: CooccurrenceArrays,
    idx: np.ndarray,
) -> np.ndarray:
    cols = edge_arrays.cols
    totals = np.asarray(edge_arrays.meta["totals"], dtype=np.int32)
    n_total = int(edge_arrays.meta["N_total"])

    source_idx = np.asarray(cols["source_taxon_index"], dtype=np.int64)[idx]
    target_idx = np.asarray(cols["target_taxon_index"], dtype=np.int64)[idx]
    shared = np.asarray(cols["shared_sample_count"], dtype=np.float64)[idx]
    source_total = totals[source_idx].astype(np.float64, copy=False)
    target_total = totals[target_idx].astype(np.float64, copy=False)

    _, phi, _ = _chi2_phi_from_counts(
        shared,
        source_total - shared,
        target_total - shared,
        n_total - source_total - target_total + shared,
    )
    return phi


def _build_cooccurrence_summary(
    edge_arrays: CooccurrenceArrays,
    taxa_universe: List[str],
    summary_n: int,
    detailed_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    idx = _edge_summary_indices(
        edge_arrays,
        taxa_universe=taxa_universe,
        summary_n=summary_n,
    )
    if detailed_df is None:
        detailed = _reconstruct_edge_frame(
            edge_arrays=edge_arrays,
            taxa_universe=taxa_universe,
            idx=idx,
            compute_fisher=False,
        )
    else:
        detailed = detailed_df.iloc[idx].copy()

    score_column = (
        "jaccard_null_q_value_bh"
        if "jaccard_null_q_value_bh" in detailed
        else "chi2_q_value_bh"
    )
    detailed["_summary_abs_phi"] = detailed["phi_coefficient"].abs()
    detailed.sort_values(
        by=[
            score_column,
            "_summary_abs_phi",
            "shared_sample_count",
            "source_taxon",
            "target_taxon",
        ],
        ascending=[True, False, False, True, True],
        na_position="last",
        kind="mergesort",
        inplace=True,
    )
    detailed.drop(columns="_summary_abs_phi", inplace=True)
    detailed.reset_index(drop=True, inplace=True)

    columns = [
        column for column in ("focal_query", "focal_taxon") if column in detailed
    ]
    columns.extend(COOCCURRENCE_SUMMARY_BASE_COLUMNS)
    columns.extend(
        column
        for column in COOCCURRENCE_SUMMARY_NULL_COLUMNS
        if column in detailed
    )
    return detailed.loc[:, columns].copy()


def _finalize_cooccurrence_edge_frame(df: pd.DataFrame, *, compute_fisher: bool) -> pd.DataFrame:
    columns = [col for col in ("focal_query", "focal_taxon") if col in df.columns]
    columns.extend(COOCCURRENCE_BASE_COLUMNS)
    if compute_fisher:
        columns.extend(COOCCURRENCE_FISHER_COLUMNS)
    columns.extend(col for col in COOCCURRENCE_NULL_COLUMNS if col in df.columns)
    columns.extend(col for col in df.columns if col not in columns)
    return df.loc[:, columns].copy()


def _write_full_edges_parquet_from_arrays(
    edge_arrays: CooccurrenceArrays,
    taxa_universe: List[str],
    output_path: str,
    *,
    null_model: str = "FE",
    batch_size: int = _COOCCURRENCE_PARQUET_BATCH_SIZE,
) -> None:
    """
    Write every detailed co-occurrence metric using bounded-memory batches.
    """
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as e:
        raise ImportError(
            "Parquet export requested for a large co-occurrence edge table, but "
            "no parquet engine is installed. Install 'pyarrow'."
        ) from e

    totals = np.asarray(edge_arrays.meta["totals"], dtype=np.int32)
    taxa_arr = np.asarray(taxa_universe, dtype=object)
    compute_fisher = bool(edge_arrays.meta.get("compute_fisher", False))
    n_rows = edge_arrays.n_rows
    metadata = edge_arrays.meta.get("null_metadata")
    fallback_null_model = "FE" if str(null_model).upper() == "FE" else None

    if batch_size <= 0:
        raise ValueError("batch_size must be positive")
    if compute_fisher:
        print(
            "WARNING: --compute_fisher will calculate Fisher's exact test for "
            f"all {n_rows:,} edges while writing the detailed Parquet output. "
            "This can substantially increase export time."
        )

    edges_path = f"{output_path}.parquet"
    taxa_path = f"{output_path}_taxa.parquet"

    n_batches = (n_rows + batch_size - 1) // batch_size
    writer = None
    try:
        for batch_number, start in enumerate(range(0, n_rows, batch_size), start=1):
            stop = min(start + batch_size, n_rows)
            if n_batches > 1:
                print(
                    "INFO: Writing detailed cooccurrence batch "
                    f"{batch_number:,}/{n_batches:,} ({start:,}:{stop:,})."
                )

            batch = _reconstruct_edge_frame(
                edge_arrays=edge_arrays,
                taxa_universe=taxa_universe,
                idx=slice(start, stop),
                compute_fisher=compute_fisher,
            )
            batch = with_compact_null_metadata(
                batch,
                metadata,
                fallback_null_model=fallback_null_model,
            )
            source_idx = np.asarray(
                edge_arrays.cols["source_taxon_index"][start:stop],
                dtype=np.int32,
            )
            target_idx = np.asarray(
                edge_arrays.cols["target_taxon_index"][start:stop],
                dtype=np.int32,
            )
            batch.drop(columns=["source_taxon", "target_taxon"], inplace=True)
            id_position = sum(
                column in batch.columns for column in ("focal_query", "focal_taxon")
            )
            batch.insert(id_position, "source_taxon_index", source_idx)
            batch.insert(id_position + 1, "target_taxon_index", target_idx)

            table = pa.Table.from_pandas(batch, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(
                    edges_path,
                    table.schema,
                    compression="snappy",
                )
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()

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
    np.minimum(log_q, 0.0, out=log_q)

    if nan_mask.any():
        log_q[nan_mask] = np.nan

    return log_q


def _valid_edge_mask_from_counts(
    shared_sample_count: np.ndarray,
    source_taxon_sample_count: np.ndarray,
    target_taxon_sample_count: np.ndarray,
    N_total: int,
) -> np.ndarray:
    """
    Cheap validity mask from integer counts only.
    """
    neither_source_nor_target_count = (
        N_total - source_taxon_sample_count - target_taxon_sample_count + shared_sample_count
    )

    return (
        (shared_sample_count >= 0) &
        (source_taxon_sample_count > 0) & (source_taxon_sample_count < N_total) &
        (target_taxon_sample_count > 0) & (target_taxon_sample_count < N_total) &
        (neither_source_nor_target_count >= 0)
    )


def _compute_core_edge_metrics(
    shared_sample_count: np.ndarray,
    source_taxon_sample_count: np.ndarray,
    target_taxon_sample_count: np.ndarray,
    N_total: int,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute only the core per-edge metrics we actually need to keep in memory.
    """
    shared = shared_sample_count.astype(np.float64, copy=False)
    source_total = source_taxon_sample_count.astype(np.float64, copy=False)
    target_total = target_taxon_sample_count.astype(np.float64, copy=False)
    n = float(N_total)

    cross = shared * n - source_total * target_total

    denom = source_total * (n - source_total)
    denom *= target_total
    denom *= (n - target_total)

    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = n * (cross * cross) / denom
        chi2_log_p_value = np.log(special.gammaincc(0.5, 0.5 * chi2))

    union = source_total + target_total - shared
    jaccard_taxon_pair = np.divide(
        shared,
        union,
        out=np.zeros(shared.shape[0], dtype=np.float32),
        where=union > 0,
    )

    return chi2_log_p_value, jaccard_taxon_pair


def _reconstruct_edge_frame(
    edge_arrays: CooccurrenceArrays,
    taxa_universe: List[str],
    idx: np.ndarray | slice,
    compute_fisher: bool | None = None,
) -> pd.DataFrame:
    cols = edge_arrays.cols
    totals = np.asarray(edge_arrays.meta["totals"], dtype=np.int32)
    N_total = int(edge_arrays.meta["N_total"])
    if compute_fisher is None:
        compute_fisher = bool(edge_arrays.meta.get("compute_fisher", False))

    source_taxon_index = cols["source_taxon_index"][idx]
    target_taxon_index = cols["target_taxon_index"][idx]
    shared_sample_count = cols["shared_sample_count"][idx]

    source_taxon_sample_count = totals[source_taxon_index]
    target_taxon_sample_count = totals[target_taxon_index]

    source_only_sample_count = source_taxon_sample_count - shared_sample_count
    target_only_sample_count = target_taxon_sample_count - shared_sample_count
    neither_source_nor_target_sample_count = (
        N_total - source_taxon_sample_count - target_taxon_sample_count + shared_sample_count
    )

    mets = table_metrics_arrays(
        shared_sample_count.astype(np.float64, copy=False),
        source_only_sample_count.astype(np.float64, copy=False),
        target_only_sample_count.astype(np.float64, copy=False),
        neither_source_nor_target_sample_count.astype(np.float64, copy=False),
        compute_fisher=compute_fisher,
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        source_taxon_prevalence = source_taxon_sample_count / float(N_total)
        target_taxon_prevalence = target_taxon_sample_count / float(N_total)
        cooccurrence_prevalence = shared_sample_count / float(N_total)

        denom_lift = source_taxon_prevalence * target_taxon_prevalence
        lift_taxon_pair = np.divide(
            cooccurrence_prevalence,
            denom_lift,
            out=np.zeros_like(cooccurrence_prevalence, dtype=float),
            where=denom_lift != 0,
        )

        p_target_given_source = np.divide(
            shared_sample_count,
            source_taxon_sample_count,
            out=np.zeros_like(shared_sample_count, dtype=float),
            where=source_taxon_sample_count > 0,
        )

        p_source_given_target = np.divide(
            shared_sample_count,
            target_taxon_sample_count,
            out=np.zeros_like(shared_sample_count, dtype=float),
            where=target_taxon_sample_count > 0,
        )

    taxa_arr = np.asarray(taxa_universe, dtype=object)

    out = {
        "source_taxon": taxa_arr[source_taxon_index],
        "target_taxon": taxa_arr[target_taxon_index],
        "p_target_given_source": p_target_given_source,
        "p_source_given_target": p_source_given_target,
        "shared_sample_count": shared_sample_count,
        "source_taxon_sample_count": source_taxon_sample_count,
        "target_taxon_sample_count": target_taxon_sample_count,
        "source_only_sample_count": source_only_sample_count,
        "target_only_sample_count": target_only_sample_count,
        "neither_source_nor_target_sample_count": neither_source_nor_target_sample_count,
        "background_sample_count": np.full(shared_sample_count.shape[0], N_total, dtype=np.int32),
        "source_taxon_prevalence": source_taxon_prevalence,
        "target_taxon_prevalence": target_taxon_prevalence,
        "cooccurrence_prevalence": cooccurrence_prevalence,
        "lift_taxon_pair": lift_taxon_pair,
        "jaccard_taxon_pair": cols["jaccard_taxon_pair"][idx],
        "phi_coefficient": mets["phi"],
        "chi2_statistic": mets["chi2"],
        "chi2_p_value": np.exp(cols["chi2_log_p_value"][idx]),
        "chi2_q_value_bh": np.exp(cols["chi2_log_q_value_bh"][idx]),
        "chi2_log_p_value": cols["chi2_log_p_value"][idx],
        "chi2_log_q_value_bh": cols["chi2_log_q_value_bh"][idx],
        "rr_target_given_source_vs_without_source": mets["RR_A_to_B"],
        "rr_source_given_target_vs_without_target": mets["RR_B_to_A"],
        "ln_rr_target_given_source_vs_without_source": mets["logRR_A_to_B"],
        "ln_rr_source_given_target_vs_without_target": mets["logRR_B_to_A"],
    }

    if compute_fisher:
        out.update(
            {
                "fisher_odds_ratio": mets["fisher_odds"],
                "fisher_p_value": mets["fisher_p"],
                "fisher_log_p_value": mets["log_fisher_p"],
            }
        )

    for key, arr in cols.items():
        if key in {
            "source_taxon_index",
            "target_taxon_index",
            "shared_sample_count",
            "jaccard_taxon_pair",
            "chi2_log_p_value",
            "chi2_log_q_value_bh",
        }:
            continue

        if key == "jaccard_null_log_q_value_bh":
            out["jaccard_null_q_value_bh"] = np.exp(arr[idx])
            out[key] = arr[idx]
        else:
            out[key] = arr[idx]

    return _finalize_cooccurrence_edge_frame(
        pd.DataFrame(out),
        compute_fisher=compute_fisher,
    )
