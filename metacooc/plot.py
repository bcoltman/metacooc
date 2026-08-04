#!/usr/bin/env python3

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype


_DEFAULT_Q_THRESHOLD = 0.10
_DEFAULT_LABEL_TOP_N = 10


def _validate_plot_options(
    df: pd.DataFrame,
    *,
    q_metric: str,
    q_thresh: float,
    label_top_n: int,
    x_metric: str | None,
    y_metric: str | None,
) -> None:
    if not 0.0 <= float(q_thresh) <= 1.0:
        raise ValueError("q_thresh must be between 0 and 1")
    if int(label_top_n) < 0:
        raise ValueError("label_top_n must be non-negative")
    if (x_metric is None) != (y_metric is None):
        raise ValueError("x_metric and y_metric must be provided together")

    required = [q_metric, "phi_coefficient"]
    if x_metric is not None:
        required.extend([x_metric, y_metric])
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise ValueError(
            "Plot input is missing required column(s): " + ", ".join(missing)
        )

    numeric_columns = [q_metric, "phi_coefficient"]
    if x_metric is not None:
        numeric_columns.extend([x_metric, y_metric])
    nonnumeric = [
        column for column in dict.fromkeys(numeric_columns)
        if not is_numeric_dtype(df[column])
    ]
    if nonnumeric:
        raise ValueError(
            "Plot metric column(s) must be numeric: " + ", ".join(nonnumeric)
        )

    finite_q = df[q_metric].replace([np.inf, -np.inf], np.nan).dropna()
    if ((finite_q < 0.0) | (finite_q > 1.0)).any():
        raise ValueError(f"Significance metric '{q_metric}' must contain values in [0, 1]")


def _short_taxon_label(value: object) -> str:
    parts = [part.strip() for part in str(value).split(";") if part.strip()]
    return parts[-1] if parts else str(value)


def _draw_points(
    ax,
    df: pd.DataFrame,
    *,
    x: str,
    y: str,
    significant: pd.Series,
    plot_all: bool,
) -> None:
    finite = np.isfinite(df[x].to_numpy(dtype=float, copy=False)) & np.isfinite(
        df[y].to_numpy(dtype=float, copy=False)
    )
    significant_np = significant.to_numpy(dtype=bool, copy=False) & finite

    if plot_all:
        context = finite & ~significant_np
        if context.any():
            ax.scatter(
                df.loc[context, x],
                df.loc[context, y],
                color="#bdbdbd",
                s=12,
                alpha=0.18,
                linewidths=0,
            )

    if significant_np.any():
        ax.scatter(
            df.loc[significant_np, x],
            df.loc[significant_np, y],
            c=df.loc[significant_np, "phi_coefficient"],
            cmap="coolwarm",
            vmin=-1.0,
            vmax=1.0,
            s=28,
            alpha=0.85,
            linewidths=0,
        )


def _top_association_rows(
    df: pd.DataFrame,
    *,
    significant: pd.Series,
    q_metric: str,
    label_top_n: int,
    x_metric: str,
    y_metric: str,
) -> pd.DataFrame:
    if label_top_n <= 0:
        return df.iloc[:0]

    candidates = df.loc[significant].copy()
    candidates = candidates[
        np.isfinite(candidates[x_metric].to_numpy(dtype=float, copy=False))
        & np.isfinite(candidates[y_metric].to_numpy(dtype=float, copy=False))
    ]
    specificity = candidates["p_cohort_given_taxon"].to_numpy(dtype=float)
    sensitivity = candidates["p_taxon_given_cohort"].to_numpy(dtype=float)
    denominator = specificity + sensitivity
    balanced_score = np.zeros_like(denominator)
    np.divide(
        2 * specificity * sensitivity,
        denominator,
        out=balanced_score,
        where=denominator > 0,
    )
    candidates["_balanced_score"] = balanced_score
    candidates.sort_values(
        ["_balanced_score", q_metric, "phi_coefficient", "taxon"],
        ascending=[False, True, False, True],
        kind="mergesort",
        inplace=True,
    )
    return candidates.head(int(label_top_n))


def _annotate_association_rows(
    ax,
    rows: pd.DataFrame,
    *,
    x_metric: str,
    y_metric: str,
) -> None:
    for offset, (_, row) in enumerate(rows.iterrows()):
        ax.annotate(
            _short_taxon_label(row["taxon"]),
            (row[x_metric], row[y_metric]),
            xytext=(4, 5 + (offset % 3) * 4),
            textcoords="offset points",
            fontsize=7,
        )


def _save_empty_plot(out_file: str, message: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axis("off")
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_association(
    df: pd.DataFrame,
    out_file: str,
    *,
    q_thresh: float,
    q_metric: str | None,
    plot_all: bool,
    label_top_n: int,
    x_metric: str | None,
    y_metric: str | None,
) -> None:
    q_metric = q_metric or "chi2_q_value_bh"
    _validate_plot_options(
        df,
        q_metric=q_metric,
        q_thresh=q_thresh,
        label_top_n=label_top_n,
        x_metric=x_metric,
        y_metric=y_metric,
    )

    required_defaults = set()
    if x_metric is None or label_top_n > 0:
        required_defaults.update(
            {"taxon", "p_cohort_given_taxon", "p_taxon_given_cohort"}
        )
    if required_defaults:
        missing = sorted(required_defaults.difference(df.columns))
        if missing:
            raise ValueError(
                "Association plot input is missing required column(s): "
                + ", ".join(missing)
            )

    clean = df.replace([np.inf, -np.inf], np.nan).copy()
    significant = (
        clean[q_metric].notna()
        & (clean[q_metric] <= q_thresh)
        & (clean["phi_coefficient"] > 0)
    )
    if not significant.any() and not plot_all:
        _save_empty_plot(
            out_file,
            f"No positive association results at {q_metric} ≤ {q_thresh:g}",
        )
        return

    if x_metric is not None:
        fig, ax = plt.subplots(figsize=(8, 6))
        _draw_points(
            ax,
            clean,
            x=x_metric,
            y=y_metric,
            significant=significant,
            plot_all=plot_all,
        )
        labels = _top_association_rows(
            clean,
            significant=significant,
            q_metric=q_metric,
            label_top_n=label_top_n,
            x_metric=x_metric,
            y_metric=y_metric,
        )
        _annotate_association_rows(
            ax,
            labels,
            x_metric=x_metric,
            y_metric=y_metric,
        )
        ax.set_xlabel(x_metric)
        ax.set_ylabel(y_metric)
        ax.set_title(f"Association: {x_metric} vs {y_metric}")
        axes = [ax]
    else:
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        primary_x = "p_cohort_given_taxon"
        primary_y = "phi_coefficient"
        _draw_points(
            axes[0],
            clean,
            x=primary_x,
            y=primary_y,
            significant=significant,
            plot_all=plot_all,
        )
        labels = _top_association_rows(
            clean,
            significant=significant,
            q_metric=q_metric,
            label_top_n=label_top_n,
            x_metric="p_cohort_given_taxon",
            y_metric="p_taxon_given_cohort",
        )
        _annotate_association_rows(
            axes[1],
            labels,
            x_metric="p_cohort_given_taxon",
            y_metric="p_taxon_given_cohort",
        )
        axes[0].set_xlabel("P(cohort | taxon)")
        axes[0].set_ylabel("Phi coefficient")
        axes[0].set_title("Specificity vs association strength")

        _draw_points(
            axes[1],
            clean,
            x="p_cohort_given_taxon",
            y="p_taxon_given_cohort",
            significant=significant,
            plot_all=plot_all,
        )
        axes[1].set_xlabel("P(cohort | taxon)")
        axes[1].set_ylabel("P(taxon | cohort)")
        axes[1].set_title("Specificity vs sensitivity")

        ranked = clean.sort_values(
            "p_cohort_given_taxon",
            ascending=False,
            na_position="last",
            kind="mergesort",
        ).copy()
        ranked["specificity_rank"] = np.arange(1, len(ranked) + 1)
        ranked_significant = significant.reindex(ranked.index)
        _draw_points(
            axes[2],
            ranked,
            x="specificity_rank",
            y="phi_coefficient",
            significant=ranked_significant,
            plot_all=plot_all,
        )
        axes[2].set_xlabel("Taxon rank by P(cohort | taxon)")
        axes[2].set_ylabel("Phi coefficient")
        axes[2].set_title("Association strength by specificity rank")

    for ax in axes:
        ax.axhline(0, color="black", linewidth=0.8, alpha=0.7)
        ax.grid(True, linestyle="--", alpha=0.25)

    fig.suptitle(
        f"Association results: {int(significant.sum()):,} positive associations of "
        f"{len(clean):,} (q ≤ {q_thresh:g})"
    )
    fig.tight_layout()
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_analysis_obj(
    df: pd.DataFrame,
    out_file: str,
    q_thresh: float = _DEFAULT_Q_THRESHOLD,
    *,
    analysis_type: str = "association",
    q_metric: str | None = None,
    plot_all: bool = False,
    label_top_n: int = _DEFAULT_LABEL_TOP_N,
    x_metric: str | None = None,
    y_metric: str | None = None,
) -> None:
    """Plot an in-memory association or co-occurrence result."""
    if df.empty:
        raise ValueError("DataFrame is empty — nothing to plot.")
    if analysis_type != "association":
        raise ValueError(f"Unsupported analysis_type: {analysis_type!r}")

    _plot_association(
        df,
        out_file,
        q_thresh=q_thresh,
        q_metric=q_metric,
        plot_all=plot_all,
        label_top_n=label_top_n,
        x_metric=x_metric,
        y_metric=y_metric,
    )
    print(f"[plot_analysis] Saved: {out_file}")


def plot_analysis(
    df_file: str,
    output_dir: str,
    tag: str = "",
    q_thresh: float = _DEFAULT_Q_THRESHOLD,
    *,
    analysis_type: str = "association",
    q_metric: str | None = None,
    plot_all: bool = False,
    label_top_n: int = _DEFAULT_LABEL_TOP_N,
    x_metric: str | None = None,
    y_metric: str | None = None,
) -> str:
    if not os.path.exists(df_file):
        raise FileNotFoundError(df_file)

    os.makedirs(output_dir, exist_ok=True)
    df = pd.read_csv(df_file, sep="\t")
    out_file = os.path.join(output_dir, f"{tag}{analysis_type}_plot.png")

    plot_analysis_obj(
        df,
        out_file,
        q_thresh=q_thresh,
        analysis_type=analysis_type,
        q_metric=q_metric,
        plot_all=plot_all,
        label_top_n=label_top_n,
        x_metric=x_metric,
        y_metric=y_metric,
    )
    return out_file
