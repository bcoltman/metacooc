#!/usr/bin/env python3

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from pandas.api.types import is_numeric_dtype

_DEFAULT_Q_THRESHOLD = 0.10
_DEFAULT_LABEL_TOP_N = 10
_PLOT_CHUNK_SIZE = 100_000


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
        column
        for column in dict.fromkeys(numeric_columns)
        if not is_numeric_dtype(df[column])
    ]
    if nonnumeric:
        raise ValueError(
            "Plot metric column(s) must be numeric: " + ", ".join(nonnumeric)
        )

    finite_q = df[q_metric].replace([np.inf, -np.inf], np.nan).dropna()
    if ((finite_q < 0.0) | (finite_q > 1.0)).any():
        raise ValueError(
            f"Significance metric '{q_metric}' must contain values in [0, 1]"
        )


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


def _annotate_points(ax, labels: list[tuple[str, float, float]]) -> None:
    if not labels:
        return

    texts = [ax.text(x, y, label, fontsize=7) for label, x, y in labels]
    adjust_text(
        texts,
        target_x=[x for _, x, _ in labels],
        target_y=[y for _, _, y in labels],
        ax=ax,
        expand=(1.2, 1.4),
        force_text=(0.2, 0.4),
        force_static=(0.2, 0.4),
        arrowprops={"arrowstyle": "-", "color": "#666666", "lw": 0.5},
    )


def _association_labels(
    rows: pd.DataFrame,
    *,
    x_metric: str,
    y_metric: str,
) -> list[tuple[str, float, float]]:
    labels = []
    for _, row in rows.iterrows():
        labels.append(
            (
                _short_taxon_label(row["taxon"]),
                float(row[x_metric]),
                float(row[y_metric]),
            )
        )
    return labels


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
        label_sets = [
            (
                ax,
                _association_labels(
                    labels,
                    x_metric=x_metric,
                    y_metric=y_metric,
                ),
            )
        ]
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

        ranked_by_taxon = ranked.set_index("taxon")["specificity_rank"]
        labels_with_rank = labels.copy()
        labels_with_rank["specificity_rank"] = labels_with_rank["taxon"].map(
            ranked_by_taxon
        )
        label_sets = [
            (
                axes[0],
                _association_labels(
                    labels,
                    x_metric="p_cohort_given_taxon",
                    y_metric="phi_coefficient",
                ),
            ),
            (
                axes[1],
                _association_labels(
                    labels,
                    x_metric="p_cohort_given_taxon",
                    y_metric="p_taxon_given_cohort",
                ),
            ),
            (
                axes[2],
                _association_labels(
                    labels_with_rank,
                    x_metric="specificity_rank",
                    y_metric="phi_coefficient",
                ),
            ),
        ]

    for ax in axes:
        ax.axhline(0, color="black", linewidth=0.8, alpha=0.7)
        ax.grid(True, linestyle="--", alpha=0.25)

    fig.suptitle(
        f"Association results: {int(significant.sum()):,} positive associations of "
        f"{len(clean):,} (q ≤ {q_thresh:g})"
    )
    fig.tight_layout()
    for label_axis, labels in label_sets:
        _annotate_points(label_axis, labels)
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _cooccurrence_q_metric(columns, requested: str | None) -> str:
    if requested is not None:
        return requested
    if "jaccard_null_q_value_bh" in columns:
        return "jaccard_null_q_value_bh"
    return "chi2_q_value_bh"


def _cooccurrence_identifier_columns(columns) -> tuple[str, str]:
    if {"source_taxon", "target_taxon"}.issubset(columns):
        return "source_taxon", "target_taxon"
    if {"source_taxon_index", "target_taxon_index"}.issubset(columns):
        return "source_taxon_index", "target_taxon_index"
    raise ValueError(
        "Co-occurrence plot input must contain source/target taxon names or indices"
    )


def _top_cooccurrence_rows(
    df: pd.DataFrame,
    *,
    significant: pd.Series,
    q_metric: str,
    label_top_n: int,
    x_metric: str,
    y_metric: str,
    source_column: str,
    target_column: str,
) -> pd.DataFrame:
    if label_top_n <= 0:
        return df.iloc[:0]

    candidates = df.loc[significant].copy()
    candidates = candidates[
        np.isfinite(candidates[x_metric].to_numpy(dtype=float, copy=False))
        & np.isfinite(candidates[y_metric].to_numpy(dtype=float, copy=False))
    ]
    candidates["_abs_phi"] = candidates["phi_coefficient"].abs()
    candidates["_support"] = (
        candidates["shared_sample_count"] if "shared_sample_count" in candidates else 0
    )
    candidates["_source_sort"] = candidates[source_column].astype(str)
    candidates["_target_sort"] = candidates[target_column].astype(str)
    candidates.sort_values(
        [q_metric, "_abs_phi", "_support", "_source_sort", "_target_sort"],
        ascending=[True, False, False, True, True],
        kind="mergesort",
        inplace=True,
    )
    return candidates.head(int(label_top_n))


def _edge_label(
    row: pd.Series,
    *,
    source_column: str,
    target_column: str,
    taxa_by_id: dict[int, str] | None,
) -> str:
    source = row[source_column]
    target = row[target_column]
    if taxa_by_id is not None:
        source = taxa_by_id.get(int(source), source)
        target = taxa_by_id.get(int(target), target)
    label = f"{_short_taxon_label(source)} → {_short_taxon_label(target)}"
    if "focal_query" in row and pd.notna(row["focal_query"]):
        label = f"{row['focal_query']}: {label}"
    return label


def _plot_cooccurrence_chunks(
    chunks,
    columns,
    out_file: str,
    *,
    q_thresh: float,
    q_metric: str | None,
    plot_all: bool,
    label_top_n: int,
    x_metric: str | None,
    y_metric: str | None,
    taxa_by_id: dict[int, str] | None = None,
    report_progress: bool = False,
) -> None:
    q_metric = _cooccurrence_q_metric(columns, q_metric)
    x_metric = x_metric or "p_target_given_source"
    y_metric = y_metric or "phi_coefficient"
    source_column, target_column = _cooccurrence_identifier_columns(set(columns))

    required = {q_metric, "phi_coefficient", x_metric, y_metric}
    missing = sorted(required.difference(columns))
    if missing:
        raise ValueError(
            "Co-occurrence plot input is missing required column(s): "
            + ", ".join(missing)
        )

    fig, ax = plt.subplots(figsize=(8, 6))
    total_rows = 0
    significant_rows = 0
    label_candidates = []
    try:
        for chunk_number, chunk in enumerate(chunks, start=1):
            if chunk.empty:
                continue
            if report_progress:
                print(
                    "INFO: Plotting co-occurrence chunk "
                    f"{chunk_number:,} ({len(chunk):,} rows)."
                )
            _validate_plot_options(
                chunk,
                q_metric=q_metric,
                q_thresh=q_thresh,
                label_top_n=label_top_n,
                x_metric=x_metric,
                y_metric=y_metric,
            )
            clean = chunk.replace([np.inf, -np.inf], np.nan)
            significant = (
                clean[q_metric].notna()
                & (clean[q_metric] <= q_thresh)
                & (clean["phi_coefficient"] > 0)
            )
            _draw_points(
                ax,
                clean,
                x=x_metric,
                y=y_metric,
                significant=significant,
                plot_all=plot_all,
            )
            total_rows += len(clean)
            significant_rows += int(significant.sum())
            if label_top_n > 0:
                label_candidates.append(
                    _top_cooccurrence_rows(
                        clean,
                        significant=significant,
                        q_metric=q_metric,
                        label_top_n=label_top_n,
                        x_metric=x_metric,
                        y_metric=y_metric,
                        source_column=source_column,
                        target_column=target_column,
                    )
                )
    except Exception:
        plt.close(fig)
        raise

    if total_rows == 0:
        plt.close(fig)
        raise ValueError("Co-occurrence result is empty — nothing to plot.")
    if significant_rows == 0 and not plot_all:
        plt.close(fig)
        _save_empty_plot(
            out_file,
            f"No positive co-occurrence results at {q_metric} ≤ {q_thresh:g}",
        )
        return

    if label_candidates:
        candidates = pd.concat(label_candidates, ignore_index=True)
        candidates["_abs_phi"] = candidates["phi_coefficient"].abs()
        candidates.sort_values(
            [q_metric, "_abs_phi", "_support", "_source_sort", "_target_sort"],
            ascending=[True, False, False, True, True],
            kind="mergesort",
            inplace=True,
        )
        labels = []
        for _, row in candidates.head(label_top_n).iterrows():
            labels.append(
                (
                    _edge_label(
                        row,
                        source_column=source_column,
                        target_column=target_column,
                        taxa_by_id=taxa_by_id,
                    ),
                    float(row[x_metric]),
                    float(row[y_metric]),
                )
            )
    ax.axhline(0, color="black", linewidth=0.8, alpha=0.7)
    ax.grid(True, linestyle="--", alpha=0.25)
    ax.set_xlabel(
        "P(target | source)" if x_metric == "p_target_given_source" else x_metric
    )
    ax.set_ylabel("Phi coefficient" if y_metric == "phi_coefficient" else y_metric)
    ax.set_title(f"Co-occurrence: {x_metric} vs {y_metric}")
    fig.suptitle(
        f"Co-occurrence results: {significant_rows:,} positive co-occurrences of "
        f"{total_rows:,} (q ≤ {q_thresh:g})"
    )
    fig.tight_layout()
    if label_candidates:
        _annotate_points(ax, labels)
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _cooccurrence_file_columns(df_file: str) -> list[str]:
    if df_file.lower().endswith(".parquet"):
        import pyarrow.parquet as pq

        return pq.ParquetFile(df_file).schema_arrow.names
    if df_file.lower().endswith(".tsv"):
        return pd.read_csv(df_file, sep="\t", nrows=0).columns.tolist()
    raise ValueError("Plot input must be a .tsv or .parquet file")


def _cooccurrence_file_chunks(df_file: str, columns: list[str]):
    if df_file.lower().endswith(".parquet"):
        import pyarrow.parquet as pq

        parquet = pq.ParquetFile(df_file)
        for batch in parquet.iter_batches(
            batch_size=_PLOT_CHUNK_SIZE,
            columns=columns,
        ):
            yield batch.to_pandas()
        return

    yield from pd.read_csv(
        df_file,
        sep="\t",
        usecols=columns,
        chunksize=_PLOT_CHUNK_SIZE,
    )


def _cooccurrence_taxa_mapping(df_file: str) -> dict[int, str] | None:
    if not df_file.lower().endswith(".parquet"):
        return None
    taxa_path = f"{os.path.splitext(df_file)[0]}_taxa.parquet"
    if not os.path.exists(taxa_path):
        return None
    taxa = pd.read_parquet(taxa_path, columns=["taxon_id", "taxon"])
    return dict(zip(taxa["taxon_id"].astype(int), taxa["taxon"].astype(str)))


def _plot_cooccurrence_file(
    df_file: str,
    out_file: str,
    *,
    q_thresh: float,
    q_metric: str | None,
    plot_all: bool,
    label_top_n: int,
    x_metric: str | None,
    y_metric: str | None,
) -> None:
    available = _cooccurrence_file_columns(df_file)
    selected_q = _cooccurrence_q_metric(available, q_metric)
    selected_x = x_metric or "p_target_given_source"
    selected_y = y_metric or "phi_coefficient"
    source_column, target_column = _cooccurrence_identifier_columns(set(available))
    selected = {
        selected_q,
        selected_x,
        selected_y,
        "phi_coefficient",
        source_column,
        target_column,
    }
    for optional in ("shared_sample_count", "focal_query"):
        if optional in available:
            selected.add(optional)

    _plot_cooccurrence_chunks(
        _cooccurrence_file_chunks(df_file, sorted(selected)),
        available,
        out_file,
        q_thresh=q_thresh,
        q_metric=selected_q,
        plot_all=plot_all,
        label_top_n=label_top_n,
        x_metric=selected_x,
        y_metric=selected_y,
        taxa_by_id=_cooccurrence_taxa_mapping(df_file),
        report_progress=True,
    )


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
    if analysis_type == "association":
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
    elif analysis_type == "cooccurrence":
        _plot_cooccurrence_chunks(
            [df],
            df.columns,
            out_file,
            q_thresh=q_thresh,
            q_metric=q_metric,
            plot_all=plot_all,
            label_top_n=label_top_n,
            x_metric=x_metric,
            y_metric=y_metric,
        )
    else:
        raise ValueError(f"Unsupported analysis_type: {analysis_type!r}")
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
    out_file = os.path.join(output_dir, f"{tag}{analysis_type}_plot.png")
    if analysis_type == "association":
        if df_file.lower().endswith(".parquet"):
            df = pd.read_parquet(df_file)
        elif df_file.lower().endswith(".tsv"):
            df = pd.read_csv(df_file, sep="\t")
        else:
            raise ValueError("Plot input must be a .tsv or .parquet file")
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
    elif analysis_type == "cooccurrence":
        _plot_cooccurrence_file(
            df_file,
            out_file,
            q_thresh=q_thresh,
            q_metric=q_metric,
            plot_all=plot_all,
            label_top_n=label_top_n,
            x_metric=x_metric,
            y_metric=y_metric,
        )
        print(f"[plot_analysis] Saved: {out_file}")
    else:
        raise ValueError(f"Unsupported analysis_type: {analysis_type!r}")
    return out_file
