from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import pandas as pd
import pytest

import metacooc.plot as plot_module
from metacooc.plot import plot_analysis, plot_analysis_obj


@pytest.fixture
def association_plot_df():
    return pd.DataFrame(
        {
            "taxon": [
                "d__Bacteria; g__Alpha; s__one",
                "d__Bacteria; g__Alpha; s__two",
                "d__Bacteria; g__Beta; s__three",
                "d__Bacteria; g__Beta; s__four",
            ],
            "p_cohort_given_taxon": [0.9, 0.7, 0.4, 0.2],
            "p_taxon_given_cohort": [0.3, 0.6, 0.8, 0.1],
            "phi_coefficient": [0.8, -0.5, 0.2, -0.1],
            "chi2_q_value_bh": [0.01, 0.08, 0.2, float("nan")],
            "lift_taxon_cohort": [2.0, 1.5, 0.9, 0.5],
        }
    )


@pytest.fixture
def cooccurrence_plot_df():
    return pd.DataFrame(
        {
            "source_taxon": [
                "g__Alpha; s__one",
                "g__Alpha; s__two",
                "g__Beta; s__three",
                "g__Beta; s__four",
            ],
            "target_taxon": [
                "g__Gamma; s__a",
                "g__Gamma; s__b",
                "g__Delta; s__c",
                "g__Delta; s__d",
            ],
            "p_target_given_source": [0.9, 0.7, 0.5, 0.3],
            "p_source_given_target": [0.4, 0.5, 0.6, 0.7],
            "phi_coefficient": [0.8, -0.5, 0.3, -0.2],
            "chi2_q_value_bh": [0.01, 0.02, 0.03, 0.04],
            "jaccard_null_q_value_bh": [0.02, 0.2, 0.08, 0.9],
            "shared_sample_count": [20, 15, 10, 5],
            "lift_taxon_pair": [2.5, 1.4, 1.1, 0.8],
        }
    )


def _capture_closed_figures(monkeypatch):
    figures = []
    monkeypatch.setattr(plot_module.plt, "close", lambda figure: figures.append(figure))
    return figures


def test_association_default_plot_has_three_focused_panels(
    tmp_path,
    monkeypatch,
    association_plot_df,
):
    figures = _capture_closed_figures(monkeypatch)
    output = tmp_path / "association.png"

    plot_analysis_obj(
        association_plot_df,
        str(output),
        analysis_type="association",
        label_top_n=1,
    )

    assert output.exists()
    assert len(figures) == 1
    axes = figures[0].axes
    assert len(axes) == 3
    assert axes[0].get_xlabel() == "P(cohort | taxon)"
    assert axes[0].get_ylabel() == "Phi coefficient"
    assert axes[1].get_ylabel() == "P(taxon | cohort)"
    assert axes[2].get_xlabel() == "Taxon rank by P(cohort | taxon)"
    assert sum(len(collection.get_offsets()) for collection in axes[0].collections) == 1
    assert [text.get_text() for text in axes[1].texts] == ["s__one"]


def test_association_custom_metrics_switch_to_one_panel(
    tmp_path,
    monkeypatch,
    association_plot_df,
):
    figures = _capture_closed_figures(monkeypatch)

    plot_analysis_obj(
        association_plot_df,
        str(tmp_path / "custom.png"),
        x_metric="lift_taxon_cohort",
        y_metric="p_taxon_given_cohort",
        label_top_n=0,
    )

    assert len(figures[0].axes) == 1
    assert figures[0].axes[0].get_xlabel() == "lift_taxon_cohort"
    assert figures[0].axes[0].get_ylabel() == "p_taxon_given_cohort"


def test_association_plot_all_adds_insignificant_context(
    tmp_path,
    monkeypatch,
    association_plot_df,
):
    figures = _capture_closed_figures(monkeypatch)

    plot_analysis_obj(
        association_plot_df,
        str(tmp_path / "all.png"),
        plot_all=True,
        label_top_n=0,
    )

    primary = figures[0].axes[0]
    assert sum(len(collection.get_offsets()) for collection in primary.collections) == 4


def test_association_no_significant_rows_writes_message_plot(
    tmp_path,
    monkeypatch,
    association_plot_df,
):
    figures = _capture_closed_figures(monkeypatch)

    plot_analysis_obj(
        association_plot_df,
        str(tmp_path / "empty.png"),
        q_thresh=0.0,
    )

    assert len(figures[0].axes) == 1
    assert "No positive association results" in figures[0].axes[0].texts[0].get_text()


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"x_metric": "lift_taxon_cohort"}, "provided together"),
        ({"x_metric": "taxon", "y_metric": "phi_coefficient"}, "must be numeric"),
        ({"q_metric": "missing_q"}, "missing required"),
    ],
)
def test_association_plot_rejects_invalid_metric_options(
    tmp_path,
    association_plot_df,
    kwargs,
    message,
):
    with pytest.raises(ValueError, match=message):
        plot_analysis_obj(
            association_plot_df,
            str(tmp_path / "invalid.png"),
            **kwargs,
        )


def test_plot_analysis_uses_typed_association_filename(
    tmp_path,
    association_plot_df,
):
    input_path = tmp_path / "association.tsv"
    association_plot_df.to_csv(input_path, sep="\t", index=False)

    output_path = plot_analysis(
        str(input_path),
        str(tmp_path),
        tag="soil_",
        analysis_type="association",
        label_top_n=0,
    )

    assert output_path == str(tmp_path / "soil_association_plot.png")
    assert (tmp_path / "soil_association_plot.png").exists()


def test_cooccurrence_default_uses_empirical_q_and_one_panel(
    tmp_path,
    monkeypatch,
    cooccurrence_plot_df,
):
    figures = _capture_closed_figures(monkeypatch)

    plot_analysis_obj(
        cooccurrence_plot_df,
        str(tmp_path / "cooccurrence.png"),
        analysis_type="cooccurrence",
        label_top_n=1,
    )

    assert len(figures[0].axes) == 1
    axis = figures[0].axes[0]
    assert axis.get_xlabel() == "P(target | source)"
    assert axis.get_ylabel() == "Phi coefficient"
    assert sum(len(collection.get_offsets()) for collection in axis.collections) == 2
    assert [text.get_text() for text in axis.texts] == ["s__one → s__a"]


def test_cooccurrence_custom_metrics_and_q_override(
    tmp_path,
    monkeypatch,
    cooccurrence_plot_df,
):
    figures = _capture_closed_figures(monkeypatch)

    plot_analysis_obj(
        cooccurrence_plot_df,
        str(tmp_path / "cooccurrence_custom.png"),
        analysis_type="cooccurrence",
        q_metric="chi2_q_value_bh",
        x_metric="lift_taxon_pair",
        y_metric="p_source_given_target",
        label_top_n=0,
    )

    axis = figures[0].axes[0]
    assert axis.get_xlabel() == "lift_taxon_pair"
    assert axis.get_ylabel() == "p_source_given_target"
    assert sum(len(collection.get_offsets()) for collection in axis.collections) == 2


def test_cooccurrence_parquet_streams_every_edge_and_maps_labels(
    tmp_path,
    monkeypatch,
    capsys,
    cooccurrence_plot_df,
):
    import pyarrow as pa
    import pyarrow.parquet as pq

    figures = _capture_closed_figures(monkeypatch)
    monkeypatch.setattr(plot_module, "_PLOT_CHUNK_SIZE", 2)

    edges = cooccurrence_plot_df.drop(columns=["source_taxon", "target_taxon"]).copy()
    edges.insert(0, "source_taxon_index", [0, 1, 2, 3])
    edges.insert(1, "target_taxon_index", [4, 5, 6, 7])
    edges_path = tmp_path / "global_edges.parquet"
    pq.write_table(pa.Table.from_pandas(edges), edges_path, row_group_size=2)

    taxa = pd.DataFrame(
        {
            "taxon_id": range(8),
            "taxon": [
                "g__Alpha; s__one",
                "g__Alpha; s__two",
                "g__Beta; s__three",
                "g__Beta; s__four",
                "g__Gamma; s__a",
                "g__Gamma; s__b",
                "g__Delta; s__c",
                "g__Delta; s__d",
            ],
        }
    )
    pq.write_table(
        pa.Table.from_pandas(taxa),
        tmp_path / "global_edges_taxa.parquet",
    )

    output_path = plot_analysis(
        str(edges_path),
        str(tmp_path),
        analysis_type="cooccurrence",
        plot_all=True,
        label_top_n=1,
    )

    assert output_path == str(tmp_path / "cooccurrence_plot.png")
    assert sum(
        len(collection.get_offsets()) for collection in figures[0].axes[0].collections
    ) == len(edges)
    assert figures[0].axes[0].texts[0].get_text() == "s__one → s__a"
    output = capsys.readouterr().out
    assert "chunk 1" in output
    assert "chunk 2" in output
