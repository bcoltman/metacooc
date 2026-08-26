# Plotting

Association and co-occurrence pipelines create analysis-specific plots by default. The default view uses the selected q-value metric, includes significant positive-phi results up to `q <= 0.10`, and labels the strongest results. Labels are placed with `adjustText` where available.

Use `--plot_all` to add nonsignificant rows as faint context, `--label_top_n` to change the number of labels, and paired `--x_metric`/`--y_metric` options to choose axes. `--q_metric` selects which probability column controls the significance filter. `--no_plot` disables automatic plotting.

The standalone `metacooc plot` command requires `--analysis_type association` or `--analysis_type cooccurrence` and can replot TSV or Parquet outputs. Structure output is intentionally left as a table because its metrics have different scales.
