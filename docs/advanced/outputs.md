# Outputs and scaling

Association writes a reduced summary TSV, a detailed TSV, and (unless disabled) an association plot. Co-occurrence writes a reduced edge summary, detailed edges, nodes, and a co-occurrence plot. Structure writes one focused TSV and is intentionally not plotted. Biome distribution writes a taxon-by-biome TSV and does not run a null model.

The common association summary columns are `taxon`, `phi_coefficient`, `p_cohort_given_taxon`, `p_taxon_given_cohort`, support counts, and `chi2_q_value_bh`, with compact empirical-null columns when applicable. Co-occurrence summaries contain the directed source and target, both conditional probabilities, support counts, phi, and adjusted significance. Detailed files contain the full metric set, optional Fisher fields, and row-level null statistics.

Association, co-occurrence, and structure results repeat `null_model`,
`null_replicates`, `null_replicates_failed`, and `null_seed`. Biome-distribution
output does not run a null model and does not include these fields. Separate
metadata sidecars are not generated.

When a co-occurrence edge table is large, the detailed output is written in bounded-memory batches as `global_edges.parquet`, with `global_edges_taxa.parquet` mapping integer taxon identifiers back to names. The summary remains a TSV intended for quick inspection. Use `--large` only when the pair count and output size are understood.
