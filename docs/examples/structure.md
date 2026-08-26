# Community structure example

This recipe asks whether a soil cohort has more or less segregation,
similarity, or nestedness than matrices generated under an FE null model.

```bash
metacooc structure \
  --data_release R226_globdb_rev1 \
  --search_mode biome \
  --search_string soil \
  --null_model FE \
  --nm_n_reps 1000 \
  --output_dir results/soil-structure
```

The focused output is `global_structure.tsv`. It reports observed C-score,
mean Jaccard, and NODF values together with empirical null summaries. Interpret
each metric against its own null distribution: the three statistics measure
different aspects of community structure and are not directly comparable in
magnitude. Structure is not automatically plotted.
