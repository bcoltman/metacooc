# Biome distribution example

This recipe creates a broad view of where taxa have been observed across the
biomes stored in a published Ingredients release.

```bash
metacooc biome_distribution \
  --data_release R226_globdb_rev1 \
  --output_dir results/biome-distribution
```

The pipeline writes `taxa_biome_distribution.tsv`, with taxa as rows and biome
labels as columns. Values describe occurrence across the archived samples; they
do not correct for differences in sampling effort among biomes. This pipeline
does not run a null model.
