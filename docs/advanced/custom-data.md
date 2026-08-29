# Custom data

Custom Ingredients are useful for a pilot analysis or for a taxonomic profile set that is not in the published registry. Supply the Ingredients directory explicitly with `--custom_ingredients` and keep a local metadata TSV if a metadata search is needed.

The profile table must contain `sample`, `taxonomy`, and `coverage`. A biome mapping must contain `accession`, `level_1`, and `level_2`. Validate a custom input with the `search` and `filter` commands before running an analysis.
