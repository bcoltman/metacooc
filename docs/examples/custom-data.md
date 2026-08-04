# Custom data example

This recipe starts from an existing custom Ingredients directory and a matching
metadata TSV. Custom profiles use the columns `sample`, `taxonomy`, and
`coverage`; creating the serialized Ingredients is a separate preparation step.
Provide the resulting directory explicitly rather than selecting a published
data release:

```bash
metacooc association \
  --custom_ingredients path/to/ingredients \
  --metadata_file path/to/metadata.tsv \
  --search_mode metadata \
  --search_string "marine" \
  --output_dir results/custom-association
```

For a metadata search, sample identifiers in the metadata file must correspond
to those in the custom Ingredients. The command produces the same association
summary, detailed table, and plot as a published-data analysis. The scientific
interpretation still depends on how the custom samples, taxonomy, and metadata
were assembled.
