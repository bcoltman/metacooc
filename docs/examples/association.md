# Association examples

These examples ask which species are associated with two metadata-defined
cohorts. Both use the same pinned GlobDB release and require its optional
metadata bundle:

```bash
metacooc download \
  --data_release R226_globdb_rev1 \
  --include_metadata
```

## Bee gut

The first cohort contains samples whose metadata includes `bee gut`. The FE
baseline compares the occurrence of each species inside and outside that
cohort. The command retains all positive conditional probabilities in the
tables and labels the six strongest significant results in the plot.

```bash
metacooc association \
  --data_release R226_globdb_rev1 \
  --search_mode metadata \
  --search_string "bee gut" \
  --filter_rank species \
  --min_conditional_probability 0.0 \
  --null_model FE \
  --q_threshold 0.10 \
  --label_top_n 6 \
  --output_dir results/examples/bee-gut-association
```

```{figure} ../../examples/poster_plots/bee_gut_association.png
:alt: Bee-gut association plot
:width: 90%

Species associated with the bee-gut metadata cohort.
```

## Human gut

The second cohort uses the metadata term `human gut metagenome`. Additional
count filters require at least 50 taxa in a sample and at least 20 samples for
a taxon before association is tested.

```bash
metacooc association \
  --data_release R226_globdb_rev1 \
  --search_mode metadata \
  --search_string "human gut metagenome" \
  --filter_rank species \
  --min_taxa_count 50 \
  --min_sample_count 20 \
  --min_conditional_probability 0.0 \
  --null_model FE \
  --q_threshold 0.10 \
  --label_top_n 6 \
  --output_dir results/examples/human-gut-association
```

```{figure} ../../examples/poster_plots/human_gut_association.png
:alt: Human-gut association plot
:width: 90%

Species associated with the human-gut metadata cohort after count filtering.
```

## Reading the results

Each output directory contains:

```text
global_association_summary.tsv
global_association.tsv
global_association_plot.png
```

Start with the summary to compare phi, conditional probabilities, support, and
adjusted significance. Use the detailed table when the additional metrics are
needed. A metadata query defines a cohort from matching archive records; it
does not establish that every sample has the same experimental design or host
context. These are occurrence associations and do not demonstrate a direct
interaction or mechanism.
