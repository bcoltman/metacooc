# *Methanocatella smithii* co-occurrence

GTDB R226 places the organism historically known as *Methanobrevibacter
smithii* under the name *Methanocatella smithii*. This workflow first retrieves
its species-level co-occurrence neighbourhood in a human-gut background, then
examines its association with *Christensenella* at species and aggregated-genus
levels.

Download the pinned release and its metadata bundle once:

```bash
metacooc download \
  --data_release R226_globdb_rev1 \
  --include_metadata
```

## 1. Retrieve the neighbourhood

This broad search reports species co-occurring with *M. smithii*. It uses a
0.10 conditional-probability cutoff to limit the neighbourhood and its output.

```bash
metacooc cooccurrence \
  --data_release R226_globdb_rev1 \
  --search_mode focal_taxa \
  --search_string "s__Methanocatella smithii" \
  --min_taxa_count 10 \
  --taxa_count_rank species \
  --min_sample_count 50 \
  --filter_rank species \
  --min_conditional_probability 0.10 \
  --null_model FE \
  --q_threshold 0.10 \
  --null_scope metadata \
  --null_metadata_query "human gut metagenome" \
  --max_pairs 150000 \
  --label_top_n 5 \
  --output_dir results/examples/smithii-cooccurrence/neighbourhood
```

```{figure} ../../examples/poster_plots/smithii_human_gut_cooccurrence.png
:alt: Methanocatella smithii co-occurrence plot
:width: 90%

The species-level *M. smithii* neighbourhood in the human-gut background.
```

## 2. Test *Christensenella* species

The `LHS -> RHS` form uses *Christensenella* species to define the focal cohort
and restricts the reported targets to *M. smithii*. A zero conditional-
probability cutoff retains the targeted result if it passes the other filters.

```bash
metacooc cooccurrence \
  --data_release R226_globdb_rev1 \
  --search_mode focal_taxa \
  --search_string "g__Christensenella -> s__Methanocatella smithii" \
  --min_taxa_count 10 \
  --taxa_count_rank species \
  --min_sample_count 50 \
  --filter_rank species \
  --min_conditional_probability 0.0 \
  --null_model FE \
  --q_threshold 0.10 \
  --null_scope metadata \
  --null_metadata_query "human gut metagenome" \
  --max_pairs 5000000 \
  --label_top_n 3 \
  --output_dir results/examples/smithii-cooccurrence/christensenella-species
```

## 3. Test the aggregated genus

The last search treats the aggregated *Christensenella* genus as the focal
taxon rather than testing its species separately.

```bash
metacooc cooccurrence \
  --data_release R226_globdb_rev1 \
  --aggregated \
  --search_mode focal_taxa \
  --search_string "g__Christensenella AGGREGATED -> s__Methanocatella smithii" \
  --min_taxa_count 10 \
  --taxa_count_rank species \
  --min_sample_count 50 \
  --min_conditional_probability 0.0 \
  --null_model FE \
  --q_threshold 0.10 \
  --null_scope metadata \
  --null_metadata_query "human gut metagenome" \
  --max_pairs 500000 \
  --label_top_n 1 \
  --output_dir results/examples/smithii-cooccurrence/christensenella-aggregated
```

## Reading the results

Each stage writes a reduced edge summary, detailed edges, nodes, and a plot.
Start with `global_edges_summary.tsv`, then use `global_edges.tsv` for the full
metrics. Very large detailed edge tables are written as Parquet instead. The
two conditional probabilities are directional, so read them in relation to
the source and target columns rather than as interchangeable values.

These occurrence patterns do not by themselves demonstrate physical
association, metabolite transfer, predation, or sample independence. The
motivating experimental study is Ruaud et al. (2020), [mBio
11:e03235-19](https://doi.org/10.1128/mBio.03235-19).
