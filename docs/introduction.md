# Introduction

## What MetaCoOc does

MetaCoOc converts taxonomic profiles into a binary taxon-by-sample matrix and asks questions such as:

- Which taxa are associated with a metadata term or biome?
- Which taxa tend to occur together?
- Does a community show checkerboard, overlap, or nestedness patterns?
- How often is each taxon observed in a set of biomes?

The package provides both a command-line interface and Python functions. The command line is the easiest place to begin and is also the reference for the supported options.

## Cohorts and backgrounds

A **cohort** is the set of samples selected by your query. The **background** is the set of samples against which the cohort is compared. These are separate choices. For example, a soil cohort can be compared with all available samples, or with a narrower background selected by `--null_scope`.

This distinction matters most for association analysis. A taxon can look specific to a cohort simply because the background contains very different environments. Record both queries when you report a result.

## What the main analyses return

- **Association** reports taxon-to-cohort relationships. Specificity is `p_cohort_given_taxon`; sensitivity is `p_taxon_given_cohort`. Effect sizes include phi and relative-risk measures.
- **Co-occurrence** reports directed taxon-to-taxon edges. The two conditional probabilities answer different questions, so keep the source and target labels.
- **Structure** reports C-score, mean Jaccard, and NODF, with empirical null summaries when null simulation is enabled.
- **Biome distribution** reports taxon occurrence across the supplied biome mapping and does not run a null model.

The summary files are intended for first inspection. The detailed files retain the metric-level results for downstream analysis.
