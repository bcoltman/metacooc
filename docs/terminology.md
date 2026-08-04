# CLI terminology

MetaCoOc uses the following terms consistently across its commands, result
files, and Python functions.

## Data and samples

**Ingredients** are MetaCoOc's prepared input data: sparse taxonomic profiles,
sample identifiers, and biome mappings. A **data release** is an immutable
published snapshot of these files, such as `R226_globdb_rev1`. It is separate
from the MetaCoOc package version and the Ingredients `format_version`.

SingleM profiles report **unfilled coverage** at the rank where a lineage was
assigned. Coverage assigned to a species is not automatically counted in the
coverage of its genus. **Raw Ingredients** retain those original rows.
**Aggregated Ingredients** retain the raw rows and add marked higher-rank rows
that sum descendant coverage. For example, an aggregated genus combines its
species and lineages assigned only to that genus. Use `--aggregated` when the
whole clade, rather than the original assignment rank, matches the question.
See [Data releases and Ingredients](advanced/data.md) for the consequences of
that choice.

A **sample** is one metagenome. A **taxon** is a named microbial lineage in a
taxonomic profile. A **biome** is an environment label from the supplied biome
ontology. **Metadata** is descriptive information submitted with a sample,
such as its environment, host, or collection context.

## Cohort selection

`--search_mode` determines how `--search_string` selects samples and, for a
focal search, which taxon pairs are reported.

| Search mode | What it selects | Analysis workflows |
| --- | --- | --- |
| `biome` | Samples carrying a matching biome label | Association, co-occurrence, structure |
| `metadata` | Samples containing a word or phrase in their metadata | Association, co-occurrence, structure |
| `taxa_context` | Samples containing the requested taxon combination | Association, co-occurrence, structure |
| `focal_taxa` | Samples containing a focal taxon, with results restricted to pairs involving it | Co-occurrence only |

Biome, metadata, and taxa-context searches define a cohort whose complete taxon
set is passed to the selected analysis. Focal-taxa searches define a cohort but
restrict co-occurrence reporting to edges involving the focal taxon. The
optional `LHS -> RHS` syntax further limits reported targets and is specific to
focal co-occurrence.

The **cohort** starts as the samples selected by the search. The **background**,
also called the **null background**, is the comparison population selected by
`--null_scope`. The final analysis can use only cohort samples present in that
background. Taxa-based scopes can also restrict the taxon universe.

## Filtering and presence

`--min_sample_count` keeps taxa observed in enough samples.
`--min_taxa_count` keeps samples containing enough taxa. `--filter_rank` keeps
taxa whose terminal rank exactly matches the selected rank.

MetaCoOc analyses use presence or absence. `--min_coverage` and
`--min_relative_abundance` set the minimum coverage or relative-abundance value
that counts as present. Increasing a threshold removes low-abundance
observations; whether that is more appropriate depends on the profiling noise
and biological question.

## Analysis and null-model terms

An **association** compares the prevalence of a taxon in a cohort with its
prevalence in the background. A **co-occurrence** compares the prevalence and
overlap of two taxa across the same samples. A **null model** describes the
pattern expected if the relationship of interest were absent while specified
properties of the taxa-by-sample matrix were retained.

Association and co-occurrence report analytical contingency-table statistics,
including chi-square by default and optional Fisher testing. FF, EF, and EE can
also generate empirical Jaccard comparisons from randomized matrices. FE is the
fast analytical baseline for these two workflows, while structure uses FE
empirically. See [Null models and statistical choices](advanced/null-models.md)
before interpreting these labels across pipelines.

## Reading results

The **summary TSV** is the smaller table containing the metrics intended for
initial inspection. The **detailed TSV** or Parquet file retains the full set of
calculated metrics. The **plot** is an automatic visual summary controlled by
options such as `--plot_all` and `--label_top_n`.

The two association probabilities are:

- `p_cohort_given_taxon`: among samples containing the taxon, the proportion
  belonging to the cohort. This is the specificity direction.
- `p_taxon_given_cohort`: among cohort samples, the proportion containing the
  taxon. This is the sensitivity direction.

**Phi** is a signed measure of association. Positive values indicate that the
taxon and cohort, or the two taxa, occur together more often than expected from
their separate frequencies. Negative values indicate avoidance or segregation.
A **q-value** is a multiple-testing-adjusted significance value. It does not
replace checking effect size, support counts, and biological context.
