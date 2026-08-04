# Null models and statistical choices

## Why use a null model?

An observed association can arise because two taxa share an environment, one
taxon is common almost everywhere, or samples differ greatly in richness. A
null model asks which occurrence pattern would be expected if the relationship
of interest were absent while selected properties of the taxa-by-sample matrix
were retained. Choosing those retained properties defines the ecological null
hypothesis; a more constrained model is not automatically more appropriate.

## Analytical evidence

Association and co-occurrence always calculate contingency-table quantities,
including conditional probabilities, phi, and a chi-square test. Add
`--compute_fisher` when an exact 2-by-2 table test is useful, particularly for
sparse or small tables. Fisher's exact test conditions on the margins of that
table. It is not the analytical equivalent of shuffling the full matrix under
MetaCoOc's FE constraint.

For association and co-occurrence, `--null_model FE` is the fast analytical
baseline. MetaCoOc fixes taxon prevalence conceptually but does not generate
empirical Jaccard replicates for these two pipelines when FE is selected. This
is useful for checking the cohort, filters, and main effect sizes before a more
computationally demanding sensitivity analysis.

## Empirical matrix comparisons

FF, EF, and EE generate randomized matrices and compare the observed Jaccard
value with its distribution across those replicates:

- **FF** fixes taxon prevalence and sample richness;
- **EF** fixes sample richness while taxon prevalence can vary;
- **EE** fixes only the total number of presences.

FF is often a useful follow-up when both variation in taxon prevalence and
variation in sample richness should be explained by the null. EF and EE answer
different questions and should be used when those weaker constraints match the
intended comparison.

Start with 1,000 replicates and set `--nm_seed` so the analysis can be repeated.
With 1,000 successful replicates, the smallest attainable empirical p-value is
approximately `1 / 1001`. Increase the replicate count when tail probabilities
need finer resolution or when conclusions near a threshold are unstable. The
worker options control parallel execution; FF also exposes burn-in and
steps-per-replicate settings.

## Choosing the background

`--null_scope` restricts the samples forming the background before analysis:

- omitted: all samples, subject to the matrix thresholds;
- `biome`: samples matching `--null_biome_query`;
- `metadata`: samples matching `--null_metadata_query`;
- `taxa`: a neighbourhood around `--null_taxa_query`;
- `biome_taxa` and `metadata_taxa`: both restrictions applied in sequence.

Taxon-neighbourhood expansion is controlled by `--taxa_degree` and
`--min_shared_samples_between_taxa`. Background selection and matrix
randomization are separate decisions: the first defines which observations are
compared, while the second defines which matrix properties are preserved.

## Structure analysis is different

The structure pipeline reports C-score, mean Jaccard, and NODF for the whole
matrix. It generates empirical matrix-level summaries for FE as well as FF,
EF, and EE. Do not transfer the association/co-occurrence interpretation of FE
directly to structure.
