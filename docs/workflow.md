# Recommended workflow

The sequence below is intended to keep the biological comparison visible while
the command-line options are chosen. The [terminology guide](terminology.md)
defines the CLI terms.

## 1. Define the comparison

Write down the biological question before selecting a search mode. Identify the
cohort, the background it should be compared with, the taxonomic rank, and what
would count as a meaningful effect. A changed cohort or background is a changed
analysis, even if every other option stays the same.

Use an exact `--data_release` for analyses that will be compared or published.
Omitting it is convenient for exploration because the CLI selects the registry
default, but that default can move to a later scientific snapshot.

## 2. Check the selected data

Use `metacooc search` or a small pilot run to inspect the samples and taxa
selected by the query. This catches spelling, taxonomy, overly broad metadata
terms, and thresholds that remove most of the cohort before computationally
expensive analysis begins.

Decide whether raw or aggregated Ingredients match the question. Aggregation
is useful for whole-clade questions but changes what one taxon row represents.
Likewise, presence thresholds define which profile observations contribute to
every downstream metric.

## 3. Run a fast first analysis

For association and co-occurrence, begin with FE, the default analytical
baseline. It quickly reports contingency-table evidence and effect sizes,
making it useful for checking the cohort, background, and filters. FE is not a
universally correct ecological model. Add `--compute_fisher` when an exact
2-by-2 test is useful for sparse or small tables.

Structure analysis differs: it generates empirical matrix-level null summaries
when null simulation is enabled, including with `--null_model FE`.

## 4. Inspect effect size and support

Read the summary TSV and automatic plot before opening the full table. Check
phi, both directional probabilities, support counts, and adjusted significance
together. A very specific association supported by few samples and a widespread
association with modest specificity may lead to different follow-up questions.

The default plots show significant positive-phi results at `q <= 0.10` and
label the strongest points. Use `--plot_all` when the nonsignificant context is
important, or change `--label_top_n` and the paired axis metrics for a focused
view.

## 5. Test sensitivity to ecological structure

Rerun important association or co-occurrence results with FF when preserving
both taxon prevalence and sample richness is a defensible null hypothesis.
Start with 1,000 replicates and an explicit `--nm_seed`. This gives a minimum
empirical p-value of approximately `1 / 1001`; use more replicates when the
tail probability needs finer resolution.

EF and EE are available when their different constraints match the question.
They are alternative null hypotheses, not automatic upgrades. Confirm that the
selected background still contains the intended cohort before comparing null
models.

## 6. Preserve the analysis context

Keep the package version, data release, complete command, cohort and background
queries, filters, null model, seed, replicate count, and output directory with
any reported result. Result tables repeat the compact run-level null metadata,
but they cannot record the biological reasoning behind the selected comparison.
