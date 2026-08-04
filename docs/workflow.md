# Recommended workflow

## 1. Pin the data and define the question

Use an exact `--data_release` for analyses that will be compared or published. Write down the cohort query, the filters, and the null/background query. A changed cohort or background is a changed analysis, even if the package version is unchanged.

Start with a small search or a custom Ingredients directory to check that the selected samples and taxa are sensible. The `search` command and the analysis progress messages are useful for this check.

## 2. Run a fast FE baseline

For association and co-occurrence, FE is the default analytical baseline. It preserves taxon row totals and uses the analytical chi-square calculation, so it is quick and useful for finding problems with the query or thresholds. It does not mean that FE is the universally correct ecological null.

For sparse or small contingency tables, add `--compute_fisher` and compare the Fisher result with the chi-square result.

Structure analysis is different: it computes empirical matrix-level null summaries when null simulation is enabled, including when `--null_model FE` is selected.

## 3. Inspect the compact outputs

Read the summary TSV and the automatic plot before opening the full table. Association summaries contain the taxon, phi, the two directional probabilities, support counts, and adjusted analytical significance. Co-occurrence summaries contain the directed source/target pair, both conditional probabilities, shared counts, phi, and adjusted significance.

The default plots select significant positive-phi results at `q <= 0.10` and label the strongest points. Use `--plot_all` when you need the nonsignificant context, or `--label_top_n` and paired `--x_metric`/`--y_metric` options for a focused view.

## 4. Check a structure-aware null

Rerun the key analysis with FF when fixing both taxon prevalence and sample richness is a defensible representation of the sampling process. Set `--nm_seed` and choose an explicit `--nm_n_reps` so that the run can be repeated. EF and EE are available when their weaker or different constraints match the biological question; they are not automatic upgrades over FE.

For association, keep the cohort query separate from the null scope. `--null_scope` changes the background used for the null calculation; it does not redefine the cohort itself.

## 5. Preserve the context

Keep the release identifier, command, seed, null model, replicate count, and output directory alongside any result you report. The result files repeat four compact run-level fields (`null_model`, `null_replicates`, `null_replicates_failed`, and `null_seed`) so that a table remains self-describing.
