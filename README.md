# MetaCoOc

MetaCoOc is a Python package for co-occurrence, association, and community-structure analysis of microorganisms in shotgun metagenomes.

It can:

- define cohorts from taxa, metadata terms, or biome labels;
- compare a cohort with an explicit null/background population;
- calculate directed taxon–taxon co-occurrence and taxon–cohort association;
- quantify C-score, mean Jaccard, and NODF structure; and
- write compact summary tables alongside detailed results and plots.

## Install

MetaCoOc requires Python 3.10 or newer:

```bash
python -m pip install metacooc
```

## Quick start

Download a pinned published Ingredients release and run a first biome association:

```bash
metacooc download --data_release R226_globdb_rev1
metacooc association \
  --data_release R226_globdb_rev1 \
  --search_mode biome \
  --search_string soil \
  --output_dir results/soil-association
```

Start with `global_association_summary.tsv`, then use the detailed TSV and plot when you need the full result. Metadata searches require the optional metadata download; the command-line guide explains the difference between cohort and null/background selection.

## Documentation

Read the full guide at [metacooc.readthedocs.io](https://metacooc.readthedocs.io), including:

- [installation and quick setup](https://metacooc.readthedocs.io/en/latest/quickstart.html);
- the [recommended FE-first workflow](https://metacooc.readthedocs.io/en/latest/workflow.html);
- advanced guidance on releases, search grammar, null models, outputs, and plotting; and
- worked examples for association, co-occurrence, structure, biome distribution, and custom data.

The generated [CLI reference](https://metacooc.readthedocs.io/en/latest/cli.html) lists the current options and defaults.

## Development

```bash
uv sync --extra dev
uv run pytest -q
```

Build the documentation with `uv sync --extra docs` followed by the command in `docs/local_build.md`.

## License

MetaCoOc is released under the GNU General Public License v3.0 or later.
