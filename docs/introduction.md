# Introduction

## Why analyse occurrence at this scale?

Individual metagenomic studies provide detailed views of particular sites,
hosts, or experiments. Public archives provide a complementary view: repeated
observations across many environments and studies can reveal which microbial
lineages are consistently associated with an environment and which taxa are
often found together.

MetaCoOc makes more than 700,000 public shotgun metagenomes and over 300,000
species available for this type of screening. The aim is to identify robust
environmental preferences, community patterns, and candidate ecological
relationships that can be followed up with targeted datasets or experiments.
An archive-scale association is useful evidence of occurrence, but it does not
establish direct interaction, causality, or independence among samples.

## What MetaCoOc does

MetaCoOc works with taxonomic profiles generated from shotgun metagenomes.
Each profile records which microorganisms were found in a sample and their
coverage. MetaCoOc turns these profiles into a sparse taxa-by-sample matrix and
uses presence or absence to ask questions such as:

- Which taxa are associated with a biome or metadata-defined environment?
- Which taxa tend to occur with a focal microorganism?
- Does a community show segregation, overlap, or nestedness?
- Across which biomes has a taxon been observed?

The command-line workflows connect sample selection, filtering, analysis, and
plotting. Python functions expose the same main analyses for integration into
other pipelines. The [terminology guide](terminology.md) defines the terms used
by the CLI, and the [command reference](cli.md) lists every option.

## Cohorts and backgrounds

A **cohort** is the set of samples selected by the biological query. A
**background** is the population against which that cohort or its taxa are
compared. These are separate choices because the appropriate comparison
depends on the question.

For example, a bee-gut cohort could be compared with the complete archive or
with a narrower host-associated background. The first asks whether a taxon is
specific relative to all sampled environments; the second asks whether it is
specific within a more similar ecological setting. A broad background may
highlight large environmental differences, while a narrow background may
remove useful contrast.

The cohort must be represented within the chosen background. If the background
excludes cohort samples, the usable cohort can shrink or become empty. Record
both selections when reporting an analysis.

## What the analyses return

- **Association** compares the prevalence of each taxon in a cohort with its
  prevalence in the background. Its directional probabilities describe
  specificity and sensitivity.
- **Co-occurrence** reports taxon pairs observed in the same samples. Its
  conditional probabilities are directional: the probability of the target
  given the source need not equal the reverse.
- **Structure** reports C-score, mean Jaccard, and NODF for the community
  matrix, describing segregation, overlap, and nestedness.
- **Biome distribution** reports how often taxa occur across supplied biome
  labels and does not run a null model.

Summary files and automatic plots support initial inspection. Detailed files
retain the full metric set and counts for downstream analysis.
