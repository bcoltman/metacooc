# MetaCoOc

MetaCoOc helps microbial ecologists investigate which microorganisms are
associated with an environment and which tend to occur together. It searches
taxonomic profiles from more than 700,000 public shotgun metagenomes, covering
more than 300,000 species. Samples can be selected through a biome ontology, a
shared metadata term, or the presence of a focal taxon. MetaCoOc then compares
prevalence within that cohort with a defined background and reports ecological
effect sizes and statistical evidence.

These archive-scale patterns can recover known environmental preferences and
prioritise candidate ecological relationships for closer study. They provide
evidence of repeated occurrence, not proof of interaction or mechanism.

This guide is written for microbial ecologists and computational biologists.
The quick setup introduces the complete workflow and its main decisions. The
advanced guide explains data aggregation, filtering, null models, plotting,
and larger analyses in more detail.

```{toctree}
:maxdepth: 2
:caption: Start here

introduction
terminology
installation
quickstart
workflow
```

```{toctree}
:maxdepth: 2
:caption: Advanced

advanced/index
```

```{toctree}
:maxdepth: 2
:caption: Worked examples

examples/index
```

```{toctree}
:maxdepth: 2
:caption: Reference

cli
api
```

```{toctree}
:hidden:

local_build
```
