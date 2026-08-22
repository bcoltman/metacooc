# API reference

The Python API is useful when an existing workflow already holds data in
memory, or when you want to combine MetaCoOc with other Python analysis. The
CLI remains the most direct interface for a complete search/filter/analyse
workflow. Python callers should pass an explicit release or Ingredients path;
the API does not silently select a published release.

The generated pages below document the high-level public entry points. Private
helpers and internal metric/array functions are intentionally omitted.

## Data and search

```{eval-rst}
.. autofunction:: metacooc.download.download_data
.. autofunction:: metacooc.pantry.load_ingredients
.. autoclass:: metacooc.pantry.Ingredients
   :no-members:
.. autofunction:: metacooc.search.search_data_obj
.. autofunction:: metacooc.search.search_data
```

## Analysis pipelines

```{eval-rst}
.. autofunction:: metacooc.analysis.association_obj
.. autofunction:: metacooc.analysis.association
.. autofunction:: metacooc.analysis.cooccurrence_obj
.. autofunction:: metacooc.analysis.cooccurrence
.. autofunction:: metacooc.structure.structure_obj
.. autofunction:: metacooc.structure.structure
.. automethod:: metacooc.pantry.Ingredients.biome_distribution
```

## Plotting

```{eval-rst}
.. autofunction:: metacooc.plot.plot_analysis_obj
.. autofunction:: metacooc.plot.plot_analysis
```

See the [recommended workflow](workflow.md) and [output guide](advanced/outputs.md)
for how these functions map onto the CLI pipelines and result files.
