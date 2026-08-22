# Python API

The high-level Python entry points expose the same association, co-occurrence,
structure, and biome-distribution operations as the CLI. Data-loading and
search functions require an explicit published release or custom source;
analysis functions receive an `Ingredients` object that you load yourself.
Python APIs do not silently resolve a registry default.

Use the CLI reference for current argument names and the package docstrings for function signatures. The CLI is the more stable interface for reproducible workflows; lower-level functions are useful when an existing Python pipeline already holds an `Ingredients` object in memory.
