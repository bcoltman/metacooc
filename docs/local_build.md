# Build the documentation locally

From a checkout of MetaCoOc, install the package with its documentation extra:

```bash
uv sync --extra docs
uv run sphinx-build -W -b html docs docs/_build/html
```

Open `docs/_build/html/index.html` in a browser. `-W` treats warnings as errors, which is also how Read the Docs builds the site.
