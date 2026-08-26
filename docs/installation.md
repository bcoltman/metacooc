# Installation

MetaCoOc requires Python 3.10 or newer. PyPI and Bioconda provide alternative
ways to install the same command-line tool.

## Bioconda

Create a dedicated environment, install MetaCoOc from Bioconda, and activate
the environment:

```bash
conda create --name metacooc \
  --channel conda-forge \
  --channel bioconda \
  metacooc
conda activate metacooc
```

## PyPI

To install MetaCoOc into an existing Python environment:

```bash
python -m pip install metacooc
```

## Source checkout

For an editable installation from the source repository:

```bash
git clone https://github.com/bcoltman/metacooc.git
cd metacooc
python -m pip install -e .
```

Whichever route you use, check that the command-line interface is available:

```bash
metacooc --help
```

For development or documentation work, use an isolated environment. The
documentation dependencies are installed with the `docs` extra:

```bash
uv venv /tmp/metacooc-docs
UV_CACHE_DIR=/tmp/metacooc-uv-cache \
  uv pip install --python /tmp/metacooc-docs/bin/python -e '.[docs]'
```

Continue with the [quick setup](quickstart.md), or build the site locally with
the instructions in [Local build](local_build.md).
