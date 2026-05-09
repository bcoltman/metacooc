from __future__ import annotations

import os
from pathlib import Path

from platformdirs import user_data_dir


DATA_DIR_ENV_VAR = "METACOOC_DATA_DIR"


def default_data_dir() -> Path:
    configured = os.environ.get(DATA_DIR_ENV_VAR)
    if configured:
        return Path(configured).expanduser()
    return Path(user_data_dir("metacooc", appauthor=False)) / "data"
