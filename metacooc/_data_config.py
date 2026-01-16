# _data_config.py

from __future__ import annotations

VARIANTS = ("gtdb", "globdb")
DEFAULT_VARIANT = "gtdb"

# Map *base* release -> Zenodo record id
RELEASES = {
    "1.1.0": "****",
    # "1.0.1": "16030719",
    # ...
}

NO_LONGER_COMPATIBLE = ["1.0.1", "1.0.0", "0.3.0"]


def _semver_key(v: str) -> tuple[int, int, int]:
    return tuple(map(int, v.split(".")))


LATEST_BASE_VERSION = max(RELEASES.keys(), key=_semver_key)
LATEST_VERSION = f"{LATEST_BASE_VERSION}_{DEFAULT_VARIANT}"


def available_versions() -> list[str]:
    out = []
    for base in RELEASES.keys():
        for variant in VARIANTS:
            out.append(f"{base}_{variant}")
    return sorted(out, key=lambda s: (_semver_key(s.split("_", 1)[0]), s))


def _parse_version_strict(version: str) -> tuple[str, str, str]:
    """
    Require a variant suffix: <base>_<variant>, e.g. 1.1.0_gtdb.
    Returns (base_version, variant, full_version).
    """
    if "_" not in version:
        raise ValueError(
            f"Version '{version}' must include a variant suffix: "
            f"{', '.join(VARIANTS)}. Example: '{LATEST_VERSION}'."
        )
        
    base, variant = version.split("_", 1)
    
    if variant not in VARIANTS:
        raise ValueError(
            f"Unknown variant '{variant}'. Allowed variants: {', '.join(VARIANTS)}. "
            f"Example: '{base}_{DEFAULT_VARIANT}'."
        )
        
    return base, variant, f"{base}_{variant}"


def get_file_info(version: str):
    base_version, variant, full_version = _parse_version_strict(version)
    
    if base_version in NO_LONGER_COMPATIBLE:
        raise ValueError(
            f"Version '{base_version}' is no longer compatible. "
            f"Available versions: {', '.join(available_versions())}"
        )
    if base_version not in RELEASES:
        raise ValueError(
            f"Version '{base_version}' is not available. "
            f"Available versions: {', '.join(available_versions())}"
        )
        
    record_id = RELEASES[base_version]
    base_url = f"https://zenodo.org/records/{record_id}/files"
    
    filenames = {
        # variant-specific
        "ingredients_raw": f"ingredients_raw_{full_version}.pkl",
        "ingredients_aggregated": f"ingredients_aggregated_{full_version}.pkl",
        # common (base only)
        "sra_metadata": f"sra_metadata_{base_version}.tsv",
        "sample_to_biome": f"sample_to_biome_{base_version}.csv",
    }
    
    download_urls = {
        filenames["ingredients_raw"]: f"{base_url}/{filenames['ingredients_raw']}.gz?download=1",
        filenames["ingredients_aggregated"]: f"{base_url}/{filenames['ingredients_aggregated']}.gz?download=1",
        filenames["sra_metadata"] + ".gz": f"{base_url}/{filenames['sra_metadata']}.gz?download=1",
        filenames["sample_to_biome"] + ".gz": f"{base_url}/{filenames['sample_to_biome']}.gz?download=1",
    }

    return filenames, download_urls
