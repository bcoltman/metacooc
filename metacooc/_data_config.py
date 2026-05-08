# _data_config.py

from __future__ import annotations

import os

VARIANTS = ("gtdb", "globdb")
DEFAULT_VARIANT = "gtdb"

# Map *base* release -> Zenodo record id
RELEASES = {
    "1.1.0": "18245268",
}

NO_LONGER_COMPATIBLE = ["1.0.1", "1.0.0", "0.3.0"]


def _semver_key(v: str) -> tuple[int, int, int]:
    return tuple(map(int, v.split(".")))


LATEST_BASE_VERSION = max(RELEASES.keys(), key=_semver_key)
LATEST_VERSION = f"{LATEST_BASE_VERSION}_{DEFAULT_VARIANT}"


class DataVersionError(FileNotFoundError):
    """Raised when a configured data version is missing locally."""


def available_versions() -> list[str]:
    out = []
    for base in RELEASES.keys():
        for variant in VARIANTS:
            out.append(f"{base}_{variant}")
    return sorted(out, key=lambda s: (_semver_key(s.split("_", 1)[0]), s))


def describe_data_version(version: str) -> str:
    try:
        base, variant, full_version = _parse_version_strict(version)
    except ValueError:
        return str(version)

    if variant == "globdb" and base == "1.1.0":
        return f"{full_version} (GlobDB r226)"
    if variant == "globdb":
        return f"{full_version} (GlobDB)"
    if variant == "gtdb":
        return f"{full_version} (GTDB)"
    return full_version


def local_data_versions(data_dir: str | os.PathLike | None) -> list[str]:
    if not data_dir or not os.path.isdir(data_dir):
        return []

    versions = set()
    prefixes = ("ingredients_raw_", "ingredients_aggregated_")
    for name in os.listdir(data_dir):
        if not name.endswith(".pkl"):
            continue
        for prefix in prefixes:
            if name.startswith(prefix):
                version = name[len(prefix):-4]
                try:
                    _parse_version_strict(version)
                except ValueError:
                    continue
                versions.add(version)

    return sorted(versions, key=lambda s: (_semver_key(s.split("_", 1)[0]), s))


def format_local_data_versions(data_dir: str | os.PathLike | None) -> str:
    versions = local_data_versions(data_dir)
    if not versions:
        return "none"
    return ", ".join(describe_data_version(v) for v in versions)


def missing_data_version_message(
    *,
    data_dir: str | os.PathLike | None,
    data_version: str,
    missing_path: str,
    file_kind: str,
    defaulted: bool,
) -> str:
    if defaulted:
        first = (
            f"No --data_version was specified, so metacooc defaulted to "
            f"{describe_data_version(data_version)}."
        )
    else:
        first = f"Data version {describe_data_version(data_version)} was requested."

    local_versions = local_data_versions(data_dir)
    local_display = (
        ", ".join(describe_data_version(v) for v in local_versions)
        if local_versions
        else "none"
    )
    suggestion = ""
    if local_versions:
        suggestion = (
            f"\nTo use an available local data source, rerun with "
            f"--data_version {local_versions[0]}."
        )

    return (
        f"{first}\n"
        f"Missing required {file_kind} file: {missing_path}\n"
        f"Available local data sources in {data_dir}: {local_display}."
        f"{suggestion}\n"
        f"To fetch the requested data source, run: "
        f"metacooc download --data_dir {data_dir} --data_version {data_version}"
    )


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
        "sample_to_biome": f"sample_to_biome_{base_version}.tsv",
    }
    
    download_urls = {
        filenames["ingredients_raw"]: f"{base_url}/{filenames['ingredients_raw']}.gz?download=1",
        filenames["ingredients_aggregated"]: f"{base_url}/{filenames['ingredients_aggregated']}.gz?download=1",
        filenames["sra_metadata"] + ".gz": f"{base_url}/{filenames['sra_metadata']}.gz?download=1",
        filenames["sample_to_biome"] + ".gz": f"{base_url}/{filenames['sample_to_biome']}.gz?download=1",
    }

    return filenames, download_urls
