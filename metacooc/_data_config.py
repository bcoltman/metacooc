"""Configuration and validation for downloadable MetaCoOc data releases."""

from __future__ import annotations

from dataclasses import dataclass
import os
import re


VARIANTS = ("gtdb", "globdb")
DEFAULT_VARIANT = "gtdb"
_RELEASE_RE = re.compile(r"^R(?P<number>[0-9]+)_(?P<variant>gtdb|globdb)$")


@dataclass(frozen=True)
class ReleaseSpec:
    """Published artifacts belonging to one Sandpiper release."""

    zenodo_record: str | None
    ingredients_raw: str = "ingredients_raw_{release}_{variant}"
    ingredients_aggregated: str = "ingredients_aggregated_{release}_{variant}"
    sra_metadata: str = "sra_metadata_{release}.tsv"


# Keep this registry explicit: release identifiers describe the underlying
# database release, while artifact names and Zenodo records are publication
# details. The temporary ``None`` value only needs replacing when the R226
# archive is published; releases do not otherwise have a separate publication
# state.
RELEASES: dict[str, ReleaseSpec] = {
    "R226": ReleaseSpec(zenodo_record=None),
}


class DataReleaseError(FileNotFoundError):
    """Raised when a requested data release is invalid or unavailable."""


def _release_key(release: str) -> int:
    return int(release[1:])


def _parse_release_strict(data_release: str) -> tuple[str, str, str]:
    """Return ``(base_release, variant, canonical_release)``."""
    if not isinstance(data_release, str):
        raise ValueError("Data release must be a string such as 'R226_gtdb'.")

    match = _RELEASE_RE.fullmatch(data_release)
    if match is None:
        raise ValueError(
            f"Data release '{data_release}' must match R<number>_<variant>, "
            f"for example 'R226_{DEFAULT_VARIANT}'. Allowed variants: "
            f"{', '.join(VARIANTS)}."
        )

    base = f"R{match.group('number')}"
    variant = match.group("variant")
    return base, variant, f"{base}_{variant}"


def available_releases() -> list[str]:
    return [
        f"{base}_{variant}"
        for base in sorted(RELEASES, key=_release_key)
        for variant in VARIANTS
    ]


LATEST_BASE_RELEASE = max(RELEASES, key=_release_key)
LATEST_DATA_RELEASE = f"{LATEST_BASE_RELEASE}_{DEFAULT_VARIANT}"


def describe_data_release(data_release: str) -> str:
    try:
        base, variant, canonical = _parse_release_strict(data_release)
    except ValueError:
        return str(data_release)
    return f"{canonical} ({variant.upper()}, database release {base})"


def local_data_releases(data_dir: str | os.PathLike | None) -> list[str]:
    if not data_dir or not os.path.isdir(data_dir):
        return []

    releases = set()
    prefixes = ("ingredients_raw_", "ingredients_aggregated_")
    for name in os.listdir(data_dir):
        path = os.path.join(data_dir, name)
        if not os.path.isdir(path):
            continue
        for prefix in prefixes:
            if name.startswith(prefix):
                release = name[len(prefix):]
                try:
                    _parse_release_strict(release)
                except ValueError:
                    continue
                releases.add(release)

    return sorted(
        releases,
        key=lambda value: (_release_key(value.split("_", 1)[0]), value),
    )


def format_local_data_releases(data_dir: str | os.PathLike | None) -> str:
    releases = local_data_releases(data_dir)
    if not releases:
        return "none"
    return ", ".join(describe_data_release(value) for value in releases)


def missing_data_release_message(
    *,
    data_dir: str | os.PathLike | None,
    data_release: str,
    missing_path: str,
    file_kind: str,
    defaulted: bool,
) -> str:
    if defaulted:
        first = (
            f"No --data-release was specified, so metacooc defaulted to "
            f"{describe_data_release(data_release)}."
        )
    else:
        first = f"Data release {describe_data_release(data_release)} was requested."

    local_releases = local_data_releases(data_dir)
    local_display = (
        ", ".join(describe_data_release(value) for value in local_releases)
        if local_releases
        else "none"
    )
    suggestion = ""
    if local_releases:
        suggestion = (
            f"\nTo use an available local data source, rerun with "
            f"--data-release {local_releases[0]}."
        )

    return (
        f"{first}\n"
        f"Missing required {file_kind} file: {missing_path}\n"
        f"Available local data sources in {data_dir}: {local_display}."
        f"{suggestion}\n"
        f"To fetch the requested data source, run: "
        f"metacooc download --data_dir {data_dir} --data-release {data_release}"
    )


def get_file_info(data_release: str):
    try:
        base_release, variant, canonical = _parse_release_strict(data_release)
    except ValueError as exc:
        raise DataReleaseError(str(exc)) from exc
    if base_release not in RELEASES:
        raise DataReleaseError(
            f"Data release '{canonical}' is not available. "
            f"Available releases: {', '.join(available_releases()) or 'none'}"
        )

    spec = RELEASES[base_release]

    filenames = {
        "ingredients_raw": spec.ingredients_raw.format(
            release=base_release, variant=variant
        ),
        "ingredients_aggregated": spec.ingredients_aggregated.format(
            release=base_release, variant=variant
        ),
        "sra_metadata": spec.sra_metadata.format(
            release=base_release, variant=variant
        ),
    }
    base_url = f"https://zenodo.org/records/{spec.zenodo_record}/files"
    download_urls = {
        filenames["ingredients_raw"] + ".tar.gz": (
            f"{base_url}/{filenames['ingredients_raw']}.tar.gz?download=1"
        ),
        filenames["ingredients_aggregated"] + ".tar.gz": (
            f"{base_url}/{filenames['ingredients_aggregated']}.tar.gz?download=1"
        ),
        filenames["sra_metadata"] + ".gz": (
            f"{base_url}/{filenames['sra_metadata']}.gz?download=1"
        ),
    }
    return filenames, download_urls
