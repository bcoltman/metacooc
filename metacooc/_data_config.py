"""Configuration and resolution for downloadable MetaCoOc data releases."""

from __future__ import annotations

from dataclasses import dataclass
from importlib import resources
import json
import os
from pathlib import Path
import re
import tempfile
import time
import warnings

import requests

from metacooc._paths import default_cache_dir


INGREDIENTS_FORMAT_VERSION = 1
REGISTRY_FORMAT_VERSION = 1
REGISTRY_FILENAME = "data_registry.v1.json"
REGISTRY_URL = (
    "https://raw.githubusercontent.com/bcoltman/metacooc/"
    f"main/metacooc/{REGISTRY_FILENAME}"
)
REGISTRY_CACHE_SECONDS = 24 * 60 * 60
REGISTRY_TIMEOUT = (2, 5)
VARIANTS = ("gtdb", "globdb")
DEFAULT_VARIANT = "gtdb"

_RELEASE_RE = re.compile(
    r"^R(?P<number>[0-9]+)_(?P<variant>gtdb|globdb)_rev(?P<revision>[1-9][0-9]*)$"
)
_LOCAL_RELEASE_RE = re.compile(
    r"^ingredients_(?:raw|aggregated)_"
    r"(?P<release>R[0-9]+_(?:gtdb|globdb)_rev[1-9][0-9]*)_"
    r"format(?P<format>[1-9][0-9]*)$"
)
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")

_memory_registry: dict | None = None
_memory_registry_expires = 0.0
_refresh_warning_emitted = False


class DataReleaseError(FileNotFoundError):
    """Raised when a requested data release is invalid or unavailable."""


class IngredientsFormatError(ValueError):
    """Raised when an Ingredients directory uses an unsupported format."""


@dataclass(frozen=True)
class ParsedDataRelease:
    base: str
    variant: str
    revision: int
    canonical: str


@dataclass(frozen=True)
class Artifact:
    """One immutable artifact published on Zenodo."""

    filename: str
    zenodo_record: str
    sha256: str

    @property
    def url(self) -> str:
        return (
            f"https://zenodo.org/records/{self.zenodo_record}/files/"
            f"{self.filename}?download=1"
        )


@dataclass(frozen=True)
class ResolvedRelease:
    """Artifacts selected for one exact release and Ingredients format."""

    data_release: str
    format_version: int
    current: bool
    artifacts: dict[str, Artifact]

    @property
    def filenames(self) -> dict[str, str]:
        return {
            key: _installed_name(artifact.filename)
            for key, artifact in self.artifacts.items()
        }


def parse_data_release(data_release: str) -> ParsedDataRelease:
    """Parse an exact published-data selector."""
    if not isinstance(data_release, str):
        raise ValueError(
            "Data release must be a string such as 'R226_gtdb_rev1'."
        )

    match = _RELEASE_RE.fullmatch(data_release)
    if match is None:
        raise ValueError(
            f"Data release '{data_release}' must match "
            "R<number>_<variant>_rev<number>, for example "
            f"'R226_{DEFAULT_VARIANT}_rev1'. Allowed variants: "
            f"{', '.join(VARIANTS)}."
        )

    base = f"R{match.group('number')}"
    variant = match.group("variant")
    revision = int(match.group("revision"))
    canonical = f"{base}_{variant}_rev{revision}"
    return ParsedDataRelease(base, variant, revision, canonical)


def is_canonical_data_release(value: object) -> bool:
    if not isinstance(value, str):
        return False
    return _RELEASE_RE.fullmatch(value) is not None


def _release_number(base: str) -> int:
    return int(base[1:])


def _installed_name(filename: str) -> str:
    if filename.endswith(".tar.gz"):
        return filename[:-7]
    if filename.endswith(".gz"):
        return filename[:-3]
    return filename


def _validate_artifact(
    value: object,
    *,
    context: str,
    zenodo_record: str,
    expected_filename: str,
) -> None:
    if not isinstance(value, dict):
        raise ValueError(f"{context} must be an object.")
    filename = value.get("filename")
    checksum = value.get("sha256")
    if filename != expected_filename:
        raise ValueError(
            f"{context}.filename must be {expected_filename!r}, got {filename!r}."
        )
    if Path(filename).name != filename:
        raise ValueError(f"{context}.filename must be a basename.")
    if not isinstance(checksum, str) or _SHA256_RE.fullmatch(checksum) is None:
        raise ValueError(f"{context}.sha256 must be 64 lowercase hexadecimal characters.")
    if not zenodo_record.isdigit() or int(zenodo_record) <= 0:
        raise ValueError(f"{context} has an invalid Zenodo record ID.")


def validate_registry(registry: object) -> dict:
    """Validate registry structure and deterministic official filenames."""
    if not isinstance(registry, dict):
        raise ValueError("Data registry must be a JSON object.")
    if registry.get("registry_format_version") != REGISTRY_FORMAT_VERSION:
        raise ValueError(
            "Unsupported data registry format version "
            f"{registry.get('registry_format_version')!r}; expected "
            f"{REGISTRY_FORMAT_VERSION}."
        )
    releases = registry.get("releases")
    if not isinstance(releases, dict):
        raise ValueError("Data registry 'releases' must be an object.")
    if not isinstance(registry.get("updated"), str):
        raise ValueError("Data registry 'updated' must be a string.")

    for base, release in releases.items():
        if not isinstance(base, str) or re.fullmatch(r"R[0-9]+", base) is None:
            raise ValueError(f"Invalid base release {base!r}.")
        if not isinstance(release, dict):
            raise ValueError(f"Release {base} must be an object.")
        revisions = release.get("revisions")
        current = release.get("current_revision")
        if not isinstance(revisions, dict) or not revisions:
            raise ValueError(f"Release {base} must contain published revisions.")
        if type(current) is not int or str(current) not in revisions:
            raise ValueError(f"Release {base} has an invalid current_revision.")

        for revision_key, revision in revisions.items():
            if re.fullmatch(r"[1-9][0-9]*", revision_key) is None:
                raise ValueError(f"Release {base} has invalid revision {revision_key!r}.")
            if not isinstance(revision, dict):
                raise ValueError(f"Release {base} revision {revision_key} must be an object.")
            revision_number = int(revision_key)

            support = revision.get("support")
            if not isinstance(support, dict):
                raise ValueError(f"Release {base} rev{revision_key} lacks support files.")
            support_record = support.get("zenodo_record")
            if not isinstance(support_record, str):
                raise ValueError("Zenodo record IDs must be strings.")
            for key, stem in (
                ("sra_metadata", "sra_metadata"),
                ("sample_to_biome", "sample_to_biome"),
            ):
                _validate_artifact(
                    support.get(key),
                    context=f"{base} rev{revision_key} {key}",
                    zenodo_record=support_record,
                    expected_filename=f"{stem}_{base}_rev{revision_number}.tsv.gz",
                )

            formats = revision.get("ingredients_formats")
            if not isinstance(formats, dict) or not formats:
                raise ValueError(f"Release {base} rev{revision_key} has no Ingredients formats.")
            for format_key, format_entry in formats.items():
                if re.fullmatch(r"[1-9][0-9]*", format_key) is None:
                    raise ValueError(
                        f"Release {base} rev{revision_key} has invalid Ingredients "
                        f"format {format_key!r}."
                    )
                if not isinstance(format_entry, dict):
                    raise ValueError("Ingredients format entries must be objects.")
                format_record = format_entry.get("zenodo_record")
                if not isinstance(format_record, str):
                    raise ValueError("Zenodo record IDs must be strings.")
                for variant in VARIANTS:
                    variant_entry = format_entry.get(variant)
                    if not isinstance(variant_entry, dict):
                        raise ValueError(
                            f"Release {base} rev{revision_key} format {format_key} "
                            f"lacks {variant}."
                        )
                    for kind in ("raw", "aggregated"):
                        expected = (
                            f"ingredients_{kind}_{base}_{variant}_rev{revision_number}_"
                            f"format{int(format_key)}.tar.gz"
                        )
                        _validate_artifact(
                            variant_entry.get(kind),
                            context=(
                                f"{base} rev{revision_key} format {format_key} "
                                f"{variant} {kind}"
                            ),
                            zenodo_record=format_record,
                            expected_filename=expected,
                        )
    return registry


def _cache_paths() -> tuple[Path, Path]:
    cache_dir = default_cache_dir()
    return cache_dir / REGISTRY_FILENAME, cache_dir / f"{REGISTRY_FILENAME}.meta"


def _read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json_atomic(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, path)


def _read_valid_cache(path: Path) -> dict | None:
    try:
        return validate_registry(_read_json(path))
    except (OSError, json.JSONDecodeError, ValueError):
        return None


def _read_cache_metadata(path: Path) -> dict:
    try:
        value = _read_json(path)
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _read_bundled_registry() -> dict:
    try:
        resource = resources.files("metacooc").joinpath(REGISTRY_FILENAME)
        registry = json.loads(resource.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise DataReleaseError("The packaged data registry is unavailable or invalid.") from exc
    try:
        return validate_registry(registry)
    except ValueError as exc:
        raise DataReleaseError(f"The packaged data registry is invalid: {exc}") from exc


def _warn_refresh_failure(exc: Exception) -> None:
    global _refresh_warning_emitted
    if _refresh_warning_emitted:
        return
    warnings.warn(
        "Could not refresh the MetaCoOc data registry; using the last valid "
        f"local registry ({exc}).",
        RuntimeWarning,
        stacklevel=3,
    )
    _refresh_warning_emitted = True


def _warn_cache_failure(exc: OSError) -> None:
    global _refresh_warning_emitted
    if _refresh_warning_emitted:
        return
    warnings.warn(
        "The MetaCoOc data registry was refreshed but could not be cached; "
        f"the validated remote registry will be used for this process ({exc}).",
        RuntimeWarning,
        stacklevel=3,
    )
    _refresh_warning_emitted = True


def load_registry(*, now: float | None = None) -> dict:
    """Return the current validated registry, refreshing its daily cache lazily."""
    global _memory_registry, _memory_registry_expires
    current_time = time.time() if now is None else now
    if _memory_registry is not None and current_time < _memory_registry_expires:
        return _memory_registry

    registry_path, metadata_path = _cache_paths()
    cached = _read_valid_cache(registry_path)
    metadata = _read_cache_metadata(metadata_path)
    checked_at = metadata.get("checked_at", 0)
    if cached is not None and isinstance(checked_at, (int, float)):
        if current_time - checked_at < REGISTRY_CACHE_SECONDS:
            _memory_registry = cached
            _memory_registry_expires = checked_at + REGISTRY_CACHE_SECONDS
            return cached

    headers = {}
    if cached is not None and isinstance(metadata.get("etag"), str):
        headers["If-None-Match"] = metadata["etag"]

    try:
        response = requests.get(REGISTRY_URL, headers=headers, timeout=REGISTRY_TIMEOUT)
        if response.status_code == 304:
            if cached is None:
                raise ValueError("registry server returned 304 without a cached registry")
            registry = cached
            etag = metadata.get("etag")
        else:
            response.raise_for_status()
            registry = validate_registry(response.json())
            etag = response.headers.get("ETag")
    except (requests.RequestException, OSError, ValueError, json.JSONDecodeError) as exc:
        _warn_refresh_failure(exc)
        registry = cached if cached is not None else _read_bundled_registry()
    else:
        try:
            if response.status_code != 304:
                _write_json_atomic(registry_path, registry)
            cache_metadata = {"checked_at": current_time}
            if isinstance(etag, str):
                cache_metadata["etag"] = etag
            _write_json_atomic(metadata_path, cache_metadata)
        except OSError as exc:
            _warn_cache_failure(exc)

    _memory_registry = registry
    _memory_registry_expires = current_time + REGISTRY_CACHE_SECONDS
    return registry


def available_releases(registry: dict | None = None) -> list[str]:
    registry = registry or load_registry()
    result = []
    for base in sorted(registry["releases"], key=_release_number):
        revisions = registry["releases"][base]["revisions"]
        for revision in sorted((int(value) for value in revisions)):
            for variant in VARIANTS:
                result.append(f"{base}_{variant}_rev{revision}")
    return result


def is_current_release(data_release: str, registry: dict | None = None) -> bool:
    parsed = parse_data_release(data_release)
    registry = registry or load_registry()
    release = registry.get("releases", {}).get(parsed.base)
    return bool(release and release.get("current_revision") == parsed.revision)


def describe_data_release(data_release: str, registry: dict | None = None) -> str:
    try:
        parsed = parse_data_release(data_release)
    except ValueError:
        return str(data_release)
    current = "current; " if registry and is_current_release(data_release, registry) else ""
    return (
        f"{parsed.canonical} ({current}{parsed.variant.upper()}, database release "
        f"{parsed.base}, revision {parsed.revision})"
    )


def local_data_releases(data_dir: str | os.PathLike | None) -> list[str]:
    if not data_dir or not os.path.isdir(data_dir):
        return []

    releases = set()
    for name in os.listdir(data_dir):
        path = os.path.join(data_dir, name)
        if not os.path.isdir(path):
            continue
        match = _LOCAL_RELEASE_RE.fullmatch(name)
        if match and int(match.group("format")) == INGREDIENTS_FORMAT_VERSION:
            releases.add(match.group("release"))

    def sort_key(value: str) -> tuple[int, int, str]:
        parsed = parse_data_release(value)
        return _release_number(parsed.base), parsed.revision, parsed.variant

    return sorted(releases, key=sort_key)


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
    defaulted: bool = False,
) -> str:
    local_releases = local_data_releases(data_dir)
    local_display = (
        ", ".join(describe_data_release(value) for value in local_releases)
        if local_releases
        else "none"
    )
    suggestion = ""
    if local_releases:
        suggestion = (
            "\nTo use an available local data source, rerun with "
            f"--data-release {local_releases[0]}."
        )
    metadata_option = " --include-metadata" if file_kind == "metadata" else ""
    return (
        f"Data release {describe_data_release(data_release)} was requested.\n"
        f"Missing required {file_kind} file: {missing_path}\n"
        f"Available local data sources in {data_dir}: {local_display}."
        f"{suggestion}\n"
        "To fetch the requested data source, run: "
        f"metacooc download --data_dir {data_dir} --data-release {data_release}"
        f"{metadata_option}"
    )


def resolve_release(
    data_release: str,
    *,
    format_version: int = INGREDIENTS_FORMAT_VERSION,
    registry: dict | None = None,
) -> ResolvedRelease:
    try:
        parsed = parse_data_release(data_release)
    except ValueError as exc:
        raise DataReleaseError(str(exc)) from exc
    registry = registry or load_registry()
    release = registry.get("releases", {}).get(parsed.base)
    revision = (
        release.get("revisions", {}).get(str(parsed.revision)) if release else None
    )
    if revision is None:
        available = ", ".join(available_releases(registry)) or "none"
        raise DataReleaseError(
            f"Data release '{parsed.canonical}' is not available. "
            f"Available releases: {available}"
        )

    format_entry = revision["ingredients_formats"].get(str(format_version))
    if format_entry is None:
        formats = ", ".join(sorted(revision["ingredients_formats"], key=int))
        raise DataReleaseError(
            f"Data release '{parsed.canonical}' is published for Ingredients "
            f"format(s) {formats}, but this MetaCoOc version requires format "
            f"{format_version}. Update MetaCoOc to use this release."
        )

    def artifact(value: dict, record: str) -> Artifact:
        return Artifact(value["filename"], record, value["sha256"])

    variant = format_entry[parsed.variant]
    support = revision["support"]
    artifacts = {
        "ingredients_raw": artifact(variant["raw"], format_entry["zenodo_record"]),
        "ingredients_aggregated": artifact(
            variant["aggregated"], format_entry["zenodo_record"]
        ),
        "sra_metadata": artifact(support["sra_metadata"], support["zenodo_record"]),
        "sample_to_biome": artifact(
            support["sample_to_biome"], support["zenodo_record"]
        ),
    }
    return ResolvedRelease(
        parsed.canonical,
        format_version,
        release["current_revision"] == parsed.revision,
        artifacts,
    )


def get_file_info(data_release: str) -> tuple[dict[str, str], dict[str, str]]:
    """Compatibility view of resolved filenames and download URLs."""
    resolved = resolve_release(data_release)
    urls = {
        artifact.filename: artifact.url
        for artifact in resolved.artifacts.values()
    }
    return resolved.filenames, urls
