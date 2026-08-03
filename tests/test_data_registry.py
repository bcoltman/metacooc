from __future__ import annotations

import json

import pytest
import requests

import metacooc._data_config as config


class RegistryResponse:
    def __init__(self, *, status_code=200, payload=None, etag=None):
        self.status_code = status_code
        self._payload = payload
        self.headers = {"ETag": etag} if etag else {}

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(str(self.status_code))

    def json(self):
        return self._payload


@pytest.fixture(autouse=True)
def reset_registry_memory(monkeypatch, tmp_path):
    monkeypatch.setattr(config, "_memory_registry", None)
    monkeypatch.setattr(config, "_memory_registry_expires", 0.0)
    monkeypatch.setattr(config, "_refresh_warning_emitted", False)
    monkeypatch.setattr(config, "default_cache_dir", lambda: tmp_path / "cache")


def _write_cache(tmp_path, registry, *, checked_at, etag=None):
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    (cache_dir / config.REGISTRY_FILENAME).write_text(
        json.dumps(registry), encoding="utf-8"
    )
    metadata = {"checked_at": checked_at}
    if etag:
        metadata["etag"] = etag
    (cache_dir / f"{config.REGISTRY_FILENAME}.meta").write_text(
        json.dumps(metadata), encoding="utf-8"
    )


def test_fresh_registry_cache_avoids_network(monkeypatch, tmp_path, registry_factory):
    registry = registry_factory()
    _write_cache(tmp_path, registry, checked_at=100)
    monkeypatch.setattr(
        config.requests,
        "get",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("network used")),
    )

    assert config.load_registry(now=101) == registry


def test_remote_registry_is_validated_and_cached(monkeypatch, tmp_path, registry_factory):
    registry = registry_factory()

    def fake_get(url, headers, timeout):
        assert url == config.REGISTRY_URL
        assert headers == {}
        assert timeout == config.REGISTRY_TIMEOUT
        return RegistryResponse(payload=registry, etag='"registry-1"')

    monkeypatch.setattr(config.requests, "get", fake_get)
    assert config.load_registry(now=100) == registry

    cache_dir = tmp_path / "cache"
    assert json.loads((cache_dir / config.REGISTRY_FILENAME).read_text()) == registry
    metadata = json.loads(
        (cache_dir / f"{config.REGISTRY_FILENAME}.meta").read_text()
    )
    assert metadata == {"checked_at": 100, "etag": '"registry-1"'}


def test_stale_cache_uses_etag_and_handles_not_modified(
    monkeypatch, tmp_path, registry_factory
):
    registry = registry_factory()
    _write_cache(tmp_path, registry, checked_at=0, etag='"registry-1"')

    def fake_get(url, headers, timeout):
        assert headers == {"If-None-Match": '"registry-1"'}
        return RegistryResponse(status_code=304)

    monkeypatch.setattr(config.requests, "get", fake_get)
    assert config.load_registry(now=config.REGISTRY_CACHE_SECONDS + 1) == registry

    metadata_path = tmp_path / "cache" / f"{config.REGISTRY_FILENAME}.meta"
    assert json.loads(metadata_path.read_text())["checked_at"] == (
        config.REGISTRY_CACHE_SECONDS + 1
    )


def test_failed_refresh_warns_once_and_preserves_stale_cache(
    monkeypatch, tmp_path, registry_factory
):
    registry = registry_factory()
    _write_cache(tmp_path, registry, checked_at=0)
    original = (tmp_path / "cache" / config.REGISTRY_FILENAME).read_text()
    monkeypatch.setattr(
        config.requests,
        "get",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            requests.ConnectionError("offline")
        ),
    )

    with pytest.warns(RuntimeWarning, match="last valid"):
        assert config.load_registry(now=config.REGISTRY_CACHE_SECONDS + 1) == registry
    assert config.load_registry(now=config.REGISTRY_CACHE_SECONDS + 2) == registry
    assert (tmp_path / "cache" / config.REGISTRY_FILENAME).read_text() == original


def test_invalid_remote_never_replaces_valid_cache(
    monkeypatch, tmp_path, registry_factory
):
    registry = registry_factory()
    _write_cache(tmp_path, registry, checked_at=0)
    invalid = registry_factory()
    invalid["registry_format_version"] = 99
    monkeypatch.setattr(
        config.requests,
        "get",
        lambda *args, **kwargs: RegistryResponse(payload=invalid),
    )

    with pytest.warns(RuntimeWarning):
        assert config.load_registry(now=config.REGISTRY_CACHE_SECONDS + 1) == registry
    stored = json.loads(
        (tmp_path / "cache" / config.REGISTRY_FILENAME).read_text()
    )
    assert stored == registry


def test_valid_remote_registry_is_used_when_cache_is_read_only(
    monkeypatch, registry_factory
):
    registry = registry_factory()
    monkeypatch.setattr(
        config.requests,
        "get",
        lambda *args, **kwargs: RegistryResponse(payload=registry),
    )
    monkeypatch.setattr(
        config,
        "_write_json_atomic",
        lambda *args, **kwargs: (_ for _ in ()).throw(OSError("read-only")),
    )

    with pytest.warns(RuntimeWarning, match="could not be cached"):
        assert config.load_registry(now=100) == registry


def test_packaged_registry_is_valid_and_contains_only_published_snapshots():
    registry = config._read_bundled_registry()
    assert config.validate_registry(registry) is registry
    assert registry["registry_format_version"] == 1
