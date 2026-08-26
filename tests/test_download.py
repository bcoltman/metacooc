from __future__ import annotations

import builtins
import gzip
import hashlib
import io
import tarfile

import pytest

from metacooc import download
from metacooc._data_config import Artifact, DataReleaseError, ResolvedRelease


def _resolved_release(*, raw_sha="c" * 64, metadata_sha="a" * 64):
    return ResolvedRelease(
        data_release="R226_gtdb_rev1",
        format_version=1,
        current=True,
        artifacts={
            "ingredients_raw": Artifact(
                "ingredients_raw_R226_gtdb_rev1_format1.tar.gz", "1001", raw_sha
            ),
            "ingredients_aggregated": Artifact(
                "ingredients_aggregated_R226_gtdb_rev1_format1.tar.gz",
                "1001",
                "d" * 64,
            ),
            "sra_metadata": Artifact(
                "sra_metadata_R226_rev1.tsv.gz", "1001", metadata_sha
            ),
            "sample_to_biome": Artifact(
                "sample_to_biome_R226_rev1.tsv.gz", "1001", "b" * 64
            ),
        },
    )


class FakeResponse:
    def __init__(self, chunks, headers=None):
        self._chunks = chunks
        self.headers = headers or {}

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def raise_for_status(self):
        return None

    def iter_content(self, chunk_size):
        assert chunk_size == download.CHUNK_SIZE
        yield from self._chunks


class FailingResponse(FakeResponse):
    def iter_content(self, chunk_size):
        yield b"partial"
        raise OSError("connection lost")


def test_download_stream_hashes_chunks_without_buffering(monkeypatch, tmp_path):
    chunks = [b"abc", b"", b"defg", b"h"]
    response = FakeResponse(chunks, headers={"content-length": "8"})
    progress_updates = []

    def fake_get(url, stream=False, timeout=None):
        assert url == "https://example.test/file.gz"
        assert stream is True
        assert timeout == (5, 60)
        return response

    class FakeTqdm:
        def __init__(self, **kwargs):
            assert kwargs["total"] == 8
            assert kwargs["desc"] == "file.tmp.gz"

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def update(self, n):
            progress_updates.append(n)

    monkeypatch.setattr(download.requests, "get", fake_get)
    monkeypatch.setattr(download, "tqdm", FakeTqdm)

    out = tmp_path / "file.tmp.gz"
    expected = hashlib.sha256(b"abcdefgh").hexdigest()
    assert download._download_stream("https://example.test/file.gz", out, expected) == expected
    assert out.read_bytes() == b"abcdefgh"
    assert progress_updates == [3, 4, 1]


def test_download_stream_removes_checksum_mismatch(monkeypatch, tmp_path):
    monkeypatch.setattr(
        download.requests,
        "get",
        lambda *args, **kwargs: FakeResponse([b"wrong"]),
    )
    monkeypatch.setattr(download, "tqdm", lambda **kwargs: _NullProgress())
    out = tmp_path / "bad.tmp.gz"
    with pytest.raises(DataReleaseError, match="Checksum mismatch"):
        download._download_stream("https://example.test/file.gz", out, "0" * 64)
    assert not out.exists()


def test_download_stream_removes_partial_file_after_transport_failure(
    monkeypatch, tmp_path
):
    monkeypatch.setattr(
        download.requests,
        "get",
        lambda *args, **kwargs: FailingResponse([]),
    )
    monkeypatch.setattr(download, "tqdm", lambda **kwargs: _NullProgress())
    out = tmp_path / "partial.tmp.gz"

    with pytest.raises(OSError, match="connection lost"):
        download._download_stream("https://example.test/file.gz", out, "0" * 64)
    assert not out.exists()


class _NullProgress:
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def update(self, _n):
        return None


def test_download_data_extracts_verified_ingredients(monkeypatch, tmp_path):
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w:gz") as tar:
        data = b'{"format_version": 1}\n'
        info = tarfile.TarInfo(
            "ingredients_raw_R226_gtdb_rev1_format1/manifest.json"
        )
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    archive_bytes = payload.getvalue()
    resolved = _resolved_release(raw_sha=hashlib.sha256(archive_bytes).hexdigest())
    (tmp_path / resolved.filenames["ingredients_aggregated"]).mkdir()
    monkeypatch.setattr(download, "resolve_release", lambda release: resolved)
    monkeypatch.setattr(builtins, "input", lambda prompt: "y")

    def fake_download_stream(url, temp_path, expected_sha256):
        assert expected_sha256 == resolved.artifacts["ingredients_raw"].sha256
        with open(temp_path, "wb") as handle:
            handle.write(archive_bytes)

    monkeypatch.setattr(download, "_download_stream", fake_download_stream)
    download.download_data(tmp_path, data_release="R226_gtdb_rev1")

    target = tmp_path / "ingredients_raw_R226_gtdb_rev1_format1"
    assert (target / "manifest.json").exists()
    assert not (tmp_path / f"{target.name}.tmp.tar.gz").exists()


def test_download_excludes_all_support_files_by_default(monkeypatch, tmp_path):
    resolved = _resolved_release()
    for key in ("ingredients_raw", "ingredients_aggregated"):
        (tmp_path / resolved.filenames[key]).mkdir()
    monkeypatch.setattr(download, "resolve_release", lambda release: resolved)
    monkeypatch.setattr(
        builtins,
        "input",
        lambda prompt: (_ for _ in ()).throw(AssertionError("unexpected prompt")),
    )

    download.download_data(tmp_path, data_release="R226_gtdb_rev1")
    assert not (tmp_path / resolved.filenames["sra_metadata"]).exists()
    assert not (tmp_path / resolved.filenames["sample_to_biome"]).exists()


def test_download_include_metadata_adds_only_sra_metadata(monkeypatch, tmp_path):
    content = b"accession\tbiome\nS1\tsoil\n"
    compressed = gzip.compress(content)
    resolved = _resolved_release(metadata_sha=hashlib.sha256(compressed).hexdigest())
    for key in ("ingredients_raw", "ingredients_aggregated"):
        (tmp_path / resolved.filenames[key]).mkdir()
    monkeypatch.setattr(download, "resolve_release", lambda release: resolved)
    monkeypatch.setattr(builtins, "input", lambda prompt: "y")

    def fake_download_stream(url, temp_path, expected_sha256):
        assert "sra_metadata" in url
        assert expected_sha256 == hashlib.sha256(compressed).hexdigest()
        with open(temp_path, "wb") as handle:
            handle.write(compressed)

    monkeypatch.setattr(download, "_download_stream", fake_download_stream)
    download.download_data(
        tmp_path,
        data_release="R226_gtdb_rev1",
        include_metadata=True,
    )

    assert (tmp_path / resolved.filenames["sra_metadata"]).read_bytes() == content
    assert not (tmp_path / resolved.filenames["sample_to_biome"]).exists()


def test_download_requires_exact_release(tmp_path):
    with pytest.raises(DataReleaseError, match="exact --data_release"):
        download.download_data(tmp_path)


def test_download_listing_shows_all_snapshots_and_marks_current(
    monkeypatch, capsys, registry_factory
):
    registry = registry_factory(revisions=(1, 2), current=2)
    monkeypatch.setattr(download, "load_registry", lambda: registry)

    download.download_data(list_data_releases=True)

    output = capsys.readouterr().out
    assert "R226_gtdb_rev1" in output
    assert "R226_globdb_rev1" in output
    assert "R226_gtdb_rev2" in output
    assert "R226_globdb_rev2" in output
    assert output.count("[current") == 2
    assert "R226_globdb_rev2 (GLOBDB" in output
    assert "[current, default]" in output
