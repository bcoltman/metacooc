from __future__ import annotations

import builtins
import io
import tarfile

from metacooc import download


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


def test_download_stream_writes_chunks_without_buffering_response(monkeypatch, tmp_path):
    chunks = [b"abc", b"", b"defg", b"h"]
    response = FakeResponse(chunks, headers={"content-length": "8"})
    progress_updates = []

    def fake_get(url, stream=False):
        assert url == "https://example.test/file.gz"
        assert stream is True
        return response

    class FakeTqdm:
        def __init__(self, **kwargs):
            assert kwargs["total"] == 8
            assert kwargs["unit"] == "B"
            assert kwargs["unit_scale"] is True
            assert kwargs["unit_divisor"] == 1024
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
    download._download_stream("https://example.test/file.gz", out)

    assert out.read_bytes() == b"abcdefgh"
    assert progress_updates == [3, 4, 1]


def test_download_data_extracts_local_ingredients_tarball(monkeypatch, tmp_path):
    payload = io.BytesIO()
    with tarfile.open(fileobj=payload, mode="w:gz") as tar:
        data = b'{"format_version": 1}\n'
        info = tarfile.TarInfo("ingredients_raw_1.1.0_gtdb/manifest.json")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    archive_bytes = payload.getvalue()

    monkeypatch.setattr(
        download,
        "get_file_info",
        lambda version: (
            {"ingredients_raw": "ingredients_raw_1.1.0_gtdb"},
            {"ingredients_raw_1.1.0_gtdb.tar.gz": "https://example.test/raw.tar.gz"},
        ),
    )
    monkeypatch.setattr(builtins, "input", lambda prompt: "y")

    def fake_download_stream(url, temp_path):
        assert url == "https://example.test/raw.tar.gz"
        with open(temp_path, "wb") as f:
            f.write(archive_bytes)

    monkeypatch.setattr(download, "_download_stream", fake_download_stream)

    download.download_data(tmp_path, data_version="1.1.0_gtdb")

    assert (tmp_path / "ingredients_raw_1.1.0_gtdb" / "manifest.json").exists()
