#!/usr/bin/env python3
"""
download.py

Download initial data files for metacooc.

This script downloads the following default files into the specified data directory:
    - ingredients_raw_<data_release>_format<format_version>/
    - ingredients_aggregated_<data_release>_format<format_version>/

The shared sra_metadata_<base_release>_rev<revision>.tsv file is downloaded only when
--include-metadata is specified.

Ingredients are downloaded as .tar.gz archives, extracted to Ingredients
directories, and the temporary archives are removed. Metadata files are
downloaded as .gz files, decompressed, and the temporary files are removed.
Biome annotations are stored inside each Ingredients directory manifest rather
than downloaded as a separate file.
Use the --force flag to re-download files even if they already exist.

Usage (CLI):
    metacooc download [--data-release R226_gtdb_rev1] [--data_dir /path/to/data] [--force]
"""

import os
import hashlib
import requests
import gzip
import shutil
import tarfile
from tqdm import tqdm


from metacooc._data_config import (
    DataReleaseError,
    available_releases,
    describe_data_release,
    get_default_data_release,
    is_current_release,
    load_registry,
    resolve_release,
)
from metacooc._paths import default_data_dir

CHUNK_SIZE = 8 * 1024 * 1024


def _download_stream(url, temp_path, expected_sha256=None):
    """
    Stream a remote file to disk without materializing the full response in RAM.
    """
    digest = hashlib.sha256()
    try:
        with requests.get(url, stream=True, timeout=(5, 60)) as response:
            response.raise_for_status()
            total = int(response.headers.get("content-length") or 0)

            with open(temp_path, "wb") as f, tqdm(
                total=total or None,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
                desc=os.path.basename(temp_path),
            ) as progress:
                for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                    if not chunk:
                        continue
                    f.write(chunk)
                    digest.update(chunk)
                    progress.update(len(chunk))
    except Exception:
        try:
            os.remove(temp_path)
        except FileNotFoundError:
            pass
        raise
    actual_sha256 = digest.hexdigest()
    if expected_sha256 is not None and actual_sha256 != expected_sha256:
        try:
            os.remove(temp_path)
        except FileNotFoundError:
            pass
        raise DataReleaseError(
            f"Checksum mismatch for {os.path.basename(temp_path)}: expected "
            f"{expected_sha256}, got {actual_sha256}. The temporary download "
            "was removed."
        )
    return actual_sha256


def download_data(
    data_dir=None,
    list_data_releases=False,
    data_release=None,
    force=False,
    include_metadata=False,
):
    """
    Download data files for a specific Sandpiper data_release into data_dir.
    
    Parameters:
        data_dir (str): Directory where data files will be saved.
        force (bool): If True, force re-download even if the file exists.
        data_release (str): exact data release to download.
        include_metadata (bool): Also download the shared SRA metadata table.
    """
    if list_data_releases:
        registry = load_registry()
        releases = available_releases(registry)
        if not releases:
            print("Available: none")
            return
        print("Published data releases:")
        for release in releases:
            labels = []
            if is_current_release(release, registry):
                labels.append("current")
            if release == get_default_data_release(registry):
                labels.append("default")
            marker = f" [{', '.join(labels)}]" if labels else ""
            print(f"  {describe_data_release(release)}{marker}")
        return

    if data_release is None:
        raise DataReleaseError(
            "An exact --data-release is required for downloads, for example "
            "'R226_gtdb_rev1'. Use --list-data-releases to list published releases."
        )
    data_dir = os.fspath(data_dir or default_data_dir())

    resolved = resolve_release(data_release)
    selected_keys = ["ingredients_raw", "ingredients_aggregated"]
    if include_metadata:
        selected_keys.append("sra_metadata")
    artifacts = {key: resolved.artifacts[key] for key in selected_keys}
        
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        print(f"Created data directory: {data_dir}")
        
    # Determine how many files need to be downloaded
    missing_files = 0
    for artifact in artifacts.values():
        final_name = artifact.filename
        if final_name.endswith(".tar.gz"):
            target_name = final_name[:-7]
        elif final_name.endswith(".gz"):
            target_name = final_name[:-3]
        else:
            target_name = final_name
            
        target_path = os.path.join(data_dir, target_name)
        
        if not os.path.exists(target_path) or force:
            missing_files += 1
            
    if missing_files == 0:
        print("All files already exist; skipping download.")
        return
    
    print(f"This script is looking for the download files of {data_release}. If you want an alternative release, "
          "please specify it with --data-release. To see which releases are available, please use --list-data-releases")
    # Prompt user for confirmation
    user_input = input(f"Do you want to download {missing_files} missing files to {data_dir}? (y/n): ").strip().lower()
    if user_input != 'y':
        print("Download cancelled by user.")
        return
        
    for artifact in artifacts.values():
        final_name = artifact.filename
        url = artifact.url
        if final_name.endswith(".tar.gz"):
            target_name = final_name[:-7]
            temp_path = os.path.join(data_dir, target_name + ".tmp.tar.gz")
        elif final_name.endswith(".gz"):
            target_name = final_name[:-3]
            temp_path = os.path.join(data_dir, target_name + ".tmp.gz")
        else:
            target_name = final_name
            temp_path = os.path.join(data_dir, target_name + ".tmp")
            
        target_path = os.path.join(data_dir, target_name)
        
        if os.path.exists(target_path) and not force:
            print(f"{target_path} already exists; skipping download.")
            continue
            
        print(f"Downloading {url} to {temp_path} ...")
        _download_stream(url, temp_path, artifact.sha256)
        print(f"Downloaded {temp_path}")
        
        try:
            if final_name.endswith(".tar.gz"):
                print(f"Extracting {temp_path} to {data_dir} ...")
                with tarfile.open(temp_path, "r:gz") as tar:
                    try:
                        tar.extractall(path=data_dir, filter="data")
                    except TypeError:
                        tar.extractall(path=data_dir)
                print(f"Extracted to {target_path}")
            elif final_name.endswith(".gz"):
                print(f"Unzipping {temp_path} to {target_path} ...")
                with gzip.open(temp_path, 'rb') as f_in, open(target_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                print(f"Unzipped to {target_path}")
            else:
                shutil.move(temp_path, target_path)
                print(f"Moved to {target_path}")
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)
                print(f"Removed temporary file {temp_path}")
