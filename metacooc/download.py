#!/usr/bin/env python3
"""
download.py

Download initial data files for metacooc.

This script downloads the following default files into the specified data directory:
    - ingredients_raw_<data_release>/
    - ingredients_aggregated_<data_release>/
    - sra_metadata_<base_release>.tsv

Ingredients are downloaded as .tar.gz archives, extracted to Ingredients
directories, and the temporary archives are removed. Metadata files are
downloaded as .gz files, decompressed, and the temporary files are removed.
Biome annotations are stored inside each Ingredients directory manifest rather
than downloaded as a separate file.
Use the --force flag to re-download files even if they already exist.

Usage (CLI):
    metacooc download --data_dir /path/to/data [--force]
"""

import os
import argparse
import requests
import gzip
import shutil
import tarfile
from tqdm import tqdm


from metacooc._data_config import *
from metacooc._paths import default_data_dir

CHUNK_SIZE = 8 * 1024 * 1024


def _download_stream(url, temp_path):
    """
    Stream a remote file to disk without materializing the full response in RAM.
    """
    with requests.get(url, stream=True) as response:
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
                progress.update(len(chunk))


def download_data(data_dir=None, list_data_releases=False, data_release=None, force=False):
    """
    Download data files for a specific Sandpiper data_release into data_dir.
    
    Parameters:
        data_dir (str): Directory where data files will be saved.
        force (bool): If True, force re-download even if the file exists.
        data_release (str): data release to download (default: latest available).
    """
    if list_data_releases:
        available = ", ".join(available_releases()) or "none"
        print(f"Available: {available}")
        return

    data_release = data_release or LATEST_DATA_RELEASE
    data_dir = os.fspath(data_dir or default_data_dir())

    filenames, download_urls = get_file_info(data_release)
        
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        print(f"Created data directory: {data_dir}")
        
    # Determine how many files need to be downloaded
    missing_files = 0
    for final_name, url in download_urls.items():
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
        
    for final_name, url in download_urls.items():
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
        _download_stream(url, temp_path)
        print(f"Downloaded {temp_path}")
        
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
        
        os.remove(temp_path)
        print(f"Removed temporary file {temp_path}")


def main():
    parser = argparse.ArgumentParser(description="Download metacooc data files")
    parser.add_argument(
        "--data_dir",
        type=str,
        default=str(default_data_dir()),
        help="Target directory for data files (default: %(default)s)",
    )
    parser.add_argument("--force", action="store_true", help="Force re-download even if files exist")
    parser.add_argument("--data-release", default=None, help="Specify which data release to download (default: latest)")
    parser.add_argument("--list-data-releases", action="store_true", help="List available data releases")
    args = parser.parse_args()
    
    download_data(
        args.data_dir,
        list_data_releases=args.list_data_releases,
        data_release=args.data_release,
        force=args.force,
    )

if __name__ == "__main__":
    main()
