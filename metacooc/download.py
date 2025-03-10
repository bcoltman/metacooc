#!/usr/bin/env python3
"""
download.py

Download initial data files for metacooc.

This script downloads the following default files into the specified data directory:
    - ingredients_raw.pkl
    - ingredients_aggregated.pkl
    - broad_metadata_index.pkl
    - broad_strict_metadata_index.pkl

These files are used by default by the other commands unless alternative file paths are provided.
Use the --force flag to re-download files even if they already exist.

Usage (CLI):
    metacooc download --data_dir /path/to/data [--force]
"""

import os
import argparse
import requests

# Dictionary mapping default filenames to their download URLs.
DOWNLOAD_FILES = {
    "ingredients_raw.pkl": "https://example.com/ingredients_raw.pkl",
    "ingredients_aggregated.pkl": "https://example.com/ingredients_aggregated.pkl",
    "broad_metadata_index.pkl": "https://example.com/broad_metadata_index.pkl",
    "broad_strict_metadata_index.pkl": "https://example.com/broad_strict_metadata_index.pkl",
    # Add additional files (e.g., exact indices) as needed.
}

def download_file(url, dest_path):
    """Download a file from the given URL to dest_path."""
    print(f"Downloading {url} to {dest_path} ...")
    response = requests.get(url)
    response.raise_for_status()  # Raise an error on bad status.
    with open(dest_path, "wb") as f:
        f.write(response.content)
    print(f"Downloaded {dest_path}")

def download_data(data_dir, force=False):
    """
    Download default data files into data_dir.

    Parameters:
        data_dir (str): Directory where data files will be saved.
        force (bool): If True, force re-download even if file exists.
    """
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        print(f"Created data directory: {data_dir}")

    for filename, url in DOWNLOAD_FILES.items():
        dest_path = os.path.join(data_dir, filename)
        if os.path.exists(dest_path) and not force:
            print(f"{filename} already exists; skipping download.")
        else:
            try:
                download_file(url, dest_path)
            except Exception as e:
                print(f"Error downloading {filename} from {url}: {e}")