#!/usr/bin/env python3
"""
download.py

Download initial data files for metacooc.

This script downloads the following default files into the specified data directory:
    - ingredients_raw.pkl
    - ingredients_aggregated_genus.pkl
    - sra_metadata.tsv (unzipped from sra_metadata.tsv.gz)

All download URLs point to gzip-compressed files. This script downloads each file
to a temporary .gz file, unzips it, and then removes the temporary file.
Use the --force flag to re-download files even if they already exist.

Usage (CLI):
    metacooc download --data_dir /path/to/data [--force]
"""

import os
import argparse
import requests
import gzip
import shutil

# Mapping of final filenames to their download URLs.
# Note: The URL always points to a .gz file.
DOWNLOAD_FILES = {
    "ingredients_raw.pkl": "https://zenodo.org/records/15025528/files/ingredients_raw.pkl.gz?download=1",
    "ingredients_aggregated_genus.pkl": "https://zenodo.org/records/15025528/files/ingredients_aggregated_genus.pkl.gz?download=1",
    "sra_metadata.tsv.gz": "https://zenodo.org/records/15025528/files/sra_metadata.tsv.gz?download=1",
}

def download_data(data_dir, force=False):
    """
    Download default data files into data_dir.
    
    Parameters:
        data_dir (str): Directory where data files will be saved.
        force (bool): If True, force re-download even if the file exists.
    """
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        print(f"Created data directory: {data_dir}")
        
    for final_name, url in DOWNLOAD_FILES.items():
        # Determine target file: if final_name ends with .gz, remove the extension after unzipping.
        if final_name.endswith(".gz"):
            target_name = final_name[:-3]
        else:
            target_name = final_name
            
        target_path = os.path.join(data_dir, target_name)
        temp_path = os.path.join(data_dir, target_name + ".tmp.gz")
        
        if os.path.exists(target_path) and not force:
            print(f"{target_path} already exists; skipping download.")
            continue
            
        # Download the file to a temporary path.
        print(f"Downloading {url} to {temp_path} ...")
        response = requests.get(url)
        response.raise_for_status()
        with open(temp_path, "wb") as f:
            f.write(response.content)
        print(f"Downloaded {temp_path}")
        
        # Unzip the downloaded file.
        print(f"Unzipping {temp_path} to {target_path} ...")
        with gzip.open(temp_path, 'rb') as f_in, open(target_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"Unzipped to {target_path}")
        
        # Remove the temporary gz file.
        os.remove(temp_path)
        print(f"Removed temporary file {temp_path}")