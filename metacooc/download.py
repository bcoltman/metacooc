#!/usr/bin/env python3
"""
download.py

Download initial data files for metacooc.

This script downloads the following default files into the specified data directory:
    - ingredients_raw_<version>.pkl
    - ingredients_aggregated_genus_<version>.pkl
    - sra_metadata_<version>.tsv

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
from metacooc._data_config import DOWNLOAD_URLS

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
        
    for final_name, url in DOWNLOAD_URLS.items():
        # If the filename includes ".tsv.gz", we keep the .gz and decompress to .tsv
        if final_name.endswith(".gz"):
            target_name = final_name[:-3]
        else:
            target_name = final_name
            
        target_path = os.path.join(data_dir, target_name)
        temp_path = os.path.join(data_dir, target_name + ".tmp.gz")
        
        if os.path.exists(target_path) and not force:
            print(f"{target_path} already exists; skipping download.")
            continue
            
        print(f"Downloading {url} to {temp_path} ...")
        response = requests.get(url)
        response.raise_for_status()
        with open(temp_path, "wb") as f:
            f.write(response.content)
        print(f"Downloaded {temp_path}")
        
        print(f"Unzipping {temp_path} to {target_path} ...")
        with gzip.open(temp_path, 'rb') as f_in, open(target_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        print(f"Unzipped to {target_path}")
        
        os.remove(temp_path)
        print(f"Removed temporary file {temp_path}")

def main():
    parser = argparse.ArgumentParser(description="Download metacooc data files")
    parser.add_argument("--data_dir", type=str, required=True, help="Target directory for data files")
    parser.add_argument("--force", action="store_true", help="Force re-download even if files exist")
    args = parser.parse_args()
    
    download_data(args.data_dir, force=args.force)

if __name__ == "__main__":
    main()
