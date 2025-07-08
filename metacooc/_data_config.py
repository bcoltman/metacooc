# _data_config.py

RELEASES = {
    "1.0.0": "15283587",
    "0.3.0": "15025528"
            }

LATEST_VERSION = max(RELEASES.keys(), key=lambda v: tuple(map(int, v.split("."))))

def get_file_info(version):
    if version not in RELEASES:
        raise ValueError(
            f"Version '{version}' is not available. "
            f"Available versions: {', '.join(RELEASES.keys())}"
        )
    
    record_id = RELEASES[version]
    base_url = f"https://zenodo.org/records/{record_id}/files"
    
    filenames = {
        "ingredients_raw": f"ingredients_raw_{version}.pkl",
        "ingredients_aggregated": f"ingredients_aggregated_genus_{version}.pkl",
        "sra_metadata": f"sra_metadata_{version}.tsv",
    }
    
    download_urls = {
        filenames["ingredients_raw"]: f"{base_url}/{filenames['ingredients_raw']}.gz?download=1",
        filenames["ingredients_aggregated"]: f"{base_url}/{filenames['ingredients_aggregated']}.gz?download=1",
        filenames["sra_metadata"] + ".gz": f"{base_url}/{filenames['sra_metadata']}.gz?download=1",
    }
    
    return filenames, download_urls



# DATA_VERSION = 
# ZENODO_RECORD_ID = 
# ZENODO_BASE_URL = f"https://zenodo.org/records/{ZENODO_RECORD_ID}/files"

# FILENAMES = {
    # "ingredients_raw": f"ingredients_raw_{DATA_VERSION}.pkl",
    # "ingredients_aggregated": f"ingredients_aggregated_genus_{DATA_VERSION}.pkl",
    # "sra_metadata": f"sra_metadata_{DATA_VERSION}.tsv",
# }

# # URLs point to gzip-compressed files (even if file doesn't have .gz extension in the name)
# DOWNLOAD_URLS = {
    # FILENAMES["ingredients_raw"]: f"{ZENODO_BASE_URL}/{FILENAMES['ingredients_raw']}.gz?download=1",
    # FILENAMES["ingredients_aggregated"]: f"{ZENODO_BASE_URL}/{FILENAMES['ingredients_aggregated']}.gz?download=1",
    # FILENAMES["sra_metadata"] + ".gz": f"{ZENODO_BASE_URL}/{FILENAMES['sra_metadata']}.gz?download=1",
# }
