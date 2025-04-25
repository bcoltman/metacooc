# _data_config.py

DATA_VERSION = "1.0.0"
ZENODO_RECORD_ID = "15283587"
ZENODO_BASE_URL = f"https://zenodo.org/records/{ZENODO_RECORD_ID}/files"

FILENAMES = {
    "ingredients_raw": f"ingredients_raw_{DATA_VERSION}.pkl",
    "ingredients_aggregated": f"ingredients_aggregated_genus_{DATA_VERSION}.pkl",
    "sra_metadata": f"sra_metadata_{DATA_VERSION}.tsv",
}

# URLs point to gzip-compressed files (even if file doesn't have .gz extension in the name)
DOWNLOAD_URLS = {
    FILENAMES["ingredients_raw"]: f"{ZENODO_BASE_URL}/{FILENAMES['ingredients_raw']}.gz?download=1",
    FILENAMES["ingredients_aggregated"]: f"{ZENODO_BASE_URL}/{FILENAMES['ingredients_aggregated']}.gz?download=1",
    FILENAMES["sra_metadata"] + ".gz": f"{ZENODO_BASE_URL}/{FILENAMES['sra_metadata']}.gz?download=1",
}
