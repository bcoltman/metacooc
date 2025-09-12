# _data_config.py

RELEASES = {
    "1.0.1": "16030719",
    # "1.0.0": "15283587",
    # "0.3.0": "15025528"
            }

NO_LONGER_COMPATIBLE = ["1.0.0", "0.3.0"]

LATEST_VERSION = max(RELEASES.keys(), key=lambda v: tuple(map(int, v.split("."))))

def get_file_info(version):
    if version in NO_LONGER_COMPATIBLE:
        raise ValueError(
            f"Version '{version}' is no longer compatible. "
            f"Available versions: {', '.join(RELEASES.keys())}"
        )
    if version not in RELEASES:
        raise ValueError(
            f"Version '{version}' is not available. "
            f"Available versions: {', '.join(RELEASES.keys())}"
        )
    
    record_id = RELEASES[version]
    base_url = f"https://zenodo.org/records/{record_id}/files"
    
    filenames = {
        "ingredients_raw": f"ingredients_raw_{version}.pkl",
        "ingredients_aggregated": f"ingredients_aggregated_{version}.pkl",
        "sra_metadata": f"sra_metadata_{version}.tsv",
        "sample_to_biome": f"sample_to_biome_{version}.csv",
    }
    
    download_urls = {
        filenames["ingredients_raw"]: f"{base_url}/{filenames['ingredients_raw']}.gz?download=1",
        filenames["ingredients_aggregated"]: f"{base_url}/{filenames['ingredients_aggregated']}.gz?download=1",
        filenames["sra_metadata"] + ".gz": f"{base_url}/{filenames['sra_metadata']}.gz?download=1",
        filenames["sample_to_biome"] + ".gz": f"{base_url}/{filenames['sample_to_biome']}.gz?download=1",
    }
    
    return filenames, download_urls