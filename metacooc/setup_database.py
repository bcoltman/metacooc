import logging
import pathlib

import pooch


logger = logging.getLogger(__name__)



DEFAULT_SANDPIPER = {
	doi:"10.5281/zenodo.11516218",
    filename:"sandpiper0.3.0.condensed.csv.gz",
    md5:"5cc1a638b7bb0a984814162e84feee12",
    url:"https://zenodo.org/api/",
)


def get_cached_filepath(settings) -> pathlib.Path:
    """Get data from a file stored in Zenodo (or from cache, if available).

    Parameters
    ----------
    settings : ZenodoSettings
        Configuration for the interaction with Zenodo.org.

    Returns
    -------
    pathlib.Path
        The path to the locally cached file.

    """
    cache_directory = pathlib.Path(appdirs.user_cache_dir(appname="equilibrator"))
    cache_directory.mkdir(parents=True, exist_ok=True)

    cache_fname = pooch.retrieve(
        path=cache_directory,
        fname=settings["filename"],
        url="doi:" + settings["doi"] + "/" + settings["filename"],
        known_hash="md5:" + settings["md5"],
        progressbar=True,
    )
    return pathlib.Path(cache_fname)


def _download_file(self, file_url, out_file, progress_bar=False):
        """Download a file to disk
        Streams a file from URL to disk.
        Can optionally use tqdm for a visual download bar
        Arguments:
            file_url (str): URL of file to download
            out_file (str): Target file path
            progress_bar (bool): Display graphical progresss bar
        """
        if progress_bar:
            logging.info('Downloading {} to {}.'.format(file_url, out_file))
            response = requests.get(file_url, stream=True)
            total_size_in_bytes = int(response.headers.get('content-length', 0))
            block_size = 1024
            progress_bar = tqdm(total=total_size_in_bytes, unit='iB', unit_scale=True)
            with open(out_file, 'wb') as file:
                for data in response.iter_content(block_size):
                    progress_bar.update(len(data))
                    file.write(data)
            progress_bar.close()

        else:
            with requests.get(file_url, stream=True) as r:
                with open(out_file, 'wb') as f:
                    shutil.copyfileobj(r.raw, f)