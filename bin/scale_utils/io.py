"""Utilities related to file I/O and parsing"""

import json
from collections import OrderedDict
from pathlib import Path
from typing import Dict


def readJSON(file: Path, preserveDictOrder: bool = False):
    """
    Function to read in JSON file and return it as a dictionary

    Args:
        file: Path to .json file
        preserveDictOrder: Flag to indicate whether to
            read the file while preserving the order of entries

    Returns:
        Dictionary with contents of json file
    """
    with open(file) as f:
        str = f.read()
        strStripped = str.rstrip()
        pairs_hook = OrderedDict if preserveDictOrder else None
        parsedJSON = json.loads(strStripped, object_pairs_hook=pairs_hook)
    return parsedJSON

def ensurePathsExist(filePaths: Dict[str, Path]):
    """
    Function to ensure all paths mentioned in the given @filePaths exist
    Raises FileNotFoundError if any file not found
    """
    for key, value in filePaths.items():
        if not value.exists():
            raise FileNotFoundError(f"{key} was assumed to be located at '{str(value)}'. It is missing")

def resolve_sample_specific_file_paths(STARsolo_out, feature_type, matrix_type):
    """
    Returns a dictionary of the paths to the necessary files in the STARsolo output directory. 

    Args:
        STARsolo_out (Path): Path to the STARsolo output directory for the sample
        feature_type (str): STARsolo feature type used (e.g. GeneFull_Ex50pAS, etc.)
        matrix_type (str): STARsolo matrix type used (i.e. uniquely mapped vs multimapped)

    Returns:
        A dictionary where each key is a custom file identifier and each value is a path to the identified file in the STARsolo output directory (Dict[str, Path])
    """
    mtx_prefix = STARsolo_out / feature_type
    files = dict(features = mtx_prefix / "raw" / "features.tsv",
                 barcodes = mtx_prefix / "raw" / "barcodes.tsv",
                 mtx = mtx_prefix / "raw" / matrix_type,
                 summary = mtx_prefix / "Summary.csv",
                 stats = mtx_prefix / "CellReads.stats")
    ensurePathsExist(files)

    return files
