#Utilities related to file I/O and parsing

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