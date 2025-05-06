"""Utilities related to file I/O and parsing"""

import json
import gzip
import duckdb
import numpy as np
import scipy.sparse as sp
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
        A dictionary where each key is a custom file identifier and each value is a path
        to the identified file in the STARsolo output directory (Dict[str, Path])
    """
    mtx_prefix = STARsolo_out / feature_type
    files = dict(
        features=mtx_prefix / "raw" / "features.tsv.gz",
        barcodes=mtx_prefix / "raw" / "barcodes.tsv.gz",
        mtx=mtx_prefix / "raw" / f"{matrix_type}.gz",
        stats=mtx_prefix / "CellReads.stats",
    )
    ensurePathsExist(files)

    return files


def load_mtx_table(mtx_path: Path, con: duckdb.DuckDBPyConnection) -> None:
    """
    Load a Matrix Market file into a duckdb table

    Args:
        mtx_path: Path to the mtx file
        con: duckdb connection
    """
    skip_lines = 0
    with gzip.open(mtx_path, "rt") as file:
        for line in file:
            skip_lines += 1
            if not line.startswith("%"):
                dims = line.strip().split()
                break
    nrows, ncols, _ = map(int, dims)
    con.sql(
        f"""
    CREATE OR REPLACE TABLE mtx_metadata (
        nrows UINTEGER,
        ncols UINTEGER
    );
    INSERT INTO mtx_metadata VALUES ({nrows}, {ncols});
    """
    )
    con.sql(
        f"""
    CREATE OR REPLACE TABLE mtx AS
    FROM read_csv(
        '{mtx_path}',
        delim=' ',
        skip={skip_lines},
        columns = {{
            'gene': 'UINTEGER',
            'barcode': 'UINTEGER',
            'count': 'FLOAT'
        }});
    """
    )


def sum_counts_by(con: duckdb.DuckDBPyConnection, cols: np.ndarray, by="gene") -> np.ndarray:
    """
    Sum rounded gene counts from a sparse matrix

    Args:
        con: duckdb connection with existing mtx and mtx_metadata tables
        cols: 0-based column indices to include in sum
        by: 'gene' or 'barcode'

    Returns:
        A numpy array of counts summed by gene or barcode
    """
    nrows, ncols = con.sql("SELECT nrows, ncols FROM mtx_metadata;").fetchall()[0]
    counts = np.zeros(nrows, dtype=int) if by == "gene" else np.zeros(ncols, dtype=int)
    # Create table with only requested columns
    bcs = cols + 1  # mtx file has 1-based index  # noqa: F841
    con.sql("CREATE OR REPLACE TABLE filter_cols (index UINTEGER);")
    con.sql("INSERT INTO filter_cols (SELECT * from bcs);")
    totals = con.sql(
        f"""
    SELECT
        {by},
        SUM(round(count)) as total
    FROM mtx
    WHERE barcode IN (SELECT index FROM filter_cols)
    GROUP BY {by};
    """
    ).fetchnumpy()
    counts[totals[by] - 1] = totals["total"]
    if by == "barcode":
        # subset to requested columns
        counts = counts[cols]

    return counts


def load_partial_mtx(con: duckdb.DuckDBPyConnection, cols: np.ndarray) -> sp.coo_array:
    """
    Load subset of columns from a sparse mtx into a coo array.

    Args:
        con: duckdb connection with existing mtx and mtx_metadata tables
        cols: 0-based column indices to read from the mtx

    Returns:
        A coo array with same dimensions as full array but only specified cols loaded
    """
    nrows, ncols = con.sql("SELECT nrows, ncols FROM mtx_metadata;").fetchall()[0]
    # Create table with only requested columns
    bcs = cols + 1  # mtx file has 1-based index  # noqa: F841
    con.sql("CREATE OR REPLACE TABLE filter_cols (index UINTEGER);")
    con.sql("INSERT INTO filter_cols (SELECT * from bcs);")
    triplets = con.sql(
        """
    SELECT
        gene,
        barcode,
        count
    FROM mtx
    WHERE barcode IN (SELECT index FROM filter_cols);
    """
    ).fetchnumpy()

    return sp.coo_array((triplets["count"], (triplets["gene"] - 1, triplets["barcode"] - 1)), shape=(nrows, ncols))
