#!/usr/bin/env python
"""
Concatenate either CSV or parquet files, in order provided,
into single output file in the format defined by --outputFile extension
"""

import duckdb
import argparse
from pathlib import Path


def concat_files(files: list[Path], cols: list[str], output_file: Path):
    """
    Concatenate files with same columns

    Args:
        files: List of files to concatenate
        cols: List of columns to include
        output_file: Output file path
    """
    select_stmt = ""
    for idx, file in enumerate(files):
        select_stmt += f"""
        SELECT {",".join(cols)}
        FROM '{file}'
        -- UNION ALL with next file if not last
        {"" if idx == len(files) - 1 else "UNION ALL BY NAME"}
        """
    duckdb.sql(
        f"""
    COPY ({select_stmt})
    TO '{output_file}';
    """
    )


def main():
    parser = argparse.ArgumentParser("Vertically stack columnar data files, e.g. csv/parquet")
    parser.add_argument("inputFiles", nargs="+", type=Path, help="Files to merge")
    parser.add_argument("--outputFile", required=True, type=Path, help="Output path for merged file")
    parser.add_argument(
        "--columns", nargs="+", default="*", type=str, help="duckdb syntax for subset of columns to include"
    )
    args = parser.parse_args()

    concat_files(files=args.inputFiles, cols=args.columns, output_file=args.outputFile)


if __name__ == "__main__":
    main()
