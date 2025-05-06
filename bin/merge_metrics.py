#!/usr/bin/env python
"""
Merge cell metrics parquet files
"""
import argparse
from pathlib import Path
import duckdb


def merge_metrics(id: str, pq_files: list[Path], rna_ids: list[str]):
    """
    Merge cell metrics

    rnaId is a column added to the samples channel in nextflow when analyzing ScalePlex,
    it refers to the matching RNA id (sample.library) for that ScalePlex sample.library.
    Here it is used in the merge samples case to add the RNA sample.libName to the merged barcode cell_id,
    so they match the RNA barcodes that are appended to in a similar manner.

    Args:
        id: Sample/lib name
        pq_files: List of cell metrics files to merge
        rna_ids: Associated rna ids
    """
    select_stmt = ""
    for idx, file in enumerate(pq_files):
        cols = (
            " * "
            if not rna_ids
            else f" * EXCLUDE(Cell_Barcode), concat_ws('_', Cell_Barcode, '{rna_ids[idx]}') AS Cell_Barcode"
        )
        select_stmt += f"""
        SELECT {cols}
        FROM read_parquet('{file}')
        {"" if idx == len(pq_files) - 1 else "UNION ALL"}
        """
    duckdb.sql(
        f"""
    COPY ({select_stmt})
    TO '{id}.cellMetrics.parquet'
    """
    )


def main():
    parser = argparse.ArgumentParser("Merge individual metrics files for sample")
    parser.add_argument("metrics", nargs="+", type=Path, help="Metrics files to merge")
    parser.add_argument("--id", required=True, help="Sample and lib name")
    parser.add_argument("--rnaId", nargs="+", type=str, help="Sample and lib name")
    args = parser.parse_args()

    merge_metrics(id=args.id, pq_files=args.metrics, rna_ids=args.rnaId)


if __name__ == "__main__":
    main()
