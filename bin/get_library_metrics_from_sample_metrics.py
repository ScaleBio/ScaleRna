#!/usr/bin/env python
"""
Generate library metrics from all sample metrics for a library
"""
import argparse
import duckdb
from pathlib import Path
from scale_utils.base_logger import logger


def concat_sample_metrics(sample_metrics, libName):
    """
    Concatenate sample metrics of one library into single metrics file

    Args:
        sample_metrics (list): List of paths to sample metrics
        libName (str): Library name
    """
    # Metrics will be written to the library_metrics folder
    metricsDir = Path(".", f"library_{libName}_metrics")
    metricsDir.mkdir()
    logger.debug(f"Writing allCellsBetweenFiles to {str(metricsDir.resolve())}")
    duckdb.sql(
        f"""
    COPY (
        SELECT *
        FROM read_parquet({sample_metrics})
    ) TO '{metricsDir}/allCellsBetweenFiles.parquet';
    """
    )


def main():
    parser = argparse.ArgumentParser(description="Generate metrics for library report")
    parser.add_argument(
        "--sample_metrics", nargs="+", type=str, required=True, help="Per sample metrics that need to be concatenated"
    )
    parser.add_argument(
        "--libName", type=str, help="libName specified in passed samplesheet for which to generate a fastq report"
    )

    args = parser.parse_args()

    concat_sample_metrics(args.sample_metrics, args.libName)


if __name__ == "__main__":
    main()
