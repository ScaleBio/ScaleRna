#!/usr/bin/env python
"""
Python script to merge reportStatistics.tsv for all samples

Args:
    samples (str): reportStatistics.tsv for all samples
"""
import argparse
import sys
import pandas as pd
from typing import List
from pathlib import Path


def merge(sampleFns: List[Path]):
    """
    Function to merge reportStatistics.csv for all samples

    Args:
        sampleFns: List of all per-sample reportStatistics.csv
    """
    merged = pd.DataFrame(columns=["Category", "Metric"])
    for fn in sorted(sampleFns):
        sampleStats = pd.read_csv(fn, names=["Category", "Metric", "Value"], header=None)
        sample_name = sampleStats["Value"][0]  # First row of Value column is the sample name
        sample_df = sampleStats.iloc[1:].copy()  # Skip the first row
        sample_df.rename(columns={"Value": sample_name}, inplace=True)
        merged = pd.merge(
            merged, sample_df, on=["Category", "Metric"], how="outer", sort=False  # Use outer join to include all rows
        )
        merged = merged.set_index(["Category", "Metric"])
        # Always use the order in sample_df
        merged = merged.reindex(sample_df.set_index(["Category", "Metric"]).index.union(merged.index, sort=False))
        merged = merged.reset_index()
    merged.to_csv(sys.stdout, index=False, header=True)


def main():
    parser = argparse.ArgumentParser(description="Merge reportStatistics.csv files for muiltiple samples into one")
    parser.add_argument("samples", nargs="+", type=str, help="reportStatistics.csv for all samples")
    args = parser.parse_args()

    merge(args.samples)


if __name__ == "__main__":
    main()
