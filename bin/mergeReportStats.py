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
    merged = {}
    for fn in sorted(sampleFns):
        sampleStats = pd.read_csv(fn, names=["Category", "Metric", "Value"], header=None)
        if 'Category' not in merged:
            merged['Category'] = sampleStats.Category[1:]
        if 'Metric' not in merged:
            merged['Metric'] = sampleStats.Metric[1:]
        else:
            if (merged['Metric'] != sampleStats.Metric[1:]).any():
                raise Exception(f"Mismatched sample metrics in {fn}")
        name = sampleStats['Value'][0]
        merged[name] = sampleStats.Value[1:]

    mergedTab = pd.DataFrame(merged)
    mergedTab.to_csv(sys.stdout, index=False, header=True)


def main():
    parser = argparse.ArgumentParser(
        description="Merge reportStatistics.csv files for muiltiple samples into one")
    parser.add_argument("samples", nargs='+', type=str, help="reportStatistics.csv for all samples")
    args = parser.parse_args()

    merge(args.samples)


if __name__ == "__main__":
    main()
