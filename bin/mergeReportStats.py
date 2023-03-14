#!/usr/bin/env python
"""
Python script to merge reportStatistics.tsv for all samples

Args:
    samples (str): reportStatistics.tsv for all samples
"""
import argparse
import sys
import pandas as pd


def merge(sampleFns):
    """
    Function to merge reportStatistics.tsv for all samples

    Args:
        sampleFns (list): List of all reportStatistics.tsv
    """
    merged = {}

    for fn in sampleFns:
        tab = pd.read_csv(fn, names=["Category", "Metric", "Value"],
                          header=None, sep="\t")

        if 'Category' not in merged:
            merged['Category'] = tab.Category[1:]

        if 'Metric' not in merged:
            merged['Metric'] = tab.Metric[1:]

        name = tab['Value'][0]
        merged[name] = tab.Value[1:]

    merged = pd.DataFrame(merged)
    merged.to_csv(sys.stdout, sep="\t", index=False, header=True)


def main():
    parser = argparse.ArgumentParser(
        description="Merge reportStatistics.tsv files "
                    "for muiltiple samples into one")
    parser.add_argument("samples", nargs='+', type=str,
                        help="reportStatistics.tsv for all samples")
    args = parser.parse_args()

    merge(args.samples)


if __name__ == "__main__":
    main()
