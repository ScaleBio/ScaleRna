#!/usr/bin/env python
"""
Generate library metrics from all sample metrics for a library
"""
import argparse
import pandas as pd
import os
import shutil
from pathlib import Path
from utils.base_logger import logger

def concat_sample_metrics(sample_metrics, libName):
    """
    Concatenate sample metrics of one library into single metrics file

    Args:
        sample_metrics (list): List of paths to sample metrics
        libName (str): Library name
    """
    list_of_dfs = []

    # Generate list of dataframes with each dataframe representing the
    # metrics of one sample
    for idx, fname in enumerate(sample_metrics):
        df = pd.read_csv(fname, index_col=0, header=0)
        df['Barcode'] = df.index
        
        list_of_dfs.append(df)

    # Concatenate list of dataframes into one dataframe
    frame = pd.concat(list_of_dfs, axis=0)
    
    # Build list that contains indices for concatenated dataframe
    frame_index = []
    for sample in frame["sample"].unique():
        sample_entries = frame[frame["sample"] == sample]
        frame_index.append(list(range(0, len(sample_entries))))

    # Flatten out list of lists that has indices for concatenated dataframe
    frame_index = [item for sublist in frame_index for item in sublist]
    frame[''] = frame_index
    frame.set_index("", inplace=True)

    # Metrics will be written to the library_metrics folder
    metricsDir = Path(".", f"library_{libName}_metrics")
    metricsDir.mkdir()
    logger.debug(f"Writing allCellsBetweenFiles to {str(metricsDir.resolve())}")
    frame.to_csv(f"{metricsDir}/allCellsBetweenFiles.csv")

def main():
    parser = argparse.ArgumentParser(description="Generate metrics for library report")
    parser.add_argument("--sample_metrics", nargs='+', type=str, required=True,
                        help="Per sample metrics that need to be concatenated")
    parser.add_argument("--libName", type=str,
                        help="libName specified in passed samplesheet for which to generate a fastq report")
    
    args = parser.parse_args()

    concat_sample_metrics(args.sample_metrics, args.libName)

if __name__ == "__main__":
    main()
