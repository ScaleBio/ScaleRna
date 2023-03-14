#!/usr/bin/env python
"""
Python script to generate library metrics from sample metrics
"""
import argparse
import pandas as pd
import os
import shutil
from pathlib import Path
from utils.base_logger import logger

def concat_sample_metrics(sample_metrics, samplesheet, libName):
    """
    Function to concatenate sample metrics of one library into
    a single metrics file representing that library

    Args:
        sample_metrics (list): List of paths to sample metrics
        samplesheet (str): Path to sample sheet
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
    metricsDir = Path(".", "library_metrics")

    # If "reports" directory exists log contents of reports directory
    # and continue execution
    try:
        os.mkdir(metricsDir)
    except OSError as error:
        logger.error(error)
        logger.debug(f"Directory reports exists with contents {os.listdir('reports')}")

    # Assumption: There will always be a sample_metrics1 folder
    # Demux jsons are same for the same library so copying from folder1
    # is equivalent to copying from folder2, etc
    if os.path.isdir("sample_metrics1"):
        shutil.copyfile(f"sample_metrics1/demuxJson.json", f"{metricsDir}/demuxJson.json")
    else:
        shutil.copyfile(f"sample_metrics/demuxJson.json", f"{metricsDir}/demuxJson.json")
    logger.debug(f"Writing allCellsBetweenFiles to {str(metricsDir.resolve())}")

    frame.to_csv(f"{metricsDir}/allCellsBetweenFiles.csv")

def main():
    parser = argparse.ArgumentParser(
        description="Generate metrics for library report")
    parser.add_argument(
        "--sample_metrics", nargs='+', type=str,
        help="Per sample metrics that need to be concatenated",
        required=True)
    parser.add_argument(
        "--samplesheet", required=True,
        help="Path to the samples.csv containing information about libName")
    parser.add_argument(
        "--libName", type=str,
        help="libName specified in passed samplesheet for "
        "which to generate a fastq report")
    
    args = parser.parse_args()

    concat_sample_metrics(args.sample_metrics, args.samplesheet, args.libName)

if __name__ == "__main__":
    main()
