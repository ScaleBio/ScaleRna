#!/usr/bin/env python
"""
Python script to merge multiple STAR outputs
"""
import argparse
import pandas as pd
import shutil
import os
from pathlib import Path


def merge_matrix_and_barcode(star_dirs, star_feature, star_matrix, sample_name):
    """
    Function to merge barcodes.tsv and matrix.mtx of each star output directory

    Args:
        star_dirs (list): List containing paths to star output directories
        star_feature (str): Name of star feature used
        star_matrix (str): Filename of the STAR matrix output to use
    
    Output:
        Writes raw matrix.mtx and raw barcodes.tsv
    Returns:
        List of folder names one level down from @star_dirs
    """
    # Generate combined list of barcodes from individual barcodes files
    list_of_barcodes = [open(f"{x}/{star_feature}/raw/barcodes.tsv").read().splitlines() for idx, x in enumerate(star_dirs)]
    combined_list_of_barcodes = [item for sublist in list_of_barcodes for item in sublist]
    combined_dict_of_barcodes = {k: v for v, k in enumerate(combined_list_of_barcodes)}

    if len(combined_list_of_barcodes) != len(set(combined_list_of_barcodes)):
        raise ValueError("Same cell in multiple star output files")

    barcode_row_number = 0
    line_count = 0

    # Open temp matrix for writing combined matrices
    tmp_matrix_fn = "tmp_matrix.mtx"
    with open(tmp_matrix_fn, "w") as f_tmp_mtx:
        # Iterate over all the star directories and merge the matrices
        for idx, star_dir in enumerate(star_dirs):
            f = open(f"{star_dir}/{star_feature}/raw/{star_matrix}", "r")
            
            header = f.readline().strip() # Datatype header
            f.readline(); f.readline() # Skip remaining header lines
            
            for line in f:
                split_line = line.split()
                feature, barcode, umi = int(split_line[0]), int(split_line[1]), float(split_line[2])
                barcode_sequence = list_of_barcodes[idx][barcode-1]
                new_barcode = combined_dict_of_barcodes[barcode_sequence]+1
                
                barcode_row_number += 1
                line_count += 1

                f_tmp_mtx.write(f"{feature} {new_barcode} {umi}\n")

    # Write header information to final merged matrix
    feature_count = len(pd.read_csv(f"{star_dirs[0]}/{star_feature}/raw/features.tsv", sep="\t", header=None).index)
    barcode_count = barcode_row_number
    with open(f"{sample_name}.star.solo/{star_feature}/raw/{star_matrix}", "w") as f_merged_mtx:
        f_merged_mtx.write(f"{header}\n%\n")
        f_merged_mtx.write(f"{feature_count} {barcode_count} {line_count}\n")
    
    # Append contents of temp matrix to final merged matrix
    os.system(f"cat {tmp_matrix_fn} >> {sample_name}.star.solo/{star_feature}/raw/{star_matrix}")

    # Write concatenated list of barcodes to barcodes.tsv
    with open(f"{sample_name}.star.solo/{star_feature}/raw/barcodes.tsv", "w") as f:
        for barcode in combined_list_of_barcodes:
            f.write(barcode+"\n")
    
    # Copy features.tsv from any of the star directories since it'll be same for all
    shutil.copyfile(f"{star_dirs[0]}/{star_feature}/raw/features.tsv",
                    f"{sample_name}.star.solo/{star_feature}/raw/features.tsv")
    
    # We just need "Estimated Number of Cells" from summary.csv to calculate cell thresholding
    # in a downstream process. The code block below gets all metrics from the individual summary.csv
    # files and just concatenates the "Estimated Number of Cells" metric
    list_of_stats = []
    for idx, star_dir in enumerate(star_dirs):
        stats = {}
        with open(f"{star_dir}/{star_feature}/Summary.csv") as f:
            for line in f:
                split_line = line.strip().split(',')
                stats[split_line[0]] = split_line[1]
            list_of_stats.append(stats)
    total = 0
    for stats in list_of_stats:
        total += int(stats["Estimated Number of Cells"])
    list_of_stats[0]["Estimated Number of Cells"] = total

    # Write "Estimated Number of Cells" to summary.csv
    with open(f"{sample_name}.star.solo/{star_feature}/Summary.csv", "w") as f:
        f.write(f"Estimated Number of Cells,{list_of_stats[0]['Estimated Number of Cells']}\n")


def merge_cell_reads(star_dirs, star_feature, sample_name):
    """
    Function to merge CellReads.stats of multiple star outputs

    Args:
        star_dirs (list): List of paths to star output directories
        star_feature (str): Name of star feature used
        sample_name (str): Name of sample
    """
    # Get header from any of the CellReads.stats files
    with open(f"{star_dirs[0]}/{star_feature}/CellReads.stats") as f:
        first_line = f.readline()
        second_line = f.readline()

    # Write the headers to the merged CellReads.stats file
    f_merged_cell_reads = open(f"{sample_name}.star.solo/{star_feature}/CellReads.stats", "w")
    f_merged_cell_reads.write(first_line+second_line)
    
    # Iterate over all the star directories
    for idx, star_dir in enumerate(star_dirs):
        f = open(f"{star_dir}/{star_feature}/CellReads.stats")
        f.readline()
        f.readline()
        # Write the contents of each CellReads.stats file to the merged CellReads.stats file
        for line in f:
            f_merged_cell_reads.write(line)
        f.close()
    f_merged_cell_reads.close()

def main():
    parser = argparse.ArgumentParser(
        description="Merge star solo outputs into one")
    parser.add_argument("--star_dirs", nargs='+', type=str,
                        help="star output directories that need to be concatenated",
                        required=True)
    parser.add_argument("--star_matrix", type=str, default='matrix.mtx')
    parser.add_argument("--star_feature", type=str)
    parser.add_argument("--sample_name", type=str)
    args = parser.parse_args()

    merged_star_path = Path(f"{args.sample_name}.star.solo/{args.star_feature}/raw")
    merged_star_path.mkdir(parents=True, exist_ok=True)
    
    merge_matrix_and_barcode(args.star_dirs, args.star_feature, args.star_matrix, args.sample_name)
    merge_cell_reads(args.star_dirs, args.star_feature, args.sample_name)


if __name__ == "__main__":
    main()
