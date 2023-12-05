#!/usr/bin/env python
"""
Python script to merge multiple STAR outputs
"""
import argparse
import pandas as pd
import shutil
import os
import sys
from pathlib import Path


def merge_matrix_and_barcode(star_dirs, star_feature, star_matrix, sample_ids, out_dir):
    """
    Function to merge barcodes.tsv and matrix.mtx of each star output directory

    Args:
        star_dirs (list): List containing paths to star output directories
        star_feature (str): Name of star feature used
        star_matrix (str): Filename of the STAR matrix output to use
        sample_ids: optional list of names for each input sample (star_dir)
    
    Output:
        Writes raw matrix.mtx and raw barcodes.tsv
    Returns:
        List of folder names one level down from @star_dirs
    """
    # Generate combined list of barcodes from individual barcodes files
    
    # List of list of barcodes across all samples
    # For each barcode we store the new name (original name + sample_id postfix to make unique if used)
    # and the new barcode id number (unique across all merged samples; 1-based)
    combined_barcodes = [] 
    all_barcodes = set() # Used to check for uniqueness
    for idx, star in enumerate(star_dirs):
        combined_barcodes.append([])
        for bc in open(star / star_feature / "raw" / "barcodes.tsv"):
            bc = bc.strip()
            if sample_ids:
                bc = f"{bc}_{sample_ids[idx]}"
            if bc in all_barcodes:
                raise ValueError("Duplicate barcode in merged STARSolo files: " + bc)
            all_barcodes.add(bc)
            # Append barcode and id (1-based) to barcode list for current sample
            combined_barcodes[-1].append( (bc, len(all_barcodes)) )

    line_count = 0
    # Write merged matrix contents (without header) to a tmp file
    # (We only know the correct size afterwards)
    tmp_matrix_fn = "tmp_matrix.mtx"
    feature_count = None # Number of genes in the annotation (same for all samples)
    with open(tmp_matrix_fn, "w") as tmp_mtx:
        # Iterate over all the star directories and merge the matrices
        for idx, star_dir in enumerate(star_dirs):
            f = open(star_dir / star_feature / "raw" / star_matrix, "r")
            header = f.readline().strip() # Datatype header
            f.readline() # Skip '%' line
            counts = f.readline().split() # Matrix size header
            if not feature_count:
                feature_count = int(counts[0])
            elif feature_count != int(counts[0]):
                raise ValueError("Matrix dimension mismatch for " + star_dir)
            barcode_count = int(counts[1])
            if barcode_count > len(combined_barcodes[idx]):
                raise ValueError("Matrix dimension mismatch for " + star_dir)
            
            for line in f:
                split_line = line.split()
                feature, barcode_id, umi = int(split_line[0]), int(split_line[1]), float(split_line[2])
                new_bc, new_id = combined_barcodes[idx][barcode_id-1]
                line_count += 1
                tmp_mtx.write(f"{feature} {new_id} {umi}\n")
    # Write merged matrix header and then append contents
    mtx_dir = out_dir / "raw"
    mtx_dir.mkdir(parents=True)
    merged_mtx_fn = mtx_dir / star_matrix
    total_barcode_count = sum(len(bcs) for bcs in combined_barcodes)
    with open(merged_mtx_fn, "w") as merged_mtx:
        merged_mtx.write(f"{header}\n%\n")
        merged_mtx.write(f"{feature_count} {total_barcode_count} {line_count}\n")
    os.system(f"cat {tmp_matrix_fn} >> {merged_mtx_fn}")

    # Write concatenated list of barcodes
    merged_barcode_fn = mtx_dir / "barcodes.tsv"
    with open(merged_barcode_fn, "w") as f:
        for barcodes in combined_barcodes:
            for bc,bc_id in barcodes:
                f.write(bc+"\n")
    
    # Copy features.tsv from any of the star directories since it'll be same for all
    merged_feature_fn = mtx_dir / "features.tsv"
    shutil.copyfile(f"{star_dirs[0]}/{star_feature}/raw/features.tsv", merged_feature_fn)
    
    # We just need "Estimated Number of Cells" from summary.csv to calculate cell thresholding
    # in a downstream process. The code block below gets all metrics from the individual summary.csv
    # files and just concatenates the "Estimated Number of Cells" metric
    list_of_stats = []
    for idx, star_dir in enumerate(star_dirs):
        stats = {}
        with open(star_dir / star_feature / "Summary.csv") as f:
            for line in f:
                split_line = line.strip().split(',')
                stats[split_line[0]] = split_line[1]
            list_of_stats.append(stats)
    total = 0
    for stats in list_of_stats:
        total += int(stats["Estimated Number of Cells"])
    list_of_stats[0]["Estimated Number of Cells"] = total

    # Write "Estimated Number of Cells" to summary.csv
    stats_fn = out_dir / "Summary.csv"
    
    with open(stats_fn, "w") as f:
        f.write(f"Estimated Number of Cells,{list_of_stats[0]['Estimated Number of Cells']}\n")


def merge_cell_reads(star_dirs, star_feature, sample_ids, out_dir):
    """
    Function to merge CellReads.stats of multiple star outputs

    Args:
        star_dirs (list): List of paths to star output directories
        star_feature (str): Name of star feature used
        sample_name (str): Name of sample
    """
    # Get header from any of the CellReads.stats files
    headers = []
    with open(star_dirs[0] / star_feature / "CellReads.stats") as f:
        headers.append(f.readline())
        headers.append(f.readline())
    # Write the merged CellReads.stats file
    with open(out_dir / "CellReads.stats", "w") as merged_cellreads:
        for h in headers:
            merged_cellreads.write(h)
        for idx, star_dir in enumerate(star_dirs):
            f = open(star_dir / star_feature / "CellReads.stats")
            f.readline(); f.readline() # skip reader
            # Append contents of each CellReads.stats file to the merged file
            for line in f:
                if sample_ids: # Append sample_id to each cell-barcode if given
                    line = line.split("\t")
                    line[0] = f"{line[0]}_{sample_ids[idx]}"
                    line = "\t".join(line)
                merged_cellreads.write(line)


def main():
    parser = argparse.ArgumentParser(description="Merge star solo outputs into one")
    parser.add_argument("--star_dirs", nargs='+', type=Path, required=True,
                        help="star output directories that need to be concatenated")
    parser.add_argument("--sample_ids", type=str, default=None, nargs="+",
                        help="Sample IDs for each STAR output directory")
    parser.add_argument("--star_matrix", type=str, default='matrix.mtx',
                        help='.mtx file name')
    parser.add_argument("--star_feature", type=str,
                        help='Feature-type (e.g. GeneFull) in the STAR output to merge')
    parser.add_argument("--out_dir", type=Path, required=True,
                        help='merged sample output directory')
    args = parser.parse_args()
    
    if args.sample_ids and len(args.star_dirs) != len(args.sample_ids):
        raise ValueError("Number of sample IDs not equal to number of star output directories")

    merged_star_path = args.out_dir / args.star_feature
    merged_star_path.mkdir(parents=True, exist_ok=True)
    merge_matrix_and_barcode(args.star_dirs, args.star_feature, args.star_matrix, args.sample_ids, merged_star_path)
    merge_cell_reads(args.star_dirs, args.star_feature, args.sample_ids, merged_star_path)


if __name__ == "__main__":
    main()
