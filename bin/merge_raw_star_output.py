#!/usr/bin/env python
"""
Python script to merge multiple STAR outputs
"""
import argparse
import shutil
import gzip
from pathlib import Path
from scale_utils import stats


def merge_matrix_and_barcode(star_dirs: list, star_feature: str, star_matrix: str, sample_ids: list, out_dir: Path):
    """
    Function to merge barcodes.tsv and matrix.mtx of each star output directory

    Args:
        star_dirs: List containing paths to star output directories
        star_feature: Name of star feature used
        star_matrix: Filename of the STAR matrix output to use
        sample_ids: List of names for each input sample (star_dir)
        out_dir: Output directory for the merged matrix.mtx and barcodes.tsv

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
    all_barcodes = set()  # Used to check for uniqueness
    for idx, star in enumerate(star_dirs):
        combined_barcodes.append([])
        with gzip.open(star / star_feature / "raw" / "barcodes.tsv.gz", "rt") as f:
            for bc in f:
                bc = bc.strip()
                if sample_ids:
                    bc = f"{bc}_{sample_ids[idx]}"
                if bc in all_barcodes:
                    raise ValueError("Duplicate barcode in merged STARSolo files: " + bc)
                all_barcodes.add(bc)
                # Append barcode and id (1-based) to barcode list for current sample
                combined_barcodes[-1].append((bc, len(all_barcodes)))

    line_count = 0
    feature_count = None  # Number of genes in the annotation (same for all samples)
    # Write merged matrix
    mtx_dir = out_dir / "raw"
    mtx_dir.mkdir(parents=True)
    merged_mtx_fn = mtx_dir / star_matrix
    total_barcode_count = sum(len(bcs) for bcs in combined_barcodes)
    with open(merged_mtx_fn, "wt") as merged_mtx:
        # Iterate over all the star directories and merge the matrices
        for idx, star_dir in enumerate(star_dirs):
            with gzip.open(star_dir / star_feature / "raw" / f"{star_matrix}.gz", "rt") as f:
                header = f.readline().strip()  # Datatype header
                f.readline()  # Skip '%' line
                counts = f.readline().split()  # Matrix size header
                if not feature_count:
                    feature_count = int(counts[0])
                elif feature_count != int(counts[0]):
                    raise ValueError("Matrix dimension mismatch for " + star_dir)
                barcode_count = int(counts[1])
                if barcode_count > len(combined_barcodes[idx]):
                    raise ValueError("Matrix dimension mismatch for " + star_dir)

                # Write matrix header
                if idx == 0:
                    merged_mtx.write(f"{header}\n%\n")
                    pos = merged_mtx.tell()
                    merged_mtx.write(f"{feature_count} {total_barcode_count} {'0'.zfill(13)}\n")

                for line in f:
                    split_line = line.split()
                    feature, barcode_id, umi = int(split_line[0]), int(split_line[1]), split_line[2]
                    _, new_id = combined_barcodes[idx][barcode_id - 1]
                    line_count += 1
                    merged_mtx.write(f"{feature} {new_id} {umi}\n")
        merged_mtx.seek(pos)
        merged_mtx.write(f"{feature_count} {total_barcode_count} {str(line_count).zfill(13)}\n")

    # Write concatenated list of barcodes
    merged_barcode_fn = mtx_dir / "barcodes.tsv.gz"
    with gzip.open(merged_barcode_fn, "wt") as f:
        for barcodes in combined_barcodes:
            for bc, _ in barcodes:
                f.write(bc + "\n")

    # Copy features.tsv from any of the star directories since it'll be same for all
    merged_feature_fn = mtx_dir / "features.tsv.gz"
    shutil.copyfile(f"{star_dirs[0]}/{star_feature}/raw/features.tsv.gz", merged_feature_fn)


def merge_cell_reads(star_dirs: list, star_feature: str, sample_ids: str, out_dir: Path):
    """
    Function to merge CellReads.stats of multiple star outputs

    Args:
        star_dirs: List of paths to star output directories
        star_feature: Name of star feature used
        sample_name: Name of sample
        out_dir: Output directory for the merged CellReads.stats file

    Output:
        Writes merged CellReads.stats file to out_dir
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
            f.readline()
            f.readline()  # skip reader
            # Append contents of each CellReads.stats file to the merged file
            for line in f:
                if sample_ids:  # Append sample_id to each cell-barcode if given
                    line = line.split("\t")
                    line[0] = f"{line[0]}_{sample_ids[idx]}"
                    line = "\t".join(line)
                merged_cellreads.write(line)


def merge_log_files(star_logs: list, log_out_dir: Path):
    """
    Function to merge Log.final.out file of multiple STAR outputs to get a single combined file

    Args:
        star_logs: Log.final.out files
        log_out_dir: Output directory for the merged log file

    Output:
        Writes a single merged Log.out file to log_out_dir
    """
    input_reads = []
    input_read_length = []
    mapped_length = []
    mismatch_percents = []
    err_short = 0
    err_mismatches = 0
    err_other = 0
    mapped_reads = 0
    mapped_too_many_loci = 0
    mapped_multiple_loci = 0
    for star_log in star_logs:
        with open(star_log) as f:
            for line in f:
                if "Number of input reads" in line:
                    input_reads.append(float(line.split("Number of input reads |")[1].strip()))
                if "Average input read length" in line:
                    input_read_length.append(float(line.split("Average input read length |")[1].strip()))
                if "Average mapped length" in line:
                    mapped_length.append(float(line.split("Average mapped length |")[1].strip()))
                if "Number of reads unmapped: too many mismatches" in line:
                    err_mismatches += float(line.split("Number of reads unmapped: too many mismatches |")[1].strip())
                if "Number of reads unmapped: too short" in line:
                    err_short += float(line.split("Number of reads unmapped: too short |")[1].strip())
                if "Number of reads unmapped: other" in line:
                    err_other += float(line.split("Number of reads unmapped: other |")[1].strip())
                if "Mismatch rate per base, % |" in line:
                    mismatch_percents.append(
                        float(line.split("Mismatch rate per base, % |")[1].strip().split("%")[0]) / 100
                    )
                if "Uniquely mapped reads number |" in line:
                    mapped_reads += int(line.split("Uniquely mapped reads number |")[1].strip())
                if "Number of reads mapped to too many loci" in line:
                    mapped_too_many_loci += int(line.split("Number of reads mapped to too many loci |")[1].strip())
                if "Number of reads mapped to multiple loci" in line:
                    mapped_multiple_loci += int(line.split("Number of reads mapped to multiple loci |")[1].strip())

    with open(log_out_dir / "Log.final.out", "w") as f:
        avg_input_read_length = stats.calculate_weighted_average(input_reads, input_read_length, 1)
        avg_mapped_length = stats.calculate_weighted_average(input_reads, mapped_length, 1)
        avg_mismatch_percent = stats.calculate_weighted_average(input_reads, mismatch_percents, 4)
        f.write(f"Number of input reads | {sum(input_reads)}\n")
        f.write(f"Average input read length | {avg_input_read_length}\n")
        f.write(f"Average mapped length | {avg_mapped_length}\n")
        f.write(f"Number of reads unmapped: too many mismatches | {round(err_mismatches, 2)}\n")
        f.write(f"Number of reads unmapped: too short | {round(err_short, 2)}\n")
        f.write(f"Number of reads unmapped: other | {round(err_other, 2)}\n")
        f.write(f"Mismatch rate per base, % | {round(avg_mismatch_percent*100, 2)}%\n")
        f.write(f"Uniquely mapped reads number | {mapped_reads}\n")
        f.write(f"Number of reads mapped to multiple loci | {mapped_multiple_loci}\n")
        f.write(f"Number of reads mapped to too many loci | {mapped_too_many_loci}\n")
        f.write(f"% of reads mapped to too many loci | {round(mapped_too_many_loci / sum(input_reads) * 100, 2)}%\n")


def main():
    parser = argparse.ArgumentParser(description="Merge star solo outputs into one")
    parser.add_argument(
        "--star_dirs", nargs="+", type=Path, required=True, help="star output directories that need to be concatenated"
    )
    parser.add_argument(
        "--star_log", nargs="+", type=Path, required=True, help="STAR log files that need to be concatenated"
    )
    parser.add_argument(
        "--sample_ids", type=str, default=None, nargs="+", help="Sample IDs for each STAR output directory"
    )
    parser.add_argument("--star_matrix", type=str, default="matrix.mtx", help=".mtx file name")
    parser.add_argument("--star_feature", type=str, help="Feature-type (e.g. GeneFull) in the STAR output to merge")
    parser.add_argument("--out_dir", type=Path, required=True, help="merged sample output directory")
    parser.add_argument("--log_out_dir", type=Path, required=True, help="merged sample log output directory")
    args = parser.parse_args()

    if args.sample_ids and len(args.star_dirs) != len(args.sample_ids):
        raise ValueError("Number of sample IDs not equal to number of star output directories")

    merged_star_path = args.out_dir / args.star_feature
    merged_star_path.mkdir(parents=True, exist_ok=True)
    args.log_out_dir.mkdir(parents=True, exist_ok=True)
    merge_matrix_and_barcode(args.star_dirs, args.star_feature, args.star_matrix, args.sample_ids, merged_star_path)
    merge_cell_reads(args.star_dirs, args.star_feature, args.sample_ids, merged_star_path)
    merge_log_files(args.star_log, args.log_out_dir)


if __name__ == "__main__":
    main()
