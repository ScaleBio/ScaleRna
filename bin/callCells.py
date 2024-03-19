#!/usr/bin/env python
"""
Filter barcodes to call cells and generate filtered expression matrix. 
"""

import numpy as np
import pandas as pd
import os
import argparse
from pathlib import Path
import shutil
import sys
from scipy.io import mmread
from scipy.stats import median_abs_deviation
from scipy.sparse import csc_matrix
from scale_utils import io
from scale_utils.stats import extrapolate_unique
from dataclasses import dataclass
from cellFinder import construct_ambient_profile, test_ambiguous_barcodes 

@dataclass(frozen=True)
class CellCallingOptions:
    """
    This class holds options for cell calling.

    Attributes:
        fixedCells: If true the top fixedCells ranked by UTC are called as cells.
        expectedCells: Expected number of cells to call.
        topCellPercent: Percentile of cells sorted descending by UMI counts to help calculate the unique transcript counts threshold.
        minCellRatio: Minimum ratio of unique transcript counts to top cell UMI counts to help calculate the unique transcript counts threshold.
        minUTC: The minimum number of unique transcript counts a barcode must be associated with to be considered a potential cell.
        UTC: The minimum number of unique transcript counts a barcode must be associated with to be considered a called cell.
        cellFinder: If true Cell Finder will be used for cell calling.
        FDR: False discovery rate to use for rescuing cells based on deviation from ambient profile.
    """
    fixedCells: bool
    expectedCells: int
    topCellPercent: int
    minCellRatio: int
    minUTC: int
    UTC: int
    cellFinder: bool
    FDR: float

@dataclass(frozen=True)
class OutlierOptions:
    """
    This class holds options for outlier filtering.

    Attributes:
        filter_outliers: Whether to filter out flagged cells.
        num_mad_genes: Flag cells with gene counts more than num_mad_genes absolute deviations of the median.
        num_mad_umis: Flag cells with umi counts more than num_mad_umis absolute deviations of the median.
        num_mad_mito: Flag cells with higher than num_mad_mito absolute deviations of the median mitochondrial read percentage.
    """
    filter_outliers: bool
    num_mad_genes: int
    num_mad_umis: int
    num_mad_mito: int

def sort_barcodes(total_counts: np.ndarray, options: CellCallingOptions) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    (Subroutine of call_cells). Return arrays of barcode indices that are ambient, ambiguous, and called as cells.

    Args:
        total_counts: Total unique transcript counts for all barcodes
        options: Options for cell calling
    
    Returns:
        An array of indices for ambient barcodes
        An array of indices for ambiguous barcodes
        An array of indices for cell barcodes
    """

    ambiguous_barcode_indices = []
    ambient_barcode_indices = np.where(total_counts < options.minUTC)[0]
    threshold = options.UTC or calculate_UTC_threshold(total_counts, ambient_barcode_indices, options=options)

    #CellFinder Cell Thresholding
    if options.cellFinder:
        cell_barcode_indices = np.where(total_counts >= threshold)[0]
        ambiguous_barcode_indices = np.setdiff1d(np.setdiff1d(np.arange(len(total_counts)), ambient_barcode_indices), cell_barcode_indices)
        return ambient_barcode_indices, ambiguous_barcode_indices, cell_barcode_indices
    #Fixed Cell Thresholding
    elif options.fixedCells:
        if options.expectedCells < 0:
            raise ValueError('--expectedCells must be a positive integer.')
        cell_barcode_indices = np.argsort(total_counts)[::-1][:options.expectedCells]
        ambient_barcode_indices = np.setdiff1d(np.arange(len(total_counts)), cell_barcode_indices)
        return ambient_barcode_indices, ambiguous_barcode_indices, cell_barcode_indices
    #UTC Cell Thresholding
    else:
        cell_barcode_indices = np.where(total_counts >= threshold)[0]
        ambient_barcode_indices = np.setdiff1d(np.arange(len(total_counts)), cell_barcode_indices)
        return ambient_barcode_indices, ambiguous_barcode_indices, cell_barcode_indices

def calculate_UTC_threshold(total_counts: np.ndarray, ambient_barcode_indices: np.ndarray, options: CellCallingOptions) -> int:
    """
    (Subroutine of sort_barcodes). Calculates the unique transcripts counts
    threshold above which all barcodes will be called as cells.

    Args:
        total_counts: Total unique transcript counts for all barcodes
        ambient_barcode_indices: Array of indices for barcodes below minUTC
        options: Options for cell calling

    Returns: 
        The unique transcript counts threshold
    """
    
    if options.topCellPercent < 1 or options.topCellPercent > 99:
        raise ValueError('--topCellPercent must fall in the range 1 <= x <= 99.')
    if options.minCellRatio < 1:
        raise ValueError('--minCellRatio must be >= 1.')     
    if options.expectedCells > 0 and options.expectedCells < len(total_counts):
        # when --expectedCells is set, use the top --expectedCells barcodes as preliminary cells
        preliminary_cell_unique_transcript_counts = total_counts[np.argsort(total_counts)[::-1]][:options.expectedCells]
    else:
        preliminary_cell_unique_transcript_counts = total_counts[np.setdiff1d(np.arange(len(total_counts)), ambient_barcode_indices)]

    # Calculate unique transcripts counts threshold
    threshold = 0
    if len(preliminary_cell_unique_transcript_counts) != 0:
        # Calculate the UTC threshold as the UMI count associated with the --topCellPercent of the preliminary cells divided by --minCellRatio.
        threshold = np.round(np.percentile(preliminary_cell_unique_transcript_counts, options.topCellPercent) / options.minCellRatio)
    
    # Don't return threshold below minUTC
    return max(threshold, options.minUTC)

def classify_barcodes(sample_specific_file_paths: dict[str, Path], sampleMetrics: Path, options: CellCallingOptions) -> pd.DataFrame:
    """
    Call cells in this sample using thresholding options, including optionally
    rescuing ambiguous barcodes as cells with the CellFinder algorithm

    Args:
        sample_specific_file_paths: Dictionary where each key is a custom file identifier and each value is a path to the identified file in the STARsolo output directory
        sampleMetrics: Path to the allCells.csv containing metrics computed across all barcodes in this sample
        options: Options for cell calling
        
    Returns:
        The allCells.csv for this sample, with updated columns to show FDR (optionally) and classification (i.e. filters failed, optionally)
    """
    all_cells = pd.read_csv(sampleMetrics, index_col=0)
    mtx = csc_matrix(mmread(sample_specific_file_paths['mtx']))
    total_counts = all_cells['umis'].to_numpy()

    ambient_barcode_indices, ambiguous_barcode_indices, cell_barcode_indices = sort_barcodes(total_counts, options=options)

    # If there are ambiguous cells that can be rescued by CellFinder, run CellFinder
    if options.cellFinder:
        all_cells['FDR'] = 0.0

    if len(ambiguous_barcode_indices) > 0:
        discrete_mtx = mtx[:]

        # Make gene counts discrete so Good-Turing algorithm can work
        discrete_mtx.data = np.round(discrete_mtx.data)

        # Disregard irrelevant features; i.e. genes with zero counts
        discrete_mtx = discrete_mtx[np.where(discrete_mtx.sum(axis = 1).A > 0)[0], :]

        # Construct ambient profile to test barcodes against
        ambient_profile = construct_ambient_profile(discrete_mtx, ambient_barcode_indices)

        # Add FDR of each barcode to allCells.csv; calling barcodes with FDR <= --FDR as cells
        FDRs = test_ambiguous_barcodes(discrete_mtx, ambient_profile, ambient_barcode_indices, ambiguous_barcode_indices)
        FDR_column = mtx.shape[1] * [1.0]
        for i in range(len(ambiguous_barcode_indices)):
            FDR_column[ambiguous_barcode_indices[i]] = FDRs[i]
        for i in range(len(cell_barcode_indices)):
            FDR_column[cell_barcode_indices[i]] = 0.0
        all_cells['FDR'] = FDR_column
        all_cells['classification'] = ['cell' if barcode == True else 'ambient barcode' for barcode in list(all_cells['FDR'] <= options.FDR)]

    else:
        classification_column = mtx.shape[1] * ['ambient barcode']
        for cell_barcode_index in cell_barcode_indices:
            classification_column[cell_barcode_index] = 'cell'
        all_cells['classification'] = classification_column

    # Call cells according to the specified filters
    all_cells['pass'] = all_cells['classification'] != 'ambient barcode'
    return all_cells
    
def filter_cells(all_cells: pd.DataFrame, sample: str, options: OutlierOptions) -> pd.DataFrame:
    """
    Filter cells based on specified median absolute deviation options
    for gene counts, UMI counts, and mitochondrial read percentage.
    all_cells (classification and pass columns) is updated in place.

    Args:
        all_cells: A DataFrame of metrics for all barcodes
        sample: Unique string to identify this sample
        options: Options for outlier filtering
        
    Returns:
        MAD statistics used for flagging cells
    """

    # Fetch group of cells over which to identify outliers (according to median absolute deviation) to filter out
    cells = all_cells[all_cells['classification'] == 'cell']

    # Flag cells with low or high gene counts relative to other cells
    gene_counts = cells['genes']
    minimum_gene_count = np.median(gene_counts) - options.num_mad_genes * median_abs_deviation(gene_counts)
    maximum_gene_count = np.median(gene_counts) + options.num_mad_genes * median_abs_deviation(gene_counts)
    all_cells['classification'] += np.where((all_cells['classification'] == 'cell') & (all_cells['genes'] < minimum_gene_count), '; low genes', '')
    all_cells['classification'] += np.where((all_cells['classification'] == 'cell') & (all_cells['genes'] > maximum_gene_count), '; high genes', '')

    # Flag cells with low or high UMI counts relative to other cells
    umi_counts = cells['umis']
    minimum_umi_count = np.median(umi_counts) - options.num_mad_umis * median_abs_deviation(umi_counts)
    maximum_umi_count = np.median(umi_counts) + options.num_mad_umis * median_abs_deviation(umi_counts)
    all_cells['classification'] += np.where((all_cells['classification'] == 'cell') & (all_cells['umis'] > maximum_umi_count), '; high umis', '')
    all_cells['classification'] += np.where((all_cells['classification'] == 'cell') & (all_cells['umis'] < minimum_umi_count), '; low umis', '')

    # Flag cells with high percentage of Mitochondrial reads relative to other cells.
    mito_pcts = cells['mitoProp']
    median_mito_pct = np.median(mito_pcts)
    median_mito_abs_dev = median_abs_deviation(mito_pcts)
    maximum_mito_pct = median_mito_pct + options.num_mad_mito * median_mito_abs_dev
    
    # Flag cells with high percentages of mitochondrial reads relative to other cells
    all_cells['classification'] += np.where((all_cells['classification'] == 'cell') & (all_cells['mitoProp'] > maximum_mito_pct), '; mitochondrial read percentage', '')
    
    if options.filter_outliers:
        all_cells['pass'] = all_cells['classification'] == 'cell'

    stats = [
        ('MAD', 'minimum_gene_count', minimum_gene_count),
        ('MAD', 'maximum_gene_count', maximum_gene_count),
        ('MAD', 'minimum_umi_count', minimum_umi_count),
        ('MAD', 'maximum_umi_count', maximum_umi_count),
        ('MAD', 'maximum_mito_pct', maximum_mito_pct)
    ]
    return pd.DataFrame.from_records(stats, columns=['Category', 'Metric', 'Value'])


def calculate_sample_stats(all_cells: pd.DataFrame, internal_report: bool) -> pd.DataFrame:
    """
    Calculate sample-level statistics from the all_cells dataframe

    Args:
        all_cells: A dataframe containing metrics for all cells in sample
        internal_report: Whether to generate stats consistent with internal report

    Returns:
        Sample-level statistics
    """
    sample_stats = {'Reads': {}, 'Cells': {}, 'Extrapolated Complexity': {}}

    reads = sample_stats['Reads']
    reads['total_reads'] = all_cells['reads'].sum()
    reads['passing_reads'] = all_cells['passingReads'].sum()
    reads['mapped_reads'] = all_cells['mappedReads'].sum()
    reads['mapped_reads_perc'] = reads['mapped_reads'] / reads['total_reads']
    reads['gene_reads'] = all_cells['geneReads'].sum()
    reads['gene_reads_perc'] = reads['gene_reads'] / reads['mapped_reads']
    reads['exon_reads'] = all_cells['exonReads'].sum()
    reads['exon_reads_perc'] = reads['exon_reads'] / reads['mapped_reads']
    reads['antisense_reads'] = all_cells['antisenseReads'].sum()
    reads['antisense_reads_perc'] = reads['antisense_reads'] / reads['mapped_reads']
    reads['mito_reads'] = all_cells['mitoReads'].sum()
    reads['mito_reads_perc'] = reads['mito_reads'] / reads['mapped_reads']
    reads['umis'] = all_cells['umis'].sum()
    reads['saturation'] = 1 - reads['umis'] / reads['passing_reads']
    
    cells = sample_stats['Cells']
    passing_cells = all_cells.loc[all_cells['pass']]
    cells['cells_called'] = passing_cells.shape[0]
    cells['utc_threshold'] = passing_cells['umis'].min()
    cells['mean_passing_reads'] = passing_cells['reads'].mean()
    cells['median_reads'] = passing_cells['passingReads'].median()
    cells['median_utc'] = passing_cells['umis'].median()
    cells['median_genes'] = passing_cells['genes'].median()
    cells['reads_in_cells'] = passing_cells['passingReads'].sum() / all_cells['passingReads'].sum()
    cells['mito_reads'] = passing_cells['mitoReads'].sum()
    
    complexity = sample_stats['Extrapolated Complexity']
    unique_reads = []
    target_mean_reads = [100, 500, 1000, 5000, 10000, 20000]
    for d in target_mean_reads:
        if d < cells['mean_passing_reads'] or internal_report:
            # Target reads refers to the mean reads per cell, but we are estimating
            # UMIs / complexity in the median cell. Approx. scaling median reads by the same
            # factor as mean reads.
            target_median_reads = (d/cells['mean_passing_reads']) * cells['median_reads']
            unique_reads = extrapolate_unique(cells['median_reads'], cells['median_utc'], target_median_reads)
            complexity[f'target_reads_{d}'] = unique_reads
        else:
            complexity[f'target_reads_{d}'] = np.nan
    
    sample_stats_df = pd.DataFrame([], columns=["Category", "Metric", "Value"])
    for category, metrics in sample_stats.items():
        sample_stats_df = pd.concat([sample_stats_df, pd.DataFrame([{"Category": category, "Metric": metric, "Value": value} for metric, value in metrics.items()])])
    return sample_stats_df


def generate_filtered_matrix(sample_specific_file_paths: dict[str, Path], all_cells: pd.DataFrame, sample: str) -> None:
    """
    Generates and writes a filtered cell-by-gene expression matrix. 

    Args:
        sample_specific_file_paths: Dictionary where each key is a custom file identifier and each value is a path to the identified file in the STARsolo output directory
        all_cells: An allCells.csv containing metrics computed across all barcodes in this sample 
        sample: Unique string to identify this sample
    
    Returns:
        None. A filtered matrix is written to the {sample}_filtered_star_output directory
    """
    cells_passing = all_cells['pass'].to_dict()
    filtered_path = f"{sample}_filtered_star_output"

    # Create /filtered directory under current directory
    Path(".", filtered_path).mkdir()

    # Features.tsv is the same; we move it to that /filtered directory
    shutil.copyfile(sample_specific_file_paths['features'], f"{filtered_path}/features.tsv")

    # File pointer to raw matrix.mtx and barcodes.tsv
    f_raw_mtx = open(sample_specific_file_paths['mtx'])
    f_raw_barcodes = open(sample_specific_file_paths['barcodes'])

    # File pointers to newly created filtered matrix.mtx and barcodes.tsv
    with (open("tmp_matrix.mtx", "w") as f_filtered_mtx_tmp, open(f"{filtered_path}/barcodes.tsv", "w") as f_filtered_barcodes):

        # Read in barcodes from raw barcodes file
        raw_barcodes = f_raw_barcodes.readlines()
        raw_barcodes = [line.rstrip() for line in raw_barcodes]

        # Data type header
        header = f_raw_mtx.readline().strip()

        # Skip extra header lines
        f_raw_mtx.readline(); f_raw_mtx.readline() 

        # Set current barcode to 0 to ensure that the first time we compare current_barcode to the barcode received from iterating over the raw matrix, we set current_barcode to that barcode
        current_barcode = 0

        # Row number in barcodes.tsv
        barcode_row_number = 1

        # Count number of lines for writing as header later on
        line_count = 0

        # Iterate over every line in raw matrix.mtx
        for line in f_raw_mtx:
            split_line = line.split()
            feature, barcode, count = int(split_line[0]), int(split_line[1]), float(split_line[2])
            barcode_sequence = raw_barcodes[barcode - 1]

            # If barcode has passing cells in the all_cells dataframe, write to file
            if cells_passing[barcode_sequence]:
                if barcode != current_barcode:
                    current_barcode = barcode
                    barcode_row_number += 1
                    f_filtered_barcodes.write(f"{raw_barcodes[current_barcode - 1]}\n")
                
                f_filtered_mtx_tmp.write(f"{feature} {barcode_row_number-1} {count}\n")
                line_count += 1

    f_raw_mtx.close()
    f_raw_barcodes.close()
    f_filtered_barcodes.close()

    # Compute header information for the matrix; first entry is length of filtered features.tsv
    header1 = len(pd.read_csv(sample_specific_file_paths['features'], sep = "\t", header = None).index)

    # Second entry is length of filtered barcodes.tsv
    header2 = barcode_row_number - 1

    # Third entry is length of filtered matrix.mtx
    header3 = line_count

    with open(f"{filtered_path}/matrix.mtx", "w") as f_filtered_mtx:
        f_filtered_mtx.write(f"{header}\n%\n")
        f_filtered_mtx.write(f"{header1} {header2} {header3}\n")
    if line_count > 0:
        os.system(f"cat tmp_matrix.mtx >> {filtered_path}/matrix.mtx")
    try:
        if os.path.isfile("tmp_matrix.mtx"):
            os.remove("tmp_matrix.mtx")
    except Exception as e:
        print(e, file = sys.stderr)

def main():
    parser = argparse.ArgumentParser()
    
    # Required argument for specifying the path to the sampleMetrics CSV file
    parser.add_argument("--sampleMetrics", type = Path, required = True, help = "Path to the sampleMetrics CSV file.") 

    # Required and optional arguments for specifying the STARsolo outputs for this sample
    parser.add_argument("--STARsolo_out", type = Path, required = True, help = "Path to the STARsolo outputs for this sample.")
    parser.add_argument("--feature_type", type = str, required = False, default = 'GeneFull_Ex50pAS', help = "STARsolo feature type used.")
    parser.add_argument("--matrix_type", type = str, required = False, default = 'UniqueAndMult-PropUnique.mtx', help = "STARsolo matrix type used.")

    # Optional argument to specify the name of the sample for which cells are being called
    parser.add_argument("--sample", type = str, required = False, default = "example", help = "Unique string to identify this sample.")

    # Optional argument to set hard thresholds for cell calling 
    parser.add_argument("--fixedCells", required = False, action="store_true", default=False, help = "Fixed number of barcodes to call as cells.")

    # Optional arguments for use of CellFinder algorithm
    parser.add_argument("--expectedCells", type = int, required = False, default = 0, help = "Expected number of cells to call (if specified in samples.csv; see algorithm description).")
    parser.add_argument("--topCellPercent", type = int, required = False, default = 99, help = "Cell thresholding parameter (see algorithm description).")
    parser.add_argument("--minCellRatio", type = int, required = False, default = 10, help = "Cell thresholding parameter (see algorithm description).")
    parser.add_argument("--minUTC", type = int, required = False, default = 100, help = "The minimum number of unique transcript counts a barcode must be associated with to be considered a potential cell.")
    parser.add_argument("--UTC", type = int, required = False, default = 0, help = "The minimum number of unique transcript counts a barcode must be associated with to be considered a called cell.")
    parser.add_argument("--cellFinder", required = False, action = "store_true", help = "If set, use CellFinder to call cells from among barcodes with unique transript counts falling between --minUTC and the UTC (set with --UTC or calculated algorithmically when --UTC <= 0).")
    parser.add_argument("--FDR", type = float, required = False, default = 0.01, help = "False discovery rate at which to call cells (see algorithm description).")

    # Optional arguments for MAD outlier filtering
    parser.add_argument("--filter_outliers", required = False, action="store_true", help = "Number of median absolute deviations in gene count/UMI count/mitochondrial read percentage above/below which a cell will be flagged as an outlier.")
    parser.add_argument("--num_mad_genes", type=int, required = False, action = "store", default= 5, help = "If set, cells with gene counts below --MADs absolute deviations of the median will be filtered out.")
    parser.add_argument("--num_mad_umis", type=int, required = False, action = "store", default= 5, help = "If set, cells with UMI counts above --MADs absolute deviations of the median will be filtered out.")
    parser.add_argument("--num_mad_mito", type=int, required = False, action = "store", default= 3, help = "If set, cells with mitochondrial read percentages above --MADs absolute deviations of the median will be filtered out.")
    
    # Optional argument to specify whether to generate statistics for internal report
    parser.add_argument("--internalReport", action="store_true", default=False)

    args = parser.parse_args()

    sample_specific_file_paths = io.resolve_sample_specific_file_paths(args.STARsolo_out, args.feature_type, args.matrix_type)
    
    call_cells_options = CellCallingOptions(
        fixedCells=args.fixedCells,
        expectedCells=args.expectedCells,
        topCellPercent=args.topCellPercent,
        minCellRatio=args.minCellRatio,
        minUTC=args.minUTC,
        UTC=args.UTC,
        cellFinder=args.cellFinder,
        FDR=args.FDR
    )
    
    outlier_options = OutlierOptions(
        filter_outliers=args.filter_outliers,
        num_mad_genes=args.num_mad_genes,
        num_mad_umis=args.num_mad_umis,
        num_mad_mito=args.num_mad_mito
    )
    all_cells = classify_barcodes(sample_specific_file_paths, args.sampleMetrics, options=call_cells_options)
    mad_stats = filter_cells(all_cells, args.sample, options=outlier_options)
    sample_stats = calculate_sample_stats(all_cells, internal_report=args.internalReport)
    sample_stats = pd.concat([sample_stats, mad_stats])

    # Write filtered matrix for this sample
    generate_filtered_matrix(sample_specific_file_paths, all_cells, args.sample)

    metrics_dir = Path(".", f"{args.sample}_metrics")
    metrics_dir.mkdir(parents = True)

    # Write allCells.csv for this sample
    all_cells.to_csv(metrics_dir / f"{args.sample}_allCells.csv")
    sample_stats.to_csv(metrics_dir / f"{args.sample}_sample_stats.csv", index=False)

if __name__ == "__main__":
    main()
