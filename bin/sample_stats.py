#!/usr/bin/env python
"""
Create sample-level statistics from barcode metrics and STAR log
"""

import argparse
from pathlib import Path
import duckdb
import numpy as np
import pandas as pd
from scale_utils.stats import extrapolate_unique


def parse_star_log(reads: dict[str, float], star_log: Path) -> dict[str, float]:
    """
    Parse STAR alignment log file to extract additional read statistics

    Args:
        reads: Dictionary containing read statistics
        star_log: Path to the STARsolo log file

    Returns:
        Updated reads dictionary with additional statistics
    """
    with open(star_log, "r") as f:
        for line in f:
            if "Average input read length" in line:
                reads["avg_trimmed_read_len"] = float(line.split("Average input read length |")[1].strip())
            if "Average mapped length" in line:
                reads["avg_mapped_len"] = float(line.split("Average mapped length |")[1].strip())
            if "Mismatch rate per base, % |" in line:
                reads["mismatch_rate_per_base_perc"] = float(
                    line.split("Mismatch rate per base, % |")[1].strip().split("%")[0]
                )
            if "% of reads mapped to too many loci |" in line:
                reads["mapped_to_too_many_loci_perc"] = float(
                    line.split("% of reads mapped to too many loci |")[1].strip().split("%")[0]
                )
    return reads


def calculate_sample_stats(
    cell_metrics: list[Path],
    internal_report: bool,
    is_cellFinder: bool,
    filter_outliers: bool,
    is_barnyard: bool,
    total_sample_reads: int,
    star_log: Path,
) -> pd.DataFrame:
    """
    Calculate sample-level statistics from the all_cells dataframe

    Args:
        con: duckdb connection with all_barcodes table
        internal_report: Whether to generate stats consistent with internal report
        is_cellFinder: Was CellFinder used to call cells?
        filter_outliers: Were QC flagged celled filtered / failed?
        total_sample_reads: Total sample reads, obtained from bcParser; set to 0 if bcParser was not run
        star_log: Path to the STARsolo log file

    Returns:
        Sample-level statistics
    """
    sample_stats = {"Reads": {}, "Cells": {}, "Extrapolated Complexity": {}}
    barnyard_cols = ""
    if is_barnyard:
        barnyard_cols = """
        min(mouse_counts) FILTER (pass AND species == 'Mouse') AS mouse_utc,
        min(human_counts) FILTER (pass AND species == 'Human') AS human_utc,
        median(greatest(mouse_counts, human_counts)) FILTER (pass) AS median_utc,
        median(greatest(mouse_genes, human_genes)) FILTER (pass) AS median_genes,
        """
    else:
        barnyard_cols = """
        median(counts) FILTER (pass) AS median_utc,
        median(genes) FILTER (pass) AS median_genes,
        """

    query_results = (
        duckdb.sql(
            f"""
        SELECT
            sum(totalReads) AS total_reads,
            sum(countedReads) AS counted_reads,
            sum(mappedReads) AS mapped_reads,
            sum(geneReads) AS gene_reads,
            sum(exonReads) AS exon_reads,
            sum(antisenseReads) AS antisense_reads,
            sum(mitoReads) AS mito_reads,
            sum(counts) AS counts,
            count(cell_id) FILTER (pass) AS cells_called,
            min(counts) FILTER (pass) AS utc_threshold,
            count(cell_id) FILTER (pass AND flags LIKE '%cellFinder%') AS cellFinder_calls,
            count(cell_id) FILTER (flags LIKE '%low_%' OR flags LIKE '%high_%') AS qc_filtered,
            count(cell_id) FILTER (pass AND (flags LIKE '%low_%' OR flags LIKE '%high_%')) AS qc_flagged,
            mean(totalReads) FILTER (pass) AS mean_passing_reads,
            median(countedReads) FILTER (pass) AS median_reads,
            sum(countedReads) FILTER (pass) / sum(countedReads) AS reads_in_cells,
            sum(mitoReads) FILTER (pass) AS mito_reads_passing,
            {barnyard_cols}
        FROM read_parquet({[str(x) for x in cell_metrics]})
        """
        )
        .df()
        .to_dict(orient="index")[0]
    )

    reads = sample_stats["Reads"]
    if total_sample_reads != 0:
        reads["total_sample_reads_post_barcode_demux"] = total_sample_reads
    reads["total_reads"] = query_results["total_reads"]
    reads["counted_reads"] = query_results["counted_reads"]
    reads["mapped_reads"] = query_results["mapped_reads"]
    reads["mapped_reads_perc"] = reads["mapped_reads"] / reads["total_reads"]
    reads["gene_reads"] = query_results["gene_reads"]
    reads["gene_reads_perc"] = reads["gene_reads"] / reads["mapped_reads"]
    reads["exon_reads"] = query_results["exon_reads"]
    reads["exon_reads_perc"] = reads["exon_reads"] / reads["gene_reads"]
    reads["antisense_reads"] = query_results["antisense_reads"]
    reads["antisense_reads_perc"] = reads["antisense_reads"] / reads["mapped_reads"]
    reads["mito_reads"] = query_results["mito_reads"]
    reads["mito_reads_perc"] = reads["mito_reads"] / reads["mapped_reads"]
    reads["counts"] = query_results["counts"]
    reads["saturation"] = 1 - reads["counts"] / reads["counted_reads"]
    reads = parse_star_log(reads, star_log)

    cells = sample_stats["Cells"]
    cells["cells_called"] = query_results["cells_called"]
    cells["mean_passing_reads"] = query_results["mean_passing_reads"]
    if is_cellFinder:
        cells["cellFinder_calls"] = query_results["cellFinder_calls"]
    elif is_barnyard:
        utc_by = [i for i in [query_results["mouse_utc"], query_results["human_utc"]] if i is not None]
        cells["utc_threshold"] = np.min(utc_by) if utc_by else query_results["utc_threshold"]
    else:
        cells["utc_threshold"] = query_results["utc_threshold"]
    if filter_outliers:
        cells["qc_filtered"] = query_results["qc_filtered"]
    else:
        cells["qc_flagged"] = query_results["qc_flagged"]
    if cells["cells_called"] == 0:
        cells["reads_per_cell"] = 0
    else:
        # For ultima we compute reads per cell relative to reads going into alignment
        if total_sample_reads == 0:
            cells["reads_per_cell"] = reads["total_reads"] / cells["cells_called"]
        else:
            cells["reads_per_cell"] = total_sample_reads / cells["cells_called"]
    cells["median_reads"] = query_results["median_reads"]
    cells["median_utc"] = query_results["median_utc"]
    cells["median_genes"] = query_results["median_genes"]
    cells["reads_in_cells"] = query_results["reads_in_cells"]
    cells["mito_reads"] = query_results["mito_reads_passing"]

    if is_barnyard:
        sample_stats["Barnyard"] = {}
        barnyard = sample_stats["Barnyard"]
        barnyard["mouse_utc"] = query_results["mouse_utc"]
        barnyard["human_utc"] = query_results["human_utc"]

    complexity = sample_stats["Extrapolated Complexity"]
    unique_reads = []
    target_mean_reads = [100, 500, 1000, 5000, 10000, 20000]
    for d in target_mean_reads:
        if cells["reads_per_cell"] == 0 or (d >= cells["reads_per_cell"] and not internal_report):
            complexity[f"target_reads_{d}"] = np.nan
        else:
            # Target reads refers to the mean reads per cell, but we are estimating
            # UMIs / complexity in the median cell. Approx. scaling median reads by the same
            # factor as mean reads.
            target_median_reads = (d / cells["reads_per_cell"]) * cells["median_reads"]
            unique_reads = extrapolate_unique(cells["median_reads"], cells["median_utc"], target_median_reads)
            complexity[f"target_reads_{d}"] = unique_reads

    stats_stats_df = pd.DataFrame(
        [
            {"Category": category, "Metric": metric, "Value": value}
            for category, metrics in sample_stats.items()
            for metric, value in metrics.items()
        ]
    )

    return stats_stats_df


def main():
    parser = argparse.ArgumentParser()

    # Required argument for specifying the path to the sampleMetrics parquet file
    parser.add_argument(
        "--cellMetrics", nargs="+", type=Path, required=True, help="Path to the barcodes parquet file(s)."
    )

    # CellFinder parameters
    parser.add_argument(
        "--cellFinder",
        action="store_true",
        help="Use CellFinder to call cells from among barcodes with counts between --minUTC and the UTC threshold",
    )
    # MAD outlier filtering
    parser.add_argument(
        "--filterOutliers",
        required=False,
        action="store_true",
    )
    parser.add_argument("--starLog", type=Path, required=True, help="Path to the STARsolo log file.")
    parser.add_argument("--sample", type=str, required=True, help="Unique string to identify this sample.")
    parser.add_argument("--totalSampleReads", type=int, help="Total sample reads, obtained from bcParser")
    # Optional argument to specify whether to generate statistics for internal report
    parser.add_argument("--internalReport", action="store_true", default=False)
    parser.add_argument("--isBarnyard", action="store_true", default=False)
    parser.add_argument("--threads", type=int, required=False, default=1, help="Number of threads for duckdb")
    parser.add_argument("--memory", type=str, required=False, default="8 GB", help="Memory allocated to task")

    args = parser.parse_args()

    mem_limit, mem_unit = args.memory.split()
    # allocate duckdb memory based on number of threads
    mem_limit = f"{float(mem_limit) / (args.threads + 1):.1f}{mem_unit}"
    duckdb.sql(
        f"""
    SET threads TO {args.threads};
    SET memory_limit TO '{mem_limit}';
    """
    )

    sample_stats = calculate_sample_stats(
        cell_metrics=args.cellMetrics,
        internal_report=args.internalReport,
        is_cellFinder=args.cellFinder,
        filter_outliers=args.filterOutliers,
        is_barnyard=args.isBarnyard,
        total_sample_reads=args.totalSampleReads,
        star_log=args.starLog,
    )

    metrics_dir = Path(".", f"{args.sample}_metrics")
    metrics_dir.mkdir(parents=True, exist_ok=True)
    sample_stats.to_csv(metrics_dir / f"{args.sample}_sample_stats.csv", index=False)


if __name__ == "__main__":
    main()
