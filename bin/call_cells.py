#!/usr/bin/env python
"""
Filter barcodes to call cells and generate filtered expression matrix.
"""

import argparse
from dataclasses import dataclass
from pathlib import Path
import tempfile
import gzip
import shutil
from filecmp import cmp

import numpy as np
import pandas as pd
import duckdb
from scipy.stats import median_abs_deviation

from cell_finder import rescue_cells
from scale_utils import io
from scale_utils.stats import OutlierOptions


@dataclass(frozen=True)
class CellCallingOptions:
    """
    This class holds options for cell calling.

    Attributes:
        fixedCells: If true the top fixedCells ranked by UTC are called as cells.
        expectedCells: Expected number of cells to call.
        topCellPercent:
            Percentile of cells sorted descending by UMI counts
            to help calculate the unique transcript counts threshold.
        minCellRatio:
            Minimum ratio of unique transcript counts to top cell UMI counts
            to help calculate the unique transcript counts threshold.
        minUTC:
            The minimum number of unique transcript counts a barcode
            must be associated with to be considered a potential cell.
        UTC:
            The minimum number of unique transcript counts a barcode
            must be associated with to be considered a called cell.
        cellFinder: If true Cell Finder will be used for cell calling.
        FDR: False discovery rate to use for rescuing cells based on deviation from ambient profile.
        alpha:
            Overdispersion parameter for Dirichlet Multinomial distribution.
            If 0, estimated from ambient barcodes.
        medianFraction:
            When running cellFinder only cells with a UTC count above medianFraction times the median
            count of cells passing the UTC threshold are considered ambiguous and hence tested
    """

    fixedCells: bool
    expectedCells: int
    topCellPercent: int
    minCellRatio: float
    minUTC: int
    UTC: int
    cellFinder: bool
    FDR: float
    alpha: float
    medianFraction: float


@dataclass(frozen=True)
class CellfinderMetrics:
    """
    Summary metrics from the cellfinder process
    """

    UTCThreshold: int
    maxAmbiguous: int
    minAmbiguous: int
    numAmbiguous: int
    alpha: float


def classify_barcodes(
    total_counts: np.ndarray, options: CellCallingOptions
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Classify cell barcodes into ambient, cell and (optionally) ambiguous,
    based on total transcript counts per barcode

    Args:
        total_counts: Total unique transcript counts for all barcodes

    Returns:
        A boolean array indicating if ambient barcode
        A boolean array indicating if ambiguous barcode, if CellFinder is enabled
        A boolean array indicating if barcode above threshold or in top --expectedCells barcodes
    """
    cell_barcodes = np.zeros(len(total_counts), dtype=bool)  # Initialize all barcodes as non-cells
    ambiguous_barcodes = np.zeros(len(total_counts), dtype=bool)  # Only used if CellFinder is enabled
    ambient_barcodes = total_counts < options.minUTC
    threshold = options.UTC or calculate_UTC_threshold(total_counts, ambient_barcodes, options=options)
    if options.cellFinder:
        cell_barcodes = total_counts >= threshold
        ambiguous_barcodes = ~ambient_barcodes & ~cell_barcodes
        ambiguous_barcodes = ambiguous_barcodes & (
            total_counts >= np.median(total_counts[cell_barcodes]) * options.medianFraction
        )
    elif options.fixedCells:
        if options.expectedCells < 0:
            raise ValueError("--expectedCells must be a positive integer.")
        # mark the top --expectedCells barcodes as cells based on total counts
        cell_barcodes[np.argsort(total_counts)[-options.expectedCells :]] = True
        ambient_barcodes = ~cell_barcodes
    # 'Top-cells" UTC Thresholding
    else:
        cell_barcodes = total_counts >= threshold
        ambient_barcodes = ~cell_barcodes
    return ambient_barcodes, ambiguous_barcodes, cell_barcodes


def calculate_UTC_threshold(total_counts: np.ndarray, ambient_barcodes: np.ndarray, options: CellCallingOptions) -> int:
    """
    Calculate the unique transcripts counts threshold above which all barcodes are called as cells

    Args:
        total_counts: Total unique transcript counts for all barcodes
        ambient_barcodes: Boolean array if barcode is ambient
        options: Options for cell calling

    Returns:
        The unique transcript counts threshold
    """

    if options.topCellPercent < 1 or options.topCellPercent > 99:
        raise ValueError("--topCellPercent must fall in the range 1 <= x <= 99.")
    if options.minCellRatio < 1:
        raise ValueError("--minCellRatio must be >= 1.")
    if options.expectedCells > 0 and options.expectedCells < len(total_counts):
        # when --expectedCells is set, use the top --expectedCells barcodes as preliminary cells
        prelim_utcs = total_counts[np.argsort(total_counts)[::-1]][: options.expectedCells]
    else:
        prelim_utcs = total_counts[~ambient_barcodes]

    # Calculate unique transcripts counts threshold
    threshold = 0
    if len(prelim_utcs) != 0:
        # Calculate the UTC threshold as the UMI count
        # associated with the --topCellPercent of the preliminary cells divided by --minCellRatio.
        threshold = np.round(np.percentile(prelim_utcs, options.topCellPercent) / options.minCellRatio)

    # Don't return threshold below minUTC
    return max(threshold, options.minUTC)


def call_cell_barcodes(
    con: duckdb.DuckDBPyConnection,
    bcs_df: pd.DataFrame,
    options: CellCallingOptions,
    subset: np.ndarray,
    counts_col: str = "counts",
) -> pd.DataFrame:
    """
    Call cells in this sample using thresholding options, including optionally
    rescuing ambiguous barcodes as cells with the CellFinder algorithm

    Args:
        con: duckdb connection with existing mtx and mtx_metadata tables
        bcs_df: dataframe containing cell_id and umi count along with pass and flags columns
        options: Options for cell calling
        subset:
            bool array to subset bcs_df or mtx.
            Used with barnyard data to subset bcs_df/mtx to mouse/human only.
        counts_col: For barnyard data, the column in bcs_df to use for cell calling.

    Returns:
        CellCalling metrics DataFrame

    Updates the pass and flags columns in bcs_df
    """
    total_counts = bcs_df.loc[subset, counts_col].to_numpy()
    ambient_bcs, ambiguous_bcs, cell_bcs = classify_barcodes(total_counts, options=options)
    flags = pd.Series("", index=bcs_df.index[subset])
    passed = pd.Series(cell_bcs, index=bcs_df.index[subset], copy=True)
    cellfinder_metrics = [
        ("CellCalling", "ambiguous_bcs", ambiguous_bcs.sum()),
        ("CellCalling", "ambient_bcs", ambient_bcs.sum()),
    ]
    if np.any(ambiguous_bcs):  # cellfinder algorithm to rescue ambiguous barcodes
        # get indices of barcodes, is subset in barnyard data
        all_indices = np.argwhere(subset).ravel()
        gene_sums = io.sum_counts_by(con, all_indices, "gene")
        # calculate gene sums for ambient barcodes to avoid loading entire matrix
        ambient_gene_sums = io.sum_counts_by(con, all_indices[ambient_bcs], "gene")
        ambient_bc_sums = io.sum_counts_by(con, all_indices[ambient_bcs], "barcode")
        # Disregard genes with zero counts
        nzgenes = gene_sums > 0
        # Match shape of mtx
        ambient_gene_sums = ambient_gene_sums[nzgenes]
        gene_sums = gene_sums[nzgenes]

        # Call barcodes with FDR <= --FDR as cells
        FDRs, stats = rescue_cells(
            con,
            nzgenes,
            all_indices,
            ambient_bc_sums,
            ambient_gene_sums,
            gene_sums,
            ambient_bcs,
            ambiguous_bcs,
            options.alpha,
            len(cell_bcs),
        )
        passed.loc[ambiguous_bcs] = FDRs <= options.FDR
        flags.loc[ambiguous_bcs] = flags.iloc[ambiguous_bcs].where(FDRs > options.FDR, other="cellFinder")
        cellfinder_metrics.extend(
            [
                ("CellCalling", "min_ambiguous", total_counts[ambiguous_bcs].min()),
                ("CellCalling", "max_ambiguous", total_counts[ambiguous_bcs].max()),
                ("CellCalling", "alpha", stats["alpha"]),
                ("CellCalling", "min_pval", stats["min_pval"]),
            ]
        )

    bcs_df.loc[subset, "pass"] = passed
    bcs_df.loc[subset, "flags"] = flags
    return pd.DataFrame.from_records(cellfinder_metrics, columns=["Category", "Metric", "Value"])


def filter_cells(bcs_df: pd.DataFrame, options: OutlierOptions) -> pd.DataFrame:
    """
    Filter cells based on median absolute deviations of various QC metrics
    'flags' and 'pass' in all_cells are updated in place.

    Args:
        bcs_df: Metrics for all cell-barcodes
        options: Options for outlier filtering

    Returns:
        Metrics about MAD cell filtering
    """
    cells = bcs_df[bcs_df["pass"]]  # Only run outlier filtering on cells (not background)
    flags = pd.Series("", index=cells.index)
    stats = []
    if options.reads_mads:  # Flag cells with low or high total read counts
        lreads = np.log(cells["totalReads"])
        min_lreads = np.median(lreads) - options.reads_mads * median_abs_deviation(lreads)
        max_lreads = np.median(lreads) + options.reads_mads * median_abs_deviation(lreads)
        flags[lreads < min_lreads] += ";low_reads"
        flags[lreads > max_lreads] += ";high_reads"
        stats.append(("MAD", "minimum_total_reads", np.round(np.exp(min_lreads))))
        stats.append(("MAD", "maximum_total_reads", np.round(np.exp(max_lreads))))
    if options.passing_mads:  # Flag cells with low fraction passing reads (counted to a gene)
        preads = cells["countedReads"] / cells["totalReads"]
        min_preads = max(0, np.median(preads) - options.passing_mads * median_abs_deviation(preads))
        flags[preads < min_preads] += ";low_passing_reads"
        stats.append(("MAD", "minimum_passing_reads", min_preads))
    if options.mito_mads:  # Flag cells with high fraction mito reads
        mito = cells["mitoProp"]
        max_mito = max(options.mito_min_thres, np.median(mito) + options.mito_mads * median_abs_deviation(mito))
        flags[mito > max_mito] += ";high_mito"
        stats.append(("MAD", "maximum_mito", max_mito))
    if options.filter_outliers:
        bcs_df.loc[flags.index[flags != ""], "pass"] = False
    bcs_df.loc[cells.index, "flags"] += flags
    bcs_df["flags"] = bcs_df["flags"].str.strip(";")
    return pd.DataFrame.from_records(stats, columns=["Category", "Metric", "Value"])


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


def validate_bcs_df_order(bcs_df: pd.DataFrame, star_barcodes_path: Path) -> None:
    """Check that the order of bcs_df matches barcodes.tsv from merged star output

    We rely on the order of the bcs_df DataFrame to match the order of the barcodes.tsv file
    It is possible for the dataframe to be inadvertently sorted or shuffled, which would cause
    data integrity issues when writing the filtered matrix.

    Args:
        bcs_df: DataFrame with cell_id as index and a column indicating if each barcode passed cell calling
        star_barcodes_path: Path to the barcodes.tsv file from the STARsolo output

    Raises:
        AssertionError: If the order of bcs_df does not match barcodes.tsv
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        bcs_df.index.to_frame().to_csv(Path(tmpdir) / "raw_barcodes.tsv", header=False, index=False)
        # create decompressed barcodes.tsv from star output
        with gzip.open(star_barcodes_path, "rb") as f_in, open(Path(tmpdir) / "star_barcodes.tsv", "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        if not cmp(Path(tmpdir) / "star_barcodes.tsv", Path(tmpdir) / "raw_barcodes.tsv", shallow=False):
            raise AssertionError("Order of bcs_df DataFrame does not match barcodes.tsv")


def generate_filtered_matrix(
    sample_specific_file_paths: dict[str, Path],
    bcs_df: pd.DataFrame,
    sample: str,
    round_counts: bool,
    con: duckdb.DuckDBPyConnection,
) -> None:
    """
    Write filtered cell-by-gene expression matrix using passing cell-barcode indices and raw mtx file

    Args:
        sample_specific_file_paths:
            Dictionary where each key is a custom file identifier and
            each value is a path to the identified file in the STARsolo output directory
        bcs_df: DataFrame with cell_id as index a column indicating if each barcode passed cell calling
        sample: Unique string to identify this sample
        round_counts: Whether to round counts to integers
        con: duckdb connection with existing mtx and mtx_metadata tables

    Returns:
        None. A filtered matrix is written to the {sample}_filtered_star_output directory
    """
    validate_bcs_df_order(bcs_df, sample_specific_file_paths["barcodes"])
    # Create /filtered directory under current directory
    filtered_path = f"{sample}_filtered_star_output"
    Path(".", filtered_path).mkdir(exist_ok=True)
    # Features.tsv is the same; we move it to that /filtered directory
    shutil.copyfile(sample_specific_file_paths["features"], f"{filtered_path}/features.tsv.gz")
    bcs_df.index[bcs_df["pass"]].to_frame().to_csv(f"{filtered_path}/barcodes.tsv.gz", header=False, index=False)
    # drop cell_id and generate numerical index
    bcs_df = bcs_df.reset_index(drop=True)
    # filter to passing cells
    bcs_df = bcs_df[bcs_df["pass"]]
    lookup = {
        # create lookup table from barcode index in raw mtx from STAR to new index in filtered mtx
        str(barcode_idx + 1): new_barcode_idx
        for new_barcode_idx, barcode_idx in enumerate(bcs_df.index, start=1)
    }

    with (
        gzip.open(sample_specific_file_paths["mtx"], "rt") as file,
        open(f"{filtered_path}/matrix.mtx", "w+") as outfile,
    ):
        header = next(file)
        if round_counts:
            # change:   %%MatrixMarket matrix coordinate real general
            # to:       %%MatrixMarket matrix coordinate integer general
            # in the filtered mtx file. https://math.nist.gov/MatrixMarket/formats.html
            header = header.replace("real", "integer")
        # write header lines
        outfile.write(header)
        outfile.write(next(file))
        nrows = next(file).split()[0]
        pos = outfile.tell()
        # write l-padded placeholder for non-zero entries
        outfile.write(f"{nrows} {len(bcs_df)} {'0'.zfill(13)}\n")
        nnz = 0
        for line in file:
            gene, barcode, count = line.split()
            if barcode in lookup:
                count = round(float(count)) if round_counts else float(count)
                if count != 0:
                    nnz += 1
                    outfile.write(f"{gene} {lookup[barcode]} {count}\n")
        outfile.seek(pos)
        # write calculated nnz
        outfile.write(f"{nrows} {len(bcs_df)} {str(nnz).zfill(13)}\n")


def label_barnyard(bcs_df: pd.DataFrame) -> None:
    """
    In barnyard datasets, based on UMI count and minor fraction, label cells in 'species' column as either:
        (1) 'Human'
        (2) 'Mouse'
        (3) 'Mixed' (above umi threshold for both species)
        (4) 'Ambiguous' (minor fraction below background)
        (5) 'None' (all non-passing barcodes)

    Args:
        bcs_df: dataframe containing selected barcode metadata
    """
    # Label barcodes that failed cell thresholding as 'None'
    bcs_df.loc[~bcs_df["pass"], "species"] = "None"

    # Calculate the minimum number of human UMIs. Used to determine if cell is "mixed".
    if (bcs_df.species == "Human").any():
        min_human = np.percentile(bcs_df.human_counts[bcs_df.species == "Human"], 10)
    else:
        min_human = bcs_df[bcs_df["pass"]].counts.min()

    # Calculate the minimum number of mouse UMIs. Used to determine if cell is "mixed".
    if (bcs_df.species == "Mouse").any():
        min_mouse = np.percentile(bcs_df.mouse_counts[bcs_df.species == "Mouse"], 10)
    else:
        min_mouse = bcs_df[bcs_df["pass"]].counts.min()

    # Label Mixed Cells
    bcs_df.loc[(bcs_df.human_counts >= min_human) & (bcs_df.mouse_counts >= min_mouse), "species"] = "Mixed"

    # Calculate Background
    human_bg_med, human_bg_std = get_background(bcs_df[(bcs_df.species == "Mouse")].minor_frac)
    mouse_bg_med, mouse_bg_std = get_background(bcs_df[(bcs_df.species == "Human")].minor_frac)

    # Labeling mixed cells with low minor species fraction as ambiguous
    bcs_df.loc[
        (bcs_df.species == "Mixed")
        & (bcs_df.human_counts > bcs_df.mouse_counts)
        & (bcs_df.minor_frac < mouse_bg_med + 2 * mouse_bg_std),
        "species",
    ] = "Ambiguous"
    bcs_df.loc[
        (bcs_df.species == "Mixed")
        & (bcs_df.mouse_counts > bcs_df.human_counts)
        & (bcs_df.minor_frac < human_bg_med + 2 * human_bg_std),
        "species",
    ] = "Ambiguous"

    # Labelling Human or Mouse as Ambiguous
    bcs_df.loc[
        (bcs_df.species == "Human") & (bcs_df.minor_frac >= max(0.1, mouse_bg_med + 3 * mouse_bg_std)), "species"
    ] = "Ambiguous"
    bcs_df.loc[
        (bcs_df.species == "Mouse") & (bcs_df.minor_frac >= max(0.1, human_bg_med + 3 * human_bg_std)), "species"
    ] = "Ambiguous"


def get_background(minor_fracs: pd.Series) -> tuple[float, float]:
    """
    Calculate median and standard deviation for barnyard cells

    Args:
    minor_fracs: List with fraction indicating percentage of human or mouse samples in a cell

    Returns:
        Median and standard deviation
    """
    if minor_fracs.size == 0:
        return np.nan, np.nan
    median = minor_fracs.median()
    std = (np.quantile(minor_fracs, 0.75) - median) / 0.675
    return median, std


def filter_beads(con, fail_ambient: bool, bead_scores: Path, min_kl_score: float) -> None:
    filter_query = ", pass = false" if fail_ambient else ""
    con.sql(
        f"""
    UPDATE all_barcodes
    SET flags = COALESCE(flags || ';', '') || 'ambient_bead'{filter_query}
    WHERE bead_bc IN (SELECT bead_bc FROM '{bead_scores}' WHERE kl_norm < {min_kl_score}) AND pass;
    """
    )


def load_barcodes_table(
    barcode_metrics: Path,
    con: duckdb.DuckDBPyConnection,
    is_barnyard: bool,
) -> None:
    """
    Load a allBarcodes.parquet file into a duckdb table

    Args:
        barcode_metrics: Path to the parquet file
        con: duckdb connection
        is_barnyard: Whether the data is from a barnyard assay
        is_quantum: Whether the data is from a QuantumScale assay
    """
    by_cols = (
        """
    'None' AS species,
    round(LEAST(human_counts, mouse_counts) / (human_counts + mouse_counts), 3) AS minor_frac,
    """
        if is_barnyard
        else ""
    )
    con.sql(
        f"""
    CREATE OR REPLACE TABLE all_barcodes
    AS SELECT * EXCLUDE (sample),
    sample,
    CAST(null AS VARCHAR) AS flags,
    {by_cols}
    false AS pass,
    FROM '{barcode_metrics}';
    """
    )


def update_barcodes_table(
    con: duckdb.DuckDBPyConnection,
    bcs_df: pd.DataFrame,
    is_barnyard: bool,
    is_quantum: bool,
) -> None:
    """
    Update the all_barcodes table with results from cell calling

    Args:
        con: duckdb connection
        bcs_df: dataframe containing cell_id to match with table and results of cell calling
        is_barnyard: Whether the data is from a barnyard assay
        is_quantum: Whether the data is from a QuantumScale assay
    """
    # flake8 complains about unused variables, but they are used in the duckdb SQL query
    passing_cells = bcs_df[bcs_df["pass"]]  # noqa: F841
    con.sql(
        """
    UPDATE all_barcodes
    SET pass = true
    FROM passing_cells
    WHERE all_barcodes.cell_id = passing_cells.cell_id;
    """
    )

    flagged_cells = bcs_df[bcs_df["flags"] != ""]  # noqa: F841
    con.sql(
        """
    UPDATE all_barcodes
    SET flags = flagged_cells.flags
    FROM flagged_cells
    WHERE all_barcodes.cell_id = flagged_cells.cell_id;
    """
    )

    if is_barnyard:
        by_cells = bcs_df[bcs_df["species"] != "None"]  # noqa: F841
        con.sql(
            """
        UPDATE all_barcodes
        SET species = by_cells.species,
        FROM by_cells
        WHERE all_barcodes.cell_id = by_cells.cell_id;
        """
        )


def main():
    parser = argparse.ArgumentParser()

    # Required argument for specifying the path to the sampleMetrics parquet file
    parser.add_argument("--sampleMetrics", type=Path, required=True, help="Path to the sampleMetrics parquet file.")

    # Required and optional arguments for specifying the STARsolo outputs for this sample
    parser.add_argument(
        "--STARsolo_out", type=Path, required=True, help="Path to the STARsolo outputs for this sample."
    )
    parser.add_argument(
        "--feature_type", type=str, required=False, default="GeneFull_Ex50pAS", help="STARsolo feature type used."
    )
    parser.add_argument(
        "--matrix_type",
        type=str,
        required=False,
        default="UniqueAndMult-PropUnique.mtx",
        help="STARsolo matrix type used.",
    )

    # Optional argument to specify the name of the sample for which cells are being called
    parser.add_argument(
        "--sample", type=str, required=False, default="example", help="Unique string to identify this sample."
    )

    # Optional argument to set hard thresholds for cell calling
    parser.add_argument(
        "--fixedCells",
        required=False,
        action="store_true",
        default=False,
        help="Fixed number of barcodes to call as cells.",
    )

    # Optional arguments for use of CellFinder algorithm
    parser.add_argument(
        "--expectedCells",
        type=int,
        required=False,
        default=0,
        help="Expected number of cells to call (if specified in samples.csv; see algorithm description).",
    )
    parser.add_argument(
        "--topCellPercent",
        type=int,
        required=False,
        default=99,
        help="Cell thresholding parameter (see algorithm description).",
    )
    parser.add_argument(
        "--minCellRatio",
        type=float,
        required=False,
        default=10,
        help="Cell thresholding parameter (see algorithm description).",
    )
    parser.add_argument(
        "--minUTC",
        type=int,
        required=False,
        default=100,
        help="The minimum number of unique transcript counts a barcode \
            must be associated with to be considered a potential cell.",
    )
    parser.add_argument(
        "--UTC",
        type=int,
        default=0,
        help="Set a fixed unique transcript count threshold for a barcode to be called cell.",
    )
    # CellFinder parameters
    parser.add_argument(
        "--cellFinder",
        action="store_true",
        help="Use CellFinder to call cells from among barcodes with counts between --minUTC and the UTC threshold",
    )
    parser.add_argument(
        "--FDR",
        type=float,
        default=0.001,
        help="False discovery rate at which to call cells (see algorithm description).",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0,
        help="Set a fixed overdispersion (Dirichlet alpha) parameter for CellFinder. \
            If 0, estimated from ambient barcodes",
    )
    parser.add_argument(
        "--medianFraction",
        type=float,
        default=0,
        help="Only cells with at least this fraction of counts relative to the median \
            of cells passing the simple UTC threshold are tested by cellFinder",
    )
    # MAD outlier filtering
    parser.add_argument(
        "--filter_outliers",
        required=False,
        action="store_true",
        help="Number of median absolute deviations in gene count/UMI count/mitochondrial \
            read percentage above/below which a cell will be flagged as an outlier.",
    )
    parser.add_argument("--madsReads", type=float, default=np.nan, help="MAD threshold for total reads (+/- X MADs)")
    parser.add_argument(
        "--madsPassingReads", type=float, default=np.nan, help="MAD threshold for passing reads fraction (- x MADs)"
    )
    parser.add_argument("--madsMito", type=float, default=np.nan, help="MAD threshold for mito. reads (+ x MADs)")
    # Optional argument to specify whether to generate statistics for internal report
    parser.add_argument("--internalReport", action="store_true", default=False)
    parser.add_argument("--isBarnyard", action="store_true", default=False)
    parser.add_argument("--roundCounts", action="store_true", default=False)

    # Argument that determines whether or not beads are filtered from matrix due to high cell count.
    parser.add_argument("--filterBeads", action="store_true", default=False)
    parser.add_argument(
        # Minimum KL divergence score to consider a bead as non-ambient
        # Normalized range is roughly 0-1, with 0 perfectly matching the library RT count distribution
        "--minDivergence",
        default=0.05,
        type=float,
        help="Minimum normalized KL divergence score to consider a bead as non-ambient",
    )
    parser.add_argument("--beadScores", type=Path, required=False, help="Path to results of BeadFiltering process.")
    # Argument which indicates if data is from QuantumScale assay.
    parser.add_argument("--isQuantum", action="store_true", default=False)
    parser.add_argument("--threads", type=int, required=False, default=1, help="Number of threads for duckdb")
    parser.add_argument("--memory", type=str, required=False, default="8 GB", help="Memory allocated to task")

    args = parser.parse_args()

    sample_specific_file_paths = io.resolve_sample_specific_file_paths(
        args.STARsolo_out, args.feature_type, args.matrix_type
    )

    call_cells_options = CellCallingOptions(
        fixedCells=args.fixedCells,
        expectedCells=args.expectedCells,
        topCellPercent=args.topCellPercent,
        minCellRatio=args.minCellRatio,
        minUTC=args.minUTC,
        UTC=args.UTC,
        cellFinder=args.cellFinder,
        FDR=args.FDR,
        alpha=args.alpha,
        medianFraction=args.medianFraction,
    )

    outlier_options = OutlierOptions(
        filter_outliers=args.filter_outliers,
        reads_mads=args.madsReads,
        passing_mads=args.madsPassingReads,
        mito_mads=args.madsMito,
    )

    mem_limit, mem_unit = args.memory.split()
    # allocate duckdb memory based on number of threads
    mem_limit = f"{float(mem_limit) / (args.threads + 1):.1f}{mem_unit}"
    con = duckdb.connect("tables.db")
    con.sql(
        f"""
    SET threads TO {args.threads};
    SET memory_limit TO '{mem_limit}';
    """
    )
    load_barcodes_table(args.sampleMetrics, con, args.isBarnyard)
    by_cols = "human_counts, mouse_counts, species, minor_frac" if args.isBarnyard else ""
    bcs_df = con.sql(
        f"""
        SELECT cell_id, counts, pass, flags, totalReads, countedReads, mitoProp, {by_cols} FROM all_barcodes
    """
    ).df()
    if args.cellFinder:
        io.load_mtx_table(sample_specific_file_paths["mtx"], con)

    metrics: list[pd.DataFrame] = []  # Collecting metrics from different steps
    if args.isBarnyard:
        bcs_df.loc[bcs_df["human_counts"] > bcs_df["mouse_counts"], "species"] = "Human"
        bcs_df.loc[bcs_df["mouse_counts"] > bcs_df["human_counts"], "species"] = "Mouse"
        mouse_cells = bcs_df["species"] != "Human"
        human_cells = bcs_df["species"] != "Mouse"
        call_cell_barcodes(con, bcs_df, call_cells_options, mouse_cells, "mouse_counts")
        call_cell_barcodes(con, bcs_df, call_cells_options, human_cells, "human_counts")
        label_barnyard(bcs_df)
    else:
        cellfinder_metrics = call_cell_barcodes(con, bcs_df, call_cells_options, np.ones(len(bcs_df), dtype=bool))
        metrics.append(cellfinder_metrics)

    mad_stats = filter_cells(bcs_df, options=outlier_options)
    metrics.append(mad_stats)
    update_barcodes_table(con, bcs_df, is_barnyard=args.isBarnyard, is_quantum=args.isQuantum)

    if args.isQuantum:
        filter_beads(con, args.filterBeads, args.beadScores, args.minDivergence)

    # Finished calling cells, write allCells.csv for this sample
    metrics_dir = Path(".", f"{args.sample}_metrics")
    metrics_dir.mkdir(parents=True, exist_ok=True)
    cell_metrics = metrics_dir / f"{args.sample}_allCells.csv"
    barcode_metrics = metrics_dir / f"{args.sample}_allBarcodes.parquet"
    con.sql(
        f"""
    COPY (
        SELECT * EXCLUDE (pass)
        FROM all_barcodes
        WHERE pass
    ) TO '{cell_metrics}';
    """
    )
    # and write entire table in parquet format
    con.sql(
        f"""
    COPY (
        SELECT *
        FROM all_barcodes
    ) TO '{barcode_metrics}';
    """
    )

    # Write cellcalling metrics
    pd.concat(metrics).to_csv(metrics_dir / f"{args.sample}_cellcalling_stats.csv", index=False)

    # Write filtered matrix for this sample
    # only need passing column for creating filtered matrix
    bcs_df = con.sql("SELECT cell_id, pass FROM all_barcodes").df().set_index("cell_id")
    generate_filtered_matrix(sample_specific_file_paths, bcs_df, args.sample, args.roundCounts, con)
    con.close()


if __name__ == "__main__":
    main()
