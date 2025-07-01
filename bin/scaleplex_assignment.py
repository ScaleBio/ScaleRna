#!/usr/bin/env python
"""Generate Filtered UMI matrices and library level metrics
for only those cell barcodes that were called a cell in RNA analysis"""
import argparse
from pathlib import Path
from dataclasses import dataclass
import csv
import gzip
import numpy as np
import pandas as pd
import polars as pl
import polars.selectors as cs
from scipy.sparse import csr_array, csc_array, hstack
from scipy.io import mmwrite, mmread
from scipy.stats import poisson
from scipy.stats import false_discovery_control
from scale_utils.lib_json_parser import LibJsonParser


@dataclass(frozen=True)
class AssignmentCodes:
    """ScalePlex assignment strings"""

    max_fail: str = "Max_Fail"
    indeterminate: str = "Indeterminate"
    enrich_fail: str = "Enrich_Fail"
    unexpected: str = "Unexpected"

    @property
    def errors(self) -> list:
        return [self.max_fail, self.indeterminate, self.enrich_fail, self.unexpected]


codes = AssignmentCodes()


def preprocessJson(libStructJson: Path):
    """Get file with scaleplex sequences and mapping of scaleplex combos to fixation plate well"""
    lib_json_obj = LibJsonParser(libStructJson)
    for bc in lib_json_obj.json_contents["barcodes"]:
        if bc.get("name") in ["scaleplex"]:
            guide_file = bc.get("sequences")
            hash_combos = bc.get("mapping")
    aliases = [
        bc["alias"]
        for bc in lib_json_obj.json_contents["barcodes"]
        if bc.get("type", None) not in ["library_index", "umi", "target"]
    ]
    bead_bc_cols = ["bead1", "bead2", "bead3"]
    if all([bc in aliases for bc in bead_bc_cols]):
        bead_bc_cols_to_drop = [
            bc.get("alias")
            for bc in lib_json_obj.json_contents["barcodes"]
            if bc.get("alias") in bead_bc_cols and not bc.get("plate")
        ]
        aliases = [bc for bc in aliases if bc not in bead_bc_cols_to_drop]

    return hash_combos, guide_file, aliases


def create_index_dict(index_values):
    index_dict = {}
    for index, value in enumerate(index_values):
        index_dict[value] = index
    return index_dict


def filter_sparse_matrix(raw_mtx, filtered_indices, cell_stats):
    """
    Rebuild sparse matrix to contain only the cell barcodes that were called a cell in the allCells.csv from RNA.
    Cell barcodes can have passed in RNA that are not present in the raw matrix generated from the enriched library,
    and so we cannot do a simple subsetting
    Args:
        raw_mtx: raw umi or read csr sparse matrix in guides x cells orientation
        filtered_indices: passing barcodes from RNA allCells.csv
        cell_stats: cell metrics file to provide raw barcodes list
    Returns:
        filtered_mtx:
            counts matrix in guides x cells orientation where cells are only barcodes that were called a cell
            in ScaleRNA, with empty columns for those that were not detected in enriched library
    """
    # get indices in raw matrix of barcodes that were called a cell in RNA
    indices = (
        filtered_indices.select(["cell_id"])
        .with_row_index("filtered_index")
        .join(
            cell_stats.select(["Cell_Barcode"]).with_row_index(), left_on="cell_id", right_on="Cell_Barcode", how="left"
        )
    )
    idx = indices.filter(pl.col("index").is_not_null())["index"].to_numpy()
    filtered_mtx = raw_mtx[:, idx]

    # Add sparse columns for barcodes that were not detected in the enriched library
    missing_idx = indices.filter(pl.col("index").is_null())["filtered_index"].to_numpy()
    for idx in missing_idx:
        filtered_mtx = hstack([filtered_mtx[:, :idx], csr_array((filtered_mtx.shape[0], 1)), filtered_mtx[:, idx:]])

    return filtered_mtx


def guide_dictionary(guide_file_path: Path):
    """Read in whitelist for dataset features and return dictionary with keys as name of"""
    with open(guide_file_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        guides = [row[0] for row in reader]
        guide_mapping = {col: idx for idx, col in enumerate(guides)}
    return guide_mapping


def wells_for_range(range_str, all_wells):
    """Return a list of wells based on column value in samples.csv"""
    wells = []
    for item in range_str.split(";"):
        input = item.replace(" ", "")
        if "-" in input:
            start, end = input.split("-")
            wells.extend(all_wells[all_wells.index(start) : all_wells.index(end) + 1])
        else:
            wells.append(input)
    return wells


def hash_to_well(hash_combos: Path, expected_combos=None):
    """Return dict of hash combos to well"""
    with open(hash_combos, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        hash_combos_dict = {row[0]: row[1] for row in reader}
    if expected_combos:
        expected_fixation_wells = wells_for_range(expected_combos, list(hash_combos_dict.values()))
        hash_combos_dict = {k: v for k, v in hash_combos_dict.items() if v in expected_fixation_wells}
    return hash_combos_dict


def write_keys_to_tsv(dictionary, tsv_file):
    with gzip.open(tsv_file, "wt", newline="") as file:
        writer = csv.writer(file, delimiter="\t")
        for key in dictionary.keys():
            writer.writerow([key])


def compute_background_series(mtx, min_cell_thresh, min_background_val, top_n=2):
    """Builds background estimation per hash,
    where background is considered to be any non-zero and non-top n values per cell per hash
    Args:
        sparse_matrix: filtered hashes x cells UMI counts matrix
        min_cell_thresh: fraction of cells in which a background count needs to be observed per oligo to be considered in estimation
        min_background_val: value set in background estimation if scaleplex oligo fails min_cell_thresh check
        top_n: number of top values within each cell (column) to ignore when estimating background.
    """
    background = np.zeros(mtx.shape[0])
    # detect 2nd highest value in column
    nth_val = np.zeros(mtx.shape[1])
    mtx = csc_array(mtx)
    for col_idx in range(mtx.shape[1]):
        arr = mtx[:, [col_idx]].toarray().ravel()  # convert sparse column to 1-D array
        # np.partition splits array into two parts, the top_n highest values and the rest
        nth_val[col_idx] = np.partition(arr, -top_n)[-top_n]  # get nth highest value

    mtx = csr_array(mtx)
    for row_idx in range(mtx.shape[0]):
        row_arr = mtx[[row_idx], :].toarray().ravel().astype(float)
        row_arr[(row_arr >= nth_val) | (row_arr == 0)] = np.nan
        # compute median of counts for this hash excluding cells where it was zero or a top n hash
        min_cells = mtx.shape[1] * min_cell_thresh # calculate min n cells a bg count needs to occur in
        bg_counts = sum(~np.isnan(row_arr)) # number of cells a bg count takes place in
        background[row_idx] = min_background_val if bg_counts < min_cells else np.nanmedian(row_arr) # set bg to zero for oligo if too few bg counts, otherwise median
    return background


def corrected_p_values_poisson(sparse_matrix, background, features_dict, threshold):
    observed_counts = sparse_matrix.toarray()
    expected_counts = 3 * background

    expected = np.tile(expected_counts, (sparse_matrix.shape[1], 1)).T  # Repeat expected counts for each sample
    p_value = 1 - poisson.cdf(observed_counts - 1, expected)
    results = pd.DataFrame(p_value)

    results_T = results.T
    results_corrected = results_T.apply(lambda col: false_discovery_control(col, method="bh"), axis=0)
    results_corrected.columns = features_dict.keys()

    # Did a hash pass the background test

    passing_hashes = results_corrected.apply(
        lambda row: (
            ";".join([col for col in results_corrected.columns[row < threshold]]) if any(row < threshold) else "No_Pass"
        ),
        axis=1,
    )
    passing_hashes_df = pd.DataFrame(passing_hashes, columns=["passing_scaleplex"])
    return passing_hashes_df


def assignment_criteria_bg(row, hash_well_dict, toptwo_frac):
    # Failed Background Test
    assigned_hash = ""
    if row["passing_scaleplex"] == "No_Pass":
        assigned_hash = codes.indeterminate
    # Passed Enrich Threshold and had at least one hash passing Background
    elif row["topTwo"] > toptwo_frac:
        max_cols = row["topTwo_scaleplex"]
        # Specific Hash Passed Background and was also Top hash in that cell

        strings_column1 = row["topTwo_scaleplex"].split(";")
        strings_column2 = row["passing_scaleplex"].split(";")
        if all(string in strings_column2 for string in strings_column1):
            assigned_hash = hash_well_dict[max_cols] if max_cols in hash_well_dict else codes.unexpected
        # Purity Threshold Passed, but Top hash was not the one that passed background
        else:
            assigned_hash = codes.max_fail
    # Had a hash that passed background, but cell did not pass topTwo threshold
    else:
        assigned_hash = codes.enrich_fail
    return assigned_hash


def assignment_criteria_fc(row, hash_well_dict, fc_threshold):
    assigned_hash = ""
    if row["second"] != 0 and (row["third"] == 0 or (row["second"] / row["third"]) > fc_threshold):
        max_cols = row["topTwo_scaleplex"]
        assigned_hash = hash_well_dict[max_cols] if max_cols in hash_well_dict else codes.unexpected
    else:
        assigned_hash = codes.enrich_fail
    return assigned_hash


def hash_assignment(
    sparse_matrix, cell_stats, features_dict, assignment_method, toptwo_frac, hash_well_dict, fc_threshold, min_cell_percent, min_background_count
):
    """Creates a list of guides (columns) that passed the minimum UMI threshold within each cell"""
    assignment_df = pd.DataFrame(index=cell_stats.index)
    if not assignment_df.index.empty:
        if assignment_method == "bg":
            background = compute_background_series(sparse_matrix, min_cell_percent, min_background_count, top_n=2)
            passing_hashes_df = corrected_p_values_poisson(sparse_matrix, background, features_dict, 0.05)
            assignment_df["passing_scaleplex"] = passing_hashes_df["passing_scaleplex"].values
            assignment_df["assigned_scaleplex"] = cell_stats.join(assignment_df).apply(
                assignment_criteria_bg, axis=1, hash_well_dict=hash_well_dict, toptwo_frac=toptwo_frac
            )
        # Iterate over rows
        if assignment_method == "fc":
            assignment_df["assigned_scaleplex"] = cell_stats.apply(
                assignment_criteria_fc, axis=1, hash_well_dict=hash_well_dict, fc_threshold=fc_threshold
            )
    else:
        # Add a column for passing_scaleplex to enable merge with non-empty allCells.csv
        if assignment_method == "bg":
            assignment_df["passing_scaleplex"] = "No_Pass"
        assignment_df["assigned_scaleplex"] = codes.indeterminate

    return pl.DataFrame(assignment_df.reset_index())


def main(
    umi_matrix: Path,
    all_cells: Path,
    cell_stats: Path,
    references: Path,
    lib_struct: Path,
    assignment_method: str,
    toptwo_frac: float,
    fc_threshold: float,
    min_cell_count_bg: float,
    min_bg_scaleplex_count: float,
    expected_combos: str,
    id: str,
    out_dir: Path,
):

    hash_combos, guide_file, aliases = preprocessJson(lib_struct)

    guide_matrix = mmread(umi_matrix)
    guide_matrix = (
        guide_matrix.tocsr()
    )  # cast to csr for formatting necessary for filtering indexing, otherwise casts as native coo matrix
    rna = pl.read_csv(all_cells)
    cell_stats_df = pl.read_parquet(cell_stats)

    # Create output directories
    matrix_dir = out_dir / f"{id}.filtered.matrix"
    raw_matrix_dir = out_dir / f"{id}.raw.matrix"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    raw_matrix_dir.mkdir(parents=True, exist_ok=True)

    # write raw barcodes
    with gzip.open(raw_matrix_dir / "barcodes.tsv.gz", "wb") as f:
        cell_stats_df.select(["Cell_Barcode"]).write_csv(f, separator="\t", include_header=False)

    guide_matrix_filtered = filter_sparse_matrix(guide_matrix, rna, cell_stats_df)
    with gzip.open(f"{matrix_dir}/matrix.mtx.gz", "wb") as f:
        mmwrite(f, guide_matrix_filtered)
    with gzip.open(f"{matrix_dir}/barcodes.tsv.gz", "wb") as f:
        rna.select(["cell_id"]).write_csv(f, separator="\t", include_header=False)

    # Create merged allCells with all rows from RNA and matching rows from ENRICH if present
    cell_stats_passing = rna.select(["cell_id"] + aliases).join(
        cell_stats_df, how="left", left_on="cell_id", right_on="Cell_Barcode", suffix="_ENRICH"
    )
    cell_stats_passing = cell_stats_passing.drop(cs.ends_with("_ENRICH"))

    guide_path = f"{references}/{guide_file}"
    guide_dict = guide_dictionary(guide_path)
    hash_well_dict = hash_to_well(references / hash_combos, expected_combos)

    assignment_df = hash_assignment(
        guide_matrix_filtered,
        cell_stats_passing.to_pandas().set_index("cell_id"),
        guide_dict,
        assignment_method,
        toptwo_frac,
        hash_well_dict,
        fc_threshold,
        min_cell_count_bg,
        min_bg_scaleplex_count,
    )  # hash assignemt of any hash that are = max UMI value for that cell
    cell_stats_passing = cell_stats_passing.join(assignment_df, how="left", on="cell_id")
    cell_stats_passing = cell_stats_passing.rename({"cell_id": "Cell_Barcode"})
    assignment_cols = assignment_df.columns
    assignment_cols.remove("cell_id")
    cell_stats_df = pl.concat(
        [
            cell_stats_passing.select(
                # match order of columns in cell_stats_df and keep new columns from cell_stats_passing
                cell_stats_df.columns
                + assignment_cols
            ).with_columns(pl.lit(True).alias("pass")),
            cell_stats_df.filter(
                # concatenate remaing rows from cell_stats_df that were not in cell_stats_passing
                ~pl.col("Cell_Barcode").is_in(cell_stats_passing["Cell_Barcode"])
            ).with_columns(
                # mark as not passing with null values for assignment columns
                [pl.lit(None).alias(i) for i in assignment_cols]
                + [pl.lit(False).alias("pass")]
            ),
        ],
        how="vertical",
        rechunk=False,
    )
    cell_stats_df = cell_stats_df.with_columns(pl.lit(id).alias("sample"))
    cell_stats_df = cell_stats_df.rename({"Cell_Barcode": "cell_id"})
    cell_stats_df.write_parquet(out_dir / f"{id}.ScalePlex.allBarcodes.parquet")

    assigned_cells = (
        rna.join(assignment_df, how="left", on="cell_id").to_pandas().set_index("cell_id")
    )  # add assignment for RNA passing cells
    # save allCells with scaleplex assignment
    assigned_cells.to_csv(out_dir / f"{id}_allCellsWithAssignment.csv")

    # write features.tsv
    write_keys_to_tsv(guide_dict, f"{matrix_dir}/features.tsv.gz")
    write_keys_to_tsv(guide_dict, f"{raw_matrix_dir}/features.tsv.gz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate Filtered Hash umi counts matrix and cellMetrics \
            that share the same cell barcodes as the passing RNA analysis per sample"
    )
    parser.add_argument(
        "--umi_matrix", metavar="MATRIX.mtx", type=Path, help="count_hash output of hash x cells raw UMI matrix"
    )
    parser.add_argument(
        "--all_cells",
        metavar="ALLCELLS.csv",
        type=Path,
        help="nf-rna per-cell summary statistics output to filter based on pass status",
    )
    parser.add_argument(
        "--cell_stats",
        metavar="CELLMETRICS.csv",
        type=Path,
        help="count_hash output of cell metadata for HASH cell guide UMI matrix",
    )
    parser.add_argument(
        "--references", metavar="REFERENCES", type=Path, help="Path to folder containing barcode whitelists"
    )
    parser.add_argument(
        "--lib_struct",
        metavar="LIB_STRUCT",
        type=Path,
        help="library structure json used for bcParser in enrich detection",
    )
    parser.add_argument(
        "--assignment_method",
        metavar="ASSIGNMENT_METHOD",
        type=str,
        default="bg",
        help="Sample background based, or fold change based assignment",
    )
    parser.add_argument(
        "--toptwo_frac",
        metavar="TOPTWO_FRAC",
        type=float,
        help="Top two unique scaleplex UMIs must be over this fraction to be assigned",
    )
    parser.add_argument(
        "--fc_threshold", metavar="FC_THRESHOLD", type=float, help="FC threshold for fold change based assignment"
    )
    parser.add_argument(
        "--min_cell_count_bg", 
        metavar="MIN_CELL_COUNT_BG", 
        type=float,
        help="Fraction of cells that an oligo needs background counts in for consideration"
    )
    parser.add_argument(
        "--min_bg_scaleplex_count", 
        metavar="MIN_BG_SCALEPLEX_COUNT", 
        type=float,
        help="Minimum value set in the background estimation for an oligo if it fails min_cell_count_bg check"
    )
    parser.add_argument(
        "--expected_combos",
        metavar="FIXATION_PLATE_WELLS",
        type=str,
        help="Range or wells from fixation plate to supervise hash assignment",
    )
    parser.add_argument("--id", required=True, help="Sample and lib name")
    parser.add_argument("--outDir", type=Path, help="Directory for outputs")
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(
        args.umi_matrix,
        args.all_cells,
        args.cell_stats,
        args.references,
        args.lib_struct,
        args.assignment_method,
        args.toptwo_frac,
        args.fc_threshold,
        args.min_cell_count_bg,
        args.min_bg_scaleplex_count,
        args.expected_combos,
        args.id,
        args.outDir,
    )
