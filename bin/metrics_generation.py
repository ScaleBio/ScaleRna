#!/usr/bin/env python
""" Generate Filtered UMI matrices and library level metrics for only those cell barcodes that were called a cell in RNA analysis """
import argparse
from pathlib import Path
from dataclasses import dataclass
import csv
import numpy as np
import pandas as pd
import json
from scipy.sparse import lil_matrix
from scipy.io import mmwrite, mmread
from scipy.stats import poisson
from scipy.stats import false_discovery_control


@dataclass
class Metrics:
    """Overall sample statistics"""
    meanReadsPerCell: int = 0 # number of reads per cell. Should be multiplied by correct for ~usable amount
    nUMIPerCell: int = 0 # median number of HASH UMI molecules per cell
    passingPercent: int = 0 # percent of cells that had at least thresh Guide UMI's detected
    maxFailPercent: int = 0 # percent of cells with maxFail assignment
    enrichFailPercent: int = 0 # percent of cells with enrichFail assignment
    indeterminatePercent: int = 0
    unexpectedPercent: int = 0
    saturation: int = 0 # median saturation of cells
    readsInCells: int = 0 # percent of hash reads that ended up in called cells

    def print(self, outFn: Path):
        with open(outFn, 'w') as out:
            print("ReadsPerCell", f"{self.meanReadsPerCell:}", sep=',', file=out)
            print("nUMIPerCell", f"{self.nUMIPerCell:}", sep=',', file=out)
            print("passingPercent", f"{self.passingPercent:.1%}", sep=',', file=out)
            print("maxFailPercent", f"{self.maxFailPercent:.1%}", sep=',', file=out)
            print("enrichFailPercent", f"{self.enrichFailPercent:.1%}", sep=',', file=out)
            print("indeterminatePercent", f"{self.indeterminatePercent:.1%}", sep=',', file=out)
            print("unexpectedPercent", f"{self.unexpectedPercent:.1%}", sep=',', file=out)
            print("saturation", f"{self.saturation:.1%}", sep=',', file=out)
            print("readsInCells", f"{self.readsInCells:.1%}", sep=',', file=out)


@dataclass(frozen=True)
class AssignmentCodes:
    """ScalePlex assignment strings"""
    max_fail: str = 'Max_Fail'
    indeterminate: str = 'Indeterminate'
    enrich_fail: str = 'Enrich_Fail'
    unexpected: str = 'Unexpected'
    
    @property
    def errors(self) -> list:
        return [self.max_fail, self.indeterminate, self.enrich_fail, self.unexpected]
codes = AssignmentCodes()


def preprocessJson(libStructJson: Path):
    """Get file with scaleplex sequences and mapping of scaleplex combos to fixation plate well"""
    libStruct = json.load(open(libStructJson))
    for bc in libStruct['barcodes']:
        if bc.get('name') in ['scaleplex']:
            guide_file = bc.get('sequences')
            hash_combos = bc.get('mapping')

    return hash_combos, guide_file


def create_index_dict(index_values):
    index_dict = {}
    for index, value in enumerate(index_values):
        index_dict[value] = index
    return index_dict


def filter_sparse_matrix(sparse_matrix, filtered_indices, cell_stats):
    """
    Rebuild sparse matrix to contain only the cell barcodes that were called a cell in the allCells.csv from RNA.
    Cell barcodes can have passed in RNA that are not present in the raw matrix generated from the enriched library, and so we cannot do a simple subsetting
    Args:
        sparse_matrix: raw umi or read csr sparse matrix in guides x cells orientation
        filtered_indices: passing barcodes from RNA allCells.csv
        cell_stats: cell metrics file to provide raw barcodes list
    Returns:
        filtered_sparse_matrix: counts matrix in guides x cells orientation where cells are only barcodes that were called a cell in ScaleRNA, with empty columns for those that were not detected in enriched library
        filtered_column_mapping: dictionary containing keys as cell barcode alia
    """
    # Load cell barcodes from the provided TSV file
    cell_mapping = list(cell_stats.index.values)
    
    # Initialize new column mapping and filtered sparse matrix
    filtered_column_mapping = {alias: idx for idx, alias in enumerate(filtered_indices)}
    
    filtered_sparse_matrix = lil_matrix((sparse_matrix.shape[0], len(filtered_indices)))

    # Populate filtered sparse matrix based on original column mapping and new filtered column mapping
    for alias, idx in filtered_column_mapping.items():
        if alias in cell_mapping:
            row_index = cell_mapping.index(alias)  # Get the index of the alias in the row_mapping
            filtered_sparse_matrix[:, idx] = sparse_matrix[:, row_index]
    
    filtered_sparse_matrix = filtered_sparse_matrix.tocsr()

    return filtered_sparse_matrix, filtered_column_mapping


def guide_dictionary(guide_file_path: Path):
    """Read in whitelist for dataset features and return dictionary with keys as name of """
    with open(guide_file_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        guides = [row[1] for row in reader]
        guide_mapping = {col: idx for idx, col in enumerate(guides)}
    return guide_mapping


def wells_for_range(range_str, all_wells):
    """Return a list of wells based on column value in samples.csv"""
    wells = []
    for item in range_str.split(';'):
        input = item.replace(' ', '')
        if '-' in input:
            start, end = input.split('-')
            wells.extend(all_wells[all_wells.index(start):all_wells.index(end)+1])
        else:
            wells.append(input)
    return wells


def hash_to_well(hash_combos: Path, expected_combos=None):
    """Return dict of hash combos to well"""
    with open(hash_combos, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        hash_combos_dict = {row[0]: row[1] for row in reader}
    if expected_combos:
        expected_fixation_wells = wells_for_range(expected_combos, list(hash_combos_dict.values()))
        hash_combos_dict = {k: v for k, v in hash_combos_dict.items() if v in expected_fixation_wells}
    return hash_combos_dict


def write_keys_to_tsv(dictionary, tsv_file):
    with open(tsv_file, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        for key in dictionary.keys():
            writer.writerow([key])

def compute_background_series(sparse_matrix, top_n=2):
    """
    Builds background estimation per hash, where background is considered to be any non-zero and non-top n values per cell per hash
    Args:
        sparse_matrix: filtered hashes x cells UMI counts matrix
        top_n: number of top values within each cell (column) to ignore when estimating background.
    """
    background_series = []
    for row_idx in range(sparse_matrix.shape[0]):
        row_data = sparse_matrix.getrow(row_idx).toarray().flatten().astype(float)  # Extract row data as dense array
        for col_idx in range(sparse_matrix.shape[1]):
            col_data = sparse_matrix[:, col_idx].toarray().flatten().astype(float)  # Extract column data as dense array
            top_values = np.partition(col_data, -top_n)[-top_n:]  # Find the top N values
            if row_data[col_idx] in top_values:
                row_data[col_idx] = 0  # Set row's value to 0 if it's in top N values of the column
        row_data[row_data == 0] = np.nan  # Replace zeros with NaN
        background_series.append(np.nanmedian(row_data))  # Compute median after excluding top N values
    background_series = pd.Series(background_series, index=range(sparse_matrix.shape[0]))
    background_series = background_series.fillna(0)
    return background_series

def corrected_p_values_poisson(sparse_matrix, background, features_dict, threshold):
    observed_counts = sparse_matrix.toarray()
    expected_counts = 3 * background.to_numpy()

    expected = np.tile(expected_counts, (sparse_matrix.shape[1], 1)).T  # Repeat expected counts for each sample
    p_value = 1 - poisson.cdf(observed_counts - 1, expected)
    results = pd.DataFrame(p_value)

    results_T = results.T
    results_corrected = results_T.apply(lambda col: false_discovery_control(col, method='bh'), axis=0)
    results_corrected.columns = features_dict.keys()

	#Did a hash pass the background test

    passing_hashes = results_corrected.apply(
    lambda row: ';'.join([col for col in results_corrected.columns[row < threshold]])
    if any(row < threshold) else "No_Pass",
    axis=1
	)
    passing_hashes_df = pd.DataFrame(passing_hashes, columns = ['passing_scaleplex'])
    return passing_hashes_df

def assignment_criteria_bg(row, hash_well_dict, toptwo_frac):
    #Failed Background Test
    assigned_hash = ''
    if (row['passing_scaleplex'] == 'No_Pass'):
        assigned_hash = codes.indeterminate
    #Passed Enrich Threshold and had at least one hash passing Background
    elif row['topTwo'] > toptwo_frac:
        max_cols = row['topTwo_scaleplex']
        #Specific Hash Passed Background and was also Top hash in that cell
        
        strings_column1 = row['topTwo_scaleplex'].split(';')
        strings_column2 = row['passing_scaleplex'].split(';')	
        if all(string in strings_column2 for string in strings_column1):
            assigned_hash = hash_well_dict[max_cols] if max_cols in hash_well_dict else codes.unexpected
        #Purity Threshold Passed, but Top hash was not the one that passed background
        else:
            assigned_hash = codes.max_fail
    #Had a hash that passed background, but cell did not pass topTwo threshold
    else:
        assigned_hash = codes.enrich_fail
    return assigned_hash
          
def assignment_criteria_fc(row, hash_well_dict, fc_threshold):
    assigned_hash = ''
    if row['second'] != 0 and (row['third'] == 0 or (row['second'] / row['third']) > fc_threshold):
        max_cols = row['topTwo_scaleplex'] 
        assigned_hash = hash_well_dict[max_cols] if max_cols in hash_well_dict else codes.unexpected
    else:
        assigned_hash = codes.enrich_fail
    return assigned_hash


def hash_assignment(sparse_matrix, cell_stats, features_dict, assignment_method, toptwo_frac, hash_well_dict, fc_threshold):
    """Creates a list of guides (columns) that passed the minimum UMI threshold within each cell"""
    assignment_df = pd.DataFrame(index=cell_stats.index)
    if (assignment_method == 'bg'):
        background_series = compute_background_series(sparse_matrix, top_n=2)
        passing_hashes_df = corrected_p_values_poisson(sparse_matrix, background_series, features_dict, 0.05)
        assignment_df['passing_scaleplex'] = passing_hashes_df['passing_scaleplex'].values
        assignment_df['assigned_scaleplex'] = cell_stats.join(assignment_df).apply(assignment_criteria_bg, axis=1,
                                                                                   hash_well_dict=hash_well_dict,
                                                                                   toptwo_frac=toptwo_frac)
    # Iterate over rows
    if (assignment_method == 'fc'):
       assignment_df['assigned_scaleplex'] = cell_stats.apply(assignment_criteria_fc, axis=1, hash_well_dict=hash_well_dict, fc_threshold=fc_threshold) 

    return assignment_df


def main(umi_matrix: Path, all_cells: Path, cell_stats: Path, references: Path, lib_struct: Path, assignment_method: str, toptwo_frac: float, fc_threshold: float, expected_combos: str, id: str, out_dir: Path):

    metrics = Metrics()

    guide_matrix = mmread(umi_matrix)
    guide_matrix = guide_matrix.tocsr() # cast to csr for formatting necessary for filtering indexing, otherwise casts as native coo matrix
    rna = pd.read_csv(all_cells, index_col=0)
    rna = rna[rna['pass']] # only join with passing cells in RNA
    cell_stats_df = pd.read_csv(cell_stats, index_col='Cell_Barcode')
    
    # Create output directories
    matrix_dir = out_dir / f"{id}.filtered.matrix"
    raw_matrix_dir = out_dir / f"{id}.raw.matrix"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    raw_matrix_dir.mkdir(parents=True, exist_ok=True)

    #write raw barcodes
    raw_barcodes = cell_stats_df.index.to_series()
    raw_barcodes.to_csv(raw_matrix_dir / "barcodes.tsv", sep="\t", header=False, index=False)

    guide_matrix_filtered, filtered_cell_mapping = filter_sparse_matrix(guide_matrix, list(rna.index), cell_stats_df)

    mmwrite(f"{matrix_dir}/matrix.mtx", guide_matrix_filtered)

    filtered_cell_barcodes = pd.DataFrame(filtered_cell_mapping.keys())
    filtered_cell_barcodes.to_csv(f"{matrix_dir}/barcodes.tsv", sep='\t', index=False, header=False)

    # Create merged allCells with all rows from RNA and ENRICH cell metrics
    cell_stats_df = cell_stats_df.join(rna, how='outer', lsuffix='_ENRICH', rsuffix='_RNA')
    cell_stats_df['pass'] = cell_stats_df['pass'].fillna(False)
    enrich_alias_df = cell_stats_df.filter(regex='alias_ENRICH', axis=1) # Keep alias columns from the enriched library to default to for barcodes not present in RNA
    cell_stats_df = cell_stats_df.drop(enrich_alias_df.columns, axis=1) # remove original alias columns to handle case where cell only observed in RNA does not have NaN or empty values for alias values
    enrich_alias_df = enrich_alias_df.rename(columns=lambda x: x.replace('_alias_ENRICH', '_alias')) # rename alias columns for fillna
    cell_stats_df = cell_stats_df.drop(rna.columns.drop('pass'), axis=1, errors='ignore') # remove columns (except pass column) found only in the RNA allCells, ignoring cases where some are not found because it was a column in both dfs and now has a suffix
    cell_stats_df = cell_stats_df.rename(columns=lambda x: x.replace('_alias_RNA', '_alias')) # rename RNA alias columns to be final alias col format
    cell_stats_df = cell_stats_df.fillna(enrich_alias_df) # fillna with enriched library alias values
    cell_stats_df = cell_stats_df.fillna(0) # fill remaining NaN values with 0
    cell_stats_df = cell_stats_df.drop(cell_stats_df.filter(regex='_RNA').columns, axis=1) # drop remaining columns from RNA table
    cell_stats_df = cell_stats_df.rename(columns=lambda x: x.replace('_ENRICH', '')) # rename columns for compatibility with allCells.csv
    
    hash_combos, guide_file = preprocessJson(lib_struct)
    guide_path = f'{references}/{guide_file}'
    guide_dict = guide_dictionary(guide_path)
    hash_well_dict = hash_to_well(references / hash_combos, expected_combos)

    assignment_df = hash_assignment(guide_matrix_filtered, cell_stats_df.loc[rna.index], guide_dict, assignment_method, toptwo_frac, hash_well_dict, fc_threshold) #hash assignemt of any hash that are = max UMI value for that cell
    cell_stats_df = cell_stats_df.join(assignment_df) # add assignment for RNA passing cells
    cell_stats_df['sample'] = id
    # put passing RNA cells first in csv, followed by the remaining barcodes from hash library
    cell_stats_df = cell_stats_df.loc[rna.index.append(cell_stats_df.index.drop(rna.index))]
    cell_stats_df.to_csv(out_dir / f"{id}.ScalePlex.allCells.csv")

    assigned_cells = rna.join(assignment_df) # add assignment for RNA passing cells
    assigned_cells = assigned_cells[~assigned_cells['assigned_scaleplex'].isin(codes.errors)] # remove cells with assignment errors
    scaleplex_stats = assigned_cells.groupby('assigned_scaleplex').agg(
        passing_reads=('reads', 'sum'),
        cells_called=('assigned_scaleplex', 'size'),
        mean_passing_reads=('reads', 'mean'),
        median_utc=('umis', 'median'),
        median_genes=('genes', 'median')
    ).sort_values('cells_called', ascending=False).astype('int')
    scaleplex_stats.to_csv(out_dir / f"{id}.scaleplex_stats.csv")

    cell_stats_passing = cell_stats_df[cell_stats_df['pass']]
    metrics.meanReadsPerCell = round(cell_stats_passing['reads'].mean(0), 1)
    metrics.nUMIPerCell = round(cell_stats_passing['umis'].median(0), 1)
    metrics.saturation = round(cell_stats_passing['Saturation'].median(0), 1)
    metrics.readsInCells = cell_stats_passing['reads'].sum() / cell_stats_df['reads'].sum()

    assignment_prop = cell_stats_passing['assigned_scaleplex'].value_counts(normalize=True)
    metrics.passingPercent = assignment_prop.drop(codes.errors, errors='ignore').sum()
    # Report out different assignment error stats
    metrics.indeterminatePercent = assignment_prop[codes.indeterminate] if codes.indeterminate in assignment_prop else 0
    metrics.maxFailPercent = assignment_prop[codes.max_fail] if codes.max_fail in assignment_prop else 0
    metrics.enrichFailPercent = assignment_prop[codes.enrich_fail] if codes.enrich_fail in assignment_prop else 0
    metrics.unexpectedPercent = assignment_prop[codes.unexpected] if codes.unexpected in assignment_prop else 0

    metrics.print(out_dir / "metrics.csv")

    #write features.tsv
    write_keys_to_tsv(guide_dict, f"{matrix_dir}/features.tsv")
    write_keys_to_tsv(guide_dict, f"{raw_matrix_dir}/features.tsv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate Filtered Hash umi counts matrix and cellMetrics that share the same cell barcodes as the passing RNA analysis per sample')
    parser.add_argument('--umi_matrix', metavar='MATRIX.mtx', type=Path,
                        help='count_hash output of hash x cells raw UMI matrix')
    parser.add_argument('--all_cells', metavar='ALLCELLS.csv', type=Path,
                        help='nf-rna per-cell summary statistics output to filter based on pass status')
    parser.add_argument('--cell_stats', metavar='CELLMETRICS.csv', type=Path,
                        help='count_hash output of cell metadata for HASH cell guide UMI matrix')
    parser.add_argument('--references', metavar='REFERENCES', type=Path,
                        help='Path to folder containing barcode whitelists')
    parser.add_argument('--lib_struct', metavar='LIB_STRUCT', type=Path,
                        help='library structure json used for bcParser in enrich detection')
    parser.add_argument('--assignment_method', metavar='ASSIGNMENT_METHOD', type=str, default='bg',
                        help='Sample background based, or fold change based assignment')
    parser.add_argument('--toptwo_frac', metavar='TOPTWO_FRAC', type=float,
                        help='Top two unique scaleplex UMIs must be over this fraction to be assigned')
    parser.add_argument('--fc_threshold', metavar='FC_THRESHOLD', type=float,
                        help='FC threshold for fold change based assignment')
    parser.add_argument('--expected_combos', metavar='FIXATION_PLATE_WELLS', type=str,
                        help='Range or wells from fixation plate to supervise hash assignment')
    parser.add_argument("--id", required=True, help="Sample and lib name")
    parser.add_argument('--outDir', type=Path,
                        help='Directory for outputs')
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(args.umi_matrix, args.all_cells, args.cell_stats, args.references, args.lib_struct, args.assignment_method, args.toptwo_frac, args.fc_threshold, args.expected_combos, args.id, args.outDir)
