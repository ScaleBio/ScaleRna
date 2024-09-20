#!/usr/bin/env python
""" Generate raw UMI counts matrices from bcParser barcodes.csv output files """
import argparse
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, Tuple
from collections import defaultdict
import csv
import os
import numpy as np
import pandas as pd
import json
from scipy.sparse import csr_matrix
from scipy.io import mmwrite

# Positive matches in bcParser output
MATCHSTATUS = ("exact", "corrected")

@dataclass
class cellReadCounts:
    """Track reads counts for cell (e.g cell X read combination)"""
    counts: Dict[Tuple[str, ...], int] # Cell Sequence -> ReadCount
    def __init__(self):
        self.counts = {}
    
    def Hash_nReads(self, read_thres: int):
        return sum(1 for i in self.counts.values() if i >= read_thres)
    
    def addCount(self, cellSeq: Tuple[str, ...]):
        self.counts[cellSeq] = self.counts.get(cellSeq, 0) + 1


class UmiReadCounts:
    """Track reads counts for each UMI within one target (e.g cell X guide combination)"""
    counts: Dict[str, int] # UMI-sequence -> ReadCount
    def __init__(self):
        self.counts = {}
    
    def nUmis(self, read_thres: int):
        return sum(1 for i in self.counts.values() if i >= read_thres)
    
    def addRead(self, umiSeq: str):
        self.counts[umiSeq] = self.counts.get(umiSeq, 0) + 1


class GuideReadCounts:
    """UMI counts for each guide in a cell"""
    guides: Dict[str, UmiReadCounts] # GuideSequence -> UMI-Counts
    def __init__(self):
        self.guides = {}
    
    def umiCounts(self, read_thres: int):
        res = {}
        for guide in self.guides:
            res[guide] = self.guides[guide].nUmis(read_thres)
        return res
    
    def guideReadCounts(self, umi: str):
        res = {}
        for guide in self.guides:
            res[guide] = self.guides[guide].counts.get(umi,0)
        return res

    def addRead(self, guide: str, umi: str):
        if not guide in self.guides:
            self.guides[guide] = UmiReadCounts()
        self.guides[guide].addRead(umi)


def preprocessJson(libStructJson: Path):
    """ read in libstuct json and parse barcode levels for generation of cell barcode tuple and the whitelist used for detection"""
    libStruct = json.load(open(libStructJson))
    aliases_list = list()
    names_list = list()
    for bc in libStruct['barcodes']:
        if bc.get('type') not in ['library_index', 'umi']:
            if bc.get('name') not in ['scaleplex']:
                aliases_list.append(bc.get('alias'))
                names_list.append(bc['name'])
        if bc.get('name') in ['scaleplex']:
            guide_file = bc.get('sequences')

    return aliases_list, names_list, guide_file


def countGuideReads(barcodesCsv: Path, namesList: list):
    """Count reads for each cell X guide X UMI combo from bcParser barcodes.csv file output"""
    # BC -> (GUIDE,UMI) -> ReadCount
    # BC is a tuple of three strings for the three cell-barcode levels (RT, Ligation, PCR)
    cellGuideReads: Dict[Tuple[str, ...], GuideReadCounts] = defaultdict(GuideReadCounts)
    cellReads = cellReadCounts()
    info_cols = [name_val + "_info" for name_val in namesList]
    for row in csv.DictReader(open(barcodesCsv),delimiter="\t"):
        valid_barcode = True
        for col_name in info_cols: # validate each value of the info col's value agains the match status, and break row processing if any fail
                if col_name in row:
                    row_col_value = row[col_name].split(",")[0].strip()
                    if row_col_value not in MATCHSTATUS:
                        valid_barcode = False
                        break
        if valid_barcode:
            cell = tuple(row[col_name] for col_name in namesList)
            cellReads.addCount(cell) #Count a read for the cell, having had a correct match to the barcode info
            if not row['scaleplex_info'].split(',')[0] in MATCHSTATUS:
                cellGuideReads[cell].addRead('noGuide',row['Read']) #Count reads that do not have a guide sequence but did have a correct cellular barcode
                continue
            cellGuideReads[cell].addRead(row['scaleplex'],row['umi'])
    return cellGuideReads, cellReads


def create_sparse_matrix(counts_data_dict, column_file):
    """
    Generate the sparse matrix for counts and UMI's from the dictionary output by countGuideReads BC -> (GUIDE,UMI) -> ReadCount
    Args:
        counts_data_dict: dictionary output by countGuideReads in the format of BC -> (GUIDE,UMI) -> ReadCount
        column_file: the guide whitelist used by bcParser for guide detection, first column is barcode sequences and second column are the aliases
    Returns:
        umi_sparse_matrix: csr sparse matrix of umi counts, in cells x guides orientation
        reads_sparse_matrix: csr sparse matrix of read counts, in cells x guides orientation
        cell_mapping: dictionary of mapping of cell barocdes to their index with tuple of cell barcode alias's as keys, and corresponding index in matrix as values (rows here)
        guide_mapping: dictionary of mapping of guide alias to their index in the matrix (columns here)
        no_guide_df: single column pd.Dataframe containing the number of reads with no guide per cell. Eventually appends to cellMetrics and converted to a percentage
    """

    # Read in guide whitelist used for detection in bcParser to parse into a list that will become our columns in the matrix being built
    with open(column_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        guides = [row[1] for row in reader]

    # Create mappings for row and column indices
    cell_mapping = {row: idx for idx, row in enumerate(counts_data_dict.keys())}
    guide_mapping = {col: idx for idx, col in enumerate(guides)}

    # Initialize lists to store row indices, column indices, and data for the sparse matrix
    rows = []
    cols = []
    umi_data = []
    no_guide_values = {row_key: 0 for row_key in counts_data_dict.keys()} 

    # Populate the lists with data from the dictionary generated by countGuideReads, BC -> (GUIDE,UMI) -> ReadCount
    for row_key in counts_data_dict.keys():
        row_idx = cell_mapping[row_key]
        for guide in counts_data_dict[row_key].guides:
            if guide in guide_mapping:
                col_idx = guide_mapping[guide]
                rows.append(row_idx)
                cols.append(col_idx)
                umi_data_val = counts_data_dict[row_key].guides[guide].nUmis(1)
                umi_data.append(umi_data_val)
            elif guide == "noGuide":
                no_guide_values[row_key] = counts_data_dict[row_key].guides[guide].nUmis(1)

    # Create the sparse matrix, cells x guides
    umi_sparse_matrix = csr_matrix((umi_data, (rows, cols)), shape=(len(cell_mapping), len(guide_mapping)))

    # Create a separate series for "noMatch" values mapped according to row_mapping
    no_guide_series = pd.Series({row_key: no_guide_values[row_key] for row_key in counts_data_dict.keys()})

    no_guide_df = no_guide_series.to_frame(name='noGuide')

    return umi_sparse_matrix, cell_mapping, guide_mapping, no_guide_df


def find_nth_highest_value_per_row(sparse_matrix, n):
    """Return value for the nth highest detected hash per cell"""
    # Initialize an array to store nth highest values
    nth_highest_values = np.zeros(sparse_matrix.shape[0])

    # Iterate over rows
    for i in range(sparse_matrix.shape[0]):
        # Get row data
        row_data = sparse_matrix.getrow(i).data

        # Find the nth highest value
        if len(row_data) >= n:
            row_data_sorted = np.sort(row_data)
            nth_highest_values[i] = row_data_sorted[-n]
        else:
            # Handle case where there are less than n elements in the row
            nth_highest_values[i] = 0

    return nth_highest_values

def top_two_features(sparse_matrix, feature_mapping):
    """Return list of top two column's per cell barcode"""
    # Convert the feature_mapping dictionary to a list for index-based lookup
    idx_to_feature = {v: k for k, v in feature_mapping.items()}
    
    # Initialize the result list
    result = []
    
    # Iterate over each row
    for row in range(sparse_matrix.shape[0]):
        # Extract the row as a dense array
        row_data = sparse_matrix[row].toarray().flatten()
        
        # Find indices of the two largest values
        if len(row_data) >= 2:
            top_indices = np.argpartition(row_data, -2)[-2:]
            top_indices = top_indices[np.argsort(-row_data[top_indices])]
        elif len(row_data) == 1:
            top_indices = [0, 0]
        else:
            top_indices = []
        
        # Map the indices to feature names
        top_features = [idx_to_feature[i] for i in top_indices]
        
        # Format the result for this row
        result.append(";".join(sorted(top_features)))
    
    return result


def tsv_to_dict(parent_directory, file_path):
    """
    Reads a TSV file with two columns and returns a dictionary
    where the first column values are keys and the second column
    values are values.
    """
    result = {}
    full_file_path = os.path.join(parent_directory, file_path)
    with open(full_file_path, 'r') as tsv_file:
        for line in tsv_file:
            barcode, name = line.strip().split('\t')
            result[name] = barcode
    return result


def create_level_dicts(json_file_path, parent_directory, selected_aliases):
    """
    Creates dictionaries for each alias specified in the JSON file
    using the provided list of selected aliases.
    """
    with open(json_file_path, 'r') as json_file:
        json_data = json.load(json_file)
    
    level_dicts = {}
    
    # Iterate over each entry in the "barcodes" section
    for entry in json_data.get('barcodes', []):
        alias = entry.get('alias')
        if alias in selected_aliases:
            file_path = entry.get('sequences')
            
            # Read the TSV file and construct the dictionary
            level_dict = tsv_to_dict(parent_directory, file_path)
            
            # Store the dictionary using the alias as the key
            level_dicts[alias + '_alias'] = level_dict
    
    return level_dicts


def cellBarcodeConcat(lib_struct_dir,lib_struct_json, valid_alias, rows_dict):
    """
    Parse alias in cellMetrics object per cell and create a series of the actual barcode sequences
    concatenated together like we do for RNA
    """
    def concatenate_row(row):
        return ''.join(str(value) for value in row)

    levels_dict = create_level_dicts(lib_struct_json,lib_struct_dir,valid_alias)
    alias_cols = [alias_val + "_alias" for alias_val in valid_alias]
    alias_df = pd.DataFrame(rows_dict.keys(), columns = alias_cols)
    for column in alias_df.columns:
        # Check if the column corresponds to a key in level_dicts
        if column in levels_dict:
        # Replace the values in the column with the corresponding values from level_dicts
            alias_df[column] = alias_df[column].replace(levels_dict[column])
    alias_df['Cell_Barcode'] = alias_df.apply(concatenate_row, axis=1)
    cell_barcode = alias_df['Cell_Barcode']
    return cell_barcode

def main(barcodesCsv: Path, references: Path, lib_struct: Path, id: str, outDir: Path):

    aliases_list, names_list, guide_file = preprocessJson(lib_struct)
    alias_cols = [alias_val + "_alias" for alias_val in aliases_list]

    cellGuideReads, cellReads = countGuideReads(barcodesCsv, names_list)
    # Create  a table with UMI counts for each cell X guide combo
    cellMetrics = pd.DataFrame.from_dict(cellReads.counts, orient='index', columns=['reads'])
    cellMetrics.index = pd.MultiIndex.from_tuples(cellMetrics.index, names = alias_cols)

    guide_path = f'{references}/{guide_file}'
    umi_sparse_counts_matrix, rows_mapped, cols_mapped, noGuides = create_sparse_matrix(cellGuideReads,guide_path)

    cell_barcodes = cellBarcodeConcat(references, lib_struct, aliases_list, rows_mapped)

    cellMetrics['Cell_Barcode'] = cell_barcodes.values

    # Continue updating cellMetrics df

    cellMetrics['noScalePlex'] = noGuides['noGuide']/cellMetrics['reads'] #Percentage of reads with no guide detected
    cellMetrics['umis'] = umi_sparse_counts_matrix.sum(axis=1) # number of unique guide UMI detected 
    cellMetrics['scaleplex'] = umi_sparse_counts_matrix.getnnz(axis=1) # number of unique guide's detected, like unique genes detected
    cellMetrics['max'] = umi_sparse_counts_matrix.max(axis=1).toarray().flatten() # max number of UMI detected across the guides within a cell
    cellMetrics['second'] = find_nth_highest_value_per_row(umi_sparse_counts_matrix, 2) # second highest number of UMI of any guide within a cell
    cellMetrics['third'] = find_nth_highest_value_per_row(umi_sparse_counts_matrix, 3) # third highest number of UMI of any guide within a cell
    cellMetrics['purity'] = (cellMetrics['max'] / cellMetrics['umis']) 
    cellMetrics['topTwo'] = ((cellMetrics['max']+cellMetrics['second']) / cellMetrics['umis'])
    cellMetrics['minorFrac'] = cellMetrics['second'] / cellMetrics['max']
    cellMetrics['Saturation'] = 1 - (cellMetrics.umis / (cellMetrics.reads - (cellMetrics.reads * cellMetrics.noScalePlex))) #Saturation 
    cellMetrics['topTwo_scaleplex'] = top_two_features(umi_sparse_counts_matrix,cols_mapped)
    
    cellMetrics = cellMetrics.reset_index()

    # handle cells with unassigned guides in the non-filtered matrix,
    # transposed to be guides x cells to be consistent with rna expression data
    umi_sparse_counts_matrix = umi_sparse_counts_matrix.transpose()
    mmwrite(f"{outDir}/{id}.raw.umi.mtx", umi_sparse_counts_matrix)

    cellMetrics = cellMetrics.set_index('Cell_Barcode')
    cellMetrics.to_csv(outDir / f"{id}.cellMetrics.csv")    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Count ScalePlex sequences matches from bcParser output')
    parser.add_argument('barcodesCsv', metavar='BARCODES.csv', type=Path,
                        help='bcParser per-sample barcode match output')
    parser.add_argument('--references', metavar='REFERENCES', type=Path,
                        help='Path to folder containing barcode whitelists')
    parser.add_argument('--lib_struct', metavar='LIB_STRUCT', type=Path,
                        help='library structure json used for bcParser in enrich detection')
    parser.add_argument('--id', metavar='ID', type=str,
                        help='Sample and lib name for output files')
    parser.add_argument('--outDir', type=Path,
                        help='Directory for outputs')
    args = parser.parse_args()
    args.outDir.mkdir(parents=True, exist_ok=True)
    main(args.barcodesCsv, args.references, args.lib_struct, args.id, args.outDir)
