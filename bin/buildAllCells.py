#!/usr/bin/env python
"""
Calculates metrics for one sample from STARsolo outputs. Takes in the path to the STARsolo output directory for one sample and returns an allCells.csv
"""


######################
###  DEPENDENCIES  ###
######################


import argparse
from dataclasses import dataclass
import itertools as it
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict
from scale_utils import io


########################
###  CORE FUNCTIONS  ###
########################

def read_genes(sample_specific_file_paths):
    """
    Reads in the features.tsv file from the STARsolo output directory, labelling each gene as either mm (mouse) or hs (human).

    Args:
        sample_specific_file_paths (Dict[str, Path]): Dictionary where each key is a custom file identifier and each value is a path to the identified file in the STARsolo output directory
    
    Returns:
        A dataframe of genes (DataFrame)
    """
    genes = pd.read_csv(sample_specific_file_paths['features'], sep = "\t", header = None)
    genes.columns = ['ID', 'name', 'type']
    genes.index = genes.index + 1
    genes['species'] = genes.ID.map(lambda x: ["mm", "hs"][x.startswith("ENSG")])

    return genes


def read_barcodes(sample_specific_file_paths):
    """
    Reads in the barcodes.tsv file from the STARsolo output directory.

    Args:
        sample_specific_file_paths (Dict[str, Path]): Dictionary where each key is a custom file identifier and each value is a path to the identified file in the STARsolo output directory

    Returns:
        A dataframe of barcodes (DataFrame)
    """
    barcodes = pd.read_csv(sample_specific_file_paths['barcodes'], sep = '\t', header = None, names = ['barcode'])
    barcodes.index = barcodes.index + 1

    return barcodes


def count_transcripts(sample_specific_file_paths, genes, barcodes, isBarnyard):
    """
    Computes barnyard-specific transcript counts from the raw STARsolo count matrix for each barcode.

    Args:
        sample_specific_file_paths (Dict[str, Path]): Dictionary where each key is a custom file identifier and each value is a path to the identified file in the STARsolo output directory
        genes (DataFrame): A dataframe of genes
        barcodes (DataFrame): A dataframe of barcodes
        isBarnyard (bool): Whether or not this is a barnyard sample (i.e. mixed mouse and human)

    Returns:
        A dataframe of unique transcript counts for mouse and human per barcode (DataFrame)
    """
    barcode_counts = [BarcodeCounts() for barcode in barcodes.barcode]
    with open(sample_specific_file_paths['mtx'], 'r') as matrix_file:
        # Skip header
        for i in range(3):
            matrix_file.readline()
            
        for line in matrix_file:
            split_line = line.split()

            # Individual expression counts can be fractions, due to multimapper resolution
            gene, barcode, count = int(split_line[0]), int(split_line[1]), float(split_line[2])
            barcode_id = barcode - 1
                
            if count == 0:
                continue
                
            if isBarnyard:               
                if genes.species[gene] == "hs":
                    barcode_counts[barcode_id].genes += 1
                    barcode_counts[barcode_id].human_genes += 1
                    barcode_counts[barcode_id].human_umis += count
                    barcode_counts[barcode_id].umis += count                   
                elif genes.species[gene] == "mm":
                    barcode_counts[barcode_id].genes += 1
                    barcode_counts[barcode_id].mouse_genes += 1
                    barcode_counts[barcode_id].mouse_umis += count
                    barcode_counts[barcode_id].umis += count          
                else:
                    print(genes.iloc[gene])                    
            else:
                barcode_counts[barcode_id].genes += 1
                barcode_counts[barcode_id].umis += count

    # For consistency with STARsolo CellReadMetrics output, we round the total expression (UMI) count
    for barcode in barcode_counts:
        barcode.umis = round(barcode.umis)

    umi_counts_by_species = pd.DataFrame(barcode_counts, index = barcodes.barcode)

    return umi_counts_by_species


def build_allCells(sample_specific_file_paths, umi_counts_by_species, isBarnyard):
    """
    Builds the allCells.csv file for this sample.

    Args:
        sample_specific_file_paths (Dict[str, Path]): Dictionary where each key is a custom file identifier and each value is a path to the identified file in the STARsolo output directory
        umi_counts_by_species (DataFrame): Dataframe of unique transcript counts for mouse and human per barcode
        isBarnyard (bool): Whether or not this is a barnyard sample (i.e. mixed mouse and human)

    Returns:
        An allCells.csv containing metrics computed across all barcodes in this sample (DataFrame)
    """
    allCells_stats = pd.read_csv(sample_specific_file_paths['stats'], sep = '\t', index_col = 'CB')[1:]
    allCells = pd.DataFrame(index = allCells_stats.index.to_list())

    # For barnyard samples, only consider genes from the primary species represented among barcodes (not the background from the other species)
    if isBarnyard:
        allCells['human_genes'] = umi_counts_by_species['human_genes']
        allCells['mouse_genes'] = umi_counts_by_species['mouse_genes']
        allCells['human_umis'] = umi_counts_by_species['human_umis']
        allCells['mouse_umis'] = umi_counts_by_species['mouse_umis']
        allCells['genes'] = umi_counts_by_species[['human_genes', 'mouse_genes']].max(1)
        allCells['umis'] = umi_counts_by_species[['human_umis', 'mouse_umis']].max(1)
    else:
        allCells['genes'] = umi_counts_by_species['genes']
        allCells['umis'] = umi_counts_by_species['umis']
    
    # Compute metrics for each barcode
    # See https://github.com/alexdobin/STAR/discussions/1826#discussioncomment-5596619 for a description of the metrics
    allCells['reads'] = allCells_stats.cbMatch
    allCells['mappedReads'] = allCells_stats.genomeU + allCells_stats.genomeM
    allCells['geneReads'] = allCells_stats.featureU + allCells_stats.featureM
    allCells['exonReads'] = allCells_stats.exonic
    allCells['intronReads'] = allCells_stats.intronic
    allCells['antisenseReads'] = allCells_stats.exonicAS
    allCells['passingReads'] = allCells_stats.countedU + allCells_stats.countedM
    allCells['uniquePassingReads'] = allCells_stats.countedU
    allCells['Saturation'] = 1 - (allCells_stats.nUMIunique / allCells_stats.countedU)
    allCells['mitoReads'] = allCells_stats.mito
    allCells['mitoProp'] = allCells['mitoReads'] / allCells['mappedReads']

    return allCells




def split_barcodes(libStruct, allCells):
    """
    Uses 'level' and 'alias' values in each element of the library json barcodes object to split barcodes, with the assumption that each element in the 'barcodes' section of the library json defining part of the 'CB' cell barcode sequence has a 'level' field and an 'alias' field (see the current library structure jsons in the /references folder for examples).

    Args:
        libStruct (Path): Path to the library structure json for this sample
        allCells (DataFrame): An allCells.csv containing metrics computed across all barcodes in this sample 
    
    Returns:
        The allCells.csv for this sample, with updated columns to show barcode splits (DataFrame)
    """
    updated_allCells = allCells

    # Directory containing the library structure json and associated sequence files
    libDir = libStruct.parent

    libJson = io.readJSON(libStruct, preserveDictOrder = True)
    barcodeInfo = libJson['barcodes']

    # We need to keep the barcodes in the original order from the library json so that it matches the order in which they were concatenated by bcParser (cell barcode / index)
    # Hence we rely on barcode in lib.json being grouped by level to make the group-by call work
    # cellBarcodesByLevel.sort(key=lambda bc: bc['level'])
    cellBarcodesByLevel = [barcodeMeta for barcodeMeta in barcodeInfo if barcodeMeta.get('type', None) != 'library_index' and barcodeMeta.get('type', None) != 'umi']

    # Overall cell barcode is the concatenation of all individual barcodes
    curBcSeqStartIndex = 0 
    for level, bcParts in it.groupby(cellBarcodesByLevel, key = lambda x: x['level']):
        # Single barcode level might have multiple sub-barcodes (e.g. bead barcode blocks); these will be concatenated and reported as one barcode
        bcParts = list(bcParts)
        bcInfo = bcParts[0]
        barcodeFullName = bcInfo['alias']
        barcodeLen = sum([x['length'] for x in bcParts])
        bcSubSeq = slice(curBcSeqStartIndex, barcodeLen + curBcSeqStartIndex)
        barcodes = [barcode[bcSubSeq] for barcode in updated_allCells.index]
        updated_allCells[barcodeFullName] = barcodes
        curBcSeqStartIndex += barcodeLen

        # Load the ids for specific barcode sequences (e.g. well coordinates) for the current barcode level
        bcIds = {}
        with open(libDir/f"{bcInfo['sequences']}") as f:
            for line in f:
                splat = line.rstrip().split('\t')
                bcIds[splat[0]] = splat[1]

        bcSeqIds = updated_allCells[barcodeFullName].apply(lambda seq: bcIds[seq])
        updated_allCells[f'{barcodeFullName}_alias'] = bcSeqIds

    return updated_allCells


@dataclass
class BarcodeCounts:
    """
    A helper class to store UMI and gene counts for a particular barcode.
    """
    # Total number of genes detected
    genes: int = 0

    # Human genes (for barnyard only)
    human_genes: int = 0 

    # Mouse genes (for barnyard only)
    mouse_genes: int = 0

    # Total number of UMIs detected
    umis: int = 0

    # Human UMIs (for barnyard only)
    human_umis: int = 0

    # Mouse UMIs (for barnyard only)
    mouse_umis: int = 0

#####################
###  MAIN METHOD  ###
#####################


def main():
    parser = argparse.ArgumentParser()

    # Required and optional arguments for specifying the STARsolo outputs for this sample
    parser.add_argument("--STARsolo_out", type = Path, required = True, help = "Path to the STARsolo outputs for this sample.")
    parser.add_argument("--feature_type", type = str, required = False, default = 'GeneFull_Ex50pAS', help = "STARsolo feature type used.")
    parser.add_argument("--matrix_type", type = str, required = False, default = 'UniqueAndMult-PropUnique.mtx', help = "STARsolo matrix type used.")
    parser.add_argument("--isBarnyard", default = False, action = "store_true", help = "If set, this sample will be interpreted as a barnyard sample (i.e. mixed mouse and human).")

    # Optional argument to specify the name of the sample for which cells are being called
    parser.add_argument("--sample", type = str, required = False, default = "example", help = "Unique string to identify this sample.")

    # Optional argument to specify the library structure for this sample
    parser.add_argument("--libStruct", type = Path, required = False, help = "Path to the library structure json for this sample.")

    # Optional argument to specify whether this is a merged workflow
    parser.add_argument("--isMerge", default = False, action = "store_true", help = "If set, workflow is merged, and sample names will be extracted from the group name.")

    args = parser.parse_args()

    sample_specific_file_paths = io.resolve_sample_specific_file_paths(args.STARsolo_out, args.feature_type, args.matrix_type)
    genes = read_genes(sample_specific_file_paths)
    barcodes = read_barcodes(sample_specific_file_paths)
    umi_counts_by_species = count_transcripts(sample_specific_file_paths, genes, barcodes, args.isBarnyard)
    allCells = build_allCells(sample_specific_file_paths, umi_counts_by_species, args.isBarnyard).reindex(list(barcodes['barcode']))

    if (args.libStruct is not None):
        allCells = split_barcodes(args.libStruct, allCells)

    # When this is a merged sample, the sample name is actually a "group" name; in order to get the correct sample name we extract it from the barcode
    if(args.isMerge):
        allCells['sample'] = [barcode.split("_")[1] for barcode in allCells.index.tolist()]
    else:
        allCells['sample'] = args.sample

    metricsDir = Path(".", f"{args.sample}_metrics")
    metricsDir.mkdir(parents = True)

    # Write allCells.csv for this sample
    allCells.to_csv(f"{metricsDir}/{args.sample}_allCells.csv")

if __name__ == "__main__":
    main()
