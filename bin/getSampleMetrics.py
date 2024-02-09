#!/usr/bin/env python
"""
Compute cell-threshold and metrics for one sample from STAR solo output
"""
import collections
from dataclasses import dataclass
import os
import itertools as it
import json
import shutil
import argparse
import sys
from typing import Dict
from pathlib import Path

import pandas as pd
import numpy as np

from scaleReportUtils.base_logger import logger
from scaleReportUtils import fileUtils, statsUtils

@dataclass
class CellStat:
    """
    Total UMI and gene-counts per cell
    """
    genes: int = 0 # Total number of genes detected
    humanGenes: int = 0 # Human genes (for barnyard only)
    mouseGenes: int = 0 # mouse genes "
    umis: int = 0 # Number of UMIs / transcripts detected
    humanUmis: int = 0 # Per species (for barnyard)
    mouseUmis: int = 0

def sample_metric(sample: str, libStruct: Path,
                  star_out: Path, featureType: str, starMatrix: str, isBarnyard: bool,
                  minReads: int, calcCellThres: bool, cells: int|None, 
                  expectedCells: int|None, topCellPercent: int|None, minCellRatio: int|None, 
                  isMerge: bool):
    """
    Generate metrics for a single sample

    Args:
        sample: Name of sample
        libStruct: Path to library json

        star_out: STAR Solo output
        featureType: STAR feature used for quantification (exon, gene, etc.)
        starMatrix: Name of STAR matrix output to use (matrix.mtx or with multimappers)
        isBarnyard: Whether sample is barnyard or not
        minReads: Parameter used for cell thresholding calculation
        calcCellThres: Use custom cell-threshold (or the one from STAR)?
        cells: Cell thresholding
        expectedCells: Parameter used for cell thresholding calculation
        topCellPercent: Parameter used for cell thresholding calculation
        minCellRatio: Parameter calcCellThresused for cell thresholding calculation

    """
    
    sampleSpecificFilePaths = resolveSampleSpecificFilePaths(star_out, sample, featureType, starMatrix)
    detectedGenes = instantiateGenesDf(sampleSpecificFilePaths)
    cellBarcodes = instantiateBarcodesDf(sampleSpecificFilePaths)
    star_stats = getStarStats(sampleSpecificFilePaths['summaryFn']) if not calcCellThres else None
    transcriptCounts = countTranscripts(sampleSpecificFilePaths, cellBarcodes, detectedGenes, isBarnyard)
    allCells = cellStats(sampleSpecificFilePaths, transcriptCounts, isBarnyard)
    callCells(allCells, star_stats, minReads, cells, expectedCells, topCellPercent, minCellRatio)

    if (libStruct is not None):
        splitBarcodes(sample, libStruct, allCells)

    generateFilteredMatrix(sampleSpecificFilePaths, allCells, sample)
    
    #When this is a merged sample the variable sample is actually a "group" name.
    #In order to get the correct sample name we extract it from the barcode.
    if(isMerge):
        allCells["sample"] = [bc.split("_")[1] for bc in allCells.index.tolist()]
    else:
        allCells["sample"] = sample

    metricsDir = Path(".", f"{sample}_metrics")
    metricsDir.mkdir(parents=True)
    allCells.to_csv(f"{metricsDir}/allCells.csv")


def splitBarcodes(sample:str, libStruct:Path, allCells:pd.DataFrame):
    """
    Function that uses `level` and `alias` values in each element
    of lib.json barcodes object to split barcodes.
    Assumption: Each element in `barcodes` section of lib.json
    defining part of the CB sequence has a `level` and `alias` field.
    See current lib.json under references folder for examples

    Results are added to the allCells dataframe

    Args:
        sample: Sample name
        libStruct: Path to library structure definition json
        allCells: Information about all cells
    """
    libDir = libStruct.parent # Directory container lib structure Json and sequence files
    libJson = fileUtils.readJSON(libStruct, preserveDictOrder=True)
    barcodeInfo = libJson['barcodes']

    cellBarcodesByLevel = [
        barcodeMeta for barcodeMeta in barcodeInfo
        if barcodeMeta.get('type', None) != 'library_index'
        and barcodeMeta.get('type', None) != 'umi']
    # We need to keep the barcodes in the original order from lib.json so that
    # it matches the order in which they were concatinated by bcParser (cell barcode / index)
    # Hence we rely on barcode in lib.json being grouped by level to make the group-by call work
    # cellBarcodesByLevel.sort(key=lambda bc: bc['level'])

    # Overall cell-barcode is the concatination of all individual cell-barcodes
    curBcSeqStartIndex = 0 
    for level, bcParts in it.groupby(cellBarcodesByLevel, key=lambda x: x['level']):
        # Single barcode level might have multiple sub-barcodes (e.g. bead barcode blocks)
        # These will be concatinated and reported as one barcode
        bcParts = list(bcParts)
        bcInfo = bcParts[0]
        barcodeFullName = bcInfo['alias']
        barcodeLen = sum([x['length'] for x in bcParts])
        bcSubSeq = slice(curBcSeqStartIndex, barcodeLen+curBcSeqStartIndex)
        barcodes = [b[bcSubSeq] for b in allCells.index]
        allCells[barcodeFullName] = barcodes
        curBcSeqStartIndex += barcodeLen
        # Load the ids for specific barcode sequences (e.g. well coordinates) for the current barcode level
        bcIds = {}
        with open(libDir/f"{bcInfo['sequences']}") as f:
            for line in f:
                splat = line.rstrip().split('\t')
                bcIds[splat[0]] = splat[1]

        bcSeqIds = allCells[barcodeFullName].apply(lambda seq: bcIds[seq])
        allCells[f'{barcodeFullName}_alias'] = bcSeqIds


def resolveSampleSpecificFilePaths(star_out: Path, sample: str, featureType: str, starMatrix: str) -> Dict[str, Path]:
    """
    Return a dictionary of needed file paths

    Args:
        star_out: Path to STAR solo output directory
        sample: Sample name
        featureType: STARSolo feature type used
        starMatrix: STARSolo mtx file to use (unique or multimappers)
    Returns:
        Dictionary where key is custom file identifier and value is
        path to file
    """
    exprMatrixPrefix = star_out / featureType
    files = dict(featuresFn=exprMatrixPrefix / "raw" / 'features.tsv',
                barcodesFn=exprMatrixPrefix / "raw" / 'barcodes.tsv',
                matrixFn=exprMatrixPrefix / "raw" / starMatrix,
                summaryFn=exprMatrixPrefix / "Summary.csv",
                cellReadsFn=exprMatrixPrefix / "CellReads.stats")
    fileUtils.ensurePathsExist(files)
    return files


def instantiateGenesDf(sampleSpecificFilePaths):
    """
    Read in sample features.tsv from cell x gene count
    matrix as Dataframe and label each gene as belonging to
    either mm (mouse) or hs (human)

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file

    Returns:
        Dataframe with detected genes
    """
    detectedGenes = pd.read_csv(sampleSpecificFilePaths['featuresFn'], sep="\t", header=None)
    detectedGenes.columns = ["Id", "Name", "Type"]
    detectedGenes.index = detectedGenes.index + 1
    detectedGenes['species'] = detectedGenes.Id.map(lambda x: ["mm", "hs"][x.startswith("ENSG")])

    return detectedGenes


def instantiateBarcodesDf(sampleSpecificFilePaths):
    """
    Read in sample barcodes.tsv from cell x gene count
    matrix as Dataframe

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
    Returns:
        Dataframe with cell barcodes
    """
    cellBarcodes = pd.read_csv(sampleSpecificFilePaths['barcodesFn'], sep='\t', header=None, names=['Barcode'])
    cellBarcodes.index = cellBarcodes.index + 1

    return cellBarcodes


def countTranscripts(sampleSpecificFilePaths, cellBarcodes, detectedGenes, isBarnyard):
    """
    Read in statistics from raw matrix.mtx to compute barnyard specific metrics

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
        cellBarcodes (pd.DataFrame): Dataframe with cell barcodes
        detectedGenes (pd.DataFrame): Dataframe with detected genes

    Returns:
        Dataframe with human and mouse counts
    """
    cells = [CellStat() for bc in cellBarcodes.Barcode]
    with open(sampleSpecificFilePaths['matrixFn'], "r") as matrixFile:
        # Read through header
        for i in range(3):
            matrixFile.readline()
            
        for line in matrixFile:
            split_line = line.split()
            # Individual expression counts can be fractions, due to multimapper resolution
            gene, cell, count = int(split_line[0]), int(split_line[1]), float(split_line[2])
            cid = cell-1
                
            if count == 0: continue
                
            if isBarnyard:
                
                if detectedGenes.species[gene] == "hs":
                    cells[cid].genes += 1
                    cells[cid].humanGenes += 1
                    cells[cid].humanUmis += count
                    cells[cid].umis += count
                    
                elif detectedGenes.species[gene] == "mm":
                    cells[cid].genes += 1
                    cells[cid].mouseGenes += 1
                    cells[cid].mouseUmis += count
                    cells[cid].umis += count
                    
                else:
                    print(detectedGenes.iloc[gene])
                    
            else:
                cells[cid].genes += 1
                cells[cid].umis += count
    # For consistency with STAR CellReadMetrics output, we round the total expression (UMI) count
    for cell in cells:
        cell.umis = round(cell.umis)
    allCells = pd.DataFrame(cells, index=cellBarcodes.Barcode)
    return allCells


def getStarStats(starSummaryFn:Path) -> Dict[str, float | str]:
    """
    Read in star stats and return a dictionary with the stats

    Args:
        sampleSpecificFilePaths: Input f
            a STAR output file name and value is path to the file

    Returns:
        StatName -> Value (numbers converted to float)
    """
    stats = {}
    for line in open(starSummaryFn):
        split_line = line.strip().split(',')
        val = split_line[1]
        try:
            val = float(val)
        except ValueError:
            pass # Not a number
        stats[split_line[0]] = val
    return stats

def calculateCustomCellThreshold(allCells:pd.DataFrame, expectedCells:int, topCellPercent:float, minCellRatio:float, minReads:int) -> int:
    """
    Calculate and set UMI threshold for cells based on our
    custom logic

    Args:
        allCells: Dataframe with per-cell information
        expectedCells: Cell thresholding parameter
        topCellPercent: Cell thresholding parameter
        minCellRatio: Cell thresholding parameter
        minReads: Cell thresholding parameter

    Returns:
        Cell threshold
    """
    field = "umis"
    allCells = allCells.sort_values(by=field, ascending=False)
    
    if expectedCells == 0:
        expectedCells = (allCells[field] >= minReads).sum()
    if expectedCells == 0:
        return minReads
    
    topCellCount = np.percentile(allCells[field][:expectedCells], (topCellPercent)) 
    threshold = max(topCellCount/minCellRatio, minReads)
    
    logger.debug(f"Custom cell threshold {threshold}")
    logger.debug(f"Expected number of cells {expectedCells}")
    logger.debug(f"Top cell UMI value {topCellCount}")
    return threshold


def generateFilteredMatrix(sampleSpecificFilePaths, allCells, sample):
    """
    Generate filtered cell gene expression matrix from STAR solo raw output

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
        allCells (pd.DataFrame): Information about all cells
        sample (str): Sample name
    """
    subset_allCells = allCells["pass"]
    subset_allCells_dict = subset_allCells.to_dict()
    filtered_path = f"{sample}_filtered_star_output"
    logger.debug(f"Writing filtered matrix to {filtered_path}")
    # Create filtered directory under metrics directory
    os.mkdir(filtered_path)

    shutil.copyfile(sampleSpecificFilePaths['featuresFn'],
                    f"{filtered_path}/features.tsv")

    # File pointer to raw matrix.mtx and barcodes.tsv
    f_raw_mtx = open(sampleSpecificFilePaths['matrixFn'])
    f_raw_barcodes = open(sampleSpecificFilePaths['barcodesFn'])

    # File pointer to newly created, filtered matrix.mtx, barcodes.tsv
    # and features.tsv
    with (open("tmp_matrix.mtx", "w") as f_filtered_mtx_tmp,
         open(f"{filtered_path}/barcodes.tsv", "w") as f_filtered_barcodes):
        # Read in barcodes from raw barcodes file
        raw_barcodes = f_raw_barcodes.readlines()
        raw_barcodes = [line.rstrip() for line in raw_barcodes]

        header = f_raw_mtx.readline().strip() # datatype header
        f_raw_mtx.readline(); f_raw_mtx.readline() # Skip extra header lines

        # Set current barcode to 0 to ensure that the first time we
        # compare current_barcode to the barcode received from
        # iterating over the raw matrix, we set current_barcode to
        # that barcode
        current_barcode = 0

        # Row number in barcodes.tsv
        barcode_row_number = 1

        # Count number of lines for writing as header later on
        line_count = 0

        # Iterate over every line in raw matrix.mtx
        for line in f_raw_mtx:
            split_line = line.split()
            feature, barcode, count = int(split_line[0]), int(split_line[1]), float(split_line[2])

            barcode_sequence = raw_barcodes[barcode-1]

            # If barcode has passing cells in the allCells dataframe, write to file
            if subset_allCells_dict[barcode_sequence]:
                if barcode != current_barcode:
                    current_barcode = barcode
                    barcode_row_number += 1
                    f_filtered_barcodes.write(f"{raw_barcodes[current_barcode-1]}\n")
                
                f_filtered_mtx_tmp.write(f"{feature} {barcode_row_number-1} {count}\n")
                line_count += 1
    f_raw_mtx.close()
    f_raw_barcodes.close()
    f_filtered_barcodes.close()

    # Compute header information for the matrix
    # First entry is length of filtered features.tsv
    # Second entry is length of filtered barcodes.tsv
    # Third entry is length of filtered matrix.mtx
    header1 = len(pd.read_csv(sampleSpecificFilePaths['featuresFn'], sep="\t", header=None).index)
    header2 = barcode_row_number - 1
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
        print(e, file=sys.stderr)
    

def cellStats(sampleSpecificFilePaths:Dict[str, Path], barnyardCounts:pd.DataFrame, isBarnyard):
    """
    Function to read in cellstats and initialize cell metrics dataframe

    Args:
        sampleSpecificFilePaths: STAR output files indexed by name
        barnyardCounts: Dataframe with per-species counts (for barnyard exp. only)

    Returns:
        dataframe with metrics per cell
    """
    cellReadStats = pd.read_csv(sampleSpecificFilePaths['cellReadsFn'], sep='\t', index_col='CB')
    cellReadStats = cellReadStats[1:]
    allCells = pd.DataFrame(index=cellReadStats.index.to_list())

    if isBarnyard:
        # For a barnyard, we only consider genes from the main species per cell-barcode
        # (not the background from the other species)
        allCells['humanGenes'] = barnyardCounts["humanGenes"]
        allCells['mouseGenes'] = barnyardCounts["mouseGenes"]
        allCells['humanUmis'] = barnyardCounts["humanUmis"]
        allCells['mouseUmis'] = barnyardCounts["mouseUmis"]
        allCells['genes'] = barnyardCounts[["humanGenes", "mouseGenes"]].max(1)
        allCells['umis'] = barnyardCounts[["humanUmis", "mouseUmis"]].max(1)
        
    else:
        # If this is not a barnyard we use values computed from from the .mtx earlier
        allCells['genes'] = barnyardCounts["genes"]
        allCells['umis'] = barnyardCounts["umis"]

    allCells['reads'] = cellReadStats.cbMatch
    allCells['mappedReads'] = cellReadStats.genomeU + cellReadStats.genomeM
    allCells['geneReads'] = cellReadStats.featureU + cellReadStats.featureM
    allCells['exonReads'] = cellReadStats.exonic
    allCells['intronReads'] = cellReadStats.intronic
    allCells['antisenseReads'] = cellReadStats.exonicAS
    allCells['passingReads'] = cellReadStats.countedU + cellReadStats.countedM
    allCells['uniquePassingReads'] = cellReadStats.countedU
    allCells['Saturation'] = 1 - (cellReadStats.nUMIunique / cellReadStats.countedU)
    allCells['mitoReads'] = cellReadStats.mito
    return allCells


def callCells(allCells:pd.DataFrame, starStats:Dict|None, minReads:int, cells:int|None, expectedCells:int, topCellPercent:float|None, 
              minCellRatio:float|None):
    """
    Calculate the UMI threshold for passing cells and add 'pass' value to cells dataframe

    Args:
        allCells: Per-cell metrics
        starStats: STAR summary metrics (if using STAR threshold)
        minReads: Minimum number of reads for any cell to possibly 'pass'
        cells: Pre-defined number of cells to pass (or 0 for dynamic threshold)
        expectedCells: Rough estimate of cells in this dataset (or 0 to use all cells >= minReads
        topCellPercent: Cell thresholding parameter
        minCellRatio: Cell thresholding parameter
    """
    
    allCells.sort_values('umis', ascending=False, inplace=True)
    cellThreshold = minReads
    if starStats and not cells:
        cells = int(starStats["Estimated Number of Cells"])
        logger.debug(f"Passing cells from STAR {cells}")
    if cells:
        k = len(allCells.umis) - cells - 1
        cellThreshold = max(cellThreshold, np.partition(allCells.umis, k)[k])
        logger.debug(f"Fixed cell threshold {cellThreshold}")
    else:
        if topCellPercent is None or minCellRatio is None:
            raise ValueError("Need cell threshold parameters to compute a dynamic threshold")
        cellThreshold = calculateCustomCellThreshold(allCells, expectedCells, topCellPercent, minCellRatio, minReads)
        logger.debug(f"Dynamic cell threshold {cellThreshold}")

    allCells['pass'] = allCells.umis >= cellThreshold
    logger.debug(f"Passing cells {allCells['pass'].sum()}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", type=str,
                        help="The name of the sample for which a report is being generated")
    parser.add_argument("--libStruct", type=Path)
    parser.add_argument("--starFeature", default="GeneFull_Ex50pAS", type=str,
                        help="STARSolo feature type used")
    parser.add_argument("--starMatrix", default="matrix.mtx", type=str,
                        help="STARSolo gene expression matrix file to use")
    parser.add_argument("--star_out", type=Path,
                        help="STAR solo output file", required=True)

    parser.add_argument("--minReads", type=int, default=100,
                        help="Minimum number of UMIs to consider a cell(-barcode)")
    parser.add_argument("--cells", type=int,
                        help="Set a fixed number of cells to pass for UMI count threshold")
    parser.add_argument("--calcCellThreshold", action="store_true",
                        help="Compute a heuristic cell threshold based on parameters below")
    parser.add_argument("--topCellPercent", type=int, help="Cell threshold parameter")
    parser.add_argument("--minCellRatio", type=int, help="Cell threshold parameter")
    parser.add_argument("--expectedCells", type=int, help="Number of cells expected in data")

    parser.add_argument("--isBarnyard", default=False, action="store_true")
    parser.add_argument("--isMerge", default=False, action="store_true")
    
    args = parser.parse_args()

    sample_metric(args.sample, args.libStruct,
                  args.star_out, args.starFeature, args.starMatrix,
                  args.isBarnyard, args.minReads, args.calcCellThreshold, args.cells, args.expectedCells,
                  args.topCellPercent, args.minCellRatio, args.isMerge)

if __name__ == "__main__":
    main()
