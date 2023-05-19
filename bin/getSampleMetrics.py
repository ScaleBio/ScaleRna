#!/usr/bin/env python
"""
Get metrics on a per sample basis from STAR solo output
"""
import pandas as pd
import collections
import numpy as np
import os
import itertools as it
import json
import shutil
import argparse
from typing import Dict
from pathlib import Path
from utils.ReportUtil import CellStat, GeneralUtils
from utils.base_logger import logger

def sample_metric(sample: str, demuxJson: dict, referencesPath: Path, libStruct,
                  star_out: Path, featureType: str, starMatrix: str,
                  useSTARthreshold: bool, isBarnyard: bool,
                  cells: int = None, expectedCells: int = None, 
                  topCellPercent: int = None, minCellRatio: int = None,
                  minReads: int = None):
    """
    Generate metrics for a single sample

    Args:
        sample: Name of sample
        demuxJson: Demuxed metrics
        useSTARthreshold: Whether to use star thresholding or not
        isBarnyard: Whether sample is barnyard or not
        star_out: STAR Solo output
        featureType: STAR feature used for quantification (exon, gene, etc.)
        starMatrix: Name of STAR matrix output to use (matrix.mtx or with multimappers)
        cells: Cell thresholding
        libStruct: Path to library json
        expectedCells: Parameter used for cell thresholding calculation
        topCellPercent: Parameter used for cell thresholding calculation
        minCellRatio: Parameter used for cell thresholding calculation
        minReads: Parameter used for cell thresholding calculation

    """
    
    sampleSpecificFilePaths = resolveSampleSpecificFilePaths(star_out, sample, featureType, starMatrix)
    
    validateInput(sampleSpecificFilePaths)
    
    metricsDir = Path(".", f"{sample}_metrics")
    
    detectedGenes = instantiateGenesDf(sampleSpecificFilePaths)
    
    cellBarcodes = instantiateBarcodesDf(sampleSpecificFilePaths)
    
    star_cells = getStarStats(sampleSpecificFilePaths)
    
    if isBarnyard:
        barnCounts = countBarn(sampleSpecificFilePaths, cellBarcodes, detectedGenes)
        allCells, cells = addStatsAndCallCells(sampleSpecificFilePaths, cells, star_cells,
                                               expectedCells, topCellPercent, minCellRatio,
                                               minReads, barnCounts, useSTARthreshold)
    else:
        allCells, cells = addStatsAndCallCells(sampleSpecificFilePaths, cells, star_cells,
                                               expectedCells, topCellPercent, minCellRatio,
                                               minReads, None, useSTARthreshold)
    
    Path(metricsDir, "sample_metrics").mkdir(parents=True)

    if (libStruct is not None):
        barcodesToPlot = splitBarcodes(demuxJson, sample, libStruct, allCells, referencesPath)

    generateFilteredMatrix(sampleSpecificFilePaths, allCells, sample)

    logger.debug(f"Writing barcodesToPlot, allCells to {str(metricsDir.resolve())}")

    with open(f"{metricsDir}/sample_metrics/barcodesToPlot.json", "w") as f:
        json.dump(barcodesToPlot, f)
    
    with open(f"{metricsDir}/sample_metrics/demuxJson.json", "w") as f:
        json.dump(demuxJson, f)

    allCells["sample"] = sample
    allCells.to_csv(f"{metricsDir}/allCells.csv")


def splitBarcodes(demuxJson, sample, libStruct, allCells, referencesPath):
    """
    Function that uses `level` and `alias` values in each element
    of lib.json barcodes object to split barcodes.
    Assumption: Each element in `barcodes` section of lib.json
    defining part of the CB sequence has a `level` and `alias` field.
    See current lib.json under references folder for examples

    Args:
        demuxJson (dict): Demuxed metrics
        sample (str): Sample name
        libStruct (str): Path to library json
        allCells (pd.DataFrame): Information about all cells
        referencesPath (Path): Reference files
    Returns:
        Dictionary of barcodes to plot in report stage
    """

    sampleBarcodeAliases = {}
    for name, bc in demuxJson["samples"][sample]["barcodes"].items():
        sampleBarcodeAliases[bc['sequence']] = name

    libDirJSON = GeneralUtils.readJSON(referencesPath/libStruct, preserveDictOrder=True)
    sampleBarcode = libDirJSON['sample_barcode']
    barcodeInfo = libDirJSON['barcodes']

    scDefiningBarcodes = [
        barcodeMeta for barcodeMeta in barcodeInfo
        if barcodeMeta.get('type', None) != 'library_index'
        and barcodeMeta.get('type', None) != 'umi']

    groupByLevel = it.groupby(scDefiningBarcodes, key=lambda x: x['level'])

    curIndex = 0
    barcodesToPlot = {}

    for level, group in groupByLevel:
        groupAsList = list(group)
        barcodeSize = sum([x['length'] for x in groupAsList])
        barcodeFullName = groupAsList[0]['alias']
        barcodeName = groupAsList[0]['name']
        barcodes = [b[curIndex:barcodeSize+curIndex] for b in allCells.index]
        allCells[barcodeFullName] = barcodes
        curIndex += barcodeSize

        if barcodeName == sampleBarcode and len(sampleBarcodeAliases) > 0:
            allCells[f"{barcodeFullName}_alias"] = allCells[barcodeFullName].apply(lambda x: sampleBarcodeAliases[x])
            barcodesToPlot[level] = {'alias': barcodeFullName,
                                     'possibleValues': list(sampleBarcodeAliases.values()),
                                     'orderAsWellAliases': True}

        else:
            barcode_mapping = {}
            with open(referencesPath/f"{groupAsList[0]['sequences']}") as f:
                for line in f:
                    line = line.rstrip()
                    barcode_mapping[line.split("\t")[0]] = line.split("\t")[1]

            allCells[f'{barcodeFullName}_alias'] = allCells.apply(lambda row: barcode_mapping[row[barcodeFullName]], axis=1)

            barcodesToPlot[level] = {'alias': barcodeFullName,
                                     'possibleValues': list(sampleBarcodeAliases.values()),
                                     'orderAsWellAliases': False}

    return barcodesToPlot


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
    return dict(featuresFn=exprMatrixPrefix / "raw" / 'features.tsv',
                barcodesFn=exprMatrixPrefix / "raw" / 'barcodes.tsv',
                matrixFn=exprMatrixPrefix / "raw" / starMatrix,
                summaryFn=exprMatrixPrefix / "Summary.csv",
                cellReadsFn=exprMatrixPrefix / "CellReads.stats")


def validateInput(sampleSpecificFilePaths):
    """
    Function to validate that input files exist

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
    """
    GeneralUtils.ensurePathsExist(sampleSpecificFilePaths)


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
    cellBarcodes = pd.read_csv(sampleSpecificFilePaths['barcodesFn'], sep='\t', header=None)
    cellBarcodes.columns = ['Barcode']
    cellBarcodes.index = cellBarcodes.index + 1

    return cellBarcodes


def countBarn(sampleSpecificFilePaths, cellBarcodes, detectedGenes):
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
            gene, cell, count = int(split_line[0]), int(split_line[1]), float(split_line[2])
            cid = cell-1
            if count == 0: continue

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
    allCells = pd.DataFrame(cells, index=cellBarcodes.Barcode)
    return allCells


def getStarStats(sampleSpecificFilePaths):
    """
    Read in star stats and return a dictionary with the stats

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file

    Returns:
        Dictionary with the star stats
    """
    stats = {}

    for line in open(sampleSpecificFilePaths['summaryFn']):
        split_line = line.strip().split(',')
        stats[split_line[0]] = split_line[1]
    
    return stats["Estimated Number of Cells"]


def getCellThreshold(cells, allCells, star_cells):
    """
    Calculate and set the UMI threshold for cells based on STAR output
    files

    Args:
        cells (int): Cell thresholding number set by user
        allCells (pd.DataFrame): Dataframe with info on all cells
        star_cells (int): Estimated number of cells computed by STAR

    Returns:
        Cell threshold and cells variable
    """
    if not cells:
        cells = star_cells

    umis = allCells.umis
    k = len(umis) - cells - 1

    cellThreshold = max(20, np.partition(umis, k)[k])

    logger.debug(f"STAR cell threshold {cellThreshold}")

    return cellThreshold, cells


def calculateCustomCellThreshold(allCells, expectedCells, topCellPercent, minCellRatio, minReads):
    """
    Calculate and set UMI threshold for cells based on our
    custom logic

    Args:
        allCells (pd.DataFrame): Dataframe with information on all cells
        expectedCells (int): Cell thresholding parameter
        topCellPercent (int): Cell thresholding parameter
        minCellRatio (int): Cell thresholding parameter
        minReads (int): Cell thresholding parameter

    Returns:
        Cell threshold
    """
    field = "umis"
    allCells = allCells.sort_values(by=field, ascending=False)
    
    if expectedCells == 0:
        expectedCells = (allCells[field] >= minReads).sum()
    
    if expectedCells == 0:
        return minReads
    
    threshold = np.percentile(allCells[field][:expectedCells], (topCellPercent))/minCellRatio

    cellThreshold = max(threshold, minReads)
    
    logger.debug(f"Custom cell threshold {cellThreshold}")
    logger.debug(f"Expected number of cells {expectedCells}")
    logger.debug(f"Top cell UMI value {np.percentile(allCells[field][:expectedCells], topCellPercent)}")
 
    return cellThreshold


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
        
        for i in range(3):
            f_raw_mtx.readline()

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
    

def addStatsAndCallCells(sampleSpecificFilePaths,
                         cells, star_cells, expectedCells,
                         topCellPercent, minCellRatio, minReads,
                         barnyardCounts, useSTARthreshold):
    """
    Function to read in cellstats and initialize dataframe

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
        cells (int): User set threshold
        star_cells (int): Estimated number of cells computed by STAR
        expectedCells (int): Cell thresholding parameter
        topCellPercent (int): Cell thresholding parameter
        minCellRatio (int): Cell thresholding parameter
        minReads (int): Cell thresholding parameter
        barnyardCounts (pd.DataFrame): Dataframe with mouse and human information
        useSTARthreshold (bool): Whether to use star thresholding or not

    Returns:
        allCells, cellThreshold, cells
    """
    cellReadStats = pd.read_csv(sampleSpecificFilePaths['cellReadsFn'], sep='\t', index_col='CB')
    cellReadStats = cellReadStats[1:]
    allCells = pd.DataFrame(index=cellReadStats.index.to_list())
    
    if barnyardCounts is None:
        # If this is not a barnyard, we are taking gene and UMI counts from the STAR Cell stats
        # otherwise we computed from from the .mtx earlier
        allCells['genes'] = cellReadStats.nGenesUnique + cellReadStats.nGenesMulti
        allCells['umis'] = cellReadStats.nUMIunique + cellReadStats.nUMImulti
    else:
        # For a barnyard, we only consider genes from the main species per cell-barcode
        # (not the background from the other species)
        allCells['humanGenes'] = barnyardCounts["humanGenes"]
        allCells['mouseGenes'] = barnyardCounts["mouseGenes"]
        allCells['humanUmis'] = barnyardCounts["humanUmis"]
        allCells['mouseUmis'] = barnyardCounts["mouseUmis"]
        allCells['genes'] = barnyardCounts[["humanGenes", "mouseGenes"]].max(1)
        allCells['umis'] = barnyardCounts[["humanUmis", "mouseUmis"]].max(1)

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


    allCells.sort_values('umis', ascending=False, inplace=True)

    if useSTARthreshold or cells:
        cellThreshold, cells = getCellThreshold(cells, allCells, star_cells)
    else:
        cellThreshold = calculateCustomCellThreshold(allCells, expectedCells, topCellPercent, minCellRatio, minReads)
    
    allCells['pass'] = allCells.umis >= cellThreshold

    return allCells, cells

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", type=str,
                        help="The name of the sample for which a report is being generated")
    parser.add_argument("--samplesheet", type=str,
                        help="Path to the samples.csv containing information about libName")
    parser.add_argument("--libJsonName", type=str)
    parser.add_argument("--starFeature", default="GeneFull_Ex50pAS", type=str,
                        help="STARSolo feature type used")
    parser.add_argument("--starMatrix", default="matrix.mtx", type=str,
                        help="STARSolo gene expression matrix file to use")
    parser.add_argument("--cells", default=None, type=int,
                        help="Set a fixed number of cells to pass for UMI count threshold")
    parser.add_argument("--star_out", type=str,
                        help="STAR solo output file", required=True)
    parser.add_argument("--libraryStructPath", type=str,
                        help="Path to folder containing library files")
    parser.add_argument("--cellThreshold", type=int)
    parser.add_argument("--topCellPercent", type=int)
    parser.add_argument("--minCellRatio", type=int)
    parser.add_argument("--minReads", type=int)
    parser.add_argument("--isBarnyard", default=False, action="store_true")
    parser.add_argument("--useSTARthreshold", action="store_true")
    
    args = parser.parse_args()

    demuxJson = GeneralUtils.readJSON("demuxMetrics.json")

    sampleSheetDf = pd.read_csv(args.samplesheet)

    if 'expectedCells' in sampleSheetDf.columns:
        expectedCells = int(sampleSheetDf.loc[sampleSheetDf['sample'] == args.sample]['expectedCells'])
    else:
        expectedCells = 0
    
    sample_metric(args.sample, demuxJson, Path(args.libraryStructPath), args.libJsonName,
                  Path(args.star_out), args.starFeature, args.starMatrix,
                  args.useSTARthreshold, args.isBarnyard, args.cells, expectedCells,
                  args.topCellPercent, args.minCellRatio, args.minReads)

if __name__ == "__main__":
    main()
