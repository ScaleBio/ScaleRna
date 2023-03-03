#!/usr/bin/env python
"""
Python script to get metrics on a per sample basis
from STAR Solo output
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

def sample_metric(sample: str, demuxJson: dict,
                  useSTARthreshold: bool, isBarnyard: bool, star_out: Path,
                  cells: int = None, libDir: str = None,
                  featureType: str = "GeneFull_Ex50pAS",
                  expectedCells: int = None, 
                  topCellPercent: int = None, minCellRatio: int = None,
                  minReads: int = None):
    """
    Function that generates metrics for a sample

    Args:
        sample (str): Name of sample
        resultsDir (str): Results directory
        demuxJson (dict): Demuxed metrics
        useSTARthreshold (bool): Whether to use star thresholding or not
        isBarnyard (bool): Whether sample is barnyard or not
        cells (int): Cell thresholding
        libDir (str): Path to library json
        featureType (str): STAR feature
        expectedCells (int): Parameter used for cell thresholding calculation
        topCellPercent (int): Parameter used for cell thresholding calculation
        minCellRatio (int): Parameter used for cell thresholding calculation
        minReads (int): Parameter used for cell thresholding calculation

    """
    
    sampleSpecificFilePaths = resolveSampleSpecificFilePaths(star_out, sample, featureType)
    
    validateInput(sampleSpecificFilePaths)
    
    metricsDir = Path(".", f"{sample}_metrics")
    
    detectedGenes = instantiateGenesDf(sampleSpecificFilePaths)
    
    cellBarcodes = instantiateBarcodesDf(sampleSpecificFilePaths)
    
    starStats = getStarStats(sampleSpecificFilePaths, featureType)
    
    if isBarnyard:
        barnCounts = countBarn(sampleSpecificFilePaths, cellBarcodes, detectedGenes)
        allCells, cellThreshold, cells = addStatsAndCallCells(sampleSpecificFilePaths, cells, starStats,
                                                              expectedCells, topCellPercent, minCellRatio,
                                                              minReads, barnCounts, useSTARthreshold)
    else:
        allCells, cellThreshold, cells = addStatsAndCallCells(sampleSpecificFilePaths, cells, starStats,
                                                              expectedCells, topCellPercent, minCellRatio,
                                                              minReads, None, useSTARthreshold)
    try:
        os.mkdir(metricsDir)
        os.mkdir(metricsDir / "sample_metrics")
    except OSError as error:
        logger.error(error)
        logger.debug(f"Directory {metricsDir} exists with contents {os.listdir(f'{sample}_metrics')}")

    if (libDir is not None):
        barcodesToPlot = splitBarcodes(
            demuxJson, sample, libDir, allCells)

    generateFilteredMetrics(sampleSpecificFilePaths, metricsDir, cellThreshold, allCells, sample)

    logger.debug("Writing barcodesToPlot, allCells to "
                 f"{str(metricsDir.resolve())}")

    with open(f"{metricsDir}/sample_metrics/barcodesToPlot.json", "w") as f:
        json.dump(barcodesToPlot, f)
    
    with open(f"{metricsDir}/sample_metrics/demuxJson.json", "w") as f:
        json.dump(demuxJson, f)

    allCells["sample"] = sample
    allCells.to_csv(f"{metricsDir}/allCells.csv")


def splitBarcodes(demuxJson, sample, libDir, allCells):
    """
    Function that uses `level` and `alias` values in each element
    of lib.json barcodes object to split barcodes.
    Assumption: Each element in `barcodes` section of lib.json
    defining part of the CB sequence has a `level` and `alias` field.
    See current lib.json under references folder for examples

    Args:
        demuxJson (dict): Demuxed metrics
        sample (str): Sample name
        libDir (str): Path to library json
        allCells (pd.DataFrame): Information about all cells
    Returns:
    """

    sampleBarcodeAliases = {}
    for name, bc in demuxJson[
            "samples"][sample]["barcodes"].items():
        sampleBarcodeAliases[bc['sequence']] = name

    libDirJSON = GeneralUtils.readJSON(libDir, preserveDictOrder=True)
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
        barcodes = [b[
            curIndex:barcodeSize+curIndex] for b in allCells.index]
        allCells[barcodeFullName] = barcodes
        curIndex += barcodeSize

        if barcodeName == sampleBarcode and len(sampleBarcodeAliases) > 0:
            allCells[f"{barcodeFullName}_alias"] = \
                allCells[barcodeFullName].apply(
                    lambda x: sampleBarcodeAliases[x])
            barcodesToPlot[level] = \
                {'alias': barcodeFullName,
                    'possibleValues': list(sampleBarcodeAliases.values()),
                    'orderAsWellAliases': True}

        else:
            barcode_mapping = {}
            with open(f"references/"
                      f"{groupAsList[0]['sequences']}") as f:
                for line in f:
                    line = line.rstrip()
                    barcode_mapping[line.split("\t")[0]] = \
                        line.split("\t")[1]

            allCells[f'{barcodeFullName}_alias'] = \
                allCells.apply(
                    lambda row: barcode_mapping[row[barcodeFullName]],
                    axis=1)

            barcodesToPlot[level] = \
                {'alias': barcodeFullName,
                    'possibleValues': list(sampleBarcodeAliases.values()),
                    'orderAsWellAliases': False}

    return barcodesToPlot


def resolveSampleSpecificFilePaths(star_out, sample, featureType) -> Dict:
    """
    Function to return a dictionary of needed file paths

    Args:
        star_out (str): Path to STAR solo output directory
        sample (str): Sample name
        featureType (str): STAR feature type

    Returns:
        Dictionary where key is custom file identifier and value is
        path to file
    """
    outputDict = {}
    
    exprMatrixPrefix = star_out / featureType

    outputDict = dict(featuresFn=exprMatrixPrefix / "raw" / 'features.tsv',
                      barcodesFn=exprMatrixPrefix / "raw" / 'barcodes.tsv',
                      matrixFn=exprMatrixPrefix / "raw" / 'matrix.mtx',
                      summaryFn=exprMatrixPrefix / "Summary.csv",
                      cellReadsFn=exprMatrixPrefix / "CellReads.stats")
    return outputDict


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
    Function to read in sample features.tsv from cell x gene count
    matrix as Dataframe and label each gene as belonging to
    either mm (mouse) or hs (human)

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file

    Returns:
        Dataframe with detected genes
    """
    detectedGenes = pd.read_csv(
        sampleSpecificFilePaths['featuresFn'], sep="\t", header=None)
    detectedGenes.columns = ["Id", "Name", "Type"]
    detectedGenes.index = detectedGenes.index + 1
    detectedGenes['species'] = \
        detectedGenes.Id.map(lambda x: ["mm", "hs"][x.startswith("ENSG")])

    return detectedGenes


def instantiateBarcodesDf(sampleSpecificFilePaths):
    """
    Function to read in sample barcodes.tsv from cell x gene count
    matrix as Dataframe

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
    Returns:
        Dataframe with cell barcodes
    """
    cellBarcodes = pd.read_csv(
        sampleSpecificFilePaths['barcodesFn'], sep='\t', header=None)
    cellBarcodes.columns = ['Barcode']
    cellBarcodes.index = cellBarcodes.index + 1

    return cellBarcodes


def countBarn(sampleSpecificFilePaths, cellBarcodes, detectedGenes):
    """
    Function to read in statistics from raw matrix.mtx

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
        cellBarcodes (pd.DataFrame): Dataframe with cell barcodes
        detectedGenes (pd.DataFrame): Dataframe with detected genes

    Returns:
        Dataframe with human and mouse counts
    """
    cells = collections.defaultdict(CellStat)
    with open(sampleSpecificFilePaths['matrixFn'], "r") as matrixFile:
        # Read through header
        for i in range(3):
            matrixFile.readline()

        for entry in matrixFile:
            line = entry.split()
            lineParsedAsInt = [int(x) for x in line]
            gene, cell, count, *_ = lineParsedAsInt

            if detectedGenes.species[gene] == "hs":
                cells[cell].genes += 1
                cells[cell].humanGenes += 1
                cells[cell].humanUmis += count
                cells[cell].umis += count

            elif detectedGenes.species[gene] == "mm":
                cells[cell].genes += 1
                cells[cell].mouseGenes += 1
                cells[cell].mouseUmis += count
                cells[cell].umis += count

            else:
                print(detectedGenes.iloc[gene])
    allCells = pd.DataFrame(
        cells.values(),
        index=[cellBarcodes.Barcode[i] for i in cells])

    return allCells


def getStarStats(sampleSpecificFilePaths, featureType):
    """
    Function to read in star stats and return a dictionary with the stats

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
        featureType (str): STAR feature type

    Returns:
        Dictionary with the star stats
    """
    return GeneralUtils.loadStarStatsFromCsv(
        sampleSpecificFilePaths['summaryFn'], featureType)


def getCellThreshold(cells, allCells, starStats):
    """
    Calculate and set the UMI threshold for cells based on STAR output
    files

    Args:
        cells (int): Cell thresholding number set by user
        allCells (pd.DataFrame): Dataframe with info on all cells
        starStats (dict): Dictionary with information from summary.csv
            of STAR output

    Returns:
        Cell threshold and cells variable
    """
    if not cells:
        cells = int(starStats['STAR Cells'])

    umis = allCells.umis
    k = len(umis) - cells - 1

    cellThreshold = max(20, np.partition(umis, k)[k])

    logger.debug(f"STAR cell threshold {cellThreshold}")

    return cellThreshold, cells


def calculateCustomCellThreshold(allCells, expectedCells,
                                 topCellPercent, minCellRatio,
                                 minReads):
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
    threshold = np.percentile(allCells[field][:expectedCells],
                              (topCellPercent))/minCellRatio

    cellThreshold = max(threshold, minReads)
    logger.debug(f"Custom cell threshold {cellThreshold}")
    logger.debug(f"Expected number of cells {expectedCells}")
    logger.debug(f"Top cell UMI value {np.percentile(allCells[field][:expectedCells], topCellPercent)}")
 
    return cellThreshold


def generateFilteredMetrics(sampleSpecificFilePaths, metricsDir,
                            cellThreshold, allCells, sample):
    """
    Function to generate filtered metrics from STAR solo raw output files

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
        metricsDir (path): Folder to write metrics to
        cellThreshold (int): Cell threshold
    """
    subset_allCells = allCells["pass"]
    subset_allCells_dict = subset_allCells.to_dict()

    logger.debug("Writing filtered matrix to "
                 f"{sample}_filtered_star_output")

    filtered_path = f"{sample}_filtered_star_output"
    # Create filtered directory under metrics directory
    os.mkdir(filtered_path)

    shutil.copyfile(sampleSpecificFilePaths['featuresFn'],
                    f"{filtered_path}/features.tsv")

    # File pointer to raw matrix.mtx and barcodes.tsv
    f_raw_mtx = open(sampleSpecificFilePaths['matrixFn'])
    f_raw_barcodes = open(sampleSpecificFilePaths['barcodesFn'])

    # File pointer to newly created, filtered matrix.mtx, barcodes.tsv
    # and features.tsv
    f_filtered_mtx_tmp = open("tmp_matrix.mtx", "w")
    f_filtered_barcodes = open(
        f"{filtered_path}/barcodes.tsv", "w")

    # Read in barcodes from raw barcodes file
    raw_barcodes = f_raw_barcodes.readlines()
    raw_barcodes = [line.rstrip() for line in raw_barcodes]

    # Skip first three lines of file
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

        # Parse line as int
        lineParsedAsInt = [int(x) for x in split_line]
        feature, barcode, umi, *_ = lineParsedAsInt

        barcode_sequence = raw_barcodes[barcode-1]

        # If barcode has passing cells in the allCells dataframe, write to file
        if subset_allCells_dict[barcode_sequence]:
            if barcode != current_barcode:
                current_barcode = barcode
                barcode_row_number += 1
                f_filtered_barcodes.write(f"{raw_barcodes[current_barcode-1]}\n")
            
            f_filtered_mtx_tmp.write(f"{feature} {barcode_row_number-1} {umi}\n")
            line_count += 1

    # Compute header information for the matrix
    # First entry is length of filtered features.tsv
    # Second entry is length of filtered barcodes.tsv
    # Third entry is length of filtered matrix.mtx
    header1 = len(pd.read_csv(sampleSpecificFilePaths['featuresFn'], sep="\t", header=None).index)
    header2 = barcode_row_number - 1
    header3 = line_count

    if line_count == 0:
        f_raw_mtx.close()
        f_raw_barcodes.close()
        f_filtered_barcodes.close()
        f_filtered_mtx_tmp.close()
        with open(f"{filtered_path}/matrix.mtx", "w") as f_filtered_mtx:
            f_filtered_mtx.write("%%MatrixMarket matrix coordinate integer general\n%\n")
            f_filtered_mtx.write(f"{header1} {header2} {header3}\n")
        return

    f_filtered_mtx_tmp.close()
    with open(f"{filtered_path}/matrix.mtx", "w") as f_filtered_mtx:
        f_filtered_mtx.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        f_filtered_mtx.write(f"{header1} {header2} {header3}\n")

    os.system(f"cat tmp_matrix.mtx >> {filtered_path}/matrix.mtx")
    
    f_raw_mtx.close()
    f_raw_barcodes.close()
    f_filtered_barcodes.close()


def addStatsAndCallCells(sampleSpecificFilePaths,
                         cells, starStats, expectedCells,
                         topCellPercent, minCellRatio, minReads,
                         barnyardCounts, useSTARthreshold):
    """
    Function to read in cellstats and initialize dataframe

    Args:
        sampleSpecificFilePaths (dict): Dictionary where key is placeholder for
            a STAR output file name and value is path to the file
        cells (int): User set threshold
        starStats (dict): Dict with summary.csv information
        expectedCells (int): Cell thresholding parameter
        topCellPercent (int): Cell thresholding parameter
        minCellRatio (int): Cell thresholding parameter
        minReads (int): Cell thresholding parameter
        barnyardCounts (pd.DataFrame): Dataframe with mouse and human information
        useSTARthreshold (bool): Whether to use star thresholding or not

    Returns:
        allCells, cellThreshold, cells
    """
    cellReadStats = pd.read_csv(
        sampleSpecificFilePaths['cellReadsFn'],
        sep='\t', index_col='CB')
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
        cellThreshold, cells = getCellThreshold(cells, allCells, starStats)
    else:
        cellThreshold = calculateCustomCellThreshold(
            allCells, expectedCells,
            topCellPercent, minCellRatio,
            minReads)
    allCells['pass'] = allCells.umis >= cellThreshold

    return allCells, cellThreshold, cells

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample", type=str,
        help="The name of the sample for which a report is being generated")
    parser.add_argument(
        "--samplesheet", type=str,
        help="Path to the samples.csv containing information about libName")
    parser.add_argument(
        "--libStruct", type=str)
    parser.add_argument(
        "--featureType", default="GeneFull_Ex50pAS",
        type=str, help="STARSolo feature type used")
    parser.add_argument(
        "--cells", default=None, type=int,
        help="Set a fixed number of cells to pass for UMI count threshold")
    parser.add_argument(
        "--star_out", type=str,
        help="STAR solo output file",
        required=True)
    parser.add_argument(
        "--cellThreshold", type=int)
    parser.add_argument(
        "--topCellPercent", type=int)
    parser.add_argument(
        "--minCellRatio", type=int)
    parser.add_argument(
        "--minReads", type=int)
    parser.add_argument(
        "--isBarnyard", default=False, action="store_true")
    parser.add_argument(
        "--useSTARthreshold", action="store_true")
    
    args = parser.parse_args()

    demuxJson = GeneralUtils.readJSON("demuxMetrics.json")

    sampleSheetDf = pd.read_csv(args.samplesheet)

    if 'expectedCells' in sampleSheetDf.columns:
        expectedCells = int(
            sampleSheetDf.loc[
                sampleSheetDf['sample'] == args.sample]['expectedCells'])
    else:
        expectedCells = 0
    
    sample_metric(args.sample, demuxJson,
                  args.useSTARthreshold, args.isBarnyard, Path(args.star_out),
                  args.cells, args.libStruct, args.featureType,
                  expectedCells, args.topCellPercent,
                  args.minCellRatio, args.minReads)

if __name__ == "__main__":
    main()
