#!/usr/bin/env python
"""
Python script to generate fastq report from metrics
"""
import pandas as pd
import json
import datapane as dp
import functools
import plotly.express as px
import os
import argparse
from pathlib import Path
import utils.myconstants as constants
from utils.base_logger import logger
from typing import Dict
from utils.ReportUtil import (
    GeneralUtils, CalculationUtils, DatapaneUtils)


def buildFastqReport(libName, libDir, libMetrics, internalReport):
    """
    Function that builds the fastq report by calling relevant functions

    Args:
        libName (str): Library name
        libDir (str): Path to library json
        internalReport (bool): Flag denoting whether report is for internal purposes
    """
    # Read in main metrics file using which we will build the report
    allCellsBetweenFiles = pd.read_csv(f"{libMetrics}/allCellsBetweenFiles.csv",
                                       index_col=0)

    # Read in demux metrics
    with open(f"{libMetrics}/demuxJson.json") as f:
        demuxJson = json.load(f)

    # Read in library json
    with open(libDir) as f:
        libJson = json.load(f)
    
    # Reports will be written to the reports folder
    writeDir = Path(".", "reports")
    
    # If "reports" directory exists log contents of reports directory
    # and continue execution
    try:
        os.mkdir(writeDir)
    except OSError as error:
        logger.error(error)
        logger.debug(f"Directory reports exists with contents {os.listdir('reports')}")

    # Call function that builds a datapane page that depicts a cellBarcodes vs umi
    # figure, a table with barcode read status and a reads per sample figure;
    # creates a dataframe that has the number of reads that passed and the number
    # of reads that have different errors;
    # creates a matplotlib plate plot for reads per rt well
    (readPerWell, overallPassingStats, readsPage) = \
        buildReadsPage(demuxJson, allCellsBetweenFiles, libName, writeDir, internalReport)

    # Call function that builds a datapane page that depicts a table with barcode
    # related information and plate plots that show reads per rt well, reads per
    # pcr well and reads per ligation well and also creates a dataframe with
    # information related to each barcoding level
    (barcodeStatsDf, cellsPage) = buildBarcodesPage(demuxJson, libName, libJson,
                                                    readPerWell, allCellsBetweenFiles, writeDir, internalReport)

    # Concat two datapane pages
    pages = [readsPage, cellsPage]

    # Build a report object from the concatenated pages
    report = dp.Report(blocks=pages)

    
    prefix = f"library_{libName}"

    # Write to log file absolute path of reports directory
    logger.debug(f"Writing {prefix}.typeLevelMatches.tsv and {prefix}.overallMatches.tsv, "
                 f"to {str(writeDir.resolve())}")

    # Write information related to each barcoding level to a tsv file
    barcodeStatsDf.to_csv(
        writeDir / f"{prefix}.typeLevelMatches.tsv",
        sep='\t', index=False)

    # Write dataframe that has the number of reads that passed and the number
    # of reads that have different errors
    overallPassingStats.to_csv(
        writeDir / f"{prefix}.overallMatches.tsv",
        sep='\t', index=False)

    report.save(writeDir / f"{prefix}.report.html")


def getCharacterIndices(lower_limit, upper_limit):
    """
    Function to generate a list of letter indices from
    a range of ascii values

    Args:
        lower_limit (int): Lower limit to generate ascii's from
        upper_limit (int): Upper limit to generate ascii's till
    """
    return [chr(i) for i in range(lower_limit, upper_limit)]

def buildBarcodesPage(demuxJson, libName, libJson, readPerWell,
                   allCellsBetweenFiles, writeDir, internalReport):
    """
    Function that builds a datapane page that depicts a table with barcode
    related information and plate plots that show reads per rt well, reads per
    pcr well and reads per ligation well. Also create a dataframe with
    information related to each barcoding level
    Args:
        demuxJson (dict): Dictionary with the demuxed metrics
        libName (str): Library name
        libJson (dict): Dictionary with library information obtained from the
            library json
        readPerWell

    Returns:
        Dataframe with stats related to barcode and dp.Page
        object
    """

    (barcodeTypeStatsDf, barcodeTypeStats) = createBarcodeTypeMetricsTables(demuxJson)

    # Dictionary to hold barcode to well info. Barcode is the key and well
    # is the value
    ligation_wells = {}
    pcr_wells = {}

    # Read in reference file that contains barcode to well mapping for ligation
    with open(f'references/{libJson["barcodes"][0]["sequences"]}') as f:
        for line in f:
            line = line.strip()
            split_line = line.split("\t")
            ligation_wells[split_line[0]] = split_line[1]
    
    # Read in reference file that contains barcode to well mapping for pcr
    with open(f'references/{libJson["barcodes"][3]["sequences"]}') as f:
        for line in f:
            line = line.strip()
            split_line = line.split("\t")
            pcr_wells[split_line[0]] = split_line[1]
    
    # Get length of ligation and pcr barcode
    length_lig_barcode = libJson["barcodes"][0]["length"]
    length_pcr_barcode = libJson["barcodes"][3]["length"]
    
    ligation_well_list = []
    pcr_well_list = []
    
    # Get list of ligation barcodes
    ligation_barcode_list = \
        [x[:length_lig_barcode] for x in allCellsBetweenFiles["Barcode"].to_list()]
    
    # Get list of pcr barcodes
    pcr_barcode_list = \
        [x[-length_pcr_barcode:] for x in allCellsBetweenFiles["Barcode"].to_list()]
    
    # Get list of ligation and pcr wells
    ligation_well_list = [ligation_wells[x] for x in ligation_barcode_list]
    pcr_well_list = [pcr_wells[x] for x in pcr_barcode_list]
   
    allCellsBetweenFiles['ligation_well'] = ligation_well_list
    allCellsBetweenFiles['pcr_well'] = pcr_well_list

    # Build ligation and pcr dataframe that will have the letters as the dataframe
    # index and the numbers as the dataframe columns
    # This dataframe will contain the per well information used to build the
    # plate plots
    ligation_well_df = pd.DataFrame(0, columns=range(1, 25),
                                    index=getCharacterIndices(65,81))
    pcr_well_df = pd.DataFrame(0, columns=range(1, 13),
                               index=getCharacterIndices(65,73))
    
    for well in set(ligation_well_list):
        letter = well[-1]
        numbers = well[:-1]
        # Get umi count for each well
        ligation_well_df.at[letter, int(numbers)] = \
            allCellsBetweenFiles.loc[allCellsBetweenFiles['ligation_well'] == well, 'umis'].sum()
    
    for well in set(pcr_well_list):
        letter = well[-1]
        numbers = well[:-1]
        # Get umi count for each well
        pcr_well_df.at[letter, int(numbers)] = \
            allCellsBetweenFiles.loc[allCellsBetweenFiles['pcr_well'] == well, 'umis'].sum()


    ligation_well_df.to_csv(writeDir / "unique_reads_ligation_well.csv")
    # Matplotlib figure that represents umi count per well for ligation
    ligationPerWell = DatapaneUtils.buildPlatePlot(ligation_well_df,
                                                   "Ligation Plate", "Unique Transcript Counts")

    pcr_well_df.to_csv(writeDir / "unique_reads_pcr_well.csv")
    # Matplotlib figure that represents umi count per well for pcr
    pcrPerWell = DatapaneUtils.buildPlatePlot(pcr_well_df,
                                              "PCR Plate", "Unique Transcript Counts")

    # Build list where each element represents an individual block
    # of a datapane page
    if internalReport:
        blocks = [dp.Text(f"## libName: {libName}"),
                  barcodeTypeStats, readPerWell, ligationPerWell, pcrPerWell]
    else:
        blocks = [dp.Text(f"## libName: {libName}"),
                  readPerWell, ligationPerWell, pcrPerWell]
    
    return (barcodeTypeStatsDf, dp.Page(blocks=blocks, title='Barcodes'))


def createBarcodeTypeMetricsTables(demuxJson):
    """
    Function to create dataframe for barcode type and create a datapane
    object for storing a table created with the statistics in the dataframe

    Args:
        demuxJson (dict): Dictionary of demuxed metrics

    Returns:
        Dataframe and dp.Group object
    """
    barcodes = demuxJson['barcodes']
    barcodesDf = buildDfFromJSONDict(barcodes, "Barcode", "dict")
    tableGroup = []
    allBarcodeTypes = list(barcodesDf['Barcode'].unique())

    for barcodeType in allBarcodeTypes:
        fullBarcodeTypeName = constants.BARCODE_SHORTHAND_TO_NAME[barcodeType]
        subset = \
            barcodesDf[barcodesDf['Barcode'] == barcodeType][['Match', 'Reads']]
        styledDf = subset.style.pipe(GeneralUtils.styleTable,
                                     title=fullBarcodeTypeName)
        table = DatapaneUtils.createTableIgnoreWarning(styledDf)
        tableGroup.append(table)

    return (barcodesDf, dp.Group(blocks=tableGroup, columns=2))


def buildReadsPage(demuxJson, allCellsBetweenFiles, libName, writeDir, internalReport):
    """
    Function to build a datapane page for reads

    Args:
        demuxJson (dict): Dictionary with demuxed metrics
        allCellsBetweenFiles (pd.DataFrame): Dataframe containing data from
            all samples
        libName (str): Library name

    Returns:
        Dataframe with barcode reads information and dp.Page object
    """
    multiSampleKneePlot, allCellsBetweenFiles = \
        makeMultisampleKneePlot(allCellsBetweenFiles)

    barcodeReadsData = demuxJson['reads']
    
    barcodeReadsPerc = buildDfFromJSONDict(barcodeReadsData, "Type", "list")
    
    barcodeReadsTotal = buildDfFromJSONDict(barcodeReadsData, "Type", "list", 0)
    barcodeReadsTotal = barcodeReadsTotal[['Type', 'Reads']]
    barcodeReadsTotal.rename(columns={'Type': 'Status'}, inplace=True)
    barcodeReadsTotal['Percent'] = barcodeReadsPerc['Reads']
    
    if internalReport:
        barcodeReadsTotalStyled = barcodeReadsTotal.style.pipe(GeneralUtils.styleTable,
                                                           title="Barcode Read Status",
                                                           numericCols=['Reads'])
    else:
        total_reads = 0
        total_percent = 0
        df_error = barcodeReadsTotal[barcodeReadsTotal["Status"].str.contains("Error")]
        for idx, row in df_error.iterrows():
            total_reads += int(row["Reads"])
            total_percent += float(row["Percent"][:-1])
        df_pass = barcodeReadsTotal[barcodeReadsTotal["Status"].str.contains("Pass")]
        df_pass.loc[len(df_pass.index)] = ['Error', total_reads, str(round(total_percent, 1))+"%"]
        barcodeReadsTotalStyled = df_pass.style.pipe(GeneralUtils.styleTable,
                                                     title="Barcode Read Status",
                                                     numericCols=['Reads'])
    
    barcodeReadStats = DatapaneUtils.createTableIgnoreWarning(barcodeReadsTotalStyled)

    (countsPerSampleDf, tgmtCountsPerSampleDf) = buildDfFromDemuxSampleMetrics(demuxJson)
    
    sortOrder = sorted(list(tgmtCountsPerSampleDf.tgmtBc.unique()),
                       key=functools.cmp_to_key(CalculationUtils.wellStringComp))
    
    tgmtCountsPerSampleDf['tgmtBc'] = pd.Categorical(tgmtCountsPerSampleDf.tgmtBc,
                                                     sortOrder)
    
    tgmtCountsPerSampleDf.sort_values(by=['tgmtBc'], inplace=True, ascending=False)

    sampleOrder = list(tgmtCountsPerSampleDf.Sample.unique())
    sampleOrder.reverse()

    countsPerSampleDf['Sample'] = pd.Categorical(countsPerSampleDf.Sample,
                                                 sampleOrder + ["Unknown"])
    countsPerSampleDf.sort_values(by=['Sample'], inplace=True)

    colorMapToUse = matchColorsToNames(list(countsPerSampleDf['Sample'].unique()))

    readsPerSample = px.bar(
        countsPerSampleDf, x='Sample', y='TotalReads', color='Sample',
        height=900, color_discrete_map=colorMapToUse,
        template=constants.DEFAULT_FIGURE_STYLE,
        title="Reads Per Sample", labels={"TotalReads": "Total Reads"})
    readsPerSample.update_layout(showlegend=False)

    wellPlateTgmtCountDf = pd.DataFrame(
        0, columns=range(1, 13), index=getCharacterIndices(65,73))
    
    for idx, row in tgmtCountsPerSampleDf.iterrows():
        letter = row['tgmtBc'][-1]
        numbers = row['tgmtBc'][:-1]
        wellPlateTgmtCountDf.at[letter, int(numbers)] = row["ReadCount"]
    
    readPerWell = DatapaneUtils.buildPlatePlot(
        wellPlateTgmtCountDf, "RT Plate", "Unique Transcript Counts"
    )
    wellPlateTgmtCountDf.to_csv(writeDir / "unique_reads_rt_well.csv")

    return (
        readPerWell,
        barcodeReadsTotal,
        dp.Page(
            blocks=[dp.Text(f"## libName: {libName}"),
                    dp.Group(multiSampleKneePlot,
                             barcodeReadStats,
                             columns=2),
                    dp.Group(readsPerSample,
                             columns=1)], title='Reads'))


def matchColorsToNames(names) -> Dict[str, str]:
    """
    Associated colors with categorical values @names
    for consistency between plots

    Args:
        names (list): List of sample names

    Returns:
        Dictionary of color map
    """
    colorMap = {"Unknown": 'rgb(179, 188, 201)'}

    if ("Unknown" in names):
        names.remove('Unknown')

    colors = px.colors.qualitative.D3
    for i, name in enumerate(sorted(names)):
        colorMap[name] = colors[i % len(colors)]

    return colorMap


def buildDfFromDemuxSampleMetrics(demuxJson):
    """
    Function to build dataframe from demuxed metrics

    Args:
        demuxJson (dict): Dictionary with demuxed metrics

    Returns:
        Two dataframes
    """
    totalCounts = []
    tgmtCounts = []

    for sampleName, sampleDict in demuxJson['samples'].items():
        readCount = sampleDict['reads'][0]
        totalCounts.append({'Sample': sampleName, 'TotalReads': readCount})
        tgmtBarcodeCounts = sampleDict['barcodes']

        if (len(tgmtBarcodeCounts) > 0):
            for tgmtId, metrics in tgmtBarcodeCounts.items():
                tgmtCounts.append({'Sample': sampleName, 'tgmtBc': tgmtId,
                                   'ReadCount': metrics['reads']})

    return pd.DataFrame(totalCounts), pd.DataFrame(tgmtCounts)


def makeMultisampleKneePlot(allCellsBetweenFiles):
    """
    Makes a kneeplot using @field in @data; drawing a vertical line at
    @threshold

    Args:
        allCellsBetweenFiles (pd.DataFrame): Dataframe containing data from
            all samples

    Returns:
        dp.Plot object and dataframe with data from all samples
    """
    indiciesToInclude = set(
        CalculationUtils.getIndicesToInclude(
            max(allCellsBetweenFiles.index)))

    allCellsBetweenFiles['rankOrder'] = allCellsBetweenFiles.index

    plottingDf = allCellsBetweenFiles[
        allCellsBetweenFiles['rankOrder'].isin(indiciesToInclude)]

    fig = px.line(
        plottingDf, x=plottingDf.index, y='umis', color='sample',
        log_x=True, log_y=True,
        template=constants.DEFAULT_FIGURE_STYLE,
        labels={"index": "Cell Barcodes", "umis": "Unique Transcript Counts"})

    return dp.Plot(fig), allCellsBetweenFiles


def buildDfFromJSONDict(jsonDict: Dict, name: str, valueDataType: str,
                        choiceIndex=1) -> pd.DataFrame:
    """
    Function to build a dataframe from a json containing metrics.
    Used as helper to process demux metrics.json

    Args:
        jsonDict (dict): Dictionary with metrics
        name (str):
        valueDataType (str): Datatype of @name
        choiceIndex (int):

    Returns:
        Computed dataframe
    """
    dictList = []

    for (keyName, valueObj) in jsonDict.items():
        if valueDataType == 'dict':

            for key, value in valueObj.items():
                newDict = {}
                newDict[name] = keyName
                newDict['Match'] = key
                newDict['Reads'] = value[choiceIndex]
                dictList.append(newDict)

        elif valueDataType == 'list':
            newDict = {"Reads": valueObj[choiceIndex]}
            newDict[name] = keyName
            dictList.append(newDict)

        elif valueDataType == 'int':
            newDict = {"Cells": valueObj}
            newDict[name] = keyName
            dictList.append(newDict)

    return pd.DataFrame(dictList)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--libName", required=False, default=None)
    parser.add_argument("--libDir", required=False, default=None)
    parser.add_argument("--libMetrics", required=True)
    parser.add_argument("--internalReport", action="store_true", default=False)
    
    args = parser.parse_args()

    buildFastqReport(args.libName, args.libDir, Path(args.libMetrics), args.internalReport)

if __name__ == "__main__":
    main()
