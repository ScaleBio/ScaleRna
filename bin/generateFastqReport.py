#!/usr/bin/env python
"""
Python script to generate fastq report from metrics
"""
import argparse
import functools
import json
from pathlib import Path
from typing import Dict
import pandas as pd
import datapane as dp
import plotly.express as px
import utils.myconstants as constants
from utils.base_logger import logger
from utils.ReportUtil import (GeneralUtils, CalculationUtils, DatapaneUtils)


def buildFastqReport(libName, libDir, libMetrics, internalReport, referencesPath):
    """
    Build the library report by calling relevant functions

    Args:
        libName (str): Library name
        libDir (str): Path to library json
        libMetrics (str): Path to library metrics
        internalReport (bool): Flag denoting whether report is for internal purposes'
        referencesPath (Path): Reference files
    """
    # Read in main metrics file using which we will build the report
    allCellsBetweenFiles = pd.read_csv(f"{libMetrics}/allCellsBetweenFiles.csv", index_col=0)

    # Read in demux metrics
    with open(f"{libMetrics}/demuxJson.json") as f:
        demuxJson = json.load(f)

    # Read in library json
    with open(referencesPath / libDir) as f:
        libJson = json.load(f)
    
    # Reports will be written to the reports folder
    writeDir = Path(".", "reports")
    Path(writeDir, "csv").mkdir(parents=True)

    # Call function that builds a datapane page that depicts a cellBarcodes vs umi
    # figure, a table with barcode read status and a reads per sample figure;
    # creates a dataframe that has the number of reads that passed and the number
    # of reads that have different errors;
    # creates a matplotlib plate plot for reads per rt well
    (overallPassingStats, readsPage, barcodeReadStatsInternal) = buildReadsPage(demuxJson, allCellsBetweenFiles, libName)

    # Call function that builds a datapane page that depicts a table with barcode
    # related information and plate plots that show reads per rt well, reads per
    # pcr well and reads per ligation well and also creates a dataframe with
    # information related to each barcoding level
    (barcodeStatsDf, cellsPage, barcodeTypeStats) = buildBarcodesPage(demuxJson, libName, libJson, allCellsBetweenFiles,
                                                                      writeDir, referencesPath)

    if internalReport:
        pages = [readsPage, cellsPage, dp.Page(blocks=[dp.Group(barcodeReadStatsInternal, barcodeTypeStats)], title='InternalReport')]
    else:
        # Concat two datapane pages
        pages = [readsPage, cellsPage]

    # Build a report object from the concatenated pages
    report = dp.Report(blocks=pages)

    
    prefix = f"library_{libName}"

    # Write to log file absolute path of reports directory
    logger.debug(f"Writing {prefix}.typeLevelMatches.csv and {prefix}.overallMatches.csv, "
                 f"to {str(writeDir.resolve())}")

    # Write information related to each barcoding level to a csv file
    barcodeStatsDf.to_csv(writeDir / "csv" / f"{prefix}.typeLevelMatches.csv", index=False)

    # Write dataframe that has the number of reads that passed and the number
    # of reads that have different errors
    overallPassingStats.to_csv(writeDir / "csv" / f"{prefix}.overallMatches.csv", index=False)

    report.save(writeDir / f"{prefix}.report.html")

def buildDfForPlatePlot(referencesPath, allCellsBetweenFiles, libJson, alias):
    """
    Construct dataframe that will be used for plotting heatmap

    Args:
        referencesPath (Path): Reference files
        allCellsBetweenFiles (pd.DataFrame): All cells information for this library
        libJson (dict): Library json information
        alias (str): Alias string corresponding to barcoding level in library json

    Returns:
        Constructed dataframe
    """
    wells = {}

    for entry in libJson["barcodes"]:
        if "alias" in entry:
            if entry["alias"] == alias:
                lib_json_dict_entry = entry
    
    with open(referencesPath / f'{lib_json_dict_entry["sequences"]}') as f:
        for line in f:
            line = line.strip()
            split_line = line.split("\t")
            wells[split_line[0]] = split_line[1]

    barcode_list = allCellsBetweenFiles[alias].to_list()

    well_list = [wells[x] for x in barcode_list]

    max_letter, max_number = DatapaneUtils.getMaxWellNumberAndLetter(referencesPath / f'{lib_json_dict_entry["sequences"]}')

    allCellsBetweenFiles[f'{alias.lower()}_well'] = well_list

    well_df = pd.DataFrame(0, columns=range(1, max_number+1), index=DatapaneUtils.getCharacterIndices(65,ord(max_letter)+1))
    
    for well in set(well_list):
        letter = well[-1]
        numbers = well[:-1]
        # Get umi count for each well
        well_df.at[letter, int(numbers)] = allCellsBetweenFiles.loc[allCellsBetweenFiles[f'{alias.lower()}_well'] == well, 'umis'].sum()
    
    return well_df

def buildBarcodesPage(demuxJson, libName, libJson, allCellsBetweenFiles, writeDir, referencesPath):
    """
    Function that builds a datapane page that depicts a table with barcode
    related information and plate plots that show reads per rt well, reads per
    pcr well and reads per ligation well. Also create a dataframe with
    information related to each barcoding level
    
    Args:
        demuxJson (dict): Dictionary with the demuxed metrics
        libName (str): Library name
        libJson (dict): Dictionary with library information obtained from the library json
        allCellsBetweenFiles (pd.DataFrame): All cell information for this library
        writeDir (Path): Write directory
        referencesPath (Path): Reference files

    Returns:
        Dataframe with stats related to barcode and dp.Pageobject
    """

    (barcodeTypeStatsDf, barcodeTypeStats) = createBarcodeTypeMetricsTables(demuxJson)

    ligation_well_df = buildDfForPlatePlot(referencesPath, allCellsBetweenFiles, libJson, "Ligation")
    ligation_well_df.to_csv(writeDir / "csv" / f"library_{libName}_unique_reads_ligation_well.csv")
    # Matplotlib figure that represents umi count per well for ligation
    ligationPerWell = DatapaneUtils.buildPlatePlot(ligation_well_df, "Ligation Plate", 100.0, "Unique Transcript Counts")

    pcr_well_df = buildDfForPlatePlot(referencesPath, allCellsBetweenFiles, libJson, "PCR")
    pcr_well_df.to_csv(writeDir / "csv" / f"library_{libName}_unique_reads_pcr_well.csv")
    # Matplotlib figure that represents umi count per well for pcr
    pcrPerWell = DatapaneUtils.buildPlatePlot(pcr_well_df, "PCR Plate", 100.0, "Unique Transcript Counts")
    
    rt_well_df = buildDfForPlatePlot(referencesPath, allCellsBetweenFiles, libJson, "RT")
    rt_well_df.to_csv(writeDir / "csv"/ f"library_{libName}_unique_reads_rt_well.csv")
    # Matplotlib figure that represents umi count per well for rt
    readPerWell = DatapaneUtils.buildPlatePlot(rt_well_df, "RT Plate", 100.0, "Unique Transcript Counts")

    blocks = [dp.Text(f"## libName: {libName}"),
              readPerWell, ligationPerWell, pcrPerWell]
    
    return (barcodeTypeStatsDf, dp.Page(blocks=blocks, title='Barcodes'), barcodeTypeStats)


def createBarcodeTypeMetricsTables(demuxJson):
    """
    Create dataframe for barcode type and create a datapane
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
        subset = barcodesDf[barcodesDf['Barcode'] == barcodeType][['Match', 'Reads']]
        styledDf = subset.style.pipe(GeneralUtils.styleTable, title=fullBarcodeTypeName)
        table = DatapaneUtils.createTableIgnoreWarning(styledDf)
        tableGroup.append(table)

    return (barcodesDf, dp.Group(blocks=tableGroup, columns=2))


def buildReadsPage(demuxJson, allCellsBetweenFiles, libName):
    """
    Function to build a datapane page for reads

    Args:
        demuxJson (dict): Dictionary with demuxed metrics
        allCellsBetweenFiles (pd.DataFrame): Dataframe containing data from all samples
        libName (str): Library name

    Returns:
        Dataframe with barcode reads information and dp.Page object
    """
    multiSampleKneePlot, allCellsBetweenFiles = makeMultisampleKneePlot(allCellsBetweenFiles)

    barcodeReadsData = demuxJson['reads']
    
    barcodeReadsPerc = buildDfFromJSONDict(barcodeReadsData, "Type", "list")
    
    barcodeReadsTotal = buildDfFromJSONDict(barcodeReadsData, "Type", "list", 0)
    barcodeReadsTotal = barcodeReadsTotal[['Type', 'Reads']]
    barcodeReadsTotal.rename(columns={'Type': 'Status'}, inplace=True)
    barcodeReadsTotal['Percent'] = barcodeReadsPerc['Reads']
    
    barcodeReadsTotalStyledInternal = barcodeReadsTotal.style.pipe(GeneralUtils.styleTable, title="Barcode Read Status", numericCols=['Reads'])

    total_reads = 0
    total_percent = 0
    df_error = barcodeReadsTotal[barcodeReadsTotal["Status"].str.contains("Error")]
    for idx, row in df_error.iterrows():
        total_reads += int(row["Reads"])
        total_percent += float(row["Percent"][:-1])
    df_pass = barcodeReadsTotal[barcodeReadsTotal["Status"].str.contains("Pass")]
    df_pass.loc[len(df_pass.index)] = ['Error', total_reads, str(round(total_percent, 1))+"%"]
    barcodeReadsTotalStyled = df_pass.style.pipe(GeneralUtils.styleTable, title="Barcode Read Status", numericCols=['Reads'])
    
    barcodeReadStats = DatapaneUtils.createTableIgnoreWarning(barcodeReadsTotalStyled)
    barcodeReadStatsInternal = DatapaneUtils.createTableIgnoreWarning(barcodeReadsTotalStyledInternal)

    (countsPerSampleDf, rtCountsPerSampleDf) = buildDfFromDemuxSampleMetrics(demuxJson)
    
    wellOrder = sorted(list(rtCountsPerSampleDf['rtWell'].unique()),
                       key=functools.cmp_to_key(CalculationUtils.wellStringComp))
    rtCountsPerSampleDf['rtWell'] = pd.Categorical(rtCountsPerSampleDf['rtWell'], wellOrder)
    rtCountsPerSampleDf.sort_values(by=['rtWell'], inplace=True, ascending=False)

    sampleOrder = list(rtCountsPerSampleDf.Sample.unique())
    sampleOrder.reverse()
    for sample in countsPerSampleDf['Sample']: # Add samples with no reads in 'rtCounts'
        if sample not in sampleOrder and sample != 'Unknown':
            sampleOrder.append(sample)
    sampleOrder.append('Unknown')
    countsPerSampleDf['Sample'] = pd.Categorical(countsPerSampleDf.Sample, sampleOrder)
    countsPerSampleDf.sort_values(by=['Sample'], inplace=True)

    colorMap = matchColorsToNames(list(countsPerSampleDf['Sample'].unique()))
    readsPerSample = px.bar(
        countsPerSampleDf, x='Sample', y='TotalReads', color='Sample',
        height=900, color_discrete_map=colorMap,
        template=constants.DEFAULT_FIGURE_STYLE,
        title="Reads Per Sample", labels={"TotalReads": "Total Reads"})
    readsPerSample.update_layout(showlegend=False)

    return (barcodeReadsTotal, dp.Page(blocks=[dp.Text(f"## libName: {libName}"),
                                               dp.Group(multiSampleKneePlot, barcodeReadStats, columns=2),
                                               dp.Group(readsPerSample, columns=1)], title='Reads'),
            barcodeReadStatsInternal)


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
    Build dataframes with per sample / RT-well bcParser metrics

    Args:
        demuxJson (dict): Dictionary with demuxed metrics

    Returns:
        One dataframe with total read counts per sample and one with read counts by RT well
        Note: The sample DF includes samples with 0 reads, the per-RT well one does not include empty wells
    """
    totalCounts = []
    rtCounts = []

    for sampleName, sampleDict in demuxJson['samples'].items():
        readCount = sampleDict['reads'][0]
        totalCounts.append({'Sample': sampleName, 'TotalReads': readCount})
        rtBarcodeCounts = sampleDict['barcodes']

        if rtBarcodeCounts:
            for rtWell, metrics in rtBarcodeCounts.items():
                rtCounts.append({'Sample': sampleName, 'rtWell': rtWell, 'ReadCount': metrics['reads']})

    return pd.DataFrame(totalCounts), pd.DataFrame(rtCounts)


def makeMultisampleKneePlot(allCellsBetweenFiles):
    """
    Makes a kneeplot using @field in @data; drawing a vertical line at @threshold

    Args:
        allCellsBetweenFiles (pd.DataFrame): Dataframe containing data from all samples

    Returns:
        dp.Plot object and dataframe with data from all samples
    """
    indiciesToInclude = set(CalculationUtils.getIndicesToInclude(max(allCellsBetweenFiles.index)))

    allCellsBetweenFiles['rankOrder'] = allCellsBetweenFiles.index

    plottingDf = allCellsBetweenFiles[allCellsBetweenFiles['rankOrder'].isin(indiciesToInclude)]

    fig = px.line(plottingDf, x=plottingDf.index, y='umis', color='sample', log_x=True, log_y=True,
                  template=constants.DEFAULT_FIGURE_STYLE, labels={"index": "Cell Barcodes", "umis": "Unique Transcript Counts"})

    return dp.Plot(fig), allCellsBetweenFiles


def buildDfFromJSONDict(jsonDict: Dict, name: str, valueDataType: str, choiceIndex=1) -> pd.DataFrame:
    """
    Build dataframe from json
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
    parser.add_argument("--libJsonName", required=False, default=None)
    parser.add_argument("--libMetrics", required=True)
    parser.add_argument("--libraryStructPath", type=str, help="Relative path to references folder")
    parser.add_argument("--internalReport", action="store_true", default=False)
    
    args = parser.parse_args()

    buildFastqReport(args.libName, args.libJsonName, Path(args.libMetrics), args.internalReport, Path(args.libraryStructPath))

if __name__ == "__main__":
    main()
