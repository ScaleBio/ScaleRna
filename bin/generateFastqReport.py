#!/usr/bin/env python
"""
Python script to generate fastq report from metrics
"""
import argparse
import functools
import json
from pathlib import Path
from typing import Dict

import datapane as dp
import pandas as pd
import plotly.express as px

from utils import fileUtils, reportUtil, statsUtils
from utils.base_logger import logger

BARCODE_SHORTHAND_TO_NAME = {
    'drop': 'Droplet Barcodes', 'P7': 'P7 Barcodes',
    'lig': 'Ligation Barcodes', 'rt': 'Reverse Transcription Barcodes',
    'umi': 'UMI'}

def buildFastqReport(libName:str, libJson:Path, demuxJson:Path, libMetrics:Path, internalReport:bool):
    """
    Build the library report by calling relevant functions

    Args:
        libName: Library name
        libJson: Library structure json
        libMetrics: Path to library metrics
        internalReport: Flag denoting whether report is for internal purposes'
    """
    allCellsBetweenFiles = pd.read_csv(libMetrics / "allCellsBetweenFiles.csv", index_col=0)
    demuxMetrics = json.load(open(demuxJson))
    libStruct = json.load(open(libJson))
    
    # Reports will be written to the reports folder
    writeDir = Path(".", "reports")
    Path(writeDir, "csv").mkdir(parents=True, exist_ok=True)

    # Call function that builds a datapane page that depicts a cellBarcodes vs umi
    # figure, a table with barcode read status and a reads per sample figure;
    # creates a dataframe that has the number of reads that passed and the number
    # of reads that have different errors;
    # creates a matplotlib plate plot for reads per rt well
    (overallPassingStats, readsPage, barcodeReadStatsInternal) = buildReadsPage(demuxMetrics, allCellsBetweenFiles, libName)

    # Call function that builds a datapane page that depicts a table with barcode
    # related information and plate plots that show reads per rt well, reads per
    # pcr well and reads per ligation well and also creates a dataframe with
    # information related to each barcoding level
    (barcodeStatsDf, cellsPage, barcodeTypeStats) = buildBarcodesPage(demuxMetrics, libName, libJson, allCellsBetweenFiles, writeDir)

    pages = [readsPage, cellsPage]
    if internalReport:
        pages.append(dp.Page(blocks=[dp.Group(barcodeReadStatsInternal, barcodeTypeStats)], title='InternalReport'))
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

def buildDfForPlatePlot(allCellsBetweenFiles:pd.DataFrame, libJson:Path, bcName:str):
    """
    Construct dataframe that will be used for plotting heatmap

    Args:
        allCellsBetweenFiles: All cells information for this library
        libJson: Library structure information
        bcName: name of the barcode in lib.json

    Returns:
        Constructed dataframe
    """
    libStruct = json.load(open(libJson))
    libStructDir = libJson.parent # Directory containing libStruct.json and associated sequence files

    for bc in libStruct["barcodes"]:
        if bc["name"] == bcName:
            bcInfo = bc
            break
    else:
        raise ValueError(f"Unknown barcode {bcName}")
    
    wells = {}
    with open(libStructDir / f'{bcInfo["sequences"]}') as f:
        for line in f:
            line = line.strip()
            split_line = line.split("\t")
            wells[split_line[0]] = split_line[1]

    alias = bcInfo.get('alias') or bcInfo['name']
    barcode_list = allCellsBetweenFiles[alias].to_list()
    well_list = [wells[x] for x in barcode_list]
    max_letter, max_number = reportUtil.getMaxWellNumberAndLetter(libStructDir / f'{bcInfo["sequences"]}')
    allCellsBetweenFiles[f'{alias.lower()}_well'] = well_list
    well_df = pd.DataFrame(0, columns=range(1, max_number+1), index=reportUtil.getCharacterIndices(65, ord(max_letter)+1))
    
    for well in set(well_list):
        letter = well[-1]
        numbers = well[:-1]
        # Get umi count for each well
        well_df.at[letter, int(numbers)] = allCellsBetweenFiles.loc[allCellsBetweenFiles[f'{alias.lower()}_well'] == well, 'umis'].sum()
    
    return well_df

def buildBarcodesPage(demuxJson:Dict, libName:str, libJson:Path, allCellsBetweenFiles:pd.DataFrame, writeDir:Path):
    """
    Function that builds a datapane page that depicts a table with barcode
    related information and plate plots that show reads per rt well, reads per
    pcr well and reads per ligation well. Also create a dataframe with
    information related to each barcoding level
    
    Args:
        demuxJson: Dictionary with the demuxed metrics
        libName: Library name
        libJson: Library Structure Json
        allCellsBetweenFiles: All cell information for this library
        writeDir: Write directory

    Returns:
        Dataframe with stats related to barcode and dp.Pageobject
    """

    (barcodeTypeStatsDf, barcodeTypeStats) = createBarcodeTypeMetricsTables(demuxJson, libJson)

    ligation_well_df = buildDfForPlatePlot(allCellsBetweenFiles, libJson, "lig")
    ligation_well_df.to_csv(writeDir / "csv" / f"library_{libName}_unique_reads_ligation_well.csv")
    # Matplotlib figure that represents umi count per well for ligation
    ligationPerWell = reportUtil.buildPlatePlot(ligation_well_df, "Ligation Plate", 10000.0, "Unique Transcript Counts")

    pcr_well_df = buildDfForPlatePlot(allCellsBetweenFiles, libJson, "pcr")
    pcr_well_df.to_csv(writeDir / "csv" / f"library_{libName}_unique_reads_pcr_well.csv")
    # Matplotlib figure that represents umi count per well for pcr
    pcrPerWell = reportUtil.buildPlatePlot(pcr_well_df, "PCR Plate", 10000.0, "Unique Transcript Counts")
    
    rt_well_df = buildDfForPlatePlot(allCellsBetweenFiles, libJson, "rt")
    rt_well_df.to_csv(writeDir / "csv"/ f"library_{libName}_unique_reads_rt_well.csv")
    # Matplotlib figure that represents umi count per well for rt
    readPerWell = reportUtil.buildPlatePlot(rt_well_df, "RT Plate", 10000.0, "Unique Transcript Counts")

    blocks = [dp.Text(f"## libName: {libName}"),
              readPerWell, ligationPerWell, pcrPerWell]
    
    return (barcodeTypeStatsDf, dp.Page(blocks=blocks, title='Barcodes'), barcodeTypeStats)


def getBarcodeAlias(libStruct: Dict, name: str):
    """Get the full name (alias) for a barcode from the library structure Json"""
    for bc in libStruct['barcodes']:
        if bc['name'] == name:
            return bc.get('alias', name)
    return None


def createBarcodeTypeMetricsTables(demuxMetrics: Dict, libJson: Dict):
    """
    Create dataframe for barcode type and create a datapane
    object for storing a table created with the statistics in the dataframe

    Args:
        demuxMetrics: bcParser metrics (from metrics.json)
        libStruct: Library structure (lib.json)

    Returns:
        Dataframe and dp.Group object
    """
    libStruct = json.load(open(libJson))
    barcodesDf = buildDfFromJSONDict(demuxMetrics['barcodes'], "Barcode", "dict")
    tableGroup = []
    allBarcodes = list(barcodesDf['Barcode'].unique())
    for bc in allBarcodes:
        subset = barcodesDf[barcodesDf['Barcode'] == bc][['Match', 'Reads']]
        styledDf = subset.style.pipe(reportUtil.styleTable, title=f"{getBarcodeAlias(libStruct, bc)} Barcodes")
        table = reportUtil.mkTable(styledDf)
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
    
    barcodeReadsTotalStyledInternal = barcodeReadsTotal.style.pipe(reportUtil.styleTable, title="Barcode Read Status", numericCols=['Reads'])

    total_reads = 0
    total_percent = 0
    df_error = barcodeReadsTotal[barcodeReadsTotal["Status"].str.contains("Error")]
    for idx, row in df_error.iterrows():
        total_reads += int(row["Reads"])
        total_percent += float(row["Percent"][:-1])
    df_pass = barcodeReadsTotal[barcodeReadsTotal["Status"].str.contains("Pass")]
    df_pass.loc[len(df_pass.index)] = ['Error', total_reads, str(round(total_percent, 1))+"%"]
    barcodeReadsTotalStyled = df_pass.style.pipe(reportUtil.styleTable, title="Barcode Read Status", numericCols=['Reads'])
    
    barcodeReadStats = reportUtil.mkTable(barcodeReadsTotalStyled)
    barcodeReadStatsInternal = reportUtil.mkTable(barcodeReadsTotalStyledInternal)

    (countsPerSampleDf, rtCountsPerSampleDf) = buildDfFromDemuxSampleMetrics(demuxJson)
    
    wellOrder = sorted(list(rtCountsPerSampleDf['rtWell'].unique()),
                       key=functools.cmp_to_key(reportUtil.wellStringComp))
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
        template=reportUtil.DEFAULT_FIGURE_STYLE,
        title="Reads Per Sample", labels={"TotalReads": "Total Reads"})
    readsPerSample.update_layout(showlegend=False)

    readsPage = dp.Page(blocks=[dp.Text(f"## libName: {libName}"),
                       dp.Group(multiSampleKneePlot, barcodeReadStats, columns=2),
                       dp.Group(readsPerSample, columns=1)], title='Reads')
    return (barcodeReadsTotal, readsPage, barcodeReadStatsInternal)


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
    maxIndex = max(allCellsBetweenFiles.index) if len(allCellsBetweenFiles.index) else 0
    indices = set(reportUtil.sparseLogCoords(maxIndex))
    plottingDf = allCellsBetweenFiles[allCellsBetweenFiles.index.isin(indices)]
    fig = px.line(plottingDf, x=plottingDf.index, y='umis', color='sample', log_x=True, log_y=True,
                  template=reportUtil.DEFAULT_FIGURE_STYLE, labels={"index": "Cell Barcodes", "umis": "Unique Transcript Counts"})
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
    parser.add_argument("--libName", required=True)
    parser.add_argument("--libStruct", required=True, type=Path)
    parser.add_argument("--libMetrics", required=True, type=Path)
    parser.add_argument("--demuxMetrics", required=True, type=Path, help="bcParser demux metrics json")
    parser.add_argument("--internalReport", action="store_true", default=False)
    
    args = parser.parse_args()

    buildFastqReport(args.libName, args.libStruct, args.demuxMetrics, args.libMetrics, args.internalReport)

if __name__ == "__main__":
    main()
