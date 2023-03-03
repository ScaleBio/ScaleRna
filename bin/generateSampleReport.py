#!/usr/bin/env python
"""
Python script to generate sample report from metrics
"""
import pandas as pd
import json
import datapane as dp
import numpy as np
import math
import os
import argparse
import plotly.express as px
from pathlib import Path
import utils.myconstants as constants
from typing import Dict
from utils.ReportUtil import (
    GeneralUtils, CalculationUtils, DatapaneUtils)
from utils.base_logger import logger


def buildSampleReport(sample_metrics, sample, libDir, internalReport, isBarnyard):
    """
    Entry function for building sample report

    Args:
        sample (str): Sample name
        libDir (str): Path to library structure json
        internalReport (bool): Flag denoting whether report is
            for internal purposes
        isBarnyard (bool): Flag denoting whether sample is barnyard
            or not
    """
    allCells = pd.read_csv(f"{sample_metrics}/allCells.csv", index_col=0)

    with open(f"{sample_metrics}/sample_metrics/barcodesToPlot.json") as f:
        barcodesToPlot = json.load(f)

    writeDir = Path(".", "reports")

    os.mkdir(writeDir)
    os.mkdir(writeDir / f"{sample}_figures")
    
    pages = []

    readsPage, statsDf, complexity_df = buildReadsPage(
        sample, allCells, writeDir, internalReport)
    pages.append(readsPage)

    if libDir is not None:
        barcodesPage = buildBarcodesPage(
            allCells, barcodesToPlot, writeDir, sample, internalReport)
        pages.append(barcodesPage)

    if isBarnyard:
        scoreBarn(allCells)
        (barnyardPage, barnyardStat) = buildBarnyardPage(allCells, writeDir, sample)
        statsDf = pd.concat([statsDf, barnyardStat])
        pages.append(barnyardPage)

    report = dp.Report(blocks=pages)

    logger.debug("Writing reportStatistics, "
                 f"report to {str(writeDir.resolve())}")

    new_complexity_df = pd.DataFrame(columns=['Metric','Value','Category'])
    for index, row in complexity_df.iterrows():
        new_complexity_df.loc[len(new_complexity_df.index)] = \
            [f"Median unique transcript counts at {row['x']} total reads per cell",
             row["unique_read"], "Extrapolated Complexity"]
    
    statsDf = pd.concat([statsDf, new_complexity_df], axis=0)

    statsDf.to_csv(
        writeDir / f"{sample}.reportStatistics.tsv",
        index=False, sep='\t', header=False,
        columns=["Category", "Metric", "Value"])

    report.save(writeDir / f"{sample}.report.html")


def checkIfBarnyard(allCells):
    """
    Function to check if sample includes both human and mouse cells

    Args:
        allCells (pd.DataFrame): Dataframe with all the data for this sample

    Returns:
        True if both human and mouse are present. False if not
    """
    scoreBarn(allCells)
    #includedSpecies = list(allCells.species.unique())
    return False
    return "Human" in includedSpecies and "Mouse" in includedSpecies


def getBackground(minorFracs):
    """
    Function to calculate median and standard deviation for barnyard
    cells

    Args:
        minorFracs (list): List with fraction indicating percentage
            of human or mouse samples in a cell

    Returns:
        Median and standard deviation
    """
    if len(minorFracs) == 0:
        return np.nan, np.nan
    median = minorFracs.median()
    std = (np.quantile(minorFracs, 0.75) - median) / 0.675
    return median, std


def scoreBarn(allCells):
    """
    Function to label each cell as belonging to either human or mouse

    Args:
        allCells (pd.DataFrame): Dataframe with all the data for this sample
    """
    cells = allCells
    cells['species'] = "None"
    cells['minorFrac'] = cells[
        ['humanUmis', 'mouseUmis']].min(1) / cells[
            ['humanUmis', 'mouseUmis']].sum(1)
    cells.loc[
        cells['pass'] & (
            cells.humanUmis > cells.mouseUmis), 'species'] = 'Human'
    cells.loc[
        cells['pass'] & (
            cells.mouseUmis > cells.humanUmis), 'species'] = 'Mouse'

    minHuman = minMouse = cells[cells['pass']].umis.min()
    if (cells.species == 'Human').any() and (cells.species == 'Mouse').any():
        minHuman = np.percentile(
            cells.humanUmis[cells.species == 'Human'], 10)
        minMouse = np.percentile(
            cells.mouseUmis[cells.species == 'Mouse'], 10)

    # Labelling low and doublet cells
    cells.loc[
        (cells.humanUmis < minHuman) & (cells.mouseUmis < minMouse),
        'species'] = "None"
    cells.loc[
        (cells.humanUmis >= minHuman) & (cells.mouseUmis >= minMouse),
        'species'] = "Mixed"

    # Labelling as Ambiguous
    humanBackMed, humanBackStd = getBackground(
        cells[(
            cells.species == 'Mouse')].minorFrac)
    mouseBackMed, mouseBackStd = getBackground(
        cells[(
            cells.species == 'Mouse')].minorFrac)

    cells.loc[
        (cells.species == 'Mixed') &
        (cells.humanUmis > cells.mouseUmis) &
        (cells.minorFrac < mouseBackMed + 2 * mouseBackStd),
        'species'] = 'Ambiguous'

    cells.loc[
        (cells.species == 'Mixed') &
        (cells.mouseUmis > cells.humanUmis) &
        (cells.minorFrac < humanBackMed + 2 * humanBackStd),
        'species'] = 'Ambiguous'

    cells.loc[
        (cells.species == 'Human') &
        (cells.minorFrac >= max(0.1, mouseBackMed + 3 * mouseBackStd)),
        'species'] = 'Ambiguous'

    cells.loc[
        (cells.species == 'Mouse') &
        (cells.minorFrac >= max(0.1, humanBackMed + 3 * humanBackStd)),
        'species'] = 'Ambiguous'


def makeReadStatsDf(allCells):
    stats = []
    stats.append(['Sample Reads', format(allCells['reads'].sum(), ",.0f")])
    stats.append(['Reads Mapped to Genome',
                  format(allCells['mappedReads'].sum() / allCells['reads'].sum(), ".1%")])
    stats.append(['Reads Mapped to Transcriptome',
                  format(allCells['geneReads'].sum() / allCells['mappedReads'].sum(), ".1%")])
    stats.append(['Exonic Reads',
                  format(allCells['exonReads'].sum() / allCells['mappedReads'].sum(), ".1%")])
    stats.append(['Antisense Reads',
                  format(allCells['antisenseReads'].sum() / allCells['mappedReads'].sum(), ".1%")])
    stats.append(['Mitochondrial Reads',
                  format(allCells['mitoReads'].sum() / allCells['mappedReads'].sum(), ".1%")])
    stats.append(['Saturation',
                  format(1 - (allCells.umis.sum() / allCells.passingReads.sum()), ".2f")])

    statsAsDictList = [
        GeneralUtils.makeTableDict(
            ['Metric', 'Value'], valuePair) for valuePair in stats]

    statsDf = pd.DataFrame(statsAsDictList)
    statsDf['Category'] = 'Mapping'
    return statsDf

def buildReadsPage(sample, allCells, writeDir, internalReport):
    """
    Function to build a datapane page that mainly shows plots related
    to reads

    Args:
        sample (str): Sample name
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        writeDir (path): Path object pointing to write directory
        internalReport (bool): Flag indicating whether report is to be
            generated for internal consumption

    Returns:
        dp.Page object and concatenated dataframe of various statistics
    """
    statsDf = pd.DataFrame({'Metric': ['SampleName'], 'Value': [sample]})
    statsDf['Category'] = 'Sample'

    readStatsDf = makeReadStatsDf(allCells)
    cellStatsDf = makeCellStatsDf(allCells)

    uniqueReadsFig, complexity_df = makeUniqueReadsFig(allCells)

    summaryStatsTable = DatapaneUtils.createTableIgnoreWarning(
        cellStatsDf[
            constants.DISPLAY_COLUMNS].style.pipe(
                GeneralUtils.styleTable,
                title="Cell Metrics",
                hideColumnHeaders=True,
                boldColumn=['Metric'],
                numericCols=['Value']))

    # We are currently not using any summary metrics computed by STAR directly
    # starStatsDf = makeStarDfFromDict(starStats)

    mappingTable = DatapaneUtils.createTableIgnoreWarning(
        readStatsDf[
            constants.DISPLAY_COLUMNS].style.pipe(
                GeneralUtils.styleTable,
                title="Mapping Metrics",
                hideColumnHeaders=True,
                boldColumn=['Metric']))

    statsDf = pd.concat([statsDf, readStatsDf, cellStatsDf])
    rankPlot = makeRankPlot(allCells, writeDir, sample, internalReport)
    genesVersusUMIsScatter = buildGeneXReadScatter(allCells)
    saturationScatter = buildSaturationScatter(allCells)

    if internalReport:
        geneReadsDist = buildDist(
            "Distribution of Reads: Gene Matched",
            'geneReads',
            'Proportion of Genic Reads',
            allCells)
        antisenseReadsDist = buildDist(
            "Distribution of Reads: Antisense",
            "antisenseReads",
            "Proportion of Antisense Reads",
            allCells)
        exonicReadsDist = buildDist(
            "Distribution of Reads: Exonic",
            'exonReads', "Proportion of Exonic Reads", allCells)

        groupBlocks = [rankPlot, uniqueReadsFig, summaryStatsTable,
                       genesVersusUMIsScatter, saturationScatter,
                       geneReadsDist,  exonicReadsDist, antisenseReadsDist]

        anyNonZeroMitoReads = (allCells[allCells['pass']].mitoReads > 0).any()

        if anyNonZeroMitoReads:
            mitoReadsDist = buildDist(
                "Distribution of Reads: Mitochondrial", 'mitoReads',
                "Proportion of Mitochondrial Reads", allCells)
            groupBlocks.append(mitoReadsDist)
    else:
        groupBlocks = [rankPlot, uniqueReadsFig, mappingTable, summaryStatsTable, genesVersusUMIsScatter, saturationScatter]
    page = dp.Page(dp.Group(blocks=groupBlocks, columns=2), title="Cells")

    return page, statsDf, complexity_df


def buildBarcodesPage(allCells, barcodesToPlot, writeDir, sampleName, internalReport):
    """
    Function to build datapane page for barcodes

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        barcodesToPlot (dict): Dictionary containing the barcodes to plot
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name

    Returns:
        dp.Page object which holds the relevant figures
    """
    blocksToRender = []

    blocksToRender.append(createBarcodeLevelFigures(
        allCells, barcodesToPlot, writeDir, sampleName, internalReport))

    return dp.Page(blocks=blocksToRender, title="Barcodes")


def createBarcodeLevelFigures(allCells, barcodesToPlot, writeDir,
                              sampleName, internalReport):
    """
    Function to create plots for the barcodes page

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        barcodesToPlot (dict): Dictionary of barcodes to plot
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name

    Returns:
        dp.Group object with the relevant figures
    """
    passingCells = allCells[allCells['pass']]
    allIndexPlots = []
    key_to_remove = []
    if not internalReport:
        for key in barcodesToPlot:
            if barcodesToPlot[key]['alias'] == "Ligation" or barcodesToPlot[key]['alias'] == "PCR":
                key_to_remove.append(key)

    for key in key_to_remove:
        barcodesToPlot.pop(key)

    for key in sorted(barcodesToPlot):
        data = barcodesToPlot[key]
        bcColName = data['alias']
        possibleValues = list(data['possibleValues'])

        indexPlots = DatapaneUtils.barcodeLevelPlots(
            sampleName, passingCells, possibleValues, f"{bcColName}_alias",
            f"{bcColName} Index", internalReport,
            wellAliases=data['orderAsWellAliases'],
            writeDir=writeDir / f"{sampleName}_figures")

        allIndexPlots.append(indexPlots)

    return dp.Group(blocks=allIndexPlots)


def buildDist(title, field, replacementStr, allCells):
    """
    Function to build histogram

    Args:
        title (str): Title of plot
        field (str): Column name in @allCells to plot on x axis
        replacementStr (str): String that will be the x axis label
            instead of @field
        allCells (pd.DataFrame): Dataframe containing all data for this sample

    Returns:
        dp.Plot object with the figure
    """
    passingCells = allCells[allCells['pass']]

    fig = px.histogram(
        passingCells, x=field, title=title, histnorm='percent',
        template=constants.DEFAULT_FIGURE_STYLE,
        labels={field: replacementStr})

    fig.update_layout(yaxis_title="% Cells")

    return dp.Plot(fig)


def buildGeneXReadScatter(allCells):
    """
    Function to make a scatter plot of reads vs umi

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample

    Returns:
        dp.Plot object with the figure
    """
    cellsToPlot = allCells[:allCells['pass'].sum()*2]

    if len(cellsToPlot.index) > 5000:
        cellsToPlot = cellsToPlot.sample(5000)

    fig = px.scatter(
        cellsToPlot, x="reads", y="genes", color='pass',
        color_discrete_map=constants.QC_COLORMAP,
        template=constants.DEFAULT_FIGURE_STYLE,
        labels={"reads": "Total reads", "genes": "Genes detected"},
        title="Genes Detected Per Cell")

    return dp.Plot(fig)


def buildSaturationScatter(allCells):
    """
    Function to build scatter plot for saturation

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample

    Returns:
        dp.Plot object with the figure
    """
    cellsToPlot = allCells[:allCells['pass'].sum()*2]
    if len(cellsToPlot.index) > 5000:
        cellsToPlot = cellsToPlot.sample(5000)

    fig = px.scatter(
        cellsToPlot, x="reads", y="Saturation",
        labels={"reads": "Total reads"},
        template=constants.DEFAULT_FIGURE_STYLE, color='pass',
        color_discrete_map=constants.QC_COLORMAP, title="Saturation Per Cell",
        log_x=False)

    return dp.Plot(fig)


def makeRankPlot(allCells, writeDir, sampleName, internalReport):
    """
    Function to make a kneeplot using @field in @data; drawing a
    vertical line at @threshold

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name
        internalReport (bool): Whether report is being generated for internal purposes

    Returns:
        dp.Plot object with the figure
    """
    # Make copy to prevent overwriting of barcode indexes
    allCellsCopy = allCells.copy(deep=True)
    allCellsCopy.reset_index(drop=True, inplace=True)

    indices = CalculationUtils.getIndicesToInclude(len(allCellsCopy.index))
    rowsToInclude = allCellsCopy.iloc[indices]

    fig = px.line(
        rowsToInclude, x=rowsToInclude.index, y="umis",
        labels={"index": "Cell barcodes", "umis": "Unique transcript counts"},
        log_x=True, log_y=True, template=constants.DEFAULT_FIGURE_STYLE,
        title="Barcode Rank plot")

    fig.add_vline(
        x=allCellsCopy['pass'].sum(), line_dash="dash", line_color="green")
    fig.update_layout(xaxis_range=[1, math.log10(rowsToInclude.index.max())])

    if internalReport:
        saveFigureAsPng(writeDir, sampleName, fig, "BarcodeRankPlot.png")
    return dp.Plot(fig)


def makeStarDfFromDict(starStats: Dict, category="STARsolo") -> pd.DataFrame:
    """
    Function to make dataframe representing STAR output from a dictionary

    Args:
        starStats (dict): Dictionary containing the output of STAR
        category (str):

    Returns:
        Dataframe with the STAR output
    """
    stats = []

    for (key, value) in starStats.items():
        if key != 'STAR Cells' and key != "Saturation":
            stats.append([key, value])

    statsAsDictList = [
        GeneralUtils.makeTableDict(
            ['Metric', 'Value'], valuePair) for valuePair in stats]

    statsDataFrame = pd.DataFrame(statsAsDictList)
    statsDataFrame['Category'] = category

    return statsDataFrame


def saveFigureAsPng(writeDir, sampleName, fig, pngName: str):
    """
    Function to save a plotly figure object as png

    Args:
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name
        fig (obj): Plotly object containing the figure
        pngName (str): Name of png file to be saved
    """
    fig.write_image(writeDir / f"{sampleName}_figures" / pngName)


def buildBarnyardPage(allCells, writeDir, sampleName) -> dp.Page:
    """
    Function to build a datapane page for barnyard samples

    Args:
       allCells (pd.DataFrame): Dataframe containing all data for this sample

    Returns:
        dp.Page object with the necessary figures
    """
    barnyardStatsDf = makeBarnyardStatsDf(allCells)

    barnyardStatsTable = DatapaneUtils.createTableIgnoreWarning(
        barnyardStatsDf[constants.DISPLAY_COLUMNS].style.pipe(
            GeneralUtils.styleTable, hideColumnHeaders=True,
            boldColumn=['Metric'], title="", numericCols=['Value']))

    barnyardPlot = makeBarnyardScatterPlot(allCells, writeDir, sampleName)
    plotGroup = dp.Group(barnyardStatsTable, barnyardPlot, columns=2)

    return (dp.Page(plotGroup, title="Barnyard"), barnyardStatsDf)


def makeBarnyardScatterPlot(allCells, writeDir, sampleName):
    """
    Function to make a scatter plot for barnyard samples

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name

    Returns:
        dp.Plot object with the figure
    """
    cellsToPlot = allCells[:allCells['pass'].sum()*2]
    if len(cellsToPlot.index) > 5000:
        cellsToPlot = cellsToPlot.sample(5000)

    fig = px.scatter(
        cellsToPlot, x="humanUmis", y="mouseUmis", color="species", log_x=True,
        log_y=True, color_discrete_map=constants.BARNYARD_COLORMAP,
        labels={"humanUmis": "Human UMIs", "mouseUmis": "Mouse UMIs"},
        template=constants.DEFAULT_FIGURE_STYLE)

    fig.update_traces(selector=dict(name='None'), visible="legendonly")
    saveFigureAsPng(writeDir, sampleName, fig, "Barnyard_ScatterPlot.png")

    return dp.Plot(fig)


def makeBarnyardStatsDf(allCells):
    """
    Make dataframe for barnyard samples

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample

    Returns:
        Dataframe with the barnyard stats
    """
    passingBarnyardCells = allCells[
        (allCells['pass']) & (~allCells.species.isin(['Ambiguous', 'None']))]

    ncells = len(passingBarnyardCells.index)
    stats = []

    propHuman = (passingBarnyardCells.species == 'Human').sum() / ncells
    stats.append(['Human Cells', f"{propHuman:.1%}"])

    propMouse = (passingBarnyardCells.species == 'Mouse').sum() / ncells
    stats.append(['Mouse Cells', f"{propMouse:.1%}"])

    propMixed = (passingBarnyardCells.species == 'Mixed').sum() / ncells
    stats.append(['Mixed Cells', f"{propMixed:.1%}"])

    doublets = propMixed/(2*propHuman*propMouse) if propHuman > 0 and propMouse > 0 else 0
    stats.append(
        ['Estimated Doublets',
         f"{min(1, doublets):.1%}"])

    # pep8 workaround
    temp = passingBarnyardCells
    stats.append(
        ['Background',
         f"{temp[temp.species != 'Mixed'].minorFrac.median():.2%}"])

    statsAsDictList = [
        GeneralUtils.makeTableDict(
            ['Metric', 'Value'], valuePair) for valuePair in stats]

    result = pd.DataFrame(statsAsDictList)
    result['Category'] = 'Barnyard'

    return result


def makeUniqueReadsFig(allCells):
    """
    Function to make figure representing unique reads

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample

    Returns:
        dp.Plot object with the figure
    """
    passingCells = allCells[allCells['pass']]
    meanReads = passingCells.reads.mean()
    medReads = passingCells.passingReads.median()
    medUmis = passingCells.umis.median()
    
    unique_reads = []
    read_depths = [100, 500, 1000, 5000, 10000, 20000]

    for d in read_depths:
        target = CalculationUtils.downsampleTargetReadCount(meanReads, medReads, d)
        unique_reads.append(
            CalculationUtils.extrapolate_umis(medReads, medUmis, target))
    df = pd.DataFrame(dict(x=read_depths, unique_read=unique_reads))

    threshold = intOrNan(passingCells.reads.mean())
    plot_points = df.loc[df['x'] < threshold]
    fig = px.line(plot_points, x="x", y="unique_read",
                  labels={"x": "Total reads per cell",
                          "unique_read": "Median unique transcript counts (extrapolated)"},
                  template=constants.DEFAULT_FIGURE_STYLE, markers=True,
                  title="Complexity")

    return dp.Plot(fig), df


def makeCellStatsDf(allCells):
                              
    """
    Function to make dataframe of summary stats

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample
    """
    passingCells = allCells[allCells['pass']]
    stats = []

    stats.append(['Unique Transcript Counts Threshold', passingCells.umis.min()])
    stats.append(['Cells above Threshold', len(passingCells)])
    stats.append(["Mean Reads per cell", intOrNan(passingCells.reads.mean())])
    stats.append(["Median Unique Transcript Counts per cell",
                  intOrNan(passingCells.umis.median())])
    stats.append(["Median Genes per cell",
                  intOrNan(passingCells.genes.median())])
    stats.append(
        ["Reads in Cells",
         f"{passingCells.passingReads.sum()/allCells.passingReads.sum():.1%}"])

    statsAsDictList = [
        GeneralUtils.makeTableDict(
            ['Metric', 'Value'], valuePair) for valuePair in stats]

    statsDataFrame = pd.DataFrame(statsAsDictList)
    statsDataFrame['Category'] = 'Cell Calling'

    return statsDataFrame


def intOrNan(x):
    """
    Function to check whether the passed variable is NaN. If
    it isn't typecast variable to an integer and return

    Args:
        x (obj): Variable to check

    Returns:
        @x if x is NaN, else returns x typecasted as an integer
    """
    if np.isnan(x):
        return x

    return int(x)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sampleName', required=False, default=None)
    parser.add_argument("--libDir", required=False, default=None)
    parser.add_argument("--sample_metrics", required=True)
    parser.add_argument("--internalReport", action="store_true", default=False)
    parser.add_argument("--isBarnyard", action="store_true", default=False)

    args = parser.parse_args()
    buildSampleReport(Path(args.sample_metrics), args.sampleName, args.libDir, args.internalReport, args.isBarnyard)

if __name__ == "__main__":
    main()
