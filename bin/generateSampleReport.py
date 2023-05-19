#!/usr/bin/env python
"""
Generate the sample report (HTML and .csv) from per-sample metrics
"""

import argparse
import json
from pathlib import Path
from typing import Dict
import datapane as dp
import numpy as np
import pandas as pd
import plotly.express as px
import utils.myconstants as constants
from utils.base_logger import logger
from utils.ReportUtil import CalculationUtils, DatapaneUtils, GeneralUtils


def buildSampleReport(sample_metrics, sample, libDir, internalReport, isBarnyard, referencesPath, trim_stats):
    """
    Entry function for building sample report

    Args:
        sample_metrics (str): Path to folder containing csv files defining the metrics for a single sample
        sample (str): Sample name
        libDir (str): Path to library structure json
        internalReport (bool): Flag denoting whether report is for internal purposes
        isBarnyard (bool): Flag denoting whether sample is barnyard or not
        referencesPath (Path): Reference files
        trim_stats (list): Json files containing cutadapt output
    """
    allCells = pd.read_csv(f"{sample_metrics}/allCells.csv", index_col=0)

    with open(f"{sample_metrics}/sample_metrics/barcodesToPlot.json") as f:
        barcodesToPlot = json.load(f)

    writeDir = Path(".", "reports")

    filtered_short_reads, average_trimmed_length = get_trim_stats_metrics(trim_stats)

    Path(writeDir, "csv").mkdir(parents=True)
    
    if internalReport or isBarnyard:
        Path(writeDir, f"{sample}_figures").mkdir()
    
    pages = []

    if internalReport:
        readsPage, statsDf, complexity_df, internalReportBlocks = buildReadsPage(sample, allCells, writeDir, internalReport,
                                                                                 filtered_short_reads, average_trimmed_length)
    else:
        readsPage, statsDf, complexity_df= buildReadsPage(sample, allCells, writeDir, internalReport,
                                                          filtered_short_reads, average_trimmed_length)
    pages.append(readsPage)
    
    with open(referencesPath / libDir) as f:
        libJson = json.load(f)
    if internalReport:
        allIndexPlots, barcodesPage = buildBarcodesPage(referencesPath, allCells, barcodesToPlot,
                                                        writeDir, sample, internalReport, libJson)
    else:
        barcodesPage = buildBarcodesPage(referencesPath, allCells, barcodesToPlot,
                                         writeDir, sample, internalReport, libJson)
    pages.append(barcodesPage)
    
    if internalReport:
        internalReportBlocks.extend([allIndexPlots["Ligation"], allIndexPlots["PCR"]])
        internalReportPage = dp.Page(dp.Group(blocks=internalReportBlocks, columns=2), title="InternalReport")
        pages.append(internalReportPage)
    
    if isBarnyard:
        scoreBarn(allCells)
        (barnyardPage, barnyardStat) = buildBarnyardPage(allCells, writeDir, sample)
        statsDf = pd.concat([statsDf, barnyardStat])
        pages.append(barnyardPage)

    report = dp.Report(blocks=pages)

    logger.debug(f"Writing reportStatistics csv, report to {str(writeDir.resolve())}")

    new_complexity_df = pd.DataFrame(columns=['Metric','Value','Category'])
    for index, row in complexity_df.iterrows():
        new_complexity_df.loc[len(new_complexity_df.index)] = [f"Median unique transcript counts at {row['x']} total reads per cell",
                                                               row["unique_read"], "Extrapolated Complexity"]
    
    statsDf = pd.concat([statsDf, new_complexity_df], axis=0)

    statsDf.to_csv(writeDir / "csv"/ f"{sample}.reportStatistics.csv", index=False, header=False,
                   columns=["Category", "Metric", "Value"])

    report.save(writeDir / f"{sample}.report.html")

def get_trim_stats_metrics(trim_stats):
    """
    Parse cutadapt output and aggregate metrics

    Args:
        trim_stats (list): Path to trim stats json
    
    Returns:
        Filtered short reads and average trimmed length
    """
    filtered_short_reads = 0
    output_read1_bp = 0
    output_read = 0

    for trim_stat in trim_stats:
        with open(trim_stat) as f:
            trim_stat_dict = json.load(f)
            filtered_short_reads += trim_stat_dict["read_counts"]["filtered"]["too_short"]
            output_read1_bp += trim_stat_dict["basepair_counts"]["output_read1"]
            output_read += trim_stat_dict["read_counts"]["output"]
    
    average_trimmed_length = output_read1_bp/output_read

    return filtered_short_reads, average_trimmed_length


def getBackground(minorFracs):
    """
    Calculate median and standard deviation for barnyard cells

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
    Label each cell as belonging to either human or mouse

    Args:
        allCells (pd.DataFrame): Dataframe with all the data for this sample
    """
    cells = allCells
    cells['species'] = "None"
    cells['minorFrac'] = cells[['humanUmis', 'mouseUmis']].min(1) / cells[['humanUmis', 'mouseUmis']].sum(1)
    cells.loc[cells['pass'] & (cells.humanUmis > cells.mouseUmis), 'species'] = 'Human'
    cells.loc[cells['pass'] & (cells.mouseUmis > cells.humanUmis), 'species'] = 'Mouse'

    minHuman = minMouse = cells[cells['pass']].umis.min()
    if (cells.species == 'Human').any() and (cells.species == 'Mouse').any():
        minHuman = np.percentile(cells.humanUmis[cells.species == 'Human'], 10)
        minMouse = np.percentile(cells.mouseUmis[cells.species == 'Mouse'], 10)

    # Labelling low and doublet cells
    cells.loc[(cells.humanUmis < minHuman) & (cells.mouseUmis < minMouse), 'species'] = "None"
    cells.loc[(cells.humanUmis >= minHuman) & (cells.mouseUmis >= minMouse), 'species'] = "Mixed"

    # Labelling as Ambiguous
    humanBackMed, humanBackStd = getBackground(cells[(cells.species == 'Mouse')].minorFrac)
    mouseBackMed, mouseBackStd = getBackground(cells[(cells.species == 'Mouse')].minorFrac)

    cells.loc[(cells.species == 'Mixed') & (cells.humanUmis > cells.mouseUmis) & (cells.minorFrac < mouseBackMed + 2 * mouseBackStd), 'species'] = 'Ambiguous'

    cells.loc[(cells.species == 'Mixed') & (cells.mouseUmis > cells.humanUmis) & (cells.minorFrac < humanBackMed + 2 * humanBackStd), 'species'] = 'Ambiguous'

    cells.loc[(cells.species == 'Human') & (cells.minorFrac >= max(0.1, mouseBackMed + 3 * mouseBackStd)), 'species'] = 'Ambiguous'

    cells.loc[(cells.species == 'Mouse') & (cells.minorFrac >= max(0.1, humanBackMed + 3 * humanBackStd)), 'species'] = 'Ambiguous'


def makeReadStatsDf(allCells, filtered_short_reads):
    """
    Make dataframe depicting read stats that will be displayed in the sample report

    Args:
        allCells (pd.DataFrame): All cells information
        filtered_short_reads (st)
    """
    stats = []
    stats.append({'Category':'Reads', 'Metric': 'Total Sample Reads',
                  'Value': intOrNan(allCells['reads'].sum() + filtered_short_reads)})
    stats.append({'Category':'Reads', 'Metric': 'Passing Sample Reads',
                  'Value': intOrNan(allCells['reads'].sum())})
    stats.append({'Category':'Reads', 'Metric': 'Reads Mapped to Genome',
                  'Value': format(allCells['mappedReads'].sum() / allCells['reads'].sum(), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Reads Mapped to Transcriptome',
                  'Value': format(allCells['geneReads'].sum() / allCells['mappedReads'].sum(), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Exonic Reads',
                  'Value': format(allCells['exonReads'].sum() / allCells['mappedReads'].sum(), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Antisense Reads',
                  'Value': format(allCells['antisenseReads'].sum() / allCells['mappedReads'].sum(), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Mitochondrial Reads',
                  'Value': format(allCells['mitoReads'].sum() / allCells['mappedReads'].sum(), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Saturation',
                  'Value': format(1 - (allCells.umis.sum() / allCells.passingReads.sum()), ".2f")})
    return pd.DataFrame(stats)

def makeExtraReadStatsDf(allCells, average_trimmed_length):
    """
    Create a dataframe with metrics included in .csv files, but not the HTML report tables
    """
    stats = []
    stats.append({'Category':'Reads', 'Metric': 'Average Trimmed Read Length',
                  'Value': format(average_trimmed_length, ".1f")})
    return pd.DataFrame(stats)

def buildReadsPage(sample, allCells, writeDir, internalReport, filtered_short_reads, average_trimmed_length):
    """
    Function to build a datapane page that mainly shows plots related
    to reads

    Args:
        sample (str): Sample name
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        writeDir (path): Path object pointing to write directory
        internalReport (bool): Flag indicating whether report is to be
            generated for internal r&d run
        filtered_short_reads (int): Reads filtered out after trimming
        average_trimmed_length (float): Average read length after trimming

    Returns:
        dp.Page object and concatenated dataframe of various statistics
    """
    internalReportBlocks = []
    statsDf = pd.DataFrame({'Metric': ['SampleName'], 'Value': [sample]})
    statsDf['Category'] = 'Sample'

    readStatsDf = makeReadStatsDf(allCells, filtered_short_reads)
    extraStatsDf = makeExtraReadStatsDf(allCells, average_trimmed_length)
    cellStatsDf = makeCellStatsDf(allCells)

    uniqueReadsFig, complexity_df = makeUniqueReadsFig(allCells, internalReport)

    summaryStatsTable = DatapaneUtils.createTableIgnoreWarning(
        cellStatsDf[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, title="Cell Metrics", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value']))

    mappingTable = DatapaneUtils.createTableIgnoreWarning(
        readStatsDf[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, title="Read Metrics", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value']))

    statsDf = pd.concat([statsDf, readStatsDf, extraStatsDf, cellStatsDf])
    rankPlot = makeRankPlot(allCells, writeDir, sample, internalReport)
    genesVersusUMIsScatter = buildGeneXReadScatter(allCells)
    saturationScatter = buildSaturationScatter(allCells)
    groupBlocks = [rankPlot, uniqueReadsFig, mappingTable, summaryStatsTable, genesVersusUMIsScatter, saturationScatter]
    page = dp.Page(dp.Group(blocks=groupBlocks, columns=2), title="Cells")
    if internalReport:
        allCells['geneReadsToMappedReads'] = allCells['geneReads'] / allCells['mappedReads']
        allCells['antisenseReadsToMappedReads'] = allCells['antisenseReads'] / allCells['mappedReads']
        allCells['exonReadsToMappedReads'] = allCells['exonReads'] / allCells['mappedReads']
        geneReadsDist = buildDist("Distribution of Reads: Gene Matched", 'geneReadsToMappedReads', 'Fraction of Genic Reads', allCells)
        antisenseReadsDist = buildDist("Distribution of Reads: Antisense", "antisenseReadsToMappedReads", "Fraction of Antisense Reads", allCells)
        exonicReadsDist = buildDist("Distribution of Reads: Exonic", 'exonReadsToMappedReads', "Fraction of Exonic Reads", allCells)
        internalReportBlocks.extend([geneReadsDist, exonicReadsDist, antisenseReadsDist])

        anyNonZeroMitoReads = (allCells[allCells['pass']].mitoReads > 0).any()

        if anyNonZeroMitoReads:
            allCells['mitoReadsToMappedReads'] = allCells['mitoReads'] / allCells['mappedReads']
            mitoReadsDist = buildDist("Distribution of Reads: Mitochondrial", 'mitoReadsToMappedReads', "Fraction of Mitochondrial Reads", allCells)
            internalReportBlocks.append(mitoReadsDist)
        return page, statsDf, complexity_df, internalReportBlocks
    else:
        return page, statsDf, complexity_df


def buildBarcodesPage(referencesPath, allCells, barcodesToPlot, writeDir, sampleName, internalReport, libJson):
    """
    Build datapane page for barcodes

    Args:
        referencesPath (Path): Reference files
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        barcodesToPlot (dict): Dictionary containing the barcodes to plot
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name
        internalReport (bool): Flag for internal report
        libJson (dict): Library json

    Returns:
        dp.Page object which holds the relevant figures
    """
    blocksToRender = []
    allIndexPlots = createBarcodeLevelFigures(referencesPath, allCells, barcodesToPlot, writeDir, sampleName, internalReport, libJson)
    blocksToRender.append(allIndexPlots["RT"])
    
    if internalReport:
        return allIndexPlots, dp.Page(blocks=blocksToRender, title="Barcodes")
    else:
        return dp.Page(blocks=blocksToRender, title="Barcodes")


def createBarcodeLevelFigures(referencesPath, allCells, barcodesToPlot, writeDir, sampleName, internalReport, libJson):
    """
    Create plots for the barcodes page

    Args:
        referencesPath (Path): Reference files
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        barcodesToPlot (dict): Dictionary of barcodes to plot
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name
        internalReport (bool): Flag for internal report
        libJson (dict): Library json

    Returns:
        dp.Group object with the relevant figures
    """
    passingCells = allCells[allCells['pass']]
    allIndexPlots = {}
    for key in sorted(barcodesToPlot):
        data = barcodesToPlot[key]
        bcColName = data['alias']
        possibleValues = list(data['possibleValues'])

        indexPlots = DatapaneUtils.barcodeLevelPlots(referencesPath, sampleName, passingCells, possibleValues, f"{bcColName}_alias", 
                                                     f"{bcColName} Index", internalReport, libJson, wellAliases=data['orderAsWellAliases'], writeDir=writeDir)

        allIndexPlots[barcodesToPlot[key]['alias']] = indexPlots

    return allIndexPlots


def buildDist(title, field, replacementStr, allCells):
    """
    Build histogram

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

    fig = px.histogram(passingCells, x=field, title=title, histnorm='percent',
                       template=constants.DEFAULT_FIGURE_STYLE, labels={field: replacementStr})

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

    fig = px.scatter(cellsToPlot, x="reads", y="genes", color='pass', color_discrete_map=constants.QC_COLORMAP, template=constants.DEFAULT_FIGURE_STYLE,
                     labels={"reads": "Total reads", "genes": "Genes detected"}, title="Genes Detected Per Cell")

    return dp.Plot(fig)


def buildSaturationScatter(allCells):
    """
    Build scatter plot for saturation

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample

    Returns:
        dp.Plot object with the figure
    """
    cellsToPlot = allCells[:allCells['pass'].sum()*2]
    if len(cellsToPlot.index) > 5000:
        cellsToPlot = cellsToPlot.sample(5000)

    fig = px.scatter(cellsToPlot, x="reads", y="Saturation", labels={"reads": "Total reads"}, template=constants.DEFAULT_FIGURE_STYLE, color='pass',
                     color_discrete_map=constants.QC_COLORMAP, title="Saturation Per Cell", log_x=False)

    return dp.Plot(fig)


def makeRankPlot(allCells, writeDir, sampleName, internalReport):
    """
    Make a kneeplot using @field in @data; drawing a
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

    fig = px.line(rowsToInclude, x=rowsToInclude.index, y="umis", labels={"index": "Cell barcodes", "umis": "Unique transcript counts"},
                  log_x=True, log_y=True, template=constants.DEFAULT_FIGURE_STYLE, title="Barcode Rank plot")

    fig.add_vline(x=allCellsCopy['pass'].sum(), line_dash="dash", line_color="green")
    
    fig.update_layout(xaxis_range=[1, np.log10(rowsToInclude.index.max()+1)])

    if internalReport:
        saveFigureAsPng(writeDir, sampleName, fig, "BarcodeRankPlot.png")
    
    return dp.Plot(fig)


def saveFigureAsPng(writeDir, sampleName, fig, pngName: str):
    """
    Save a plotly figure object as png

    Args:
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name
        fig (obj): Plotly object containing the figure
        pngName (str): Name of png file to be saved
    """
    fig.write_image(writeDir / f"{sampleName}_figures" / pngName)


def buildBarnyardPage(allCells, writeDir, sampleName) -> dp.Page:
    """
    Build a datapane page for barnyard samples

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

    fig = px.scatter(cellsToPlot, x="humanUmis", y="mouseUmis", color="species", log_x=True, log_y=True, template=constants.DEFAULT_FIGURE_STYLE,
                     color_discrete_map=constants.BARNYARD_COLORMAP, labels={"humanUmis": "Human UMIs", "mouseUmis": "Mouse UMIs"})

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
    passingBarnyardCells = allCells[(allCells['pass']) & (~allCells.species.isin(['Ambiguous', 'None']))]

    ncells = max(1, len(passingBarnyardCells.index))
    stats = []

    propHuman = (passingBarnyardCells.species == 'Human').sum() / ncells
    stats.append(['Human Cells', f"{propHuman:.1%}"])

    propMouse = (passingBarnyardCells.species == 'Mouse').sum() / ncells
    stats.append(['Mouse Cells', f"{propMouse:.1%}"])

    propMixed = (passingBarnyardCells.species == 'Mixed').sum() / ncells
    stats.append(['Mixed Cells', f"{propMixed:.1%}"])

    doublets = propMixed/(2*propHuman*propMouse) if propHuman > 0 and propMouse > 0 else 0
    stats.append(['Estimated Doublets', f"{min(1, doublets):.1%}"])

    stats.append(['Background', f"{passingBarnyardCells[passingBarnyardCells.species != 'Mixed'].minorFrac.median():.2%}"])

    statsAsDictList = [GeneralUtils.makeTableDict(['Metric', 'Value'], valuePair) for valuePair in stats]

    result = pd.DataFrame(statsAsDictList)
    result['Category'] = 'Barnyard'

    return result


def makeUniqueReadsFig(allCells, extrapolate_upwards=False):
    """
    Estimate number of unique transcript counts (unique reads / UMIs)
    at different sequencing depths

    The plot includes a point for the observed values in this sample (not extrapolated)
    but this is excluded from the returned dataframe (For statistics.csv output)

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        extrapolate_upwards: If true, estimate complexity at depth higher than actuallt sequenced 
            (which is generally less accurate)

    Returns:
        dp.Plot object with the figure and dataframe with estimates
    """
    passingCells = allCells[allCells['pass']]
    mean_reads = passingCells.reads.mean()
    median_reads = passingCells.passingReads.median()
    median_unique = passingCells.umis.median()
    
    unique_reads = []
    target_mean_reads = [100, 500, 1000, 5000, 10000, 20000]
    for d in target_mean_reads:
        if d < mean_reads or extrapolate_upwards:
            target_median_reads = CalculationUtils.downsampleTargetReadCount(mean_reads, median_reads, d)
            unique_reads.append(CalculationUtils.extrapolate_umis(median_reads, median_unique, target_median_reads))
        else:
            unique_reads.append(np.nan)
    
    plot_df = pd.DataFrame(dict(x=target_mean_reads + [mean_reads], unique_read=unique_reads + [median_unique]))
    
    plot_df.sort_values(by='x', inplace=True) 
    fig = px.line(plot_df, x="x", y="unique_read",
                  labels={"x": "Total reads per cell",
                          "unique_read": "Median unique transcript counts (extrapolated)"},
                  template=constants.DEFAULT_FIGURE_STYLE, markers=True,
                  title="Complexity")

    stats_df = pd.DataFrame(dict(x=target_mean_reads, unique_read=unique_reads))
    return dp.Plot(fig), stats_df


def makeCellStatsDf(allCells):
    """
    Make dataframe of summary stats

    Args:
        allCells (pd.DataFrame): Dataframe containing all data for this sample
    """
    passingCells = allCells[allCells['pass']]
    stats = []

    stats.append(['Unique Transcript Counts Threshold', passingCells.umis.min()])
    stats.append(['Cells above Threshold', len(passingCells)])
    stats.append(["Mean Passing Reads per Cell", intOrNan(passingCells.reads.mean())])
    stats.append(["Median Unique Transcript Counts per Cell", intOrNan(passingCells.umis.median())])
    stats.append(["Median Genes per Cell", intOrNan(passingCells.genes.median())])
    stats.append(["Reads in Cells", f"{passingCells.passingReads.sum()/allCells.passingReads.sum():.1%}"])

    statsAsDictList = [GeneralUtils.makeTableDict(['Metric', 'Value'], valuePair) for valuePair in stats]

    statsDataFrame = pd.DataFrame(statsAsDictList)
    statsDataFrame['Category'] = 'Cells'

    return statsDataFrame


def intOrNan(x):
    """
    Check whether the passed variable is NaN. If
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
    parser.add_argument("--libJsonName", required=False, default=None)
    parser.add_argument("--sample_metrics", required=True)
    parser.add_argument("--libraryStructPath", required=True)
    parser.add_argument("--trim_stats", nargs='+')
    parser.add_argument("--internalReport", action="store_true", default=False)
    parser.add_argument("--isBarnyard", action="store_true", default=False)

    args = parser.parse_args()
    buildSampleReport(Path(args.sample_metrics), args.sampleName, args.libJsonName, args.internalReport, args.isBarnyard, Path(args.libraryStructPath), args.trim_stats)

if __name__ == "__main__":
    main()

