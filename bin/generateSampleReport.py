#!/usr/bin/env python
"""
Generate the sample report (HTML and .csv) from per-sample metrics
"""

import argparse
import json
from pathlib import Path
from typing import List

import datapane as dp
import numpy as np
import pandas as pd
import plotly.express as px

from utils import fileUtils, reportUtil, statsUtils
from utils.base_logger import logger


def buildSampleReport(sample_metrics:str, sample:str, internalReport:bool, isBarnyard:bool, libStructJson:Path, trim_stats:List[Path]):
    """
    Entry function for building sample report

    Args:
        sample_metrics: Path to folder containing csv files defining the metrics for a single sample
        sample: Sample name
        internalReport: Flag denoting whether report is for internal purposes
        isBarnyard: Flag denoting whether sample is barnyard or not
        libStructJson: Library structure json (in directory with sequence files, etc.)
        trim_stats: json files containing cutadapt output
    """
    allCells = pd.read_csv(f"{sample_metrics}/allCells.csv", index_col=0)
    filtered_short_reads, avg_trimmed_length = get_trim_stats_metrics(trim_stats)

    writeDir = Path(".", "reports")
    Path(writeDir, "csv").mkdir(parents=True, exist_ok=True)
    if internalReport or isBarnyard:
        Path(writeDir, f"{sample}_figures").mkdir()
    
    pages = []

    readsPage, statsDf, complexity_df, internalReportBlocks = \
        buildReadsPage(sample, allCells, writeDir, internalReport, filtered_short_reads, avg_trimmed_length)
    pages.append(readsPage)
    
    libStructDir = libStructJson.parent
    libStruct = json.load(open(libStructJson))
    allIndexPlots, barcodesPage = buildBarcodesPage(libStructDir, allCells, writeDir, sample, internalReport, libStruct)
    pages.append(barcodesPage)
    
    if internalReport:
        internalReportBlocks.extend([allIndexPlots["lig"], allIndexPlots["pcr"]])
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
        new_complexity_df.loc[new_complexity_df.index.size] = [f"Median unique transcript counts at {row['x']} total reads per cell",
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
        Number of filtered short reads and average trimmed read-length
        Returns NaN if no trimStats given
    """
    filtered_short_reads = 0
    output_read1_bp = 0
    output_read = 0
    if not trim_stats:
        return np.nan, np.nan
    for trim_stat in trim_stats:
        with open(trim_stat) as f:
            trim_stat_dict = json.load(f)
            filtered_short_reads += trim_stat_dict["read_counts"]["filtered"]["too_short"]
            output_read1_bp += trim_stat_dict["basepair_counts"]["output_read1"]
            output_read += trim_stat_dict["read_counts"]["output"]
    average_trimmed_length = (output_read1_bp/output_read) if output_read > 0 else np.nan
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
    if minorFracs.size == 0:
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


def makeReadStatsDf(allCells:pd.DataFrame, filtered_short_reads:float):
    """
    Make dataframe depicting read stats that will be displayed in the sample report

    Args:
        allCells: Per-cell metrics
        filtered_short_reads: Number of filtered reads (or NaN)
    """
    stats = []
    stats.append({'Category':'Reads', 'Metric': 'Total Sample Reads',
                  'Value': intOrNan(allCells['reads'].sum() + filtered_short_reads)})
    stats.append({'Category':'Reads', 'Metric': 'Passing Sample Reads',
                  'Value': intOrNan(allCells['reads'].sum())})
    stats.append({'Category':'Reads', 'Metric': 'Reads Mapped to Genome',
                  'Value': format(np.divide(allCells['mappedReads'].sum(), allCells['reads'].sum()), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Reads Mapped to Transcriptome',
                  'Value': format(np.divide(allCells['geneReads'].sum(), allCells['mappedReads'].sum()), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Exonic Reads',
                  'Value': format(np.divide(allCells['exonReads'].sum(), allCells['mappedReads'].sum()), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Antisense Reads',
                  'Value': format(np.divide(allCells['antisenseReads'].sum(), allCells['mappedReads'].sum()), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Mitochondrial Reads',
                  'Value': format(np.divide(allCells['mitoReads'].sum(), allCells['mappedReads'].sum()), ".1%")})
    stats.append({'Category':'Reads', 'Metric': 'Saturation',
                  'Value': format(1 - (np.divide(allCells.umis.sum(), allCells.passingReads.sum())), ".2f")})
    return pd.DataFrame(stats)

def makeExtraReadStatsDf(allCells, average_trimmed_length):
    """
    Create a dataframe with metrics included in .csv files, but not the HTML report tables
    """
    stats = []
    stats.append({'Category':'Reads', 'Metric': 'Average Trimmed Read Length',
                  'Value': format(average_trimmed_length, ".1f")})
    return pd.DataFrame(stats)

def buildReadsPage(sample:str, allCells:pd.DataFrame, writeDir:Path, internalReport:bool, filtered_short_reads:float, avg_trimmed_length:float):
    """
    Function to build a datapane page that mainly shows plots related
    to reads

    Args:
        sample: Sample name
        allCells: Per-cell metrics
        writeDir: Output directory
        internalReport: Flag indicating whether report is to be
            generated for internal r&d run
        filtered_short_reads: Number of reads filtered out after trimming (or NaN)
        avg_trimmed_length: Average read length after trimming (or NaN)

    Returns:
        dp.Page object and concatenated dataframe of various statistics
        List of elements for the internal report (or empty)
    """

    statsDf = pd.DataFrame({'Metric': ['SampleName'], 'Value': [sample]})
    statsDf['Category'] = 'Sample'
    readStatsDf = makeReadStatsDf(allCells, filtered_short_reads)
    extraStatsDf = makeExtraReadStatsDf(allCells, avg_trimmed_length)
    cellStatsDf = makeCellStatsDf(allCells)
    statsDf = pd.concat([statsDf, readStatsDf, extraStatsDf, cellStatsDf])

    uniqueReadsFig, complexity_df = makeUniqueReadsFig(allCells, internalReport)
    summaryStatsTable = reportUtil.createMetricTable(cellStatsDf, title="Cell Metrics")
    mappingTable = reportUtil.createMetricTable(readStatsDf, title="Read Metrics")

    rankPlot = makeRankPlot(allCells, writeDir, sample, internalReport)
    genesVersusUMIsScatter = buildGeneXReadScatter(allCells)
    saturationScatter = buildSaturationScatter(allCells)
    groupBlocks = [rankPlot, uniqueReadsFig, mappingTable, summaryStatsTable, genesVersusUMIsScatter, saturationScatter]
    page = dp.Page(dp.Group(blocks=groupBlocks, columns=2), title="Cells")

    internalReportBlocks = []
    if internalReport:
        allCells['geneReadsToMappedReads'] = allCells['geneReads'] / allCells['mappedReads']
        allCells['antisenseReadsToMappedReads'] = allCells['antisenseReads'] / allCells['mappedReads']
        allCells['exonReadsToMappedReads'] = allCells['exonReads'] / allCells['mappedReads']
        geneReadsDist = buildDist("Distribution of Reads: Gene Matched", 'geneReadsToMappedReads', 'Fraction of Genic Reads', allCells)
        antisenseReadsDist = buildDist("Distribution of Reads: Antisense", "antisenseReadsToMappedReads", "Fraction of Antisense Reads", allCells)
        exonicReadsDist = buildDist("Distribution of Reads: Exonic", 'exonReadsToMappedReads', "Fraction of Exonic Reads", allCells)
        internalReportBlocks.extend([geneReadsDist, exonicReadsDist, antisenseReadsDist])

        anyNonZeroMitoReads = (allCells.loc[allCells['pass']].mitoReads > 0).any()
        if anyNonZeroMitoReads:
            allCells['mitoReadsToMappedReads'] = allCells['mitoReads'] / allCells['mappedReads']
            mitoReadsDist = buildDist("Distribution of Reads: Mitochondrial", 'mitoReadsToMappedReads', "Fraction of Mitochondrial Reads", allCells)
            internalReportBlocks.append(mitoReadsDist)

    return page, statsDf, complexity_df, internalReportBlocks


def buildBarcodesPage(referencesPath, allCells, writeDir, sampleName, internalReport, libJson):
    """
    Build datapane page for barcodes

    Args:
        referencesPath (Path): Reference files
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name
        internalReport (bool): Flag for internal report
        libJson (dict): Library json

    Returns:
        allIndexPlots for internalReports (empty otherwise)
        dp.Page object with shared figures figures
    """
    blocksToRender = []
    allIndexPlots = createBarcodeLevelFigures(referencesPath, allCells, writeDir, sampleName, internalReport, libJson)
    blocksToRender.append(allIndexPlots["rt"])
    return allIndexPlots, dp.Page(blocks=blocksToRender, title="Barcodes")


def createBarcodeLevelFigures(referencesPath, allCells, writeDir, sampleName, internalReport, libJson):
    """
    Create plots for the barcodes page

    Args:
        referencesPath (Path): Reference files
        allCells (pd.DataFrame): Dataframe containing all data for this sample
        writeDir (path): Path object pointing to the write directory
        sampleName (str): Sample name
        internalReport (bool): Flag for internal report
        libJson (dict): Library json

    Returns:
        dp.Group object with the relevant figures
    """
    passingCells = allCells.loc[allCells['pass']] # .loc needed here to keep the column names if zero passing cells

    allIndexPlots = {}
    for bc in libJson['barcodes']:
        if bc.get('type') not in ['library_index', 'umi']:
            alias = bc.get('alias') or bc['name']
            indexPlots = reportUtil.barcodeLevelPlots(referencesPath, sampleName, passingCells, f"{alias}_alias", 
                                                      f"{alias} Index", internalReport, libJson, writeDir=writeDir)
            allIndexPlots[bc['name']] = indexPlots
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
    passingCells = allCells.loc[allCells['pass']]

    fig = px.histogram(passingCells, x=field, title=title, histnorm='percent',
                       template=reportUtil.DEFAULT_FIGURE_STYLE, labels={field: replacementStr})

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

    if cellsToPlot.index.size > 5000:
        cellsToPlot = cellsToPlot.sample(5000)

    fig = px.scatter(cellsToPlot, x="reads", y="genes", color='pass', color_discrete_map=reportUtil.QC_COLORMAP, template=reportUtil.DEFAULT_FIGURE_STYLE,
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
    if cellsToPlot.index.size > 5000:
        cellsToPlot = cellsToPlot.sample(5000)

    fig = px.scatter(cellsToPlot, x="reads", y="Saturation", labels={"reads": "Total reads"}, template=reportUtil.DEFAULT_FIGURE_STYLE, color='pass',
                     color_discrete_map=reportUtil.QC_COLORMAP, title="Saturation Per Cell", log_x=False)

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
    indices = reportUtil.sparseLogCoords(allCells.index.size)
    fig = px.line(allCells.iloc[indices], x=indices, y="umis", labels={"x": "Cell barcodes", "umis": "Unique transcript counts"},
                  log_x=True, log_y=True, template=reportUtil.DEFAULT_FIGURE_STYLE, title="Barcode Rank Plot")

    fig.add_vline(x=max(1, allCells['pass'].sum()), line_dash="dash", line_color="green")
    
    if indices.size > 0:
        fig.update_layout(xaxis_range=[1, np.log10(max(indices)+1)])

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

    barnyardStatsTable = reportUtil.createMetricTable(barnyardStatsDf, title="")
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
    if cellsToPlot.index.size > 5000:
        cellsToPlot = cellsToPlot.sample(5000)

    fig = px.scatter(cellsToPlot, x="humanUmis", y="mouseUmis", color="species", log_x=True, log_y=True, template=reportUtil.DEFAULT_FIGURE_STYLE,
                     color_discrete_map=reportUtil.BARNYARD_COLORMAP, labels={"humanUmis": "Human UMIs", "mouseUmis": "Mouse UMIs"})

    fig.update_traces(selector=dict(name='None'), visible="legendonly")
    saveFigureAsPng(writeDir, sampleName, fig, "Barnyard_ScatterPlot.png")

    return dp.Plot(fig)


def makeBarnyardStatsDf(allCells: pd.DataFrame) -> pd.DataFrame:
    """
    Compute statistics for barnyard samples

    Args:
        allCells: Per-cell metrics

    Returns:
        Dataframe with the barnyard stats
    """
    passingBarnyardCells = allCells.loc[(allCells['pass']) & (~allCells.species.isin(['Ambiguous', 'None']))]
    ncells = max(1, passingBarnyardCells.index.size)

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

    result = pd.DataFrame(stats, columns=['Metric', 'Value'])
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
    passingCells = allCells.loc[allCells['pass']]
    mean_reads = passingCells.reads.mean()
    median_reads = passingCells.passingReads.median()
    median_unique = passingCells.umis.median()
    
    unique_reads = []
    target_mean_reads = [100, 500, 1000, 5000, 10000, 20000]
    for d in target_mean_reads:
        if d < mean_reads or extrapolate_upwards:
            # Target reads refers to the mean reads per cell, but we are estimating
            # UMIs / complexity in the median cell. Approx. scaling median reads by the same
            # factor as mean reads.
            target_median_reads = (d/mean_reads) * median_reads
            unique_reads.append(statsUtils.extrapolate_unique(median_reads, median_unique, target_median_reads))
        else:
            unique_reads.append(np.nan)
    
    plot_df = pd.DataFrame(dict(x=target_mean_reads + [mean_reads], unique_read=unique_reads + [median_unique]))
    
    plot_df.sort_values(by='x', inplace=True) 
    fig = px.line(plot_df, x="x", y="unique_read",
                  labels={"x": "Total reads per cell",
                          "unique_read": "Median unique transcript counts (extrapolated)"},
                  template=reportUtil.DEFAULT_FIGURE_STYLE, markers=True,
                  title="Complexity")

    stats_df = pd.DataFrame(dict(x=target_mean_reads, unique_read=unique_reads))
    return dp.Plot(fig), stats_df


def makeCellStatsDf(allCells: pd.DataFrame) -> pd.DataFrame:
    """
    Make dataframe of summary stats

    Args:
        allCells: Per-cell metrics
    """
    passingCells = allCells.loc[allCells['pass']]
    stats = []
    stats.append(['Unique Transcript Counts Threshold', passingCells.umis.min()])
    stats.append(['Cells above Threshold', passingCells.index.size])
    stats.append(["Mean Passing Reads per Cell", intOrNan(passingCells.reads.mean())])
    stats.append(["Median Unique Transcript Counts per Cell", intOrNan(passingCells.umis.median())])
    stats.append(["Median Genes per Cell", intOrNan(passingCells.genes.median())])
    stats.append(["Reads in Cells", f"{np.divide(passingCells.passingReads.sum(), allCells.passingReads.sum()):.1%}"])

    statsDataFrame = pd.DataFrame(stats, columns=["Metric","Value"])
    statsDataFrame['Category'] = 'Cells'
    return statsDataFrame


def intOrNan(x: int|float|None) -> int|float:
    """
    NaN if @x is NaN or None, else int(x)
    """
    if x is None or np.isnan(x):
        return np.nan
    return int(x)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sampleName', required=False, default=None)
    parser.add_argument("--sample_metrics", required=True, type=Path)
    parser.add_argument("--libraryStruct", required=True, type=Path)
    parser.add_argument("--trim_stats", nargs='+', type=Path)
    parser.add_argument("--internalReport", action="store_true", default=False)
    parser.add_argument("--isBarnyard", action="store_true", default=False)

    args = parser.parse_args()
    buildSampleReport(args.sample_metrics, args.sampleName, args.internalReport, args.isBarnyard, args.libraryStruct, args.trim_stats)

if __name__ == "__main__":
    main()

