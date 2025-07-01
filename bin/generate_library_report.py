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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
from scale_utils import reporting
from scale_utils.base_logger import logger
from scale_utils.lib_json_parser import LibJsonParser


def buildLibraryReport(
    libName: str,
    lib_json_obj: LibJsonParser,
    demuxJson: Path,
    libMetrics: Path,
    internalReport: bool,
    scalePlexLib: bool,
    outDir: str,
    ultima: bool,
    minPassingSampleReads: int,
    minDivergence: float,
):
    """
    Build the library report by calling relevant functions

    Args:
        libName: Library name
        lib_json_obj: Object containing library structure information
        libMetrics: Path to library metrics
        internalReport: Flag denoting whether report is for internal purposes
        scalePlexLib: Flag denoting whether library is a ScalePlex library
        outDir: Output directory
        ultima: Flag denoting whether sequencing data was generated on an Ultima instrument
        minPassingSampleReads: Minimum number of reads a sample must have post barcode
                               demux to be included in the report
        minDivergence: threshold KL divergence score to pass bead
    """
    # need barcode cols for plate plots
    barcode_cols = [
        barcode.get("alias") or barcode["name"]
        for barcode in lib_json_obj.json_contents["barcodes"]
        if barcode.get("plate")
    ]
    allCellsBetweenFiles = pd.read_parquet(
        libMetrics / "allCellsBetweenFiles.parquet", columns=["sample", "counts"] + barcode_cols
    )
    # Within each sample sort by counts
    allCellsBetweenFiles.sort_values(by=["sample", "counts"], ascending=False, inplace=True)
    # Set index to reflect the rank based on umi sorted order
    # cumcount() returns the number of occurrences of each value up to that point
    # So each sample will have an index range starting at zero to the number of barcodes - 1 in that sample
    allCellsBetweenFiles.set_index(allCellsBetweenFiles.groupby("sample").cumcount(), inplace=True)
    if demuxJson:
        demuxMetrics = json.load(open(demuxJson))
    else:
        demuxMetrics = {}

    # Reports will be written to the <outDir> folder
    writeDir = Path(".", outDir)
    Path(writeDir, "csv").mkdir(parents=True, exist_ok=True)
    pages = []

    # Don't make multisample knee plot when scaleplex library
    # Don't make Reads Per Sample plot when ultima=True
    if ultima and scalePlexLib:
        logger.debug("Skip making reads page for scaleplex library sequenced on ultima")
    # Make only multisample knee plot when ultima=True
    elif ultima and not scalePlexLib:
        multiSampleKneePlot = makeMultisampleKneePlot(allCellsBetweenFiles)
        readsPage = dp.Page(
            blocks=[dp.Text(f"## libName: {libName}"), dp.Group(multiSampleKneePlot, columns=1)], title="Reads"
        )
        pages.append(readsPage)
    # Not ultima
    else:
        overallPassingStats, readsPage, barcodeReadStatsInternal = buildReadsPage(
            demuxMetrics, allCellsBetweenFiles, libName, scalePlexLib, minPassingSampleReads
        )
        pages.append(readsPage)

    barcodeStatsDf, cellsPage, barcodeTypeStats = buildBarcodesPage(
        demuxMetrics, libName, lib_json_obj, scalePlexLib, allCellsBetweenFiles, writeDir, ultima
    )
    pages.append(cellsPage)

    bead_bc_group = None
    library_beads_path = libMetrics / "libraryBeads.parquet"
    if library_beads_path.exists():
        Path(writeDir, "figures_internal").mkdir(exist_ok=True)
        beads = pd.read_parquet(library_beads_path)
        bead_scores = pd.read_parquet("bead_scores.parquet")
        bead_reports = make_bead_report(beads, bead_scores, writeDir, libName, minDivergence)
        bead_bc_group = dp.Group(dp.Text("## Bead BC"), dp.Group(blocks=bead_reports, columns=2))

    if internalReport:
        internal_report_page = dp.Page(
            dp.Group(
                dp.Group(dp.Text("## Read Metrics"), barcodeReadStatsInternal) if not ultima else dp.HTML("&nbsp;"),
                dp.Group(dp.Text("## Barcode Metrics"), barcodeTypeStats) if not ultima else dp.HTML("&nbsp;"),
                # bead bc group not present in 3-level assay
                bead_bc_group if bead_bc_group else dp.HTML("&nbsp;"),
                columns=1,
            ),
            title="InternalReport",
        )
        pages.append(internal_report_page)

    # Build a report object from the concatenated pages
    if len(pages) > 1:
        report = dp.Report(blocks=pages)
    # True when ultima and scaleplex
    else:
        report = dp.Report(blocks=pages[0].blocks)

    prefix = f"library_{libName}"

    if not ultima:
        # Write information related to each barcoding level to a csv file
        barcodeStatsDf.to_csv(writeDir / "csv" / f"{prefix}.typeLevelMatches.csv", index=False)
        # Write dataframe that has the number of reads that passed and the number
        # of reads that have different errors
        overallPassingStats.to_csv(writeDir / "csv" / f"{prefix}.overallMatches.csv", index=False)

    report.save(writeDir / f"{prefix}.report.html")


def buildDfForPlatePlot(allCellsBetweenFiles: pd.DataFrame, lib_json_obj: LibJsonParser, bcName: str):
    """
    Construct dataframe that will be used for plotting heatmap

    Args:
        allCellsBetweenFiles: All cells information for this library
        lib_json_obj: Object containing library structure information
        bcName: name of the barcode in lib.json

    Returns:
        Constructed dataframe
    """

    for bc in lib_json_obj.json_contents["barcodes"]:
        if bc["name"] == bcName:
            bcInfo = bc
            break
    else:
        raise ValueError(f"Unknown barcode {bcName}")

    alias = bcInfo.get("alias") or bcInfo["name"]
    well_list = allCellsBetweenFiles[alias].to_list()
    max_letter, max_number = lib_json_obj.getMaxWellNumberAndLetter(lib_json_obj.parent_dir / f'{bcInfo["sequences"]}')
    allCellsBetweenFiles[f"{alias.lower()}_well"] = well_list
    well_df = pd.DataFrame(
        0, columns=range(1, max_number + 1), index=reporting.getCharacterIndices(65, ord(max_letter) + 1)
    )

    for well in set(well_list):
        letter = well[-1]
        numbers = well[:-1]
        # Get umi count for each well
        well_df.at[letter, int(numbers)] = allCellsBetweenFiles.loc[
            allCellsBetweenFiles[f"{alias.lower()}_well"] == well, "counts"
        ].sum()

    return well_df


def buildBarcodesPage(
    demuxJson: Dict,
    libName: str,
    lib_json_obj: LibJsonParser,
    scalePlexLib: bool,
    allCellsBetweenFiles: pd.DataFrame,
    writeDir: Path,
    ultima: bool,
):
    """
    Function that builds a datapane page that depicts a table with barcode
    related information and plate plots that show reads per rt well, reads per
    pcr well and reads per ligation well. Also create a dataframe with
    information related to each barcoding level

    Args:
        demuxJson: Dictionary with the demuxed metrics
        libName: Library name
        lib_json_obj: Object containing library structure information
        scalePlexLib: True if a ScalePlex library
        allCellsBetweenFiles: All cell information for this library
        writeDir: Write directory
        ultima: True if sequencing data was generated on an Ultima instrument

    Returns:
        Dataframe with stats related to barcode and dp.Pageobject
    """
    if ultima:
        barcodeTypeStatsDf = None
        barcodeTypeStats = None
    else:
        (barcodeTypeStatsDf, barcodeTypeStats) = createBarcodeTypeMetricsTables(demuxJson)
    blocks = [dp.Text(f"## libName: {libName}")]
    y_axis_label = "Unique ScalePlex Counts" if scalePlexLib else "Unique Transcript Counts"
    for barcode in lib_json_obj.json_contents["barcodes"]:
        if not barcode.get("plate"):
            continue
        df = buildDfForPlatePlot(allCellsBetweenFiles, lib_json_obj, barcode["name"])
        df.to_csv(writeDir / "csv" / f"library_{libName}_unique_reads_{barcode['alias']}_well.csv")
        plate_plot = reporting.buildPlatePlot(df, f"{barcode['alias']} Plate", 10000.0, y_axis_label)
        blocks.append(plate_plot)
    return (barcodeTypeStatsDf, dp.Page(blocks=blocks, title="Barcodes"), barcodeTypeStats)


def createBarcodeTypeMetricsTables(demuxMetrics: Dict):
    """
    Create dataframe for barcode type and create a datapane
    object for storing a table created with the statistics in the dataframe

    Args:
        demuxMetrics: bcParser metrics (from metrics.json)

    Returns:
        Dataframe and dp.Group object
    """
    barcodesDf = buildDfFromJSONDict(demuxMetrics["barcodes"], "Barcode", "dict")
    tableGroup = []
    allBarcodes = list(barcodesDf["Barcode"].unique())
    for bc in allBarcodes:
        subset = barcodesDf[barcodesDf["Barcode"] == bc][["Match", "Reads"]]
        styledDf = subset.style.pipe(reporting.styleTable, title=f"{bc} Barcodes")
        table = reporting.make_table(styledDf)
        tableGroup.append(table)
    return (barcodesDf, dp.Group(blocks=tableGroup, columns=2))


def buildReadsPage(demuxJson, allCellsBetweenFiles, libName, scalePlexLib, minPassingSampleReads):
    """
    Function to build a datapane page for reads

    Args:
        demuxJson (dict): Dictionary with demuxed metrics
        allCellsBetweenFiles (pd.DataFrame): Dataframe containing data from all samples
        libName (str): Library name
        scalePlexLib (bool): True if ScalePlex library
        minPassingSampleReads (int): Minimum number of reads a sample must have post barcode demux\
                                     to be included in the report

    Returns:
        Dataframe with barcode reads information and dp.Page object
    """
    if not scalePlexLib:
        multiSampleKneePlot = makeMultisampleKneePlot(allCellsBetweenFiles)

    barcodeReadsData = demuxJson["reads"]

    barcodeReadsPerc = buildDfFromJSONDict(barcodeReadsData, "Type", "list")
    barcodeReadsTotal = buildDfFromJSONDict(barcodeReadsData, "Type", "list", 0)
    barcodeReadsTotal = barcodeReadsTotal[["Type", "Reads"]]
    barcodeReadsTotal.rename(columns={"Type": "Status"}, inplace=True)
    barcodeReadsTotal["Percent"] = barcodeReadsPerc["Reads"]
    table_title = "ScalePlex Read Status" if scalePlexLib else "Read Status"
    barcodeReadsTotalStyledInternal = barcodeReadsTotal.style.pipe(
        reporting.styleTable, title=table_title, numericCols=["Reads"]
    )

    total_reads = 0
    total_percent = 0
    too_short_error_reads = 0
    too_short_error_perc = 0
    df_error = barcodeReadsTotal[barcodeReadsTotal["Status"].str.contains("Error")]
    for idx, row in df_error.iterrows():
        if row["Status"] == "TooShortError":
            too_short_error_reads += int(row["Reads"])
            too_short_error_perc += float(row["Percent"][:-1])
        else:
            total_reads += int(row["Reads"])
            total_percent += float(row["Percent"][:-1])
    df_pass = barcodeReadsTotal[barcodeReadsTotal["Status"].str.contains("Pass")]
    # Include TooShortError in BarcodePass
    df_pass_percent = str(round(float(df_pass.loc[0, "Percent"].replace("%", "")) + too_short_error_perc, 1)) + "%"
    df_pass.loc[0, "Percent"] = df_pass_percent
    df_pass.at[df_pass.index[0], "Status"] = "Barcode Pass"
    df_pass.loc[len(df_pass.index) + 1] = ["Barcode Error", total_reads, str(round(total_percent, 1)) + "%"]
    if scalePlexLib:
        barcodeHashReads = demuxJson["barcodes"]["scaleplex"]
        hash_error_reads = barcodeHashReads["Ambiguous"][0] + barcodeHashReads["NoMatch"][0]
        hash_error_perc = float(barcodeHashReads["Ambiguous"][1].strip("%")) + float(
            barcodeHashReads["NoMatch"][1].strip("%")
        )
        hash_error = {
            "Status": "ScalePlex Error",
            "Reads": hash_error_reads,
            "Percent": f"{hash_error_perc:.1f}%",
        }
        df_pass = df_pass.append(hash_error, ignore_index=True)
    barcodeReadsTotalStyled = df_pass.style.pipe(reporting.styleTable, table_title, numericCols=["Reads"])

    barcodeReadStats = reporting.make_table(barcodeReadsTotalStyled)
    barcodeReadStatsInternal = reporting.make_table(barcodeReadsTotalStyledInternal)

    (countsPerSampleDf, rtCountsPerSampleDf) = buildDfFromDemuxSampleMetrics(demuxJson)

    wellOrder = sorted(list(rtCountsPerSampleDf["rtWell"].unique()), key=functools.cmp_to_key(reporting.wellStringComp))
    rtCountsPerSampleDf["rtWell"] = pd.Categorical(rtCountsPerSampleDf["rtWell"], wellOrder)
    rtCountsPerSampleDf.sort_values(by=["rtWell"], inplace=True, ascending=False)

    sampleOrder = list(rtCountsPerSampleDf.Sample.unique())
    sampleOrder.reverse()
    for sample in countsPerSampleDf["Sample"]:  # Add samples with no reads in 'rtCounts'
        if sample not in sampleOrder and sample != "Unknown":
            sampleOrder.append(sample)
    sampleOrder.append("Unknown")
    countsPerSampleDf["Sample"] = pd.Categorical(countsPerSampleDf.Sample, sampleOrder)
    countsPerSampleDf.sort_values(by=["Sample"], inplace=True)
    # Add in html magic to ensure that sample names are read in as string due to plotly express
    # not recognizing numeric sample names as strings
    countsPerSampleDf["Sample"] = ["<span>" + elem + "</span>" for elem in countsPerSampleDf["Sample"]]
    colorMap = matchColorsToNames(list(countsPerSampleDf["Sample"].unique()))
    readsPerSample = px.bar(
        countsPerSampleDf,
        x="Sample",
        y="TotalReads",
        color="Sample",
        height=900,
        color_discrete_map=colorMap,
        template=reporting.DEFAULT_FIGURE_STYLE,
        title="Reads Per Sample",
        labels={"TotalReads": "Total Reads"},
    )
    readsPerSample.update_layout(showlegend=False)
    error_msg = ""
    sample_read_data = demuxJson["samples"]
    for sample in sample_read_data:
        if sample == "Unknown":
            continue
        if sample_read_data[sample]["passingreads"] < minPassingSampleReads:
            error_msg += (
                f"Sample {sample} has fewer than {minPassingSampleReads} passing reads post barcode demux. "
                "No results are generated for this sample. \n"
            )
    if not scalePlexLib:
        if error_msg:
            html_content = f"""
            <h2 style="color:red;">WARNING!</h2>
            <p>{error_msg}</p>
            """
            error_box = dp.Text(html_content)
            first_group = dp.Group(multiSampleKneePlot, error_box, barcodeReadStats, columns=2)
        else:
            first_group = dp.Group(multiSampleKneePlot, barcodeReadStats, columns=2)
    else:
        first_group = dp.Group(barcodeReadStats, columns=1)

    readsPage = dp.Page(
        blocks=[dp.Text(f"## libName: {libName}"), first_group, dp.Group(readsPerSample, columns=1)], title="Reads"
    )
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
    colorMap = {"Unknown": "rgb(179, 188, 201)"}

    if "Unknown" in names:
        names.remove("Unknown")

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

    for sampleName, sampleDict in demuxJson["samples"].items():
        readCount = sampleDict["reads"][0]
        totalCounts.append({"Sample": sampleName, "TotalReads": readCount})
        rtBarcodeCounts = sampleDict["barcodes"]

        if rtBarcodeCounts:
            for rtWell, metrics in rtBarcodeCounts.items():
                rtCounts.append({"Sample": sampleName, "rtWell": rtWell, "ReadCount": metrics["reads"]})

    return pd.DataFrame(totalCounts), pd.DataFrame(rtCounts)


def makeMultisampleKneePlot(allCellsBetweenFiles):
    """
    Makes a kneeplot using @field in @data; drawing a vertical line at @threshold

    Args:
        allCellsBetweenFiles (pd.DataFrame): Dataframe containing data from all samples

    Returns:
        dp.Plot object
    """
    maxIndex = max(allCellsBetweenFiles.index) if len(allCellsBetweenFiles.index) else 0
    indices = set(reporting.sparseLogCoords(maxIndex))
    plottingDf = allCellsBetweenFiles[allCellsBetweenFiles.index.isin(indices)]
    fig = px.line(
        plottingDf,
        x=plottingDf.index,
        y="counts",
        color="sample",
        log_x=True,
        log_y=True,
        template=reporting.DEFAULT_FIGURE_STYLE,
        labels={"index": "Cell Barcodes", "counts": "Unique Transcript Counts"},
    )
    return dp.Plot(fig)


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

    for keyName, valueObj in jsonDict.items():
        if valueDataType == "dict":

            for key, value in valueObj.items():
                newDict = {}
                newDict[name] = keyName
                newDict["Match"] = key
                newDict["Reads"] = value[choiceIndex]
                dictList.append(newDict)

        elif valueDataType == "list":
            newDict = {"Reads": valueObj[choiceIndex]}
            newDict[name] = keyName
            dictList.append(newDict)

        elif valueDataType == "int":
            newDict = {"Cells": valueObj}
            newDict[name] = keyName
            dictList.append(newDict)

    return pd.DataFrame(dictList)


def make_bead_report(
    all_beads: pd.DataFrame,
    bead_scores: pd.DataFrame,
    write_dir: Path,
    lib_name: str,
    min_divergence: float,
) -> list[dp.Plot]:
    """Bead rank plot, histogram and bead statistics

    Args:
        bead_df: Bead-level metrics
        bead_scores: KL divergence scores for beads with > minUTC counts for > 1 RT
        write_dir: Output directory
        lib_name: Library name for file naming
        min_divergence: Minimum KL divergence score to consider a bead as non-ambient

    Returns:
        Plot to include on internal report tab
    """
    if all_beads.empty:
        return []
    all_beads = all_beads.sort_values("counts", ignore_index=True, ascending=False)
    all_beads["cell"] = (all_beads["pass"] > 0).rolling(25, min_periods=1).mean()
    indices = reporting.sparseLogCoords(all_beads.index.size)
    fig = px.scatter(
        all_beads.iloc[indices],
        x=indices,
        y="counts",
        color="cell",
        labels={"x": "Bead barcodes", "counts": "Unique transcript counts"},
        log_x=True,
        log_y=True,
        template=reporting.DEFAULT_FIGURE_STYLE,
        title="Bead Rank Plot",
        opacity=0.5,
        color_continuous_scale=["rgb(178, 181, 187)", "rgb(39, 139, 176)"],
        hover_data={"cell": True},
    )
    if indices.size > 0:
        fig.update_layout(xaxis_range=[1, np.log10(max(indices) + 1)])
    fig.write_image(write_dir / "figures_internal" / f"{lib_name}_BeadRankPlot.png")
    bead_rank_plot = dp.Plot(fig)

    cell_beads = all_beads[all_beads["pass"] > 0]
    min_bead_frac = 0.05  # Beads with >= 5% of median cell bead counts are 'top beads' (could be empty / no cell)
    bead_count_thres = int(round(cell_beads["counts"].median() * min_bead_frac)) if not cell_beads.empty else 100
    top_beads = all_beads[all_beads["counts"] >= bead_count_thres]
    # calculate histogram with matplotlib
    max_pass = top_beads["pass"].max() if not top_beads.empty else 0
    counts, bins, _ = plt.hist(top_beads["pass"], bins=range(max_pass + 2))

    yaxis_title = "Number of beads"

    # plot as bar chart in plotly so raw dataset is not saved to HTML report
    fig = px.bar(
        x=bins[:-1],
        y=counts,
        text_auto=True,
        title=f"Number of Passing Cells per Bead<br><sup>Beads with >= {bead_count_thres} counts </sup>",
        template=reporting.DEFAULT_FIGURE_STYLE,
        labels={"x": "Number of passing cells", "y": yaxis_title},
    )
    fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
    fig.update_layout(
        yaxis_title=yaxis_title,
        bargap=0.1,
    )
    fig.write_image(write_dir / "figures_internal" / f"{lib_name}_BeadPassHist.png")
    bead_pass_hist = dp.Plot(fig)

    bead_score_dist = dp.HTML("&nbsp;")
    bead_stats = []
    if not bead_scores.empty:
        num_filtered_bead = bead_scores[bead_scores["kl_norm"] < min_divergence].shape[0]
        bead_stats.append(("Beads below minimum KL divergence score", num_filtered_bead))
        counts, bins, _ = plt.hist(bead_scores["kl_norm"], weights=bead_scores["bead_reads"], bins=100)

        # plot as bar chart in plotly so raw dataset is not saved to HTML report
        fig = px.bar(
            x=bins[:-1],
            y=counts,
            title="RT barcode uniformity across beads<br><sup>Filtered to beads with > 1 RT above minUTC</sup>",
            template=reporting.DEFAULT_FIGURE_STYLE,
            labels={"x": "Normalized KL divergence score", "y": "Bead reads"},
        )
        fig.update_traces(textfont_size=12, textangle=0, textposition="outside", cliponaxis=False)
        fig.update_layout(
            bargap=0.0,
        )
        fig.add_vline(
            x=min_divergence,
            line_dash="dash",  # dashed line
            line_color="red",
            line_width=2,
            annotation_text="Filter threshold",
        )
        fig.write_image(write_dir / "figures_internal" / f"{lib_name}_BeadScoreDist.png")
        bead_score_dist = dp.Plot(fig)

    reads_in_beads = top_beads["counts"].sum() / all_beads["counts"].sum()
    reads_in_cell_beads = cell_beads["counts"].sum() / all_beads["counts"].sum()
    bead_stats = [
        ("Beads with cells", cell_beads.shape[0]),
        ("Bead Threshold", bead_count_thres),
        ("Beads above threshold", top_beads.shape[0]),
        ("Reads in beads with cells", f"{reads_in_cell_beads:.1%}"),
        ("Reads in beads above threshold", f"{reads_in_beads:.1%}"),
    ] + bead_stats
    bead_stats = pd.DataFrame.from_records(bead_stats, columns=["Metric", "Value"])
    bead_stats_table = reporting.make_table(
        bead_stats.style.pipe(reporting.styleTable, title="Bead Statistics", numericCols=["Value"])
    )

    return [bead_rank_plot, bead_pass_hist, bead_stats_table, bead_score_dist]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--libName", required=True)
    parser.add_argument("--outDir", required=True)
    parser.add_argument("--libStruct", required=True, type=Path)
    parser.add_argument("--libMetrics", required=True, type=Path)
    parser.add_argument("--demuxMetrics", type=Path, help="bcParser demux metrics json")
    parser.add_argument("--internalReport", action="store_true", default=False)
    parser.add_argument("--scalePlexLib", action="store_true", default=False)
    parser.add_argument(
        "--ultima",
        default=False,
        action="store_true",
        help="If set, sequencing data was generated on an Ultima instrument",
    )
    parser.add_argument(
        "--minPassingSampleReads",
        type=int,
        help="Minimum number of reads a sample must have post barcode demux to be included in the report",
    )
    parser.add_argument(
        # Minimum KL divergence score to consider a bead as non-ambient
        # Normalized range is roughly 0-1, with 0 perfectly matching the library RT count distribution
        "--minDivergence",
        default=0.05,
        type=float,
        help="Minimum normalized KL divergence score to consider a bead as non-ambient",
    )

    args = parser.parse_args()

    buildLibraryReport(
        args.libName,
        LibJsonParser(args.libStruct),
        args.demuxMetrics,
        args.libMetrics,
        args.internalReport,
        args.scalePlexLib,
        args.outDir,
        args.ultima,
        args.minPassingSampleReads,
        args.minDivergence,
    )


if __name__ == "__main__":
    main()
