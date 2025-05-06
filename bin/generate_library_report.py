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
from scale_utils import reporting
from scale_utils.base_logger import logger


def buildLibraryReport(
    libName: str,
    libJson: Path,
    demuxJson: Path,
    libMetrics: Path,
    internalReport: bool,
    scalePlexLib: bool,
    outDir: str,
    ultima: bool,
    minPassingSampleReads: int,
):
    """
    Build the library report by calling relevant functions

    Args:
        libName: Library name
        libJson: Library structure json
        libMetrics: Path to library metrics
        internalReport: Flag denoting whether report is for internal purposes
        scalePlexLib: Flag denoting whether library is a ScalePlex library
        outDir: Output directory
        ultima: Flag denoting whether sequencing data was generated on an Ultima instrument
        minPassingSampleReads: Minimum number of reads a sample must have post barcode
                               demux to be included in the report
    """
    libStruct = json.load(open(libJson))
    # need barcode cols for plate plots
    barcode_cols = [
        barcode.get("alias") or barcode["name"] for barcode in libStruct["barcodes"] if barcode.get("plate")
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
        demuxMetrics, libName, libJson, scalePlexLib, allCellsBetweenFiles, writeDir, ultima
    )
    pages.append(cellsPage)

    if not ultima and internalReport:
        pages.append(dp.Page(blocks=[dp.Group(barcodeReadStatsInternal, barcodeTypeStats)], title="InternalReport"))

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


def buildDfForPlatePlot(allCellsBetweenFiles: pd.DataFrame, libJson: Path, bcName: str):
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
    libStructDir = libJson.parent  # Directory containing libStruct.json and associated sequence files

    for bc in libStruct["barcodes"]:
        if bc["name"] == bcName:
            bcInfo = bc
            break
    else:
        raise ValueError(f"Unknown barcode {bcName}")

    alias = bcInfo.get("alias") or bcInfo["name"]
    well_list = allCellsBetweenFiles[alias].to_list()
    max_letter, max_number = reporting.getMaxWellNumberAndLetter(libStructDir / f'{bcInfo["sequences"]}')
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
    libJson: Path,
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
        libJson: Library Structure Json
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
        (barcodeTypeStatsDf, barcodeTypeStats) = createBarcodeTypeMetricsTables(demuxJson, libJson)
    blocks = [dp.Text(f"## libName: {libName}")]
    y_axis_label = "Unique ScalePlex Counts" if scalePlexLib else "Unique Transcript Counts"
    libStruct = json.load(open(libJson))
    for barcode in libStruct["barcodes"]:
        if not barcode.get("plate"):
            continue
        df = buildDfForPlatePlot(allCellsBetweenFiles, libJson, barcode["name"])
        df.to_csv(writeDir / "csv" / f"library_{libName}_unique_reads_{barcode['alias']}_well.csv")
        plate_plot = reporting.buildPlatePlot(df, f"{barcode['alias']} Plate", 10000.0, y_axis_label)
        blocks.append(plate_plot)
    return (barcodeTypeStatsDf, dp.Page(blocks=blocks, title="Barcodes"), barcodeTypeStats)


def getBarcodeAlias(libStruct: Dict, name: str):
    """Get the full name (alias) for a barcode from the library structure Json"""
    for bc in libStruct["barcodes"]:
        if bc["name"] == name:
            return bc.get("alias", name)
    return None


def createBarcodeTypeMetricsTables(demuxMetrics: Dict, libJson: Path):
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
    barcodesDf = buildDfFromJSONDict(demuxMetrics["barcodes"], "Barcode", "dict")
    tableGroup = []
    allBarcodes = list(barcodesDf["Barcode"].unique())
    for bc in allBarcodes:
        subset = barcodesDf[barcodesDf["Barcode"] == bc][["Match", "Reads"]]
        styledDf = subset.style.pipe(reporting.styleTable, title=f"{getBarcodeAlias(libStruct, bc)} Barcodes")
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
    df_pass.at[df_pass.index[0], "Status"] = "Barcode Pass"
    df_pass.loc[len(df_pass.index)] = [
        "Too Short Error",
        too_short_error_reads,
        str(round(too_short_error_perc, 1)) + "%",
    ]
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
                f"Sample {sample} has no passing reads post barcode demux. "
                "Sample report was not generated for this sample. \n"
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

    args = parser.parse_args()

    buildLibraryReport(
        args.libName,
        args.libStruct,
        args.demuxMetrics,
        args.libMetrics,
        args.internalReport,
        args.scalePlexLib,
        args.outDir,
        args.ultima,
        args.minPassingSampleReads,
    )


if __name__ == "__main__":
    main()
