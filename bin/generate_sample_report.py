#!/usr/bin/env python
"""
Generate the sample report (HTML and .csv) from per-sample metrics
"""

import argparse
import datetime
from pathlib import Path
from typing import Dict
import datapane as dp
import numpy as np
import pandas as pd
import polars as pl
import polars.selectors as cs
from plotly.graph_objects import Figure
import matplotlib.pyplot as plt
import plotly.express as px
from scale_utils import reporting
from scale_utils.base_logger import logger
from scale_utils.lib_json_parser import LibJsonParser
import re


class Metric:
    """
    Store all information about a single metric that is published in allSamples csv or in user facing sample report

    Args:
        csv_category: Category of the metric (e.g. "Reads", "Cells")
        display_name: Name of the metric that will be displayed in the user facing report
        col_name: Name of the metric in the csv file
        raw_value: Value of the metric pre formatting
        formatting: Formatting string to be applied to the raw_value

    Members:
        col_name: Name of the metric in the csv file
        display_name: Name of the metric that will be displayed in the user facing report or in the csv
        csv_entry: Tuple of (csv_category, display_name, formatted_value)
        warning: Boolean indicating if the metric is below threshold
        message: Warning message to be displayed if the metric is below threshold
    """

    # message is not used for now
    all_thresholds = {
        "Cells Called": {"threshold": 100, "message": "Number of cells called is below {threshold}."},
        "Reads in Cells": {"threshold": 0.5, "message": "Percentage of reads in cells is below {threshold}."},
        "Passing Sample Reads": {"threshold": 100000, "message": "Passing sample reads is below {threshold:,}."},
        "Reads Mapped to Genome": {
            "threshold": 0.5,
            "message": "Percentage of reads mapped to genome is below {threshold}.",
        },
        "Reads Mapped to Transcriptome": {
            "threshold": 0.5,
            "message": "Percentage of reads mapped to transcriptome is below {threshold}.",
        },
        "Antisense Reads": {"threshold": 0.05, "message": "Percentage of antisense reads is above {threshold}."},
        "Reads per Cell": {"threshold": 1000, "message": "Reads per cell is below {threshold}."},
        "Median Genes per Cell": {"threshold": 100, "message": "Median genes per cell is below {threshold}."},
    }

    def __init__(
        self,
        csv_category: str,
        display_name: str,
        col_name: str,
        raw_value: str | int | float,
        formatting: str = "",
    ):
        self.col_name = col_name  # class attribute not used for now, keeping in case we need it later
        self.raw_value = raw_value
        self.display_name = display_name
        if formatting:
            formatted_value = f"{raw_value:{formatting}}"
        else:
            formatted_value = raw_value
        self.report_entry = (csv_category, display_name, reporting.formatNumericVal(formatted_value))
        self.csv_entry = (csv_category, display_name, formatted_value)
        threshold = self.all_thresholds.get(display_name)
        self.warning = False
        self.message = None
        if threshold and formatting:
            if (display_name == "Antisense Reads" and raw_value > threshold["threshold"]) or (
                display_name != "Antisense Reads" and raw_value < threshold["threshold"]
            ):
                self.warning = True
                self.message = "<li>" + threshold["message"].format(threshold=threshold["threshold"]) + "</li>"
                self.report_entry = (
                    csv_category,
                    # Add hover text to display the warning message
                    f'<span title="Metric is outside the normal range of values for a typical experiment"> '
                    f"{display_name} </span>",
                    f'<span title="Metric is outside the normal range of values for a typical experiment"> '
                    f"{self.report_entry[2]} </span>",
                )


def build_sample_report(
    write_dir: Path,
    sample_metrics: Path,
    sample_id: str,
    internal_report: bool,
    is_barnyard: bool,
    is_scale_plex: bool,
    cellfinder: bool,
    lib_struct_json: Path,
    metadata: Dict[str, str],
    hash_metrics: Path = None,
    hash_stats: Path = None,
    hash_lib_name: str = None,
    hash_all_cells: Path = None,
):
    """
    Entry function for building sample report

    Args:
        write_dir: Report Output directory
        sample_metrics: Folder containing csv files defining the metrics for a single sample
        sample_id: Unique Sample ID (e.g. `SampleNam.LibName`)
        internal_report: Create additional plots for internal QC
        is_barnyard: Add barnyard plots and metrics
        is_scale_plex: Add page with hash analysis tables and plots
        cellfinder: Flag indicating whether cellfinder was used for cell calling
        lib_struct_json: Library structure json (in directory with sequence files, etc.)
        trim_stats: json files containing cutadapt output
        metadata: Information about the sample and analysis run [FieldName -> Value]
        hash_metrics: Path to hash metrics file
        hash_stats: Path to hash stats file
        hash_lib_name: Name of the hash library
        hash_all_cells: Path to hash all cells file
    """
    lib_json_obj = LibJsonParser(lib_struct_json)

    columns_to_select = [
        pl.col("pass"),
        # for barnyard also need human_counts and mouse_counts columns
        cs.ends_with("counts"),
        cs.matches("species"),
        pl.col("totalReads"),
        cs.ends_with("genes"),
        pl.col("Saturation"),
    ]
    if internal_report and any(
        ["bead" in (bc.get("alias") or bc["name"]) for bc in lib_json_obj.json_contents["barcodes"]]
    ):
        columns_to_select.append(pl.col("bead_bc"))
    all_cells = pl.scan_parquet(sample_metrics / "allBarcodes.parquet").select(*columns_to_select).collect().to_pandas()
    passing_cells = pl.read_csv(sample_metrics / "allCells.csv").to_pandas().set_index("cell_id")
    sample_stats = pd.read_csv(sample_metrics / "sample_stats.csv", index_col=["Category", "Metric"])

    Path(write_dir, "csv").mkdir(parents=True, exist_ok=True)
    if internal_report:
        Path(write_dir, "figures_internal").mkdir(exist_ok=True)

    pages = []
    reads_page, report_stats, internal_report_blocks = build_reads_page(
        sample_id,
        all_cells,
        passing_cells,
        sample_stats,
        write_dir,
        internal_report,
        metadata,
        is_barnyard,
        cellfinder,
    )
    pages.append(reads_page)

    all_index_plots, barcodes_page = build_barcodes_page(
        lib_json_obj, passing_cells, write_dir, sample_id, internal_report
    )
    pages.append(barcodes_page)

    if internal_report:
        plate_plots = []
        for bc in lib_json_obj.json_contents["barcodes"]:
            if bc.get("plate"):
                plate_plots.append(all_index_plots[bc["name"]])
        dist_group = dp.Group(dp.Text("## Read Distributions"), dp.Group(blocks=internal_report_blocks, columns=2))
        plate_group = dp.Group(blocks=plate_plots, columns=1)
        internal_report_page = dp.Page(
            dp.Group(
                dist_group,
                plate_group,
                columns=1,
            ),
            title="InternalReport",
        )
        pages.append(internal_report_page)

    if is_barnyard:
        barnyard_stats_obj = make_barnyard_stats(passing_cells)
        barnyard_stats = pd.DataFrame.from_records(
            [metric.report_entry for metric in barnyard_stats_obj], columns=["Category", "Metric", "Value"]
        )
        barnyard_page = build_barnyard_page(
            all_cells, barnyard_stats, write_dir, sample_id, internal_report, barnyard_stats_obj
        )
        report_stats = pd.concat([report_stats, barnyard_stats])
        pages.append(barnyard_page)

    if is_scale_plex:
        # Indicates that we have passing scaleplex reads for this library
        if hash_stats:
            hash_metrics_series = pd.read_csv(hash_metrics, header=None, index_col=0).squeeze(axis=1)
            assigned_hashes = pd.read_csv(hash_stats, index_col=0)
            hash_all_cells_df = (
                pl.scan_parquet(hash_all_cells)
                .select(
                    pl.col("pass"),
                    pl.col("counts"),
                    pl.col("totalReads"),
                    pl.col("Saturation"),
                )
                .collect()
                .to_pandas()
            )
            hash_passing_cells = (
                pl.scan_parquet(hash_all_cells)
                .filter(pl.col("pass"))
                .select(
                    pl.col("topTwo"),
                    pl.col("assigned_scaleplex"),
                )
                .collect()
                .to_pandas()
            )
            scale_plex_page, scale_plex_stats = build_scale_plex_page(
                hash_lib_name, hash_metrics_series, assigned_hashes, hash_all_cells_df, hash_passing_cells
            )
            report_stats = pd.concat([report_stats, scale_plex_stats])
        else:
            scale_plex_page = dp.Page(dp.Text("Library has no passing ScalePlex reads"), title="ScalePlex")
        pages.append(scale_plex_page)

    report = dp.Report(blocks=pages)

    logger.debug(f"Writing reportStatistics csv, report to {str(write_dir.resolve())}")
    # Change 'nan' back to np.nan in Value column
    report_stats.replace(to_replace=["nan", "nan%"], value=np.nan, inplace=True)
    report_stats.to_csv(
        write_dir / "csv" / f"{sample_id}.reportStatistics.csv",
        index=False,
        header=False,
        na_rep="",
        columns=["Category", "Metric", "Value"],
    )  # can change na_rep to '<NA>' or similar if needed

    report.save(write_dir / f"{sample_id}.report.html")


def get_read_metrics(sample_stats: pd.DataFrame) -> list[Metric]:
    """
    Make list of Metric objects depicting read stats that will be displayed in the sample report

    Args:
        sample_stats: Sample-level statistics
    """
    # Need to divide by 100 because the raw value is a percentage
    mapped_reads_perc = (sample_stats.loc[("Reads", "mapped_to_too_many_loci_perc")].Value / 100) + sample_stats.loc[
        ("Reads", "mapped_reads_perc")
    ].Value
    stats_obj = []
    if ("Reads", "total_sample_reads_post_barcode_demux") in sample_stats.index:
        stats_obj.append(
            Metric(
                "Reads",
                "Total Sample Reads",
                "total_sample_reads_post_barcode_demux",
                sample_stats.loc[("Reads", "total_sample_reads_post_barcode_demux")].Value,
                ".0f",
            )
        )
    stats_obj.extend(
        [
            Metric(
                "Reads", "Passing Sample Reads", "total_reads", sample_stats.loc[("Reads", "total_reads")].Value, ".0f"
            ),
            Metric(
                "Reads",
                "Reads Mapped to Genome",
                "mapped_reads_perc",
                mapped_reads_perc,
                ".1%",
            ),
            Metric(
                "Reads",
                "Passing Read Alignments",
                "passing_aligns_perc",
                1 - (sample_stats.loc[("Reads", "mapped_to_too_many_loci_perc")].Value / 100),
                ".1%",
            ),
            Metric(
                "Reads",
                "Reads Mapped to Transcriptome",
                "gene_reads_perc",
                sample_stats.loc[("Reads", "gene_reads_perc")].Value,
                ".1%",
            ),
            Metric(
                "Reads", "Exonic Reads", "exon_reads_perc", sample_stats.loc[("Reads", "exon_reads_perc")].Value, ".1%"
            ),
            Metric(
                "Reads",
                "Antisense Reads",
                "antisense_reads_perc",
                sample_stats.loc[("Reads", "antisense_reads_perc")].Value,
                ".1%",
            ),
            Metric(
                "Reads",
                "Mitochondrial Reads",
                "mito_reads_perc",
                sample_stats.loc[("Reads", "mito_reads_perc")].Value,
                ".1%",
            ),
            Metric("Reads", "Saturation", "saturation", sample_stats.loc[("Reads", "saturation")].Value, ".2f"),
            Metric(
                "Reads",
                "Average Trimmed Read Length",
                "avg_trimmed_read_len",
                sample_stats.loc[("Reads", "avg_trimmed_read_len")].Value,
                ".0f",
            ),
            Metric(
                "Reads",
                "Average Mapped Length",
                "avg_mapped_len",
                sample_stats.loc[("Reads", "avg_mapped_len")].Value,
                ".1f",
            ),
            Metric(
                "Reads",
                "Mapped Mismatch Rate(%)",
                "mismatch_rate_per_base_perc",
                sample_stats.loc[("Reads", "mismatch_rate_per_base_perc")].Value,
                ".2f",
            ),
        ]
    )
    return stats_obj


def build_reads_page(
    sample: str,
    all_cells: pd.DataFrame,
    passing_cells: pd.DataFrame,
    sample_stats: pd.DataFrame,
    write_dir: Path,
    internal_report: bool,
    metadata: Dict[str, str],
    is_barnyard: bool,
    cellfinder: bool,
):
    """
    Build a datapane page that mainly shows plots related to reads

    Args:
        sample: Sample name
        all_cells: Per-cell metrics
        passing_cells: Per-cell metrics for passing cells
        sample_stats: Sample-level statistics
        write_dir: Output directory
        internal_report: Flag indicating whether report is to be
            generated for internal r&d run
        metadata: Information about the sample and analysis run [FieldName -> Value]
        is_barnyard: Flag indicating whether barnyard analysis was performed
        cellfinder: Flag indicating whether cellfinder was used for cell calling

    Returns:
        dp.Page object and concatenated dataframe of various statistics
        List of elements for the internal report (or empty)
    """
    report_stats = pd.DataFrame({"Category": ["Sample"], "Metric": ["SampleName"], "Value": [sample]})
    metrics_obj = Metric("Sample", "SampleName", "SampleName", sample)
    read_metrics_obj = get_read_metrics(sample_stats)
    cell_metrics_obj = get_cell_metrics(sample_stats, is_barnyard)
    complexity_metrics_obj = get_complexity_metrics(sample_stats)
    all_metrics = [metrics_obj]
    all_metrics.extend(read_metrics_obj)
    all_metrics.extend(cell_metrics_obj)
    all_metrics.append(
        Metric(
            "Cells", "Mean Reads", "mean_passing_reads", sample_stats.loc[("Cells", "mean_passing_reads")].Value, ".0f"
        )
    )
    all_metrics.extend(complexity_metrics_obj)
    all_csv_entries = [metric.csv_entry for metric in all_metrics]
    report_stats = pd.DataFrame.from_records(all_csv_entries, columns=["Category", "Metric", "Value"])
    info = pd.DataFrame.from_dict({"Metric": metadata.keys(), "Value": metadata.values()})
    info_table = reporting.create_metric_table(info, title="Sample Information")
    cell_table = reporting.create_metric_table(
        pd.DataFrame.from_records(
            [metric.report_entry for metric in cell_metrics_obj], columns=["Category", "Metric", "Value"]
        ),
        title="Cell Metrics",
        obj_list=cell_metrics_obj,
    )
    mapping_table = reporting.create_metric_table(
        pd.DataFrame.from_records(
            [metric.report_entry for metric in read_metrics_obj], columns=["Category", "Metric", "Value"]
        ),
        title="Read Metrics",
        obj_list=read_metrics_obj,
        rm_nan=True,
    )

    rank_plot = make_rank_plot(all_cells, write_dir, sample, internal_report, cellfinder, is_barnyard)
    genes_counts_scatter = build_gene_read_scatter(all_cells)
    saturation_scatter = build_saturation_scatter(all_cells)
    group_blocks = [
        info_table,
        cell_table,
        mapping_table,
        rank_plot,
        genes_counts_scatter,
        saturation_scatter,
    ]
    unique_reads_fig = make_unique_reads_fig(sample_stats)
    group_blocks.insert(4, unique_reads_fig)
    page = dp.Page(dp.Group(blocks=group_blocks, columns=2), title="Summary")

    internal_report_blocks = []
    if internal_report:
        gene_reads_dist = build_dist(
            "Distribution of Reads: Gene Matched", "geneReads", "Fraction of Genic Reads", passing_cells
        )
        antisense_reads_dist = build_dist(
            "Distribution of Reads: Antisense", "antisenseReads", "Fraction of Antisense Reads", passing_cells
        )
        exonic_reads_dist = build_dist(
            "Distribution of Reads: Exonic", "exonReads", "Fraction of Exonic Reads", passing_cells
        )
        internal_report_blocks.extend([gene_reads_dist, exonic_reads_dist, antisense_reads_dist])

        if sample_stats.loc[("Cells", "mito_reads")].Value != 0:
            mito_reads_dist = build_dist(
                "Distribution of Reads: Mitochondrial", "mitoReads", "Fraction of Mitochondrial Reads", passing_cells
            )
            internal_report_blocks.append(mito_reads_dist)

    return page, report_stats, internal_report_blocks


def build_barcodes_page(
    lib_json_obj: LibJsonParser,
    passing_cells: pd.DataFrame,
    write_dir: Path,
    sample_id: str,
    internal_report: bool,
) -> tuple[dict, dp.Page]:
    """
    Build datapane page for barcodes

    Args:
        lib_json_obj: Object containing library structure information
        passing_cells: Metrics per passing cell-barcode
        write_dir: Ouput directory
        sample_id: Unique ID of this sample
        internal_report: Extra plots for internal QC

    Returns:
        Index plots for internal reports (empty otherwise)
        Page object with shared figures
    """
    blocks_to_render = []
    all_index_plots = {}
    for bc in lib_json_obj.json_contents["barcodes"]:
        if bc.get("plate"):
            alias = bc.get("alias") or bc["name"]
            index_plots = reporting.barcodeLevelPlots(
                lib_json_obj,
                sample_id,
                passing_cells,
                alias,
                f"{alias} Index",
                internal_report,
                write_dir,
            )
            all_index_plots[bc["name"]] = index_plots
    blocks_to_render.append(all_index_plots[lib_json_obj.json_contents["sample_barcode"]])

    return all_index_plots, dp.Page(blocks=blocks_to_render, title="Barcodes")


def build_plate_plot(counts: pd.Series) -> dp.Plot:
    """
    Build a fixation plate plot for assigned ScalePlex

    Returns:
        Plot object with the figure
    """
    cols = [str(i) for i in range(1, 13)]  # columns on 96 well plate are 1-12
    rows = ["A", "B", "C", "D", "E", "F", "G", "H"]  # rows on 96 well plate are A-H
    wells = [["" for j in cols] for i in rows]  # initialize empty 96 well plate
    for col_idx, j in enumerate(cols):
        for row_idx, i in enumerate(rows):
            current_well = f"{j}{i}"
            wells[row_idx][col_idx] = counts[current_well] if current_well in counts else 0
    fig = reporting.buildPlatePlot(
        pd.DataFrame(wells, index=rows, columns=cols),
        "Fixation Plate Assigned ScalePlex",
        1.0,
        subtitle="Cells called per well",
    )
    return dp.Plot(fig)


def build_scale_plex_page(
    lib_name: str,
    hash_metrics_series: pd.Series,
    assigned_hashes: pd.Series,
    hash_all_cells_df: pd.DataFrame,
    hash_passing_cells: pd.DataFrame,
) -> tuple[dp.Page, pd.DataFrame]:
    """
    Build a datapane page displaying hash metrics

    Args:
        lib_name: Hash lib name
        hash_metrics_series: Sample-level hash metrics
        assigned_hashes: Frequency counts of assigned hashes
        hash_all_cells_df: Cell-level metrics
        hash_passing_cells: Cell-level metrics for passing cells

    Returns:
        datapane page
        DataFrame with hash metrics
    """
    stats_obj = [
        Metric("ScalePlex", "Library Name", None, lib_name),
        Metric("ScalePlex", "Mean Reads per Cell", None, hash_metrics_series["ReadsPerCell"]),
        Metric("ScalePlex", "Median Unique ScalePlex Counts per Cell", None, hash_metrics_series["nUMIPerCell"]),
        Metric("ScalePlex", "Percent of Cells with Assigned ScalePlex", None, hash_metrics_series["passingPercent"]),
        Metric("ScalePlex", "Saturation", None, hash_metrics_series["saturation"]),
        Metric("ScalePlex", "Reads in Cells", None, hash_metrics_series["readsInCells"]),
    ]
    stats_df = pd.DataFrame.from_records(
        [metric.report_entry for metric in stats_obj], columns=["Category", "Metric", "Value"]
    )
    hash_table = reporting.create_metric_table(stats_df, title="ScalePlex Metrics", obj_list=stats_obj)

    fig = px.bar(
        assigned_hashes["cells_called"],
        title="Assigned ScalePlex Cell Counts",
        labels={"value": "Cell Counts", "index": ""},
        template=reporting.DEFAULT_FIGURE_STYLE,
    )
    fig.update_traces(hovertemplate="Well: %{x}<br>Cells: %{y}<extra></extra>")
    fig.update_layout(
        showlegend=False, autosize=False, margin=dict(l=50, r=50, b=150, t=50, pad=10)
    )  # left margin  # right margin  # bottom margin  # top margin  # padding
    fig.update_xaxes(tickangle=-90, title_text="Assigned ScalePlex")

    saturation_scatter = build_saturation_scatter(hash_all_cells_df)

    toptwo = px.histogram(
        x=hash_passing_cells["topTwo"],
        title="Top ScalePlex Fraction",
        histnorm="percent",
        template=reporting.DEFAULT_FIGURE_STYLE,
    )
    toptwo.update_layout(yaxis_title="% Cells", xaxis_title="")
    toptwo.update_xaxes(range=[0, 1])

    fixation_qc_plot = build_plate_plot(hash_passing_cells["assigned_scaleplex"].value_counts())

    error_stats_obj = [
        Metric("ScalePlex", "Max Fail", None, hash_metrics_series["maxFailPercent"]),
        Metric("ScalePlex", "Enrich Fail", None, hash_metrics_series["enrichFailPercent"]),
        Metric("ScalePlex", "Unexpected", None, hash_metrics_series["unexpectedPercent"]),
        Metric("ScalePlex", "Indeterminate", None, hash_metrics_series["indeterminatePercent"]),
    ]
    error_stats_df = pd.DataFrame.from_records(
        [metric.report_entry for metric in error_stats_obj], columns=["Category", "Metric", "Value"]
    )
    assignment_error_table = reporting.create_metric_table(
        error_stats_df, title="ScalePlex Assignment Errors", obj_list=error_stats_obj
    )

    # create per-assigned scaleplex metrics table
    # first sort by well position and then reset_index to bring well into a column
    assigned_hashes = assigned_hashes.sort_index(
        key=lambda x: x.str.extract(r"(\d+)(\D+)", expand=True).apply(lambda x: (int(x[0]), x[1]), axis=1)
    )
    assigned_hashes_styled = assigned_hashes.reset_index().style.pipe(
        reporting.styleTable,
        title="Per-well ScalePlex RNA Metrics",
        numericCols=["passing_reads", "cells_called", "mean_passing_reads", "median_utc"],
        boldColumn="assigned_scaleplex",
    )
    assigned_hashes_styled = assigned_hashes_styled.relabel_index(
        [
            "Well",
            "Passing Reads",
            "Cells Called",
            "Mean Passing Reads",
            "Median Unique Transcript Counts",
            "Median Genes",
        ],
        axis=1,
    )
    scaleplex_stats_table = reporting.make_table(assigned_hashes_styled)

    group_blocks = [
        hash_table,
        dp.Plot(fig),
        saturation_scatter,
        dp.Plot(toptwo),
        fixation_qc_plot,
        assignment_error_table,
        scaleplex_stats_table,
    ]
    page = dp.Page(dp.Group(blocks=group_blocks, columns=2), title="ScalePlex")

    return page, pd.concat([stats_df, error_stats_df])


def build_dist(title: str, field: str, label: str, passing_cells: pd.DataFrame) -> dp.Plot:
    """
    Build histogram for a given field divided by mapped reads

    Args:
        title: Title of plot
        field: Column name in all_cells to plot on x axis after dividing by mapped reads
        label: String that will be the x axis label
        passing_cells: Cell-level metrics

    Returns:
        Plot object with the figure
    """
    x = passing_cells[field] / passing_cells["mappedReads"]

    # calculate histogram with matplotlib
    counts, bins, _ = plt.hist(x, bins="auto")

    # calculate percentage of cells in each bin
    total = sum(counts)
    percentages = (counts / total) * 100
    yaxis_title = "% Cells"

    # plot as bar chart in plotly so raw dataset is not saved to HTML report
    fig = px.bar(
        x=bins[:-1],
        y=percentages,
        title=title,
        template=reporting.DEFAULT_FIGURE_STYLE,
        labels={"x": label, "y": yaxis_title},
    )

    fig.update_layout(
        yaxis_title=yaxis_title,
        bargap=0.0,
    )
    fig.update_traces(hovertemplate=f"{label}: %{{x:.2f}}<br>{yaxis_title}: %{{y:.2f}}%<extra></extra>")

    return dp.Plot(fig)


def subsample_all_cells(all_cells: pd.DataFrame) -> pd.DataFrame:
    """
    Subsample allCells to a maximum of 5000 cells for plotting purposes

    Too many cells results in overplotting so in this function we:
        (1) Take all passing cells
        (2) Take an equal number of the top failing cells by counts
        (3) Return a maximum of 5000 cells by subsampling if necessary

    Args:
        all_cells: Cell-level metrics

    Returns:
        Subsampled dataframe
    """
    passing_cells = all_cells.loc[all_cells["pass"]]
    top_failed_cells = (
        all_cells[~all_cells["pass"]].sort_values(by="counts", ascending=False).head(passing_cells.shape[0])
    )
    top_cells = pd.concat([passing_cells, top_failed_cells])
    if top_cells.shape[0] > 5000:
        top_cells = top_cells.sample(5000)
    return top_cells


def build_gene_read_scatter(all_cells: pd.DataFrame) -> dp.Plot:
    """
    Make a scatter plot of reads vs counts

    Args:
        all_cells: Cell-level metrics

    Returns:
        Plot object with the figure
    """
    all_cells["plt_color"] = all_cells["pass"].map({True: "Cell", False: "Background"})
    if "mouse_genes" in all_cells.columns:
        # for barnyard data don't use sum of genes across species but max species
        all_cells["totalGenes"] = all_cells[["mouse_genes", "human_genes"]].max(axis=1)
    else:
        all_cells["totalGenes"] = all_cells["genes"]
    fig = px.scatter(
        subsample_all_cells(all_cells),
        x="totalReads",
        y="totalGenes",
        color="plt_color",
        color_discrete_map=reporting.SCATTER_COLORMAP,
        template=reporting.DEFAULT_FIGURE_STYLE,
        labels={"totalReads": "Total reads", "genes": "Genes detected", "plt_color": ""},
        title="Genes Detected Per Cell",
        category_orders={"plt_color": ["Cell", "Background"]},
        opacity=0.5,
    )

    return dp.Plot(fig)


def build_saturation_scatter(all_cells: pd.DataFrame) -> dp.Plot:
    """
    Build scatter plot for saturation

    Args:
        all_cells: Cell-level metrics

    Returns:
        Plot object with the figure
    """
    all_cells["plt_color"] = all_cells["pass"].map({True: "Cell", False: "Background"})
    sample_all_cells = subsample_all_cells(all_cells)
    fig_bg = px.scatter(
        sample_all_cells[~sample_all_cells["pass"]],
        x="totalReads",
        y="Saturation",
        color="plt_color",
        labels={"totalReads": "Total reads", "plt_color": ""},
        template=reporting.DEFAULT_FIGURE_STYLE,
        color_discrete_map=reporting.SCATTER_COLORMAP,
        title="Saturation Per Cell",
        log_x=False,
        category_orders={"plt_color": ["Cell", "Background"]},
        opacity=0.5,
    )
    fig_pass = px.scatter(
        sample_all_cells[sample_all_cells["pass"]],
        x="totalReads",
        y="Saturation",
        color="plt_color",
        labels={"totalReads": "Total reads", "plt_color": ""},
        template=reporting.DEFAULT_FIGURE_STYLE,
        color_discrete_map=reporting.SCATTER_COLORMAP,
        title="Saturation Per Cell",
        log_x=False,
        category_orders={"plt_color": ["Cell", "Background"]},
        opacity=0.5,
    )
    fig = Figure(data=fig_bg.data + fig_pass.data, layout=fig_bg.layout)
    return dp.Plot(fig)


def make_rank_plot(
    all_cells: pd.DataFrame,
    write_dir: Path,
    sample_id: str,
    internal_report: bool,
    cellfinder: bool,
    barynard: bool,
) -> dp.Plot:
    """
    Make barcode rank plot and draw vertical line at threshold if cellfinder was not used

    Args:
        all_cells: Cell-level metrics
        write_dir: Output directory
        sample_id: Unique ID of this sample
        internal_report: Extra plots for internal QC
        cellfinder: Flag indicating whether cellfinder was used for cell calling
        barnyard: Flag indicating whether barnyard analysis was performed

    Returns:
        Plot object with the figure
    """
    indices = reporting.sparseLogCoords(all_cells.index.size)
    all_cells = all_cells.sort_values(by=["counts"], ascending=False).reset_index(drop=True)
    if cellfinder or barynard:
        # Calculate proportion of cells with a moving window of 25.
        # min_periods=1 allows us to calculate proportion of windows with less than 10 barcodes.
        # This prevents nans.
        all_cells["z"] = all_cells["pass"].rolling(25, min_periods=1).mean()
        # find the index of "last" passing cell
        last_passing_index = all_cells["pass"][::-1].idxmax()
        # rolling average shouldn't extend into the background
        all_cells["z"] = all_cells["z"].where(all_cells.index <= last_passing_index, other=0)
        # Color plot using new classification and update category orders to fix overplotting.
        fig = px.scatter(
            all_cells.iloc[indices],
            x=indices,
            y="counts",
            color="z",
            labels={
                "x": "Cell barcodes",
                "counts": "Unique transcript counts",
                "cell": "Cell Barcode",
                "z": "Prop. Cells",
            },
            log_x=True,
            log_y=True,
            template=reporting.DEFAULT_FIGURE_STYLE,
            title="Barcode Rank Plot",
            opacity=0.5,
            color_continuous_scale=["rgb(178, 181, 187)", "rgb(39, 139, 176)"],
            hover_data={"z": True},
        )
        fig.update_layout(
            coloraxis_colorbar=dict(
                thickness=10, tickvals=[0, 1], ticktext=["Background", "Cell"], len=0.25, y=0.75, title_text=""
            )
        )
    else:
        all_cells["plt_color"] = all_cells["pass"].map({True: "Cell", False: "Background"})
        fig = px.scatter(
            all_cells.iloc[indices],
            x=indices,
            y="counts",
            color="plt_color",
            labels={"x": "Cell barcodes", "counts": "Unique transcript counts", "plt_color": ""},
            log_x=True,
            log_y=True,
            template=reporting.DEFAULT_FIGURE_STYLE,
            title="Barcode Rank Plot",
            color_discrete_map=reporting.SCATTER_COLORMAP,
            category_orders={"plt_color": ["Cell", "Background"]},
            opacity=0.5,
            hover_data={"plt_color": False},
        )
    # If cellfinder was not used for cell calling. Add UTC threshold to plot.
    if not cellfinder and not barynard:
        fig.add_vline(x=max(1, all_cells["pass"].sum()), line_dash="dash", line_color="green")
    if indices.size > 0:
        fig.update_layout(xaxis_range=[1, np.log10(max(indices) + 1)])
    if internal_report:
        save_figure_png(write_dir, sample_id, fig, "BarcodeRankPlot.png")
    return dp.Plot(fig)


def save_figure_png(write_dir: Path, sample_id: str, fig: Figure, fname: str):
    """
    Save a plotly figure object as png

    Args:
        write_dir: Output directory
        sample_id: Unique ID of this sample
        fig: Plotly figure object
        fname: Name of the file to save
    """
    fig.write_image(write_dir / "figures_internal" / f"{sample_id}_{fname}")


def build_barnyard_page(
    all_cells: pd.DataFrame,
    barnyard_stats: pd.DataFrame,
    write_dir: Path,
    sample_id: str,
    internal_report: bool,
    barnyard_stats_obj: list[Metric],
) -> dp.Page:
    """
    Build a datapane page for barnyard analysis

    Args:
        all_cells: Cell-level metrics
        barnyard_stats: Barnyard statistics
        write_dir: Output directory
        sample_id: Unique ID of this sample
        internal_report: Figures to be published as png for internal QC
        barnyard_stats_obj: List of instances of the Metric class for barnyard stats

    Returns:
        Page object with the relevant figures
    """
    barnyard_table = reporting.create_metric_table(barnyard_stats, title="", obj_list=barnyard_stats_obj)
    barnyard_plot = make_barnyard_scatter_plot(all_cells, write_dir, sample_id, internal_report)
    plot_group = dp.Group(barnyard_table, barnyard_plot, columns=2)
    return dp.Page(plot_group, title="Barnyard")


def make_barnyard_scatter_plot(
    all_cells: pd.DataFrame, write_dir: Path, sample_id: str, internal_report: bool
) -> dp.Plot:
    """
    Scatter plot transcript count per species for all cells

    Args:
        all_cells: Cell-level metrics
        write_dir: Output directory
        sample_id: Unique ID of this sample
        internal_report: Figures to be published as png for internal QC

    Returns:
        Plot object with the figure
    """
    fig = px.scatter(
        subsample_all_cells(all_cells),
        x="human_counts",
        y="mouse_counts",
        color="species",
        log_x=True,
        log_y=True,
        template=reporting.DEFAULT_FIGURE_STYLE,
        color_discrete_map=reporting.BARNYARD_COLORMAP,
        labels={"human_counts": "Human Unique Transcript Counts", "mouse_counts": "Mouse Unique Transcript Counts"},
    )

    fig.update_traces(selector=dict(name="None"), visible="legendonly")
    if internal_report:
        save_figure_png(write_dir, sample_id, fig, "Barnyard_ScatterPlot.png")

    return dp.Plot(fig)


def make_barnyard_stats(passing_cells: pd.DataFrame) -> list[Metric]:
    """
    Compute statistics for barnyard samples

    Args:
        passing_cells: Cell-level metrics

    Returns:
        Metric objects for all barnyard statistics
    """
    passing_barnyard_cells = passing_cells.loc[~passing_cells.species.isin(["Ambiguous", "None"])]
    ncells = max(1, passing_barnyard_cells.index.size)
    prop_human = (passing_barnyard_cells.species == "Human").sum() / ncells
    prop_mouse = (passing_barnyard_cells.species == "Mouse").sum() / ncells
    prop_mixed = (passing_barnyard_cells.species == "Mixed").sum() / ncells
    doublets = prop_mixed / (2 * prop_human * prop_mouse) if prop_human > 0 and prop_mouse > 0 else 0
    doublets = min(1, doublets)
    bg = passing_barnyard_cells[passing_barnyard_cells.species != "Mixed"].minor_frac.median()
    mouse_bg = passing_barnyard_cells[
        (passing_barnyard_cells.species != "Mixed") & (passing_barnyard_cells.species != "Human")
    ].minor_frac.median()
    human_bg = passing_barnyard_cells[
        (passing_barnyard_cells.species != "Mixed") & (passing_barnyard_cells.species != "Mouse")
    ].minor_frac.median()
    human_counts_med = passing_barnyard_cells[(passing_barnyard_cells.species == "Human")].human_counts.median()
    mouse_counts_med = passing_barnyard_cells[(passing_barnyard_cells.species == "Mouse")].mouse_counts.median()
    human_genes_med = passing_barnyard_cells[(passing_barnyard_cells.species == "Human")].human_genes.median()
    mouse_genes_med = passing_barnyard_cells[(passing_barnyard_cells.species == "Mouse")].mouse_genes.median()
    human_reads_med = passing_barnyard_cells[(passing_barnyard_cells.species == "Human")].totalReads.median()
    mouse_reads_med = passing_barnyard_cells[(passing_barnyard_cells.species == "Mouse")].totalReads.median()
    human_saturation = passing_barnyard_cells[(passing_barnyard_cells.species == "Human")].Saturation.median()
    mouse_saturation = passing_barnyard_cells[(passing_barnyard_cells.species == "Mouse")].Saturation.median()

    stats_obj = [
        Metric("Barnyard", "Human Cells", None, prop_human, ".1%"),
        Metric("Barnyard", "Mouse Cells", None, prop_mouse, ".1%"),
        Metric("Barnyard", "Mixed Cells", None, prop_mixed, ".1%"),
        Metric("Barnyard", "Estimated Doublets", None, doublets, ".1%"),
        Metric("Barnyard", "Background", None, bg, ".2%"),
        Metric("Barnyard", "Mouse Background", None, mouse_bg, ".2%"),
        Metric("Barnyard", "Human Background", None, human_bg, ".2%"),
        Metric("Barnyard", "Median UMIs per Cell (Human)", None, human_counts_med, ".0f"),
        Metric("Barnyard", "Median UMIs per Cell (Mouse)", None, mouse_counts_med, ".0f"),
        Metric("Barnyard", "Median Genes per Cell (Human)", None, human_genes_med, ".0f"),
        Metric("Barnyard", "Median Genes per Cell (Mouse)", None, mouse_genes_med, ".0f"),
        Metric("Barnyard", "Median Reads per Cell (Human)", None, human_reads_med, ".0f"),
        Metric("Barnyard", "Median Reads per Cell (Mouse)", None, mouse_reads_med, ".0f"),
        Metric("Barnyard", "Median Saturation per Cell (Human)", None, human_saturation, ".2f"),
        Metric("Barnyard", "Median Saturation per Cell (Mouse)", None, mouse_saturation, ".3f"),
    ]
    return stats_obj


def make_unique_reads_fig(sample_stats: pd.DataFrame) -> dp.Plot:
    """
    Create plot for estimate of unique transcript counts at different sequencing depths

    The plot includes a point for the observed values in this sample (not extrapolated)
    but this is excluded from the dataframe (For statistics.csv output)

    Args:
        sample_stats: Sample-level statistics

    Returns:
        Plot object with the figure
    """
    # Add estimated unique reads at different sequencing depths
    target_mean_reads = [
        (int(re.match(r"^target_reads_(\d+)", name).group(1)), row["Value"])
        for name, row in sample_stats.loc["Extrapolated Complexity"].iterrows()
    ]

    # Add the observed value to the plot
    target_mean_reads.append(
        (sample_stats.loc[("Cells", "reads_per_cell")].Value, sample_stats.loc[("Cells", "median_utc")].Value)
    )

    plot_df = pd.DataFrame.from_records(target_mean_reads, columns=["x", "unique_read"])
    plot_df.sort_values(by="x", inplace=True)
    if ("Reads", "total_sample_reads_post_barcode_demux") in sample_stats.index:
        x_axis_label = "Total reads per cell"
    else:
        x_axis_label = "Passing reads per cell"
    fig = px.line(
        plot_df,
        x="x",
        y="unique_read",
        labels={"x": x_axis_label, "unique_read": "Median unique transcript counts (extrapolated)"},
        template=reporting.DEFAULT_FIGURE_STYLE,
        markers=True,
        title="Complexity",
    )

    return dp.Plot(fig)


def get_cell_metrics(sample_stats: pd.DataFrame, is_barnyard: bool) -> list[Metric]:
    """
    Make list of Metric objects that represent summary stats

    Args:
        sample_stats: Sample-level statistics
    """
    stats_obj = [
        Metric("Cells", "Cells Called", "cells_called", sample_stats.loc[("Cells", "cells_called")].Value, ".0f"),
        Metric(
            "Cells",
            "Median Unique Transcript Counts per Cell",
            "median_utc",
            sample_stats.loc[("Cells", "median_utc")].Value,
            ".0f",
        ),
        Metric(
            "Cells", "Median Genes per Cell", "median_genes", sample_stats.loc[("Cells", "median_genes")].Value, ".0f"
        ),
        Metric("Cells", "Reads in Cells", "reads_in_cells", sample_stats.loc[("Cells", "reads_in_cells")].Value, ".1%"),
    ]
    if ("Cells", "reads_per_cell") in sample_stats.index:
        stats_obj.insert(
            1,
            Metric(
                "Cells", "Reads per Cell", "reads_per_cell", sample_stats.loc[("Cells", "reads_per_cell")].Value, ".0f"
            ),
        )
    if is_barnyard and ("Cells", "utc_threshold") in sample_stats.index:
        stats_obj.insert(
            1,
            Metric(
                "Cells",
                "Unique Transcript Counts Threshold (Mouse)",
                "mouse_utc",
                sample_stats.loc[("Barnyard", "mouse_utc")].Value,
                ".0f",
            ),
        )
        stats_obj.insert(
            1,
            Metric(
                "Cells",
                "Unique Transcript Counts Threshold (Human)",
                "human_utc",
                sample_stats.loc[("Barnyard", "human_utc")].Value,
                ".0f",
            ),
        )
    if ("Cells", "utc_threshold") in sample_stats.index:
        stats_obj.insert(
            1,
            Metric(
                "Cells",
                "Unique Transcript Counts Threshold",
                "utc_threshold",
                sample_stats.loc[("Cells", "utc_threshold")].Value,
                ".0f",
            ),
        )
    return stats_obj


def get_complexity_metrics(sample_stats: pd.DataFrame) -> list[Metric]:
    """
    Make list of Metric objects that represent extrapolated complexity stats

    Args:
        sample_stats: Sample-level statistics
    """
    pattern = r"^target_reads_(\d+)"
    obj = [
        (
            Metric(
                "Extrapolated Complexity",
                f"Median unique transcript counts at {re.match(pattern, name).group(1)} total reads per cell",
                f"target_reads_{name}",
                row["Value"],
                ".0f",
            )
        )
        for name, row in sample_stats.loc["Extrapolated Complexity"].iterrows()
    ]
    return obj


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outDir", required=True, type=Path, help="Output directory")
    parser.add_argument(
        "--sampleId", required=True, help="Unique ID of the sample during the workflow (e.g. SampleName.LibName)"
    )
    parser.add_argument("--sampleMetrics", required=True, type=Path)
    parser.add_argument("--scalePlexMetrics", required=False, type=Path)
    parser.add_argument("--scalePlexAllBarcodes", required=False, type=Path)
    parser.add_argument("--scalePlexStats", required=False, type=Path)
    parser.add_argument("--scalePlexLibName", required=False, type=str)
    parser.add_argument("--libraryStruct", required=True, type=Path)
    parser.add_argument("--internalReport", action="store_true", default=False)
    parser.add_argument("--isBarnyard", action="store_true", default=False)
    parser.add_argument("--isScalePlex", action="store_true", default=False)
    parser.add_argument("--cellFinder", action="store_true", default=False)
    parser.add_argument("--sampleName", help="Metadata for report")
    parser.add_argument("--libName", default="NA", help="Metadata for report")
    parser.add_argument("--barcodes", default="NA", help="Metadata for report")
    parser.add_argument("--addOutDir", help="Add output directory path to report")
    parser.add_argument("--workflowVersion", default="NA", help="Metadata for report")
    args = parser.parse_args()

    metadata = {
        "Sample Name": args.sampleName or args.sampleId,  # Default to ID if no name is given
        "Library Name": args.libName,
        "Sample Barcodes": args.barcodes,
        "Workflow Version": args.workflowVersion,
        "Analysis Date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        "Library Structure": args.libraryStruct.stem,
    }
    if args.addOutDir:
        output_dir_lines = reporting.split_string_into_lines(args.addOutDir)
        # Join the lines with HTML line breaks
        formatted_output_dir = "<wbr>".join(output_dir_lines)
        metadata["Output Directory"] = f"<span>{formatted_output_dir}</span>"

    build_sample_report(
        write_dir=args.outDir,
        sample_metrics=args.sampleMetrics,
        sample_id=args.sampleId,
        internal_report=args.internalReport,
        is_barnyard=args.isBarnyard,
        is_scale_plex=args.isScalePlex,
        cellfinder=args.cellFinder,
        lib_struct_json=args.libraryStruct,
        metadata=metadata,
        hash_metrics=args.scalePlexMetrics,
        hash_stats=args.scalePlexStats,
        hash_lib_name=args.scalePlexLibName,
        hash_all_cells=args.scalePlexAllBarcodes,
    )


if __name__ == "__main__":
    main()
