#!/usr/bin/env python
"""
Generate the sample report (HTML and .csv) from per-sample metrics
"""

import argparse
import datetime
import json
from pathlib import Path
from typing import Dict,List

import datapane as dp
import numpy as np
import pandas as pd
from plotly.graph_objects import Figure
import plotly.express as px

from scale_utils import reporting
from scale_utils.base_logger import logger
import re


def build_sample_report(write_dir: Path, sample_metrics: Path, sample_id: str, internal_report: bool, is_barnyard: bool,
                        lib_struct_json: Path, trim_stats: List[Path], metadata:Dict[str, str]):
    """
    Entry function for building sample report

    Args:
        write_dir: Report Output directory
        sample_metrics: Folder containing csv files defining the metrics for a single sample
        sample_id: Unique Sample ID (e.g. `SampleNam.LibName`)
        internal_report: Create additional plots for internal QC
        is_barnyard: Add barnyard plots and metrics
        lib_struct_json: Library structure json (in directory with sequence files, etc.)
        trim_stats: json files containing cutadapt output
        metadata: Information about the sample and analysis run [FieldName -> Value]
    """
    all_cells = pd.read_csv(sample_metrics / "allCells.csv", index_col=0)
    sample_stats = pd.read_csv(sample_metrics / "sample_stats.csv", index_col=['Category', 'Metric'])
    filtered_short_reads, avg_trimmed_length = get_trim_stats_metrics(trim_stats)

    Path(write_dir, "csv").mkdir(parents=True, exist_ok=True)
    if internal_report:
        Path(write_dir, "figures_internal").mkdir()
    
    lib_struct_dir = lib_struct_json.parent
    lib_struct = json.load(open(lib_struct_json))

    pages = []
    reads_page, report_stats, internal_report_blocks = \
        build_reads_page(sample_id, all_cells, sample_stats, write_dir, internal_report,
                         filtered_short_reads, avg_trimmed_length, metadata)
    pages.append(reads_page)

    all_index_plots, barcodes_page = build_barcodes_page(lib_struct, lib_struct_dir, all_cells, write_dir,
                                                         sample_id, internal_report)
    pages.append(barcodes_page)
    
    if internal_report:
        internal_report_blocks.extend([all_index_plots["lig"], all_index_plots["pcr"]])
        internal_report_page = dp.Page(dp.Group(blocks=internal_report_blocks, columns=2), title="InternalReport")
        pages.append(internal_report_page)
    
    if is_barnyard:
        score_barn(all_cells, sample_stats)
        barnyard_stats = make_barnyard_stats(all_cells)
        barnyard_page = build_barnyard_page(all_cells, barnyard_stats, write_dir, sample_id, internal_report)
        report_stats = pd.concat([report_stats, barnyard_stats])
        pages.append(barnyard_page)

    report = dp.Report(blocks=pages)

    logger.debug(f"Writing reportStatistics csv, report to {str(write_dir.resolve())}")
    # Change 'nan' back to np.nan in Value column
    report_stats.replace(to_replace=['nan', 'nan%'], value=np.nan, inplace=True)
    report_stats.to_csv(write_dir / "csv"/ f"{sample_id}.reportStatistics.csv",
                        index=False, header=False, na_rep='', # can change na_rep to '<NA>' or similar if needed
                        columns=["Category", "Metric", "Value"])

    report.save(write_dir / f"{sample_id}.report.html")


def get_trim_stats_metrics(trim_stats: list) -> tuple[float, float]:
    """
    Parse cutadapt output and aggregate metrics

    Args:
        trim_stats: List of paths to trim stats json
    
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


def get_background(minor_fracs: list) -> tuple[float, float]:
    """
    Calculate median and standard deviation for barnyard cells

    Args:
        minor_fracs: List with fraction indicating percentage
            of human or mouse samples in a cell

    Returns:
        Median and standard deviation
    """
    if minor_fracs.size == 0:
        return np.nan, np.nan
    median = minor_fracs.median()
    std = (np.quantile(minor_fracs, 0.75) - median) / 0.675
    return median, std


def score_barn(all_cells: pd.DataFrame, sample_stats: pd.DataFrame):
    """
    Label each cell as belonging to either human or mouse

    Args:
        all_cells: Barcode dataframe for this sample
        sample_stats: Sample-level statistics
    """
    all_cells['species'] = 'None'
    all_cells['minor_frac'] = all_cells[['human_umis', 'mouse_umis']].min(axis=1) / all_cells[['human_umis', 'mouse_umis']].sum(axis=1)
    all_cells.loc[all_cells['pass'] & (all_cells.human_umis > all_cells.mouse_umis), 'species'] = 'Human'
    all_cells.loc[all_cells['pass'] & (all_cells.mouse_umis > all_cells.human_umis), 'species'] = 'Mouse'

    min_human = min_mouse = all_cells.loc[all_cells['pass'], 'umis'].min()
    if (all_cells.species == 'Human').any() and (all_cells.species == 'Mouse').any():
        min_human = np.percentile(all_cells.human_umis[all_cells.species == 'Human'], 10)
        min_mouse = np.percentile(all_cells.mouse_umis[all_cells.species == 'Mouse'], 10)

    # Labelling low and doublet cells
    all_cells.loc[(all_cells.human_umis < min_human) & (all_cells.mouse_umis < min_mouse), 'species'] = "None"
    all_cells.loc[(all_cells.human_umis >= min_human) & (all_cells.mouse_umis >= min_mouse), 'species'] = "Mixed"

    # Labelling as Ambiguous
    human_bg_med, human_bg_std = get_background(all_cells[(all_cells.species == 'Mouse')].minor_frac)
    mouse_bg_med, mouse_bg_std = get_background(all_cells[(all_cells.species == 'Human')].minor_frac)

    all_cells.loc[(all_cells.species == 'Mixed') & (all_cells.human_umis > all_cells.mouse_umis) & (all_cells.minor_frac < mouse_bg_med + 2 * mouse_bg_std), 'species'] = 'Ambiguous'
    all_cells.loc[(all_cells.species == 'Mixed') & (all_cells.mouse_umis > all_cells.human_umis) & (all_cells.minor_frac < human_bg_med + 2 * human_bg_std), 'species'] = 'Ambiguous'
    all_cells.loc[(all_cells.species == 'Human') & (all_cells.minor_frac >= max(0.1, mouse_bg_med + 3 * mouse_bg_std)), 'species'] = 'Ambiguous'
    all_cells.loc[(all_cells.species == 'Mouse') & (all_cells.minor_frac >= max(0.1, human_bg_med + 3 * human_bg_std)), 'species'] = 'Ambiguous'


def make_read_stats(sample_stats: pd.DataFrame, filtered_short_reads: float):
    """
    Make dataframe depicting read stats that will be displayed in the sample report

    Args:
        sample_stats: Sample-level statistics
        filtered_short_reads: Number of filtered reads (or NaN)
    """
    totalReads = sample_stats.loc[('Reads', 'total_reads')].Value + filtered_short_reads
    stats = [
        ('Reads', 'Total Sample Reads', f"{totalReads:.0f}"),
        ('Reads', 'Passing Sample Reads', f"{sample_stats.loc[('Reads', 'total_reads')].Value:.0f}"),
        ('Reads', 'Reads Mapped to Genome', f"{sample_stats.loc[('Reads', 'mapped_reads_perc')].Value:.1%}"),
        ('Reads', 'Reads Mapped to Transcriptome', f"{sample_stats.loc[('Reads', 'gene_reads_perc')].Value:.1%}"),
        ('Reads', 'Exonic Reads', f"{sample_stats.loc[('Reads', 'exon_reads_perc')].Value:.1%}"),
        ('Reads', 'Antisense Reads', f"{sample_stats.loc[('Reads', 'antisense_reads_perc')].Value:.1%}"),
        ('Reads', 'Mitochondrial Reads', f"{sample_stats.loc[('Reads', 'mito_reads_perc')].Value:.1%}"),
        ('Reads', 'Saturation', f"{sample_stats.loc[('Reads', 'saturation')].Value:.2f}")
    ]
    return pd.DataFrame.from_records(stats, columns=['Category', 'Metric', 'Value'])


def make_extra_read_stats(average_trimmed_length: float):
    """
    Create a dataframe with metrics included in .csv files, but not the HTML report tables
    
    Args:
        average_trimmed_length: Average read length after trimming
    """
    stats = [
        ('Reads', 'Average Trimmed Read Length', f"{average_trimmed_length:.1f}"),
    ]
    return pd.DataFrame.from_records(stats, columns=['Category', 'Metric', 'Value'])


def build_reads_page(sample: str, all_cells: pd.DataFrame, sample_stats: pd.DataFrame, write_dir: Path, internal_report: bool,
                   filtered_short_reads: float, avg_trimmed_length: float, metadata:Dict[str, str]):
    """
    Build a datapane page that mainly shows plots related to reads

    Args:
        sample: Sample name
        all_cells: Per-cell metrics
        sample_stats: Sample-level statistics
        write_dir: Output directory
        internal_report: Flag indicating whether report is to be
            generated for internal r&d run
        filtered_short_reads: Number of reads filtered out after trimming (or NaN)
        avg_trimmed_length: Average read length after trimming (or NaN)
        metadata: Information about the sample and analysis run [FieldName -> Value]

    Returns:
        dp.Page object and concatenated dataframe of various statistics
        List of elements for the internal report (or empty)
    """
    report_stats = pd.DataFrame({'Category': ['Sample'], 'Metric': ['SampleName'], 'Value': [sample]})
    read_stats = make_read_stats(sample_stats, filtered_short_reads)
    extra_stats = make_extra_read_stats(avg_trimmed_length)
    cell_stats = make_cell_stats(sample_stats)
    complexity_stats = make_complexity_stats(sample_stats)
    report_stats = pd.concat([report_stats, read_stats, extra_stats, cell_stats, complexity_stats])

    info = pd.DataFrame.from_dict({'Metric':metadata.keys(), 'Value':metadata.values()})
    info_table = reporting.create_metric_table(info, title="Sample Information")
    cell_table = reporting.create_metric_table(cell_stats, title="Cell Metrics")
    mapping_table = reporting.create_metric_table(read_stats, title="Read Metrics", rm_nan=True)

    rank_plot = make_rank_plot(all_cells, write_dir, sample, internal_report)
    unique_reads_fig = make_unique_reads_fig(sample_stats)
    genes_umi_scatter = build_gene_read_scatter(all_cells)
    saturation_scatter = build_saturation_scatter(all_cells)
    group_blocks = [info_table, cell_table, mapping_table, rank_plot, unique_reads_fig, genes_umi_scatter, saturation_scatter]
    page = dp.Page(dp.Group(blocks=group_blocks, columns=2), title="Summary")

    internal_report_blocks = []
    if internal_report:
        gene_reads_dist = build_dist("Distribution of Reads: Gene Matched", 'geneReads', "Fraction of Genic Reads", all_cells)
        antisense_reads_dist = build_dist("Distribution of Reads: Antisense", 'antisenseReads', "Fraction of Antisense Reads", all_cells)
        exonic_reads_dist = build_dist("Distribution of Reads: Exonic", 'exonReads', "Fraction of Exonic Reads", all_cells)
        internal_report_blocks.extend([gene_reads_dist, exonic_reads_dist, antisense_reads_dist])

        if sample_stats.loc[('Cells', 'mito_reads')].Value != 0:
            mito_reads_dist = build_dist("Distribution of Reads: Mitochondrial", 'mitoReads', "Fraction of Mitochondrial Reads", all_cells)
            internal_report_blocks.append(mito_reads_dist)

    return page, report_stats, internal_report_blocks


def build_barcodes_page(lib_struct: dict[str, str], lib_struct_dir: Path, all_cells: pd.DataFrame,
                        write_dir: Path, sample_id: str, internal_report: bool) -> tuple[dict, dp.Page]:
    """
    Build datapane page for barcodes

    Args:
        lib_struct: Library structure definition
        lib_struct_dir: Directory with barcode sequences ec.
        all_cells: Metrics per cell-barcode
        write_dir: Ouput directory
        sample_id: Unique ID of this sample
        internal_report: Extra plots for internal QC

    Returns:
        Index plots for internal reports (empty otherwise)
        Page object with shared figures
    """
    blocks_to_render = []
    passing_cells = all_cells.loc[all_cells['pass']] # .loc needed here to keep the column names if zero passing cells
    all_index_plots = {}
    for bc in lib_struct['barcodes']:
        if bc.get('type') not in ['library_index', 'umi']:
            alias = bc.get('alias') or bc['name']
            index_plots = reporting.barcodeLevelPlots(
                lib_struct, lib_struct_dir, sample_id, passing_cells, f"{alias}_alias", 
                f"{alias} Index", internal_report, write_dir)
            all_index_plots[bc['name']] = index_plots
    blocks_to_render.append(all_index_plots["rt"])

    return all_index_plots, dp.Page(blocks=blocks_to_render, title="Barcodes")


def build_dist(title: str, field: str, label: str, all_cells: pd.DataFrame) -> dp.Plot:
    """
    Build histogram for a given field divided by mapped reads

    Args:
        title: Title of plot
        field: Column name in all_cells to plot on x axis after dividing by mapped reads
        label: String that will be the x axis label
        all_cells: Cell-level metrics

    Returns:
        Plot object with the figure
    """
    passing_cells = all_cells.loc[all_cells['pass']]
    x = passing_cells[field] / passing_cells['mappedReads']

    fig = px.histogram(x=x, title=title, histnorm='percent',
                       template=reporting.DEFAULT_FIGURE_STYLE, labels={field: label})

    fig.update_layout(yaxis_title="% Cells")

    return dp.Plot(fig)


def subsample_all_cells(all_cells: pd.DataFrame) -> pd.DataFrame:
    """
    Subsample allCells to a maximum of 5000 cells for plotting purposes

    Too many cells results in overplotting so in this function we:
        (1) Take all passing cells
        (2) Take an equal number of the top failing cells by umis
        (3) Return a maximum of 5000 cells by subsampling if necessary
    
    Args:
        all_cells: Cell-level metrics

    Returns:
        Subsampled dataframe
    """
    passing_cells = all_cells.loc[all_cells['pass']]
    top_failed_cells = all_cells[~all_cells['pass']].sort_values(by='umis', ascending=False).head(passing_cells.shape[0])
    top_cells = pd.concat([passing_cells, top_failed_cells])
    if top_cells.shape[0] > 5000:
        top_cells = top_cells.sample(5000)
    return top_cells


def build_gene_read_scatter(all_cells: pd.DataFrame)-> dp.Plot:
    """
    Make a scatter plot of reads vs umi

    Args:
        all_cells: Cell-level metrics

    Returns:
        Plot object with the figure
    """
    all_cells["plt_color"] = all_cells["pass"].map({True:"Cell", False:"Background"})
    fig = px.scatter(subsample_all_cells(all_cells), x="reads", y="genes", color="plt_color",
                     color_discrete_map=reporting.SCATTER_COLORMAP, template=reporting.DEFAULT_FIGURE_STYLE,
                     labels={"reads": "Total reads", "genes": "Genes detected", "plt_color":""}, title="Genes Detected Per Cell",
                     category_orders={"plt_color":["Cell", "Background"]}, opacity=0.5)
    

    return dp.Plot(fig)


def build_saturation_scatter(all_cells: pd.DataFrame)-> dp.Plot:
    """
    Build scatter plot for saturation

    Args:
        all_cells: Cell-level metrics

    Returns:
        Plot object with the figure
    """
    all_cells["plt_color"] = all_cells["pass"].map({True:"Cell", False:"Background"})
    fig = px.scatter(subsample_all_cells(all_cells), x="reads", y="Saturation", color = "plt_color",
                     labels={"reads": "Total reads", "plt_color":""}, template=reporting.DEFAULT_FIGURE_STYLE,
                     color_discrete_map=reporting.SCATTER_COLORMAP, title="Saturation Per Cell", log_x=False,
                     category_orders={"plt_color":["Cell", "Background"]}, opacity=0.5)

    return dp.Plot(fig)


def make_rank_plot(all_cells: pd.DataFrame, write_dir: Path, sample_id: str, internal_report: bool) -> dp.Plot:
    """
    Make barcode rank plot and draw vertical line at threshold if cellfinder was not used
    
    Args:
        all_cells: Cell-level metrics
        write_dir: Output directory
        sample_id: Unique ID of this sample
        internal_report: Extra plots for internal QC

    Returns:
        Plot object with the figure
    """
    is_cellFinder = all_cells['flags'].fillna("").str.contains('cellFinder').any()
    indices = reporting.sparseLogCoords(all_cells.index.size)
    all_cells = all_cells.sort_values(by=["umis"], ascending=False)
    if is_cellFinder:
        #Calculate proportion of cells with a moving window of 10.
        #min_periods=1 allows us to calculate proportion of windows with less than 10 barcodes.
        #This prevents nans.
        all_cells['z'] = all_cells['pass'].rolling(25, min_periods=1).mean()
        #Color plot using new classification and update category orders to fix overplotting.
        all_cells = all_cells.iloc[indices]
        all_cells["x_val"] = indices
        fig = px.scatter(all_cells, x="x_val", y="umis", color = "z", 
                      labels={"x_val":"Cell barcodes", "umis": "Unique transcript counts", "cell": "Cell Barcode", "z": "Prop. Cells"},
                      log_x=True, log_y=True, template=reporting.DEFAULT_FIGURE_STYLE, title="Barcode Rank Plot", opacity=.5, color_continuous_scale=['rgb(178, 181, 187)', 'rgb(39, 139, 176)'],
                      hover_data={"umis":True, "z":True, "x_val":False})
        fig.update_layout(coloraxis_colorbar = dict(thickness=10, tickvals=[0,1], ticktext=["Background", "Cell"], len=.25, y = .75, title_text = ""))
    else:
        all_cells["plt_color"] = all_cells["pass"].map({True:"Cell", False:"Background"})
        fig = px.scatter(all_cells.iloc[indices], x=indices, y="umis", color = "plt_color", labels={"x": "Cell barcodes", "umis": "Unique transcript counts", "plt_color":""},
                      log_x=True, log_y=True, template=reporting.DEFAULT_FIGURE_STYLE, title="Barcode Rank Plot", color_discrete_map=reporting.SCATTER_COLORMAP, category_orders={"plt_color":["Cell", "Background"]}, opacity=0.5,
                      hover_data={"plt_color":False})
    # If cellfinder was not used for cell calling. Add UTC threshold to plot.
    if not is_cellFinder:
        fig.add_vline(x=max(1, all_cells['pass'].sum()), line_dash="dash", line_color="green")
    if indices.size > 0:
        fig.update_layout(xaxis_range=[1, np.log10(max(indices)+1)])
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

def build_barnyard_page(all_cells: pd.DataFrame, barnyard_stats: pd.DataFrame, write_dir: Path, sample_id: str, internal_report: bool) -> dp.Page:
    """
    Build a datapane page for barnyard analysis
    
    Args:
        all_cells: Cell-level metrics
        barnyard_stats: Barnyard statistics
        write_dir: Output directory
        sample_id: Unique ID of this sample
        internal_report: Figures to be published as png for internal QC
        
    Returns:
        Page object with the relevant figures
    """
    barnyard_table = reporting.create_metric_table(barnyard_stats, title="")
    barnyard_plot = make_barnyard_scatter_plot(all_cells, write_dir, sample_id, internal_report)
    plot_group = dp.Group(barnyard_table, barnyard_plot, columns=2)
    return dp.Page(plot_group, title="Barnyard")

def make_barnyard_scatter_plot(all_cells: pd.DataFrame, write_dir: Path, sample_id: str, internal_report: bool) -> dp.Plot:
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
    fig = px.scatter(subsample_all_cells(all_cells),
                     x="human_umis", y="mouse_umis", color="species",
                     log_x=True, log_y=True, template=reporting.DEFAULT_FIGURE_STYLE,
                     color_discrete_map=reporting.BARNYARD_COLORMAP,
                     labels={"human_umis": "Human UMIs", "mouse_umis": "Mouse UMIs"})

    fig.update_traces(selector=dict(name='None'), visible="legendonly")
    if internal_report:
        save_figure_png(write_dir, sample_id, fig, "Barnyard_ScatterPlot.png")

    return dp.Plot(fig)

def make_barnyard_stats(all_cells: pd.DataFrame) -> pd.DataFrame:
    """
    Compute statistics for barnyard samples
    
    Args:
        all_cells: Cell-level metrics

    Returns:
        Sample-level barnyard statistics
    """
    passing_barnyard_cells = all_cells.loc[(all_cells['pass']) & (~all_cells.species.isin(['Ambiguous', 'None']))]
    ncells = max(1, passing_barnyard_cells.index.size)
    prop_human = (passing_barnyard_cells.species == 'Human').sum() / ncells
    prop_mouse = (passing_barnyard_cells.species == 'Mouse').sum() / ncells
    prop_mixed = (passing_barnyard_cells.species == 'Mixed').sum() / ncells
    doublets = prop_mixed/(2*prop_human*prop_mouse) if prop_human > 0 and prop_mouse > 0 else 0
    doublets = min(1, doublets)
    bg = passing_barnyard_cells[passing_barnyard_cells.species != 'Mixed'].minor_frac.median()

    stats = [
        ('Barnyard', 'Human Cells', f"{prop_human:.1%}"),
        ('Barnyard', 'Mouse Cells', f"{prop_mouse:.1%}"), 
        ('Barnyard', 'Mixed Cells', f"{prop_mixed:.1%}"), 
        ('Barnyard', 'Estimated Doublets', f"{doublets:.1%}"),
        ('Barnyard', 'Background', f"{bg:.2%}"),
    ]
    return pd.DataFrame.from_records(stats, columns=['Category', 'Metric', 'Value'])

def make_unique_reads_fig(sample_stats: pd.DataFrame)-> dp.Plot:
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
    target_mean_reads = [(
                int(re.match(r'^target_reads_(\d+)', name).group(1)),
                row['Value']
            ) for name, row in sample_stats.loc['Extrapolated Complexity'].iterrows()]
    
    # Add the observed value to the plot
    target_mean_reads.append((sample_stats.loc[('Cells', 'mean_passing_reads')].Value,
                              sample_stats.loc[('Cells', 'median_utc')].Value))

    plot_df = pd.DataFrame.from_records(target_mean_reads, columns=['x', 'unique_read'])
    plot_df.sort_values(by='x', inplace=True) 

    fig = px.line(plot_df, x="x", y="unique_read",
                  labels={"x": "Total reads per cell",
                          "unique_read": "Median unique transcript counts (extrapolated)"},
                  template=reporting.DEFAULT_FIGURE_STYLE, markers=True,
                  title="Complexity")

    return dp.Plot(fig)


def make_cell_stats(sample_stats: pd.DataFrame) -> pd.DataFrame:
    """
    Make dataframe of summary stats

    Args:
        sample_stats: Sample-level statistics
    """
    stats = [
        ('Cells', 'Cells Called', f"{sample_stats.loc[('Cells', 'cells_called')].Value:.0f}"),
        ('Cells', 'Mean Passing Reads per Cell', f"{sample_stats.loc[('Cells', 'mean_passing_reads')].Value:.0f}"),
        ('Cells', 'Median Unique Transcript Counts per Cell', f"{sample_stats.loc[('Cells', 'median_utc')].Value:.0f}"),
        ('Cells', 'Median Genes per Cell', f"{sample_stats.loc[('Cells', 'median_genes')].Value:.0f}"),
        ('Cells', 'Reads in Cells', f"{sample_stats.loc[('Cells', 'reads_in_cells')].Value:.1%}"),
    ]
    if ('Cells', 'utc_threshold') in sample_stats.index:
        stats.insert(1, ('Cells', 'Unique Transcript Counts Threshold', f"{sample_stats.loc[('Cells', 'utc_threshold')].Value:.0f}"))
    return pd.DataFrame.from_records(stats, columns=['Category', 'Metric', 'Value'])


def make_complexity_stats(sample_stats: pd.DataFrame) -> pd.DataFrame:
    """
    Make dataframe of extrapolated complexity stats

    Args:
        sample_stats: Sample-level statistics
    """
    pattern = r'^target_reads_(\d+)'
    stats = [(
                "Extrapolated Complexity",
                f"Median unique transcript counts at {re.match(pattern, name).group(1)} total reads per cell",
                f"{row['Value']:.0f}"
            ) for name, row in sample_stats.loc['Extrapolated Complexity'].iterrows()]
    return pd.DataFrame.from_records(stats, columns=['Category', 'Metric', 'Value'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outDir", required=True, type=Path, help="Output directory")
    parser.add_argument("--sampleId", required=True, help="Unique ID of the sample during the workflow (e.g. SampleName.LibName)")
    parser.add_argument("--sampleMetrics", required=True, type=Path)
    parser.add_argument("--libraryStruct", required=True, type=Path)
    parser.add_argument("--trim_stats", nargs='+', type=Path)
    parser.add_argument("--internalReport", action="store_true", default=False)
    parser.add_argument("--isBarnyard", action="store_true", default=False)
    parser.add_argument('--sampleName', help='Metadata for report')
    parser.add_argument("--libName", default='NA', help='Metadata for report')
    parser.add_argument("--barcodes", default='NA', help='Metadata for report' )
    parser.add_argument("--workflowVersion", default='NA', help='Metadata for report')
    args = parser.parse_args()

    metadata = {
        'Sample Name': args.sampleName or args.sampleId, # Default to ID if no name is given
        'Library Name': args.libName,
        'Sample Barcodes': args.barcodes,
        'Workflow Version': args.workflowVersion,
        'Analysis Date': datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
        'Library Structure': args.libraryStruct.stem,
    }
    build_sample_report(args.outDir, args.sampleMetrics, args.sampleId, args.internalReport, args.isBarnyard, args.libraryStruct, args.trim_stats, metadata)

if __name__ == "__main__":
    main()
