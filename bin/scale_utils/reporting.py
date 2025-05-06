"""Utilities related to QC report generation"""

import statistics
import warnings
from pathlib import Path
from typing import List, Optional
import datapane as dp
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
import plotly.graph_objects as go


# Plotly express figure style to be used for all figures
DEFAULT_FIGURE_STYLE = "none"

# Number of cells sampled for scatter plots
SAMPLING_NUMBER = 4000

# Color mapping used for qc filter categorized scatter plots
QC_COLORMAP = {True: "rgb(39, 139, 176)", False: "rgb(233, 237, 245)"}
SCATTER_COLORMAP = {"Cell": "rgb(39, 139, 176)", "Background": "rgb(178, 181, 187)"}

# Color mapping used for barnyard species categorization scatter plots
# (human and mouse hardcoded as it is the only barnyard genome used
# [could be changed])
BARNYARD_COLORMAP = {
    "None": "rgb(179, 188, 201)",
    "Ambig": "rgb(201, 147, 166)",
    "Mixed": "rgb(242, 5, 33)",
    "Human": "rgb(36, 36, 227)",
    "Mouse": "rgb(27, 99, 25)",
}


def formatNumericVal(val: int | float) -> str | float:
    """
    Pretty-print a integer value
    """
    if isinstance(val, int):
        return f"{val:,}"
    elif isinstance(val, str) and val.isdigit():
        return f"{int(val):,}"
    # negative numbers
    elif isinstance(val, str) and len(val) != 0 and val[0] == "-" and val[1:].isdigit():
        return f"{int(val):,}"
    # display nan as NA
    elif isinstance(val, str) and val.startswith("nan"):
        return "NA"
    elif isinstance(val, float):
        return round(val, 2)
    else:
        return val


def styleTable(styler, title: str, hideColumnHeaders=False, boldColumn=None, numericCols=None):
    """
    Function to modify given @pd.DataFrame.Styler

    Args:
        styler ():
        title (str):
        hideColumnHeaders (bool):
        boldColumn ():
        numericCols ():

    Returns:
        Object containing modified styler
    """
    if numericCols is not None:
        styler.format(formatNumericVal, subset=numericCols)

    styler.hide(axis="index")

    if hideColumnHeaders:
        styler.hide(axis="columns")
    else:
        styler.set_table_styles([{"selector": "th", "props": [("border", "black solid !important")]}], overwrite=False)

    if boldColumn is not None:
        styler.set_properties(subset=boldColumn, **{"font-weight": "bold"})

    if title != "":
        styler.set_caption(title)
        styler.set_table_styles(
            [{"selector": "caption", "props": [("border", "black solid !important"), ("font-weight", "bold")]}],
            overwrite=False,
        )

    styler.set_properties(**{"border-color": "black", "border-style": "solid !important"})
    return styler


def highlight_below_threshold_values(row: pd.Series, obj_list: list) -> str:
    """
    Highlight values below threshold in metrics table

    Args
        row: Dataframe row which has metric which we'll check to see if it is below threshold
        obj_list: List of instances of the Metrics class

    Returns:
        CSS style
    """
    metric = row["Metric"]
    for obj in obj_list:
        # Can't do == here because if @obj.warning is true, then @metric will contain the html string
        # needed to display the warning hover text
        if obj.display_name in metric:
            if obj.warning:
                return "color: red;"
    return "color: black;"


def split_string_into_lines(s: str, length: int = 30) -> list[str]:
    """
    Split a string into smaller strings of a certain length

    Args:
        s: String to split
        length: Maximum length of the string beyond which it'll be split
    """
    return [s[i : i + length] for i in range(0, len(s), length)]


def create_metric_table(df: pd.DataFrame, title: str, obj_list: list = [], rm_nan: bool = False) -> dp.Table:
    """
    Create a datapane tabel from a dataframe with a set of metrics (name / value)

    Args:
        df: Dataframe containing metrics with columns 'Metric' and 'Value'
        title: Title of the table
        rm_nan: Exclude rows with "nan" (string) value
        obj_list: List of instances of the Metrics class

    Returns:
        A datapane table object
    """
    if rm_nan:
        # nan is applicable when function is called on non transformed metrics
        # NA is applicable when function is called on transformed metrics
        df = df[~df["Value"].str.contains("nan|NA", na=False)]
    metrics = df[["Metric", "Value"]]
    style = metrics.style.apply(
        lambda row: highlight_below_threshold_values(row, obj_list), axis=1, result_type="broadcast"
    ).pipe(styleTable, title=title, hideColumnHeaders=True, boldColumn="Metric")
    return make_table(style)


def make_table(styler, label=None) -> dp.Table:
    """
    Function that creates wrapper around dp.Table which suppresses
    warning that arises from  code within dp.Table calling Styler.render():
    'this method is deprecated in favour of `Styler.to_html()`'

    Args:
        styler ():
        label ():

    Returns:
        dp.Table object
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        return dp.Table(styler, label=label)


def wellStringComp(x: str, y: str) -> bool:
    """
    Order well coordinates (1A, 12H) by column (number) then row (letter) or vice versa

    Args:
        x: First well coordinate
        y: Second well coordinate

    Returns:
        x < y
    """
    xNum = int(x[0:-1])
    xLetter = x[-1]
    yNum = int(y[0:-1])
    yLetter = y[-1]
    return (xNum, xLetter) < (yNum, yLetter)


def getMaxWellNumberAndLetter(fname):
    """
    Get maximum well coordinate from reference file

    Args:
        fname (str): Path to file

    Returns:
        Maximum letter and maximum number corresponding to the max well coordinate
    """
    max_letter = "A"
    max_number = 1
    with open(fname) as f:
        for line in f:
            line = line.strip()
            split_line = line.split("\t")
            letter = split_line[0][-1]
            numbers = int(split_line[0][:-1])
            max_letter = letter if letter > max_letter else max_letter
            max_number = numbers if numbers > max_number else max_number
    return max_letter, max_number


def getCharacterIndices(lower_limit: int, upper_limit: int) -> List[chr]:
    """
    Generate a list of letter indices from a range of ascii values
    """
    return [chr(i) for i in range(lower_limit, upper_limit)]


def sparseLogCoords(elements: int, pointsPerDecade: int = 250) -> np.ndarray:
    """
    Returns a subset of indicies to evenly sample a vector in log-space
    For selecting x-coords in a log-scale line-plot

    Args:
        elements: Number of elements
        pointsPerDecade: Max. number of elements from each log10 step

    Returns:
        List of selected indices in ascending order
    """
    if elements <= 0:
        return np.empty(0)
    n = int(np.log10(elements) * pointsPerDecade)
    vals = np.logspace(0, np.log10(elements), n)
    inds = np.unique([min(round(v), elements - 1) for v in vals])
    return inds


def buildPlatePlot(counts, title, threshold, what="Cells", subtitle=""):
    """
    Build plate like plot for displaying information on a per well basis

    Args:
        counts (pd.DataFrame): 96 well plate data as DataFrame
        title (str): Title of the plot
        threshold (float): Range till which the colorbar is linear, after that it's logscale
        subtitle (str): Display a subtitle for heatmap

    Returns:
        plotly.graph_objects.Figure
    """
    max_val = counts.to_numpy().max()
    # Create a custom colormap
    colors = [(0, 0, 0), (128 / 255, 170 / 255, 255 / 255)]  # RGB values for black and shade of blue
    positions = [0, 1]  # Corresponding positions for 0 and 1 in the colormap
    cmap = LinearSegmentedColormap.from_list("CustomColormap", list(zip(positions, colors)))
    # Log-transform the colormap
    norm = SymLogNorm(linthresh=threshold, vmin=0, vmax=max(10, max_val), clip=True)
    # Create plotly compatible colorscale
    # https://plotly.github.io/plotly.py-docs/generated/plotly.graph_objects.Heatmap.html
    positions = np.linspace(0, max_val, 256)
    colors = ["rgb" + str((int(r * 255), int(g * 255), int(b * 255))) for r, g, b, _ in cmap(norm(positions))]
    colorscale = list(zip(np.linspace(0, 1, 256), colors))

    # create custom tickvals
    tickvals = [0]
    ticktext = ["0"]
    if max_val != 0:
        # create tickvals for log scale up to power of 10 that is below max_val
        nztickvals = [10**i for i in range(0, int(np.ceil(np.log10(max_val + 1))))]
        nzticktext = [f"10<sup>{i}</sup>" for i, _ in enumerate(nztickvals)]
        # take last two tick labels because otherwise the lower ones bunch up together
        tickvals += nztickvals[-2:]
        ticktext += nzticktext[-2:]

    fig = go.Figure(
        data=go.Heatmap(
            z=counts[::-1],
            x=[str(col) for col in counts.columns],
            y=counts.index[::-1],
            hoverongaps=False,
            hovertemplate=f"Well: %{{x}}%{{y}}<br>{what}: %{{z}}<extra></extra>",
            colorscale=colorscale,
            colorbar=dict(title=what, tickvals=tickvals, ticktext=ticktext),
            xgap=1,
            ygap=1,
            zmin=0,
            zmax=max_val,
        )
    )
    # Update layout to change axis font size
    subtitle_html = f"<br><sub>{subtitle}</sub>" if subtitle else ""
    fig.update_layout(
        title=dict(text=f"{title}{subtitle_html}", y=0.9, x=0.5, xanchor="center", yanchor="top", font=dict(size=20)),
        xaxis=dict(tickfont=dict(size=13), side="top"),
        yaxis=dict(tickfont=dict(size=13)),
    )
    return fig


def buildDfForSamplePlatePlot(libJson, index, referencesPath):
    well_dict = {}
    lib_json_entry_dict = None
    for entry in libJson["barcodes"]:
        if "alias" in entry:
            if entry["alias"] == index.split("_")[0]:
                lib_json_entry_dict = entry
    assert lib_json_entry_dict is not None, f"Could not find barcode {index.split('_')[0]} in libJson"
    max_letter, max_number = getMaxWellNumberAndLetter(referencesPath / f'{lib_json_entry_dict["sequences"]}')
    wellPlateCellCountDf = pd.DataFrame(
        0, columns=range(1, max_number + 1), index=getCharacterIndices(65, ord(max_letter) + 1)
    )
    wellPlateNumCellDf = wellPlateCellCountDf.copy()
    for i in range(65, ord(max_letter) + 1):
        for j in range(1, max_number + 1):
            key = str(j) + chr(i)
            well_dict[key] = []
    return wellPlateCellCountDf, wellPlateNumCellDf, well_dict


def barcodeLevelPlots(
    libJson: dict,
    libStructDir: Path,
    sampleId: str,
    cells: pd.DataFrame,
    index: str,
    title: str,
    internalReport: bool,
    writeDir: Optional[Path] = None,
) -> dp.Group:
    """
    Function to compute statistics to build plate like plot for
    displaying statistics on a per well basis

    Args:
        cells: Metrics per cell-barcode
        index: Column name to sort dataframe by
        title: Title of the plot
        internalReport: Extra plots for internal QC
        writeDir: Output directory

    Returns:
        dp.Group object containing all the plots
    """
    wellPlateCellCountDf, wellPlateNumCellDf, well_dict = buildDfForSamplePlatePlot(libJson, index, libStructDir)

    num_cells_dict = cells[index].value_counts().to_dict()

    for idx, row in cells.iterrows():
        letter = row[index][-1]
        numbers = row[index][:-1]
        try:
            well_dict[numbers + letter].append(row["counts"])
        except KeyError as e:
            print(f"{e}: {letter}{numbers} does not exist")

    for well in num_cells_dict:
        letter = well[-1]
        numbers = well[:-1]
        try:
            wellPlateNumCellDf.at[letter, int(numbers)] = num_cells_dict[well]
        except KeyError as e:
            print(f"{e}: {letter}{numbers} doesn't exist")

    for key in well_dict:
        letter = key[-1]
        numbers = key[:-1]
        if len(well_dict[key]) == 0:
            wellPlateCellCountDf.at[letter, int(numbers)] = 0
        else:
            wellPlateCellCountDf.at[letter, int(numbers)] = int(statistics.median(well_dict[key]))

    readsPerIndexBox = buildPlatePlot(
        wellPlateCellCountDf, "Unique Transcript Counts Per Cell", 100.0, what="Transcripts"
    )
    cellsPerIndexBar = buildPlatePlot(wellPlateNumCellDf, "Number of cells", 1.0)
    namePrefix = title.replace(" ", "_")
    if internalReport:
        cellsPerIndexBar.write_image(
            writeDir / "figures_internal" / f"{sampleId}_CellCount_By_{namePrefix}_Heatmap.png"
        )
        readsPerIndexBox.write_image(
            writeDir / "figures_internal" / f"{sampleId}_UniqueTranscriptCount_By_{namePrefix}_Heatmap.png"
        )
    wellPlateCellCountDf.to_csv(writeDir / "csv" / f"{sampleId}_unique_transcript_counts_by_{namePrefix}_well.csv")
    wellPlateNumCellDf.to_csv(writeDir / "csv" / f"{sampleId}_num_cells_by_{namePrefix}_well.csv")
    return dp.Group(dp.Text(f"## {title}"), dp.Group(dp.Plot(cellsPerIndexBar), dp.Plot(readsPerIndexBox), columns=2))
