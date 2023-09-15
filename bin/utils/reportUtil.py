# Utilities related to QC report generation

import functools
import operator
import statistics
import warnings
from pathlib import Path
from typing import List

import datapane as dp
import matplotlib.colors as colors
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


# Plotly express figure style to be used for all figures
DEFAULT_FIGURE_STYLE = "none"

# Number of cells sampled for scatter plots
SAMPLING_NUMBER = 4000

# Color mapping used for qc filter categorized scatter plots
QC_COLORMAP = {True: 'rgb(39, 139, 176)', False: 'rgb(233, 237, 245)'}

# Color mapping used for barnyard species categorization scatter plots
# (human and mouse hardcoded as it is the only barnyard genome used
# [could be changed])
BARNYARD_COLORMAP = {
    "None": 'rgb(179, 188, 201)', "Ambig": 'rgb(201, 147, 166)',
    'Mixed': 'rgb(242, 5, 33)', 'Human': 'rgb(36, 36, 227)',
    'Mouse': 'rgb(27, 99, 25)'}


def formatNumericVal(val:int|float) -> str|float:
    """
    Pretty-print a integer value
    """
    if isinstance(val, int):
        return f"{val:,}"
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
    if (numericCols is not None):
        styler.format(formatNumericVal, subset=numericCols)

    styler.hide(axis='index')

    if hideColumnHeaders:
        styler.hide(axis='columns')
    else:
        styler.set_table_styles([{'selector': 'th', 'props': [('border', 'black solid !important')]}], overwrite=False)

    if boldColumn is not None:
        styler.set_properties(subset=boldColumn, **{'font-weight': 'bold'})

    if title != "":
        styler.set_caption(title)
        styler.set_table_styles([{'selector': 'caption', 'props': [('border', 'black solid !important'), ("font-weight", "bold")]}], overwrite=False)

    styler.set_properties(**{"border-color": 'black', "border-style": 'solid !important'})
    return styler

def createMetricTable(dataframe, title, metricCol='Metric', valueCol='Value'):
    """
    Create a datapane tabel from a dataframe with a set of metrics (name / value)
    """
    df = dataframe[[metricCol, valueCol]]
    style = df.style.pipe(styleTable, title=title, hideColumnHeaders=True, boldColumn=metricCol, numericCols=[valueCol])
    return mkTable(style)

def mkTable(styler, label=None) -> dp.Table:
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
    Order well coordinates (1A, 12H) by column (number) then row (letter)

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
            letter = split_line[1][-1]
            numbers = int(split_line[1][:-1])
            max_letter = letter if letter > max_letter else max_letter
            max_number = numbers if numbers > max_number else max_number
    return max_letter, max_number

def getCharacterIndices(lower_limit:int, upper_limit:int) -> List[chr]:
    """
    Generate a list of letter indices from a range of ascii values
    """
    return [chr(i) for i in range(lower_limit, upper_limit)]

def sparseLogCoords(elements: int, pointsPerDecade:int = 250) -> np.ndarray:
    """
    Returns a subset of indicies to evenly sample a vector in log-space
    For selecting x-coords in a log-scale line-plot

    Args:
        elements: Number of elements
        pointsPerDecade: Max. number of elements from each log10 step

    Returns:
        List of selected indices in ascending order
    """
    if elements <= 0: return np.empty(0)
    n = int(np.log10(elements) * pointsPerDecade)
    vals = np.logspace(0, np.log10(elements), n)
    inds = np.unique([min(round(v), elements-1) for v in vals])
    return inds

def buildPlatePlot(df, title, threshold, colorbar_title=None):
    """
    Build plate like plot for displaying information on a per well basis

    Args:
        df (pd.DataFrame): Data to plot
        title (str): Title of the plot
        threshold (float): Range till which the colorbar is linear, after that it's logscale
        colorbar_title (str): Colorbar title

    Returns:
        Matplotlib figure
    """
    fig = plt.figure()
    cbar_args = {}
    if colorbar_title:
        cbar_args['label'] = colorbar_title
    norm = colors.SymLogNorm(linthresh=threshold, vmin=0, vmax=max(threshold*10, df.values.max()))
    ax = sns.heatmap(df, linewidth=0.5, cmap=sns.color_palette("dark:#80AAFF", as_cmap=True), cbar_kws=cbar_args, norm=norm)
    ax.set_title(title)

    return fig

def buildDfForSamplePlatePlot(libJson, index, referencesPath):
    well_dict = {}
    for entry in libJson["barcodes"]:
        if "alias" in entry:
            if entry["alias"] == index.split("_")[0]:
                lib_json_entry_dict = entry
    max_letter, max_number = getMaxWellNumberAndLetter(referencesPath / f'{lib_json_entry_dict["sequences"]}')
    wellPlateCellCountDf = pd.DataFrame(0, columns=range(1, max_number+1), index=getCharacterIndices(65,ord(max_letter)+1))
    wellPlateNumCellDf = wellPlateCellCountDf.copy()
    for i in range(65, ord(max_letter)+1):
        for j in range(1, max_number+1):
            key = str(j)+chr(i)
            well_dict[key] = []
    return wellPlateCellCountDf, wellPlateNumCellDf, well_dict

def barcodeLevelPlots(referencesPath: Path, sampleName: str, cells: pd.DataFrame,
                        index: str, title: str, internalReport: bool, libJson: dict,
                        writeDir=None) -> dp.Group:
    """
    Function to compute statistics to build plate like plot for
    displaying statistics on a per well basis

    Args:
        cells (pd.DataFrame): Dataframe containing information to plot
        index (str): Column name to sort dataframe by
        title (str): Title of the plot
        internalReport (bool): Whether report is being generated for internal purposes
        writeDir (str): Write directory

    Returns:
        dp.Group object containing all the plots
    """
    wellPlateCellCountDf, wellPlateNumCellDf, well_dict = buildDfForSamplePlatePlot(libJson, index, referencesPath)
    
    num_cells_dict = cells[index].value_counts().to_dict() 
    
    for idx, row in cells.iterrows():
        letter = row[index][-1]
        numbers = row[index][:-1]
        try:
            well_dict[numbers+letter].append(row['umis'])
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

    readsPerIndexBox = buildPlatePlot(wellPlateCellCountDf, "Unique Transcript Counts Per Cell", 100.0)
    cellsPerIndexBar = buildPlatePlot(wellPlateNumCellDf, "Number of cells", 1.0)
    namePrefix = title.replace(" ", "_")
    if writeDir is not None and internalReport:
        cellsPerIndexBar.savefig(writeDir / f"{sampleName}_figures" / f"CellCount_By_{namePrefix}_Heatmap.png")
        readsPerIndexBox.savefig(writeDir / f"{sampleName}_figures" / f"UniqueTranscriptCount_By_{namePrefix}_Heatmap.png")
    wellPlateCellCountDf.to_csv(writeDir / "csv" / f"{sampleName}_unique_transcript_counts_by_{namePrefix}_well.csv")
    wellPlateNumCellDf.to_csv(writeDir / "csv" / f"{sampleName}_num_cells_by_{namePrefix}_well.csv")
    return dp.Group(dp.Text(f'## {title}'), dp.Group(dp.Plot(cellsPerIndexBar), dp.Plot(readsPerIndexBox), columns=2))
