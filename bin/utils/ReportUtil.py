#!/usr/bin/env python
"""
Python script to define utility functions that will be
used by other scripts
"""
import numpy as np
import json
from typing import Dict, Union, List
import pandas as pd
import seaborn as sns
from os import mkdir
from os.path import exists
from pathlib import Path
import warnings
import datapane as dp
import statistics
from math import log10
import functools
import operator
import utils.myconstants as constants
from collections import OrderedDict
import plotly.express as px
from dataclasses import dataclass
from matplotlib import pyplot as plt


def compare(this, that):
    """
    Function to compare two variables

    Args:
        this (object): First variable to compare
        that (object): Second variable to compare

    Returns:
        1 if @this greater than @that
        -1 if @this less than @that
        0 if both are equal
    """
    if this > that:
        return 1
    elif this < that:
        return -1
    else:
        return 0


@dataclass
class CellStat:
    """
    Class defining the attributes of the @cells dictionary

    Attributes:
        genes (int): Total number of genes
        humanGenes (int): Number of human genes
        mouseGenes (int): Number of mouse genes
        umis (int): Total number of UMIs
        humanUmis (int): Number of human UMIs
        mouseUmis (int): Number of mouse UMIs
    """
    genes: int = 0
    humanGenes: int = 0
    mouseGenes: int = 0
    umis: int = 0
    humanUmis: int = 0
    mouseUmis: int = 0


class GeneralUtils:
    """
    Class containing general purpose utility functions related to i/o,
    input parsing or pandas manipulations
    """
    def __init__(self):
        self = self

    @staticmethod
    def loadStatsAsDict(fn):
        """
        Function to create a dictionary from a file (@fn)
        containing tab delimited statistics

        Args:
            fn (str): Path to file

        Returns:
            Dictionary containing stats
        """
        res = {}

        for line in open(fn):
            split_line = line.strip().split("\t")

            if not len(split_line) >= 2:
                continue

            stat = split_line[0]
            val = split_line[1]

            try:
                val = int(val)

            except ValueError:
                try:
                    val = float(val)
                except ValueError:
                    pass

            res[stat] = val

        return res

    @staticmethod
    def loadStarStatsFromCsv(fn, features):
        """
        Function to create a dictionary with information
        obtained from summary.csv outputted by STAR

        Args:
            fn (str): Path to summary.csv
            features (str): Feature type

        Returns:
            Dictionary containing the statistics
        """
        stats = {}

        for line in open(fn):
            split_line = line.strip().split(',')
            stats[split_line[0]] = split_line[1]

        starStats = {}
        starStats['Total Reads'] = f"{int(stats['Number of Reads']):,}"
        starStats['Saturation'] = f"{float(stats['Sequencing Saturation']):.2}"
        starStats['Reads Mapped to Genome'] =\
            f"{float(stats['Reads Mapped to Genome: Unique+Multiple']):.1%}"
        val = stats[f"Reads Mapped to {features}: Unique+Multiple {features}"]
        starStats['Reads Mapped to Transcripts'] = f"{float(val):.1%}"
        val = stats[f"Reads Mapped to {features}: Unique {features}"]
        starStats['Reads Mapped to unique Transcript'] = f"{float(val):.1%}"
        starStats['STAR Cells'] = stats["Estimated Number of Cells"]
        return starStats

    @staticmethod
    def makeTableDict(colNames, values):
        """
        Function to create a simple dictionary from @pair values
        (for use in pd.DataFrame creation)
        Assumption: @colNames and @values are equal length where elements at
        each index correspond to each other

        Args:
            colNames (list): List of column names
            values (list): List of values

        Returns:
            Dictionary where @colNames are keys and @values are values
        """
        resultsDict = {}
        for i, col in enumerate(colNames):
            resultsDict[col] = values[i]
        return resultsDict

    @staticmethod
    def reformatIfInt(val):
        """
        Function to reformat an integer variable

        Args:
            val (obj): Variable to reformat

        Returns:
            Variable after reformatting
        """
        if isinstance(val, int):
            return f"{val:,}"
        elif isinstance(val, float):
            return round(val, 2)
        else:
            return val

    @staticmethod
    def styleTable(styler, title: str, hideColumnHeaders=False,
                   boldColumn=None, numericCols=None):
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
            styler.format(GeneralUtils.reformatIfInt, subset=numericCols)

        styler.hide(axis='index')

        if hideColumnHeaders:
            styler.hide(axis='columns')
        else:
            styler.set_table_styles(
                [{'selector': 'th',
                  'props': [('border', 'black solid !important')]}],
                overwrite=False)

        if boldColumn is not None:
            styler.set_properties(
                subset=boldColumn,
                **{'font-weight': 'bold'}
            )

        if title != "":
            styler.set_caption(title)
            styler.set_table_styles(
                [{'selector': 'caption',
                  'props': [('border', 'black solid !important'),
                            ("font-weight", "bold")]}],
                overwrite=False)

        styler.set_properties(
            **{"border-color": 'black', "border-style": 'solid !important'})
        return styler

    @staticmethod
    def readJSON(file, preserveDictOrder: bool = False):
        """
        Function to read in JSON file and return it as a dictionary

        Args:
            file (str): Path to file
            preserveDictOrder (bool): Flag to indicate whether to
                read the file while preserving the order of entries

        Returns:
            Dictionary with contents of json file
        """
        with open(file) as f:
            str = f.read()
            strStripped = str.rstrip()
            pairs_hook = OrderedDict if preserveDictOrder else None
            parsedJSON = json.loads(strStripped, object_pairs_hook=pairs_hook)
        return parsedJSON

    @staticmethod
    def ensurePathsExist(filePathDict: Dict[str, Union[str, Dict]]):
        """
        Function to ensure all paths mentioned in the given @filePathDict exist

        Args:
            filePathDict (dict): Dictionary containing custom filename as key
                and path to file as value

        Output:
            Raises error if any file not found
        """
        for key, value in filePathDict.items():
            if type(value) is dict:
                GeneralUtils.ensurePathsExist(value)
            else:
                if not exists(value):
                    raise FileNotFoundError(
                        f"{key} was assumed to be located at '{str(value)}'. "
                        "It is missing")

    @staticmethod
    def makeDir(dirName: Path):
        """
        Function to make a directory

        Args:
            dirName (path): Path object pointing to a directory that is to
                be created
        """
        if not exists(dirName):
            mkdir(dirName)

    @staticmethod
    def resolveWriteDir(resultsDir: Union[str, None]) -> Path:
        """
        Function to return path object to directory for writing reports to

        Args:
            resultsDir (str): Path to results directory

        Returns:
            Path object pointing to directory to write reports to
        """
        if (resultsDir is None):
            return Path(".", "reports")
        else:
            return Path(resultsDir, 'reports')


class DatapaneUtils:
    """
    Class that contains helper functions for creating figures in Datapane
    """
    @staticmethod
    def createTableIgnoreWarning(styler, label=None) -> dp.Table:
        '''
        Function that creates wrapper around dp.Table which suppresses
        warning that arises from  code within dp.Table calling Styler.render():
        'this method is deprecated in favour of `Styler.to_html()`'

        Args:
            styler ():
            label ():

        Returns:
            dp.Table object
        '''
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")

            return dp.Table(styler, label=label)

    @staticmethod
    def makeScatterPlot(dataFrame: pd.DataFrame, xVar: str, yVar: str,
                        colorBy: str, title: str, isLogX: bool, isLogY: bool,
                        readThresh: int) -> dp.Plot:
        """
        Function to return dp.Plot rendering of an interactive plotly express
        scatter plot

        Args:
            dataFrame (pd.DataFrame): Data being visualized
            xVar (str): x variable being plotted
            yVar (str): y var being plotted
            colorBy (str): Categorical variable the points are colored by
            isLogX (bool): Flag to indicate whether x axis is on log scale
            isLogY (bool): Flag to indicate whether y axis is on log scale
            title (str): Title of plot
            readThresh (int): Minimum number of unique reads. Cells with reads
                below @readThresh will not be considered in the plot

        Returns:
            dp.Plot object with the scatter plot
        """
        dataFrame = dataFrame[dataFrame['UniqueReads'] >= readThresh]
        cellsToPlot = dataFrame[:dataFrame['pass'].sum()*2]

        if len(cellsToPlot.index) > constants.SAMPLING_NUMBER:
            cellsToPlot = cellsToPlot.sample(constants.SAMPLING_NUMBER)

        scatterPlot = px.scatter(
            cellsToPlot, x=xVar, y=yVar, color=colorBy, title=title,
            log_x=isLogX, log_y=isLogY,
            color_discrete_map=constants.QC_COLORMAP,
            template=constants.DEFAULT_FIGURE_STYLE)

        scatterPlot.update_traces(marker=dict(size=4, opacity=0.5))

        return dp.Plot(scatterPlot)

    @staticmethod
    def buildPlatePlot(df, title, colorbar_title=None):
        """
        Function to build plate like plot for displaying information on a
        per well basis

        Args:
            df (pd.DataFrame): Data to plot
            title (str): Title of the plot

        Returns:
            Matplotlib figure
        """
        fig = plt.figure()
        if colorbar_title:
            ax = sns.heatmap(df, linewidth=0.5, cbar_kws={'label': colorbar_title})
        else:
            ax = sns.heatmap(df, linewidth=0.5)
        ax.set_title(title)

        return fig

    @staticmethod
    def barcodeLevelPlots(sampleName: str, cells: pd.DataFrame, possibleIndexValues: List,
                          index: str, title: str, internalReport: bool, wellAliases=False,
                          writeDir=None) -> dp.Group:
        """
        Function to compute statistics to build plate like plot for
        displaying statistics on a per well basis

        Args:
            cells (pd.DataFrame): Dataframe containing information to plot
            possibleIndexValues (list):
            index (str): Column name to sort dataframe by
            title (str): Title of the plot
            internalReport (bool): Whether report is being generated for internal purposes
            wellAliases (bool): Flag to indicate whether well aliases are used
            writeDir (str): Write directory

        Returns:
            dp.Group object containing all the plots
        """
        allValues = set(possibleIndexValues)
        presentValues = set(cells[index].unique())
        missingIndices = allValues.difference(presentValues)

        if len(missingIndices) > 0:
            cells = cells.copy(deep=True)
            dummyValues = []

            for missingIndex in missingIndices:
                dummyValues.append({index: missingIndex, 'umis': 0})

            cells = pd.concat([cells, pd.DataFrame(dummyValues)])

        cellCountsByIndex = cells.groupby(index).size().reset_index()

        if wellAliases:
            sortOrder = sorted(
                list(cells[index].unique()),
                key=functools.cmp_to_key(CalculationUtils.wellStringComp))
            cells[index] = pd.Categorical(cells[index], sortOrder)
            cellCountsByIndex[index] = pd.Categorical(
                cellCountsByIndex[index], sortOrder)
            cells.sort_values(index, inplace=True)
            cellCountsByIndex.sort_values(index, inplace=True)

        else:
            cells.sort_values(by=[index], inplace=True)
            cellCountsByIndex.sort_values(by=[index], inplace=True)

        readsPerIndexBox = px.box(
            cells, y=index, x='umis', log_x=True, height=700,
            color_discrete_sequence=constants.DEFAULT_COLOR_SEQUENCE,
            color=index, template=constants.DEFAULT_FIGURE_STYLE,
            boxmode="overlay",
            labels={index: title, "umis": "Unique Reads Per Cell"},
            points=False)

        cellCountsByIndex.rename(columns={0: "CellCount"}, inplace=True)
        cellCountsByIndex["Number of Cells"] = \
            cellCountsByIndex.apply(
                lambda row: 0 if row['CellCount'] == 1 else row['CellCount'],
                axis=1)
        well_dict = {}
        # Plate is 24*16 for ligation
        if index == "Ligation_alias":
            wellPlateCellCountDf = \
                pd.DataFrame(
                    0, columns=range(1, 25),
                    index=[chr(i) for i in range(65, 81)])
            wellPlateNumCellDf = \
                pd.DataFrame(
                    0, columns=range(1, 25),
                    index=[chr(i) for i in range(65, 81)])
            for i in range(65, 81):
                for j in range(1, 25):
                    key = str(j)+chr(i)
                    well_dict[key] = []

        # Plate is 12*8 for pcr and rt
        else:
            wellPlateCellCountDf = \
                pd.DataFrame(
                    0, columns=range(1, 13),
                    index=[chr(i) for i in range(65, 73)])
            wellPlateNumCellDf = \
                pd.DataFrame(
                    0, columns=range(1, 13),
                    index=[chr(i) for i in range(65, 73)])
            for i in range(65, 73):
                for j in range(1, 13):
                    key = str(j)+chr(i)
                    well_dict[key] = []
        
        for idx, row in cells.iterrows():
            letter = row[index][-1]
            numbers = row[index][:-1]
            try:
                well_dict[numbers+letter].append(row['umis'])
            except KeyError as e:
                print(f"{e}: {letter}{numbers} does not exist")
            
        for idx, row in cellCountsByIndex.iterrows():
            letter = row[index][-1]
            numbers = row[index][:-1]
            try:
                wellPlateNumCellDf.at[letter, int(numbers)] =\
                    row['Number of Cells']
            except KeyError as e:
                print(f"{e}: {letter}{numbers} does not exist")

        
        for key in well_dict:
            letter = key[-1]
            numbers = key[:-1]
            if len(well_dict[key]) == 0:
                wellPlateCellCountDf.at[letter, int(numbers)] = 0
            else:
                wellPlateCellCountDf.at[letter, int(numbers)] = int(statistics.median(well_dict[key]))

        readsPerIndexBox = DatapaneUtils.buildPlatePlot(
            wellPlateCellCountDf, "Unique Transcript Counts Per Cell")
        cellsPerIndexBar = DatapaneUtils.buildPlatePlot(
            wellPlateNumCellDf, "Number of cells")
        namePrefix = title.replace(" ", "_")
        if writeDir is not None and internalReport:
            cellsPerIndexBar.savefig(
                writeDir / f"CellCount_By_{namePrefix}_Heatmap.png")
            readsPerIndexBox.savefig(
                writeDir / f"UniqueTranscriptCount_By_{namePrefix}_Heatmap.png")
        wellPlateCellCountDf.to_csv(f"reports/{sampleName}_unique_transcript_counts_by_{namePrefix}_well.csv")
        wellPlateNumCellDf.to_csv(f"reports/{sampleName}_num_cells_by_{namePrefix}_well.csv")
        return dp.Group(
            dp.Text(f'## {title}'),
            dp.Group(dp.Plot(cellsPerIndexBar),
                     dp.Plot(readsPerIndexBox), columns=2))

    def createPerWellFigures(self, cells: pd.DataFrame) -> List[dp.Plot]:
        """
        Function that returns per tagmentation barcode (per well) figures
        for Barcode tab of report:
        1. Boxplot of unique reads per cell per well
        2. Barplot of number of cells per well

        Args:
            cells (pd.DataFrame): Data to plot

        Returns:
            List of dp.Plot objects
        """
        cells.sort_values(by='tgmtBc', inplace=True, ascending=False)

        readsPerWellBox = px.box(
            y=cells.tgmtBc, x=cells.UniqueReads, height=700, log_x=True,
            color=cells.tgmtBc, template="none",
            labels={'x': "Unique Reads Per Cell"}, notched=True,
            boxmode="overlay", points=False, title="Unique reads per well")
        readsPerWellBox.update_yaxes(visible=False, showticklabels=False)
        readsPerWellBox.update_layout(showlegend=False)

        tgmtCounts = cells.groupby('tgmtBc').size().reset_index()
        tgmtCounts.rename(columns={0: "Number of Cells"}, inplace=True)
        tgmtCounts.sort_values(by='tgmtBc', inplace=True, ascending=False)

        readsPerWellBar = px.bar(
            tgmtCounts, y="tgmtBc", x="Number of Cells", height=700,
            template="none", color='tgmtBc',
            labels={'tgmtBc': 'Tagmentation Barcodes'}, title="Cells per well")

        tgmtLabels = list(tgmtCounts['tgmtBc'])
        tgmtLabels.sort(reverse=True)

        readsPerWellBox.update_yaxes(categoryorder='array',
                                     categoryarray=tgmtLabels)
        readsPerWellBar.update_yaxes(categoryorder='array',
                                     categoryarray=tgmtLabels)
        readsPerWellBar.update_layout(showlegend=False)

        return [dp.Plot(readsPerWellBar), dp.Plot(readsPerWellBox)]

    @staticmethod
    def makeKneePlot(data: pd.DataFrame, field: str, title: str) -> dp.Plot:
        """
        Function to make a kneeplot using @field in @data; drawing a
        vertical line at @threshold

        Args:
            data (pd.DataFrame): Data to plot
            field (str): Column to sort dataframe by
            title (str): Title of plot

        Returns:
            dp.Plot object that contains the figure
        """
        indices = CalculationUtils.getIndicesToInclude(len(data.index))
        vals = pd.DataFrame()

        cellbarcodeCounts = list(range(len(data.index)))
        vals['Cell Barcodes'] = [cellbarcodeCounts[i] for i in indices]

        uniqueReads = list(data.sort_values(by=field, ascending=False)[field])
        vals['Unique Reads'] = [uniqueReads[i] for i in indices]

        fig = px.line(vals, x='Cell Barcodes', y='Unique Reads',
                      title=title, log_x=True, log_y=True,
                      template=constants.DEFAULT_FIGURE_STYLE)
        fig.add_vline(x=data['pass'].sum(), line_dash="dash",
                      line_color="green")

        # Xlimit set in powers of 10 (10, maximum)
        fig.update_layout(xaxis_range=[1, log10(vals['Cell Barcodes'].max())])

        return dp.Plot(fig)

    def showAllTicks(plotlyFig):
        """
        Function to modify plotly figure to show all ticks

        Args:
            plotlyFig (obj): Plotly figure
        """
        plotlyFig.update_layout(
            xaxis=dict(tickmode='linear', tick0=1, dtick=1)
        )

    @staticmethod
    def createAliasMap(ids):
        """
        Function to create an alias map from a list

        Args:
            ids (list): List to create alias map for

        Returns:
            Dictionary where key is an entry in @ids and value is
            corresponding alias
        """
        result = {}

        for i, id in enumerate(sorted(ids)):
            result[id] = str(i)

        return result


class CalculationUtils:
    """
    Class that houses util functions for different calculations
    """
    def __init__(self):
        self = self

    def wellStringComp(x: str, y: str):
        """
        Function to compare well coordinates

        Args:
            x (str): First well coordinate
            y (str): Second well coordinate

        Returns:
            Result of comparison
        """
        xNum = int(x[0:-1])
        xLetter = x[-1]
        yNum = int(y[0:-1])
        yLetter = y[-1]
        num = compare(xNum, yNum)
        letter = compare(xLetter, yLetter)
        return num if num != 0 else letter

    @staticmethod
    def getIndicesToInclude(elementCount: int,
                            result: List[int] = []) -> List[int]:
        '''
        Function that returns a list of indices in ascending order for
        elements to be included. Used to reduce size of line
        plots while retaining overall trend shape

        Args:
            elementCount (int): Total number of elements in a list

        Returns:
            List of indices in ascending order for elements to be included.
            Used to reduce size of line plots while retaining overall trend
            shape
        '''
        logVal = int(round(log10(elementCount)))
        nextDown = int(logVal - 1)
        lowerLimit = int(10**nextDown)
        step = ((logVal)**3) if logVal > 0 else 1
        lowerLimit = lowerLimit - 1 if (lowerLimit > 0) else lowerLimit
        result.append(list(range(lowerLimit, elementCount, step)))

        if (lowerLimit > 0):
            return CalculationUtils.getIndicesToInclude(lowerLimit, result)
        else:
            return sorted(functools.reduce(operator.iconcat, result, []))

    @staticmethod
    def coverage_equation(x, c, n):
        """
        Function to

        Args:
            x ():
            c ():
            n ():

        Returns:
        """
        return c / x - 1 + np.exp(-n / x)

    @staticmethod
    def estLibSize(reads, umis):
        """
        Function to estimate library size(number of unique molecules)
        at infinite sequencing depth

        Args:
            reads (int): Total reads
            umis (int): Observed UMIs

        Returns:
            Float containing estimated library size
        """
        if umis > reads:
            raise ValueError()

        if umis <= 0:
            raise ValueError()

        m = 1
        M = 100

        while CalculationUtils.coverage_equation(M * umis, umis, reads) > 0:
            M *= 10

        for _ in range(40):
            r = (m + M) / 2
            u = CalculationUtils.coverage_equation(r * umis, umis, reads)

            if (u > 0):
                m = r

            elif (u < 0):
                M = r

            else:
                break

        return (umis * (m + M) / 2)

    @staticmethod
    def downsampleTargetReadCount(meanReads, medianPassingReads, targetReads):
        """
        Function to downsample targetRead count based on relation of
        meanReads to medianPassingReads

        Args:
            meanReads (): Mean of reads of passing cells
            medReads (): Median of passing reads of passing cells
            targetReads (int): Number of target reads

        Returns:
            Downsampled targetRead count
        """
        return (targetReads/meanReads) * medianPassingReads

    @staticmethod
    def extrapolate_umis(reads, umis, seqDepthInReads):
        """
        Function to estimate library size at infinite sequencing depth
        and extrapolate the expected number of unique molecules at
        @seqDepthInReads

        Args:
            reads (): Median of passing reads of passing cells
            umis (): Median of total UMIs of passing cells

        Returns:
            Estimated library size
        """
        estimatedLibrarySize = CalculationUtils.estLibSize(reads, umis)
        # Total UMIs
        res = (
            estimatedLibrarySize * (
                1 - (((
                    estimatedLibrarySize - 1) /
                    (estimatedLibrarySize))**seqDepthInReads)))
        return int(res) if not np.isnan(res) else res
