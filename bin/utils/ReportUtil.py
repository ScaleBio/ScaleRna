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
import functools
import operator
import utils.myconstants as constants
from collections import OrderedDict
import plotly.express as px
from dataclasses import dataclass
from matplotlib import pyplot as plt
import matplotlib.colors as colors


def compare(this, that):
    """
    Compare two variables

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
        Reformat an integer variable

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
            styler.format(GeneralUtils.reformatIfInt, subset=numericCols)

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
                    raise FileNotFoundError(f"{key} was assumed to be located at '{str(value)}'. It is missing")


class DatapaneUtils:
    """
    Class that contains helper functions for creating figures in Datapane
    """
    @staticmethod
    def createTableIgnoreWarning(styler, label=None) -> dp.Table:
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

    @staticmethod
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
        if colorbar_title:
            ax = sns.heatmap(df, linewidth=0.5, cmap=sns.color_palette("dark:#80AAFF", as_cmap=True), norm=colors.SymLogNorm(linthresh=threshold, vmin=0), cbar_kws={'label': colorbar_title})
        else:
            ax = sns.heatmap(df, linewidth=0.5, cmap=sns.color_palette("dark:#80AAFF", as_cmap=True), norm=colors.SymLogNorm(linthresh=threshold, vmin=0))
        ax.set_title(title)

        return fig
    
    def getCharacterIndices(lower_limit, upper_limit):
        """
        Generate a list of letter indices from a range of ascii values

        Args:
            lower_limit (int): Lower limit to generate ascii's from
            upper_limit (int): Upper limit to generate ascii's till
        """
        return [chr(i) for i in range(lower_limit, upper_limit)]

    def buildDfForSamplePlatePlot(libJson, index, referencesPath):
        well_dict = {}
        for entry in libJson["barcodes"]:
            if "alias" in entry:
                if entry["alias"] == index.split("_")[0]:
                    lib_json_entry_dict = entry
        max_letter, max_number = DatapaneUtils.getMaxWellNumberAndLetter(referencesPath / f'{lib_json_entry_dict["sequences"]}')
        wellPlateCellCountDf = pd.DataFrame(0, columns=range(1, max_number+1), index=DatapaneUtils.getCharacterIndices(65,ord(max_letter)+1))
        wellPlateNumCellDf = wellPlateCellCountDf.copy()
        for i in range(65, ord(max_letter)+1):
            for j in range(1, max_number+1):
                key = str(j)+chr(i)
                well_dict[key] = []
        return wellPlateCellCountDf, wellPlateNumCellDf, well_dict


    @staticmethod
    def barcodeLevelPlots(referencesPath: Path, sampleName: str, cells: pd.DataFrame, possibleIndexValues: List,
                          index: str, title: str, internalReport: bool, libJson: dict, wellAliases=False,
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
        wellPlateCellCountDf, wellPlateNumCellDf, well_dict = DatapaneUtils.buildDfForSamplePlatePlot(libJson, index, referencesPath)
        
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

        readsPerIndexBox = DatapaneUtils.buildPlatePlot(wellPlateCellCountDf, "Unique Transcript Counts Per Cell", 100.0)
        cellsPerIndexBar = DatapaneUtils.buildPlatePlot(wellPlateNumCellDf, "Number of cells", 1.0)
        namePrefix = title.replace(" ", "_")
        if writeDir is not None and internalReport:
            cellsPerIndexBar.savefig(writeDir / f"{sampleName}_figures" / f"CellCount_By_{namePrefix}_Heatmap.png")
            readsPerIndexBox.savefig(writeDir / f"{sampleName}_figures" / f"UniqueTranscriptCount_By_{namePrefix}_Heatmap.png")
        wellPlateCellCountDf.to_csv(writeDir / "csv" / f"{sampleName}_unique_transcript_counts_by_{namePrefix}_well.csv")
        wellPlateNumCellDf.to_csv(writeDir / "csv" / f"{sampleName}_num_cells_by_{namePrefix}_well.csv")
        return dp.Group(dp.Text(f'## {title}'), dp.Group(dp.Plot(cellsPerIndexBar), dp.Plot(readsPerIndexBox), columns=2))


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
    def getIndicesToInclude(elementCount: int, result: List[int] = []) -> List[int]:
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
        logVal = int(round(np.log10(elementCount+1)))
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
        res = (estimatedLibrarySize * (1 - (((estimatedLibrarySize - 1) / (estimatedLibrarySize))**seqDepthInReads)))
        return int(res) if not np.isnan(res) else res
