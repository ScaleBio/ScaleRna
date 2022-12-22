from logging import raiseExceptions
import numpy as np 
import json
from typing import Dict, Union, Tuple, List
import pandas as pd
from os import mkdir
from os.path import exists 
from pathlib import Path
import warnings
import datapane as dp
from math import log10 
import functools 
import operator
import reportGeneration.myconstants as constants
from collections import OrderedDict
import plotly.express as px
from dataclasses import dataclass


def compare(this, that):
    if this > that:
        return 1 
    elif this < that:
        return -1
    else:
        return 0 



@dataclass
class CellStat:
    totalGenes: int = 0
    humanGenes: int = 0
    mouseGenes: int = 0
    totalUmis: int = 0
    humanUmis: int = 0
    mouseUmis: int = 0

## GENERAL PURPOSE UTILITY FUNCTIONS RELATED TO IO, INPUT PARSING, OR PANDAS MANIPULATIONS  
class GeneralUtils:
    def __init__(self):
        self=self

    @staticmethod
    def loadStatsAsDict(fn):
        '''
        Create a dictionary from a file (@fn) containing tab delimited statistics
        '''
        res = {}
        for l in open(fn):
            line = l.strip().split("\t")
            if not len(line) >= 2: continue
            stat = line[0]
            val = line[1]
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
        stats = {}
        for l in open(fn):
            line = l.strip().split(',')
            stats[line[0]] = line[1]

        starStats = {}
        starStats['Total Reads'] = f"{int(stats['Number of Reads']):,}"
        starStats['Saturation'] = f"{float(stats['Sequencing Saturation']):.2}"
        starStats['Reads Mapped to Genome'] = f"{float(stats['Reads Mapped to Genome: Unique+Multiple']):.1%}"
        val = stats[f"Reads Mapped to {features}: Unique+Multiple {features}"]
        starStats['Reads Mapped to Transcripts'] = f"{float(val):.1%}"
        val = stats[f"Reads Mapped to {features}: Unique {features}"]
        starStats['Reads Mapped to unique Transcript'] = f"{float(val):.1%}"
        starStats['STAR Cells'] = stats["Estimated Number of Cells"]
        return starStats

    @staticmethod
    def makeTableDict(colNames, values):
        '''
        Create a simple dictionary from @pair values (for use in pd.DataFrame creation)
        Assumption: 
        - @colNames and @values are equal length where the elements at each index correspond to each other
        '''
        resultsDict = {}
        for i,col in enumerate(colNames):
            resultsDict[col] = values[i]
        return resultsDict

    
    @staticmethod 
    def reformatIfInt(val):
        if isinstance(val, int):
            return f"{val:,}"
        elif isinstance(val, float):
            return round(val, 2)
        else:
            return val

    @staticmethod
    def styleTable(styler, title:str, hideColumnHeaders=False, boldColumn=None, numericCols=None):
        '''
        Modifies given @pd.DataFrame.Styler 
        '''
        if (numericCols is not None):
            styler.format(GeneralUtils.reformatIfInt, subset=numericCols)
        styler.hide(axis='index')
        if hideColumnHeaders:
            styler.hide(axis='columns')
        else:
            styler.set_table_styles([{'selector':'th', 'props':[('border', 'black solid !important')]}], overwrite=False)
        if boldColumn is not None:
            styler.set_properties(
                subset=boldColumn,
                **{'font-weight': 'bold'}
            )
        if title != "":
            styler.set_caption(title)
            styler.set_table_styles([{'selector':'caption', 'props':[('border', 'black solid !important'),  ("font-weight", "bold")]}],overwrite=False)
        styler.set_properties(**{"border-color":'black', "border-style":'solid !important'})
        return styler

    @staticmethod
    def readJSON(file, preserveDictOrder:bool=False):
        """
        Reads in JSON as a python object
        """
        with open(file) as f:
            str = f.read()
            strStripped = str.rstrip()
            pairs_hook = OrderedDict if preserveDictOrder else None
            parsedJSON = json.loads(strStripped, object_pairs_hook=pairs_hook)
        return parsedJSON


    @staticmethod
    def ensurePathsExist(filePathDict:Dict[str, Union[str,Dict]]):
        '''
        Ensures all paths mentioned in the given @filePathDict exist
        '''
        for key, value in filePathDict.items():
            if type(value) is dict:
                GeneralUtils.ensurePathsExist(value)
            else: 
                if not exists(value):
                    raise FileNotFoundError(f"{key} was assumed to be located at '{str(value)}'. It is missing")

    @staticmethod
    def makeDir(dirName:Path):
        if not exists(dirName): 
            mkdir(dirName)

    @staticmethod
    def resolveWriteDir(resultsDir: Union[str,None]) -> Path: 
        if (resultsDir is None):
            return Path(".", "reports")
        else:
            return Path(resultsDir, 'reports')


### FIGURE CREATION METHODS ### 
class DatapaneUtils:
    @staticmethod
    def createTableIgnoreWarning(styler, label=None) -> dp.Table:
        '''
        Wrapper around dp.Table which suppresses warning which arrises from 
        code within dp.Table calling Styler.render():
        'this method is deprecated in favour of `Styler.to_html()`' 
        '''
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            return dp.Table(styler, label=label)


    @staticmethod
    def makeScatterPlot(dataFrame: pd.DataFrame, xVar:str, yVar:str, colorBy:str, title:str, isLogX:bool, isLogY:bool, colorMap:Dict[str,str], readThresh:int) -> dp.Plot:
        '''
        returns dp.Plot rendering of an interactive plotly express scatter plot where: 
        - @dataFrame: is data being visualized
        - @xVar: x variable being plotted 
        - @yVar: y var being plotted 
        - @colorBy: categorical variable points are colored by 
        - @logX, logY: specify whether those axes are to be log scaled
        '''
        dataFrame = dataFrame[dataFrame['UniqueReads'] >= readThresh]
        cellsToPlot = dataFrame[:dataFrame['pass'].sum()*2]
        if len(cellsToPlot.index) > constants.SAMPLING_NUMBER:
            cellsToPlot = cellsToPlot.sample(constants.SAMPLING_NUMBER)
        scatterPlot = px.scatter(cellsToPlot, x=xVar, y=yVar, color=colorBy, title=title, log_x=isLogX, log_y=isLogY,color_discrete_map=constants.QC_COLORMAP, template=constants.DEFAULT_FIGURE_STYLE)
        scatterPlot.update_traces(marker=dict(size=4, opacity=0.5))
        return dp.Plot(scatterPlot)


    @staticmethod
    def barcodeLevelPlots(cells: pd.DataFrame, possibleIndexValues:List, index:str, title:str, showLegend:bool=True, wellAliases=False, writeDir=None) -> dp.Group:
        allValues = set(possibleIndexValues)
        presentValues = set(cells[index].unique())
        missingIndices = allValues.difference(presentValues)
        if len(missingIndices) > 0:
            cells = cells.copy(deep=True)
            dummyValues = []
            for missingIndex in missingIndices:
                dummyValues.append({index:missingIndex, 'umis':0})
            cells = pd.concat([cells, pd.DataFrame(dummyValues)])
        cellCountsByIndex = cells.groupby(index).size().reset_index()
        if wellAliases:
            sortOrder = sorted(list(cells[index].unique()), key=functools.cmp_to_key(CalculationUtils.wellStringComp))
            cells[index] = pd.Categorical(cells[index], sortOrder)
            cellCountsByIndex[index] = pd.Categorical(cellCountsByIndex[index], sortOrder)
            cells.sort_values(index, inplace=True)
            cellCountsByIndex.sort_values(index,inplace=True)
        else:
            cells.sort_values(by=[index],inplace=True)
            cellCountsByIndex.sort_values(by=[index], inplace=True)
        readsPerIndexBox = px.box(cells, y=index, x='umis', log_x=True, height=700, color_discrete_sequence= constants.DEFAULT_COLOR_SEQUENCE, color=index, template=constants.DEFAULT_FIGURE_STYLE, boxmode="overlay", labels={index:title, "umis":"Unique Reads Per Cell"}, points=False)
        cellCountsByIndex.rename(columns={0: "CellCount"},inplace=True)
        cellCountsByIndex["Number of Cells"] = cellCountsByIndex.apply(lambda row: 0 if row['CellCount'] == 1 else row['CellCount'], axis=1)
        cellsPerIndexBar = px.bar(cellCountsByIndex, y=index, x="Number of Cells", color_discrete_sequence= constants.DEFAULT_COLOR_SEQUENCE, height=700,template=constants.DEFAULT_FIGURE_STYLE, color=index, labels={index:title})
        if writeDir is not None:
            namePrefix = title.replace(" ", "_")
            cellsPerIndexBar.write_image(writeDir / f"CellCount_By_{namePrefix}_Barplot.png")
            readsPerIndexBox.write_image(writeDir / f"ReadCount_By_{namePrefix}_BoxPlot.png")
        # Modifying figures after they are written as png so each can be interpreted independently
        readsPerIndexBox.update_layout(yaxis={'showticklabels':False}, yaxis_title=None)
        if not showLegend:
            readsPerIndexBox.update_layout(showlegend=False)
        cellsPerIndexBar.update_layout(showlegend=False)
        return dp.Group(dp.Text(f'## {title}'),dp.Group(dp.Plot(cellsPerIndexBar), dp.Plot(readsPerIndexBox), columns=2))

    def createPerWellFigures(self, cells: pd.DataFrame) -> List[dp.Plot]:
        '''
        Returns per tagmentation barcode (per well) figures for Barcode tab of report: 
        1. Boxplot of unique reads per cell per well 
        2. Barplot of numver of cells per well 
        '''
        cells.sort_values(by='tgmtBc', inplace=True, ascending=False)
        readsPerWellBox = px.box(y=cells.tgmtBc, x=cells.UniqueReads,height=700, log_x=True, color=cells.tgmtBc, template="none",  labels={'x':"Unique Reads Per Cell"}, notched=True, boxmode="overlay", points=False, title="Unique reads per well")
        readsPerWellBox.update_yaxes(visible=False, showticklabels=False)
        readsPerWellBox.update_layout(showlegend=False)
        tgmtCounts = cells.groupby('tgmtBc').size().reset_index()
        tgmtCounts.rename(columns={0: "Number of Cells"},inplace=True)
        tgmtCounts.sort_values(by='tgmtBc', inplace=True, ascending=False)
        readsPerWellBar = px.bar(tgmtCounts, y="tgmtBc", x="Number of Cells", height=700,template="none", color='tgmtBc', labels={'tgmtBc': 'Tagmentation Barcodes'}, title="Cells per well")
        tgmtLabels = list(tgmtCounts['tgmtBc'])
        tgmtLabels.sort(reverse=True)
        readsPerWellBox.update_yaxes(categoryorder='array', categoryarray=tgmtLabels)
        readsPerWellBar.update_yaxes(categoryorder='array', categoryarray=tgmtLabels)
        readsPerWellBar.update_layout(showlegend=False)
        return [dp.Plot(readsPerWellBar), dp.Plot(readsPerWellBox)]

    @staticmethod
    def makeKneePlot(data: pd.DataFrame, field:str, title:str) -> dp.Plot:
        '''
        Makes a kneeplot using @field in @data; drawing a vertical line at @threshold
        '''    
        indices = CalculationUtils.getIndicesToInclude(len(data.index))
        vals = pd.DataFrame()
        cellbarcodeCounts = list(range(len(data.index)))
        vals['Cell Barcodes'] = [cellbarcodeCounts[i] for i in indices]
        uniqueReads = list(data.sort_values(by=field, ascending=False)[field])
        vals['Unique Reads'] = [uniqueReads[i] for i in indices]
        fig = px.line(vals, x='Cell Barcodes', y='Unique Reads',  title=title, log_x=True, log_y=True, template=constants.DEFAULT_FIGURE_STYLE)
        fig.add_vline(x=data['pass'].sum(), line_dash="dash", line_color="green")
        # Xlimit set in powers of 10 (10, maximum)
        fig.update_layout(xaxis_range=[1, log10(vals['Cell Barcodes'].max())])
        return dp.Plot(fig)

    def showAllTicks(plotlyFig):
        plotlyFig.update_layout(
            xaxis = dict(
                tickmode='linear',
                tick0 = 1,
                dtick=1
            )
        )

    @staticmethod 
    def createAliasMap(ids):
        result = {}
        for i, id in enumerate(sorted(ids)):
            result[id] = str(i)
        return result


class CalculationUtils:
    def __init__(self):
        self=self

    def wellStringComp(x:str, y:str):
        xNum = int(x[0:-1])
        xLetter = x[-1]
        yNum = int(y[0:-1])
        yLetter = y[-1]
        num = compare(xNum, yNum)
        letter = compare(xLetter, yLetter)
        return num if num != 0 else letter

    @staticmethod
    def getIndicesToInclude(elementCount:int, result: List[int]=[])-> List[int]:
        '''
        Where @elementCount is the total number of elements in a 
        list. Returns a list of indices in ascending order for 
        elements to be included. Used to reduce size of line 
        plots while retaining overall trend shape    
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
        return c / x - 1 + np.exp(-n / x)

    @staticmethod
    def estLibSize(reads, umis):
        '''
        given @reads total reads and @umis observed UMIS, 
        estimates library size (number of unique molecules) 
        at infinite sequencing depth
        '''
        if umis > reads: raise ValueError()
        if umis <= 0: raise ValueError()
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
        '''
        Downsample targetRead count based on relation of meanReads to medianPassingReads 
        '''
        return (targetReads/meanReads) * medianPassingReads


    @staticmethod
    def extrapolate_umis_at_nreads(reads, umis, seqDepthInReads):
        '''
        1. estimates library size at infinite sequencing depth
        2. extrapolates the expected number of unique molecules 
        at @seqDepthInReads
        '''
        estimatedLibrarySize = CalculationUtils.estLibSize(reads, umis)
        # Total UMIs 
        res = (estimatedLibrarySize * (1 - (((estimatedLibrarySize - 1)/(estimatedLibrarySize))**seqDepthInReads)))
        return int(res) if not np.isnan(res) else res

    # @staticmethod
    # def extrapolate_umis_at_50sat(reads, umis):
    #     '''
    #     1. estimates library size at infinite sequencing depth
    #     2. extrapolates the expected number of unique molecules 
    #     at @seqDepthInReads
    #     '''
    #     extrapolatedLibSize = CalculationUtils.estLibSize(reads, umis)
    #     for seqDepthInReads in range(1000,1000000, 100):
    #         totalUMIs = int(extrapolatedLibSize * (1 - ((extrapolatedLibSize - 1)/(extrapolatedLibSize))**seqDepthInReads))
    #         sat = 1 - (totalUMIs/seqDepthInReads)
    #         if abs(sat - 0.50) < 0.01:
    #             return totalUMIs
    #     # Saturation of 50% couldn't be reached w/ read depth of 1M reads 
    #     return None

    # @staticmethod 
    # def getStatsAsDataframes(filterDfPath: str, cellStatsPath:Union[str,None], threshJson: Dict) -> Tuple[pd.DataFrame,pd.DataFrame]:
    #     '''
    #     Parses in @filterDfPath (outputted by QC filtering step) & @cellStatsPath as dataframes, adding additional columns 
    #     based on existing columns
    #     '''
    #     df = pd.read_csv(filterDfPath, sep='\t', index_col="Barcode")
    #     cellStats = None
    #     if (cellStatsPath is not None):
    #         cellStats = pd.read_csv(cellStatsPath, sep="\t", index_col=0)
    #         cellStats['pass'] = cellStats['UniqueReads'] >= threshJson['UniqueReads']
    #         cellStats['Saturation'] = np.around(1-cellStats['UniqueReads'] /cellStats['PassingReads'],decimals=3)
    #         cellStats['FractionMito'] = np.around(cellStats['MitoReads'] / cellStats['TotalReads'],decimals=3)
    #     return (df, cellStats)

