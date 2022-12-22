#!/usr/bin/env python
import numpy as np
import pandas as pd 
import datapane as dp
import reportGeneration.myconstants as constants
import plotly.express as px
import collections 
from reportGeneration.ReportUtil import CellStat, GeneralUtils, CalculationUtils, DatapaneUtils
from reportGeneration.AbstractReport import AbstractReport
from pathlib import Path
from typing import Dict, Union
import math
import itertools as it

def intOrNan(x):
    if np.isnan(x): return x
    return int(x)

class SampleReport(AbstractReport):
    sampleName: str
    sampleSpecificFilePaths: Dict[str,Path]
    sampleThresholds: Dict[str, float]
    featuresType: str
    genes: pd.DataFrame
    barcodes: pd.DataFrame
    allCells: pd.DataFrame

    def __init__(self, sampleName: str, resultsDir: Path, featureType: str, cells: Union[int,None], makeAdvanced:bool, demuxJson: Dict[str,float], libDir:Path, forFastqReport:bool=False):
        super().__init__(resultsDir, makeAdvanced, demuxJson, libDir)
        self.sampleName = sampleName
        self.featuresType = featureType
        self.cells = cells
        self.sampleSpecificFilePaths = self.resolveSampleSpecificFilePaths(forFastqReport)
        self.validateInput()

    def validateInput(self):
        GeneralUtils.makeDir(self.writeDir)
        GeneralUtils.makeDir(self.writeDir / f"{self.sampleName}_figures")
        GeneralUtils.ensurePathsExist(self.sampleSpecificFilePaths)

    def saveFigureAsPng(self, fig, pngName:str):
        fig.write_image(self.writeDir / f"{self.sampleName}_figures" / pngName)

    def splitBarcodes(self):
        '''
        Uses `level` and `alias` values in each element of lib.json barcodes object to split barcodes 
        ## ASSUMPTION: Each element in `barcodes` section of lib.json defining part of the CB sequence has a `level` and `alias` field 
        ## see current lib.jsons in repo references directory for examples 
        '''
        sampleBarcodeAliases = {}
        for name,bc in self.demuxJson["samples"][self.sampleName]["barcodes"].items():
            sampleBarcodeAliases[bc['sequence']] =  name 
        libDirJSON = GeneralUtils.readJSON(self.libDir, preserveDictOrder=True)
        sampleBarcode = libDirJSON['sample_barcode']
        barcodeInfo = libDirJSON['barcodes']
        scDefiningBarcodes = [barcodeMeta for barcodeMeta in barcodeInfo if barcodeMeta.get('type', None) != 'library_index' and barcodeMeta.get('type', None) != 'umi']
        groupByLevel = it.groupby(scDefiningBarcodes, key=lambda x: x['level'])
        curIndex = 0
        barcodesToPlot = {}
        for level, group in groupByLevel:
                groupAsList = list(group)
                barcodeSize = sum([x['length'] for x in groupAsList])
                barcodeFullName = groupAsList[0]['alias']
                barcodeName = groupAsList[0]['name']
                barcodes = [b[curIndex:barcodeSize+curIndex] for b in self.allCells.index]
                self.allCells[barcodeFullName] = barcodes
                curIndex += barcodeSize
                if barcodeName == sampleBarcode and len(sampleBarcodeAliases) > 0:
                    self.allCells[f"{barcodeFullName}_alias"] = self.allCells[barcodeFullName].apply(lambda x: sampleBarcodeAliases[x])
                    barcodesToPlot[level] = {'alias': barcodeFullName, 'possibleValues':sampleBarcodeAliases.values(), 'orderAsWellAliases':True} 
                else: 
                    indexAliases = DatapaneUtils.createAliasMap(list(self.allCells[barcodeFullName].unique()))
                    self.allCells[f'{barcodeFullName}_alias'] = self.allCells.apply (lambda row: indexAliases[row[barcodeFullName]],axis=1)
                    barcodesToPlot[level] = {'alias':barcodeFullName, 'possibleValues':indexAliases.values(), 'orderAsWellAliases':False} 
        self.barcodesToPlot = barcodesToPlot

    def getCellThreshold(self):
        """
        Calculate and set the UMI threshold for cells
        Either based on explicit @cells parameter to STAR call
        """
        if not self.cells:
            self.cells = int(self.starStats['STAR Cells'])
        umis = self.allCells.totalUmis
        k = len(umis) - self.cells - 1
        self.cellThreshold = max(20, np.partition(umis, k)[k])

    def build(self):
        pages = []
        self.instantiateGenesDf()
        self.instantiateBarcodesDf()
        self.instantiateCellCountsDf()

        self.starStats = self.getStarStats()
        self.getCellThreshold()
        self.addStatsAndCallCells(self.starStats)

        (readsPage, statsDf) = self.buildReadsPage()
        pages.append(readsPage)
        if (self.libDir is not None):
            self.splitBarcodes()
            barcodesPage = self.buildBarcodesPage()
            pages.append(barcodesPage)
        isBarnyard = self.checkIfBarnyard()
        if isBarnyard:
            (barnyardPage, barnyardStat) = self.buildBarnyardPage()
            statsDf = pd.concat([statsDf, barnyardStat])
            pages.append(barnyardPage)
        report = dp.Report(
            blocks=pages
        )
        statsDf.to_csv(self.writeDir/ f"{self.sampleName}.reportStatistics.tsv", index=False, sep='\t',header=False, columns=["Category", "Metric", "Value"])
        report.save(self.writeDir / f"{self.sampleName}.report.html")


    def buildBarnyardPage(self) -> dp.Page:
        barnyardStatsDf = self.makeBarnyardStatsDf()
        barnyardStatsTable = DatapaneUtils.createTableIgnoreWarning(barnyardStatsDf[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, hideColumnHeaders=True, boldColumn=['Metric'], title="", numericCols=['Value']))
        barnyardPlot = self.makeBarnyardScatterPlot()
        plotGroup = dp.Group(barnyardStatsTable, barnyardPlot, columns=2)
        return (dp.Page(plotGroup, title="Barnyard"), barnyardStatsDf)

    def makeBarnyardScatterPlot(self):
        cellsToPlot = self.allCells[:self.allCells['pass'].sum()*2]
        if len(cellsToPlot.index) > 5000:
            cellsToPlot = cellsToPlot.sample(5000)
        fig = px.scatter(cellsToPlot, x="humanUmis", y="mouseUmis", color="species", log_x=True, log_y=True, color_discrete_map=constants.BARNYARD_COLORMAP, labels={"humanUmis":"Human UMIs", "mouseUmis":"Mouse UMIs"}, template=constants.DEFAULT_FIGURE_STYLE)        
        fig.update_traces(selector=dict(name='None'), visible="legendonly")
        self.saveFigureAsPng(fig, "Barnyard_ScatterPlot.png")
        return dp.Plot(fig)

    def buildBarcodesPage(self):
        blocksToRender=[]
        blocksToRender.append(self.createBarcodeLevelFigures())
        return dp.Page(blocks=blocksToRender, title="Barcodes")

    def createBarcodeLevelFigures(self):
        passingCells = self.allCells[self.allCells['pass']]
        allIndexPlots = []
        for key in sorted(self.barcodesToPlot):
            data = self.barcodesToPlot[key]
            bcColName = data['alias']
            possibleValues = list(data['possibleValues'])
            indexPlots = DatapaneUtils.barcodeLevelPlots(passingCells, possibleValues, f"{bcColName}_alias", f"{bcColName} Index", False, wellAliases=data['orderAsWellAliases'], writeDir=self.writeDir / f"{self.sampleName}_figures")
            allIndexPlots.append(indexPlots)
        return dp.Group(blocks=allIndexPlots)

    def getStarStats(self):
        return GeneralUtils.loadStarStatsFromCsv(self.sampleSpecificFilePaths['summaryFn'], self.featuresType)

    def buildReadsPage(self):
        statsDf =  pd.DataFrame({'Metric': ['SampleName'], 'Value': [self.sampleName]})
        statsDf['Category'] = 'Sample'
        statsToAdd = {"Genic Reads":'propGeneMatched', "Exonic Reads":'exon', "Antisense Reads": "antisense", "Mitochondrial Reads":"mito"}
        summaryStatsDf = self.makeSummaryStatsDisplayDf(statsToAdd)
        summaryStatsTable = DatapaneUtils.createTableIgnoreWarning(summaryStatsDf[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, title="Cell Calling", hideColumnHeaders=True, boldColumn=['Metric'], numericCols=['Value']))
        starStatsDf = SampleReport.makeStarDfFromDict(self.starStats)
        starStatsTable = DatapaneUtils.createTableIgnoreWarning(starStatsDf[constants.DISPLAY_COLUMNS].style.pipe(GeneralUtils.styleTable, title="STARSolo Metrics", hideColumnHeaders=True, boldColumn=['Metric']))
        statsDf = pd.concat([statsDf, summaryStatsDf,starStatsDf])
        print(statsDf)
        rankPlot = self.makeRankPlot()
        genesVersusUMIsScatter = self.buildGeneXUMIScatter()
        saturationScatter = self.buildSaturationScatter()
        geneReadsDist = self.buildDist("Distribution of Reads: Gene Matched", 'propGeneMatched', 'Proportion of Genic Reads')
        antisenseReadsDist = self.buildDist("Distribution of Reads: Antisense", "antisense", "Proportion of Antisense Reads")
        exonicReadsDist = self.buildDist("Distribution of Reads: Exonic", 'exon', "Proportion of Exonic Reads")
        groupBlocks = [rankPlot, summaryStatsTable, genesVersusUMIsScatter, saturationScatter, geneReadsDist,  exonicReadsDist, antisenseReadsDist]
        anyNonZeroMitoReads = (self.allCells[self.allCells['pass']].mito > 0).any()
        if anyNonZeroMitoReads:
            mitoReadsDist = self.buildDist("Distribution of Reads: Mitochondrial", 'mito', "Proportion of Mitochondrial Reads")
            groupBlocks.append(mitoReadsDist)
        page = dp.Page(starStatsTable,dp.Group(blocks=groupBlocks, columns=2), title="Cells")
        return (page, statsDf)

        
    def buildGeneXUMIScatter(self):
        '''
        Makes a scatter plot of  
        '''
        cellsToPlot = self.allCells[:self.allCells['pass'].sum()*2]
        if len(cellsToPlot.index) > 5000:
            cellsToPlot = cellsToPlot.sample(5000)
        fig = px.scatter(cellsToPlot, x="umis", y="genes", color='pass', color_discrete_map=constants.QC_COLORMAP, template=constants.DEFAULT_FIGURE_STYLE, labels={"umis": "UMIs detected", "genes":"Genes detected"}, title="Genes vs. UMIs detected")
        return dp.Plot(fig)

    def buildDist(self, title, field, replacementStr):
        passingCells = self.allCells[self.allCells['pass']]
        fig = px.histogram(passingCells, x=field, title=title, histnorm='percent',template=constants.DEFAULT_FIGURE_STYLE, labels={field:replacementStr})
        fig.update_layout(yaxis_title="% Cells")
        return dp.Plot(fig)

    def buildSaturationScatter(self):
        cellsToPlot = self.allCells[:self.allCells['pass'].sum()*2]
        if len(cellsToPlot.index) > 5000:
            cellsToPlot = cellsToPlot.sample(5000)
        fig = px.scatter(cellsToPlot, x="reads", y="Saturation", labels={"reads":"Total Reads"}, template=constants.DEFAULT_FIGURE_STYLE, color='pass', color_discrete_map=constants.QC_COLORMAP, title="Saturation Per Cell")
        return dp.Plot(fig)
    
    def makeSummaryStatsDisplayDf(self, statsToDeriveMediansFor:Dict[str,str]):
        passingCells = self.allCells[self.allCells['pass']]
        stats = []
        stats.append(['Cells above Threshold', len(passingCells)])
        stats.append(['Unique Reads Threshold', passingCells.totalUmis.min()])
        stats.append(["Mean Reads per cell", intOrNan(passingCells.reads.mean())])
        stats.append(["Median UMIs per cell", intOrNan(passingCells.totalUmis.median())])
        stats.append(["Median Genes per cell", intOrNan(passingCells.totalGenes.median())])
        stats.append(["Reads in Cells",f"{passingCells.passingReads.sum() / self.allCells.passingReads.sum():.1%}"]) 
        stats.append(["Median Saturation", passingCells.Saturation.median()])
        meanReads = passingCells.reads.mean()
        medReads = passingCells.passingReads.median()
        medUmis = passingCells.totalUmis.median()
        stats.append(['Unique reads @ 1000 mean reads (extrapolated)', f"{CalculationUtils.extrapolate_umis_at_nreads(medReads, medUmis, CalculationUtils.downsampleTargetReadCount(meanReads, medReads, 1000)):,}"])
        stats.append(['Unique reads @ 5000 mean reads (extrapolated)', f"{CalculationUtils.extrapolate_umis_at_nreads(medReads, medUmis, CalculationUtils.downsampleTargetReadCount(meanReads, medReads, 5000)):,}"])
        stats.append(['Unique reads @ 10k mean reads (extrapolated)', f"{CalculationUtils.extrapolate_umis_at_nreads(medReads, medUmis, CalculationUtils.downsampleTargetReadCount(meanReads, medReads, 10000)):,}"])
        stats.append(['Unique reads @ 20k mean reads (extrapolated)', f"{CalculationUtils.extrapolate_umis_at_nreads(medReads, medUmis, CalculationUtils.downsampleTargetReadCount(meanReads, medReads, 20000)):,}"])
        stats.append(['Unique reads @ 100k mean reads (extrapolated)', f"{CalculationUtils.extrapolate_umis_at_nreads(medReads, medUmis, CalculationUtils.downsampleTargetReadCount(meanReads, medReads, 100000)):,}"])
        
        for (title,colName) in statsToDeriveMediansFor.items():
            stats.append([f"Median {title}", f"{passingCells[colName].median():.1%}"])
        statsAsDictList = [GeneralUtils.makeTableDict(['Metric','Value'], valuePair) for valuePair in stats]
        statsDataFrame =  pd.DataFrame(statsAsDictList)
        statsDataFrame['Category'] = 'Cell Calling'
        return statsDataFrame

    def makeRankPlot(self):
        '''
        Makes a kneeplot using @field in @data; drawing a vertical line at @threshold
        '''
        # Make copy to prevent overwriting of barcode indexes
        allCellsCopy = self.allCells.copy(deep=True)
        allCellsCopy.reset_index(drop=True,inplace=True) 
        indices = CalculationUtils.getIndicesToInclude(len(allCellsCopy.index))
        rowsToInclude = allCellsCopy.iloc[indices]
        fig = px.line(rowsToInclude, x=rowsToInclude.index, y="totalUmis", labels={"index":"Cell Barcodes", "umis":"Unique Reads"}, log_x=True, log_y=True, template=constants.DEFAULT_FIGURE_STYLE)
        fig.add_vline(x=allCellsCopy['pass'].sum(), line_dash="dash", line_color="green")
        fig.update_layout(xaxis_range=[1, math.log10(rowsToInclude.index.max())])
        self.saveFigureAsPng(fig, "UMIRankPlot.png")
        return dp.Plot(fig)
    
    def addStatsAndCallCells(self, starStats):
        cellReadStats = pd.read_csv(self.sampleSpecificFilePaths['cellReadsFn'], sep='\t', index_col='CB')
        cellReadStats = cellReadStats[1:]
        self.allCells['reads'] = cellReadStats.cbMatch
        self.allCells['genomeMappedReads'] = cellReadStats.genomeU / self.allCells.reads
        self.allCells['geneReads'] = cellReadStats.featureU + cellReadStats.featureM
        self.allCells['exon'] = cellReadStats.exonic / self.allCells.geneReads
        self.allCells['intron'] = cellReadStats.intronic / self.allCells.geneReads
        self.allCells['antisense'] = (cellReadStats.exonicAS + cellReadStats.intronicAS) / self.allCells.geneReads
        self.allCells['uniqueMappingReads'] = (cellReadStats.countedU)
        self.allCells['passingReads'] = (cellReadStats.countedU + cellReadStats.countedM)
        self.allCells['Saturation'] = 1-(self.allCells.totalUmis / self.allCells.passingReads)
        self.allCells['mito'] = cellReadStats.mito / (cellReadStats.genomeU + cellReadStats.genomeM)
        self.allCells['genes'] = self.allCells[["humanGenes", "mouseGenes"]].max(1)
        self.allCells['umis']  = self.allCells[["humanUmis", "mouseUmis"]].max(1)
        self.allCells['propGeneMatched'] = self.allCells.geneReads / self.allCells.reads
        self.allCells.sort_values('totalUmis', ascending=False, inplace=True)
        self.allCells['pass'] = self.allCells.totalUmis >= self.cellThreshold

    def instantiateCellCountsDf(self):
        cells = collections.defaultdict(CellStat)
        with open(self.sampleSpecificFilePaths['matrixFn'], "r") as matrixFile:
            # Read through header 
            matrixFile.readline();matrixFile.readline();matrixFile.readline()
            for entry in matrixFile:
                line = entry.split()
                lineParsedAsInt = [int(x) for x in line]
                gene, cell, count, *_ = lineParsedAsInt
                if self.genes.species[gene] == "hs":
                    cells[cell].totalGenes += 1
                    cells[cell].humanGenes += 1
                    cells[cell].humanUmis += count
                    cells[cell].totalUmis += count
                elif self.genes.species[gene] == "mm":
                    cells[cell].totalGenes += 1
                    cells[cell].mouseGenes += 1
                    cells[cell].mouseUmis += count
                    cells[cell].totalUmis += count
                else:
                    print(self.genes.iloc[gene])
        self.allCells = pd.DataFrame(cells.values(), index=[self.barcodes.Barcode[i] for i in cells])

    def instantiateGenesDf(self):
        '''
        Read in sample features.tsv from cell x gene count matrix as Dataframe
        label each gene as belonging to either mm (mouse) or hs (human)
        '''
        detectedGenes = pd.read_csv(self.sampleSpecificFilePaths['featuresFn'], sep="\t", header=None)
        detectedGenes.columns = ["Id", "Name", "Type"]
        detectedGenes.index = detectedGenes.index + 1
        detectedGenes['species'] = detectedGenes.Id.map(lambda x: ["mm","hs"][x.startswith("ENSG")])
        self.genes = detectedGenes

    def instantiateBarcodesDf(self):
        '''
        Read in sample barcodes.tsv from cell x gene count matrix as Dataframe
        '''
        cellBarcodes = pd.read_csv(self.sampleSpecificFilePaths['barcodesFn'], sep='\t', header=None)
        cellBarcodes.columns=['Barcode']
        cellBarcodes.index = cellBarcodes.index + 1
        self.barcodes = cellBarcodes

    def checkIfBarnyard(self):
        self.scoreBarn()
        includedSpecies = list(self.allCells.species.unique())
        return "Human" in includedSpecies and "Mouse" in includedSpecies

    @staticmethod   
    def makeStarDfFromDict(starStats: Dict, category="STARsolo") -> pd.DataFrame:
        stats = []
        for (key, value) in starStats.items():
            if key != 'STAR Cells' and key != "Saturation":
                stats.append([key, value])
        statsAsDictList = [GeneralUtils.makeTableDict(['Metric','Value'], valuePair) for valuePair in stats]
        statsDataFrame =  pd.DataFrame(statsAsDictList)
        statsDataFrame['Category'] = category
        return statsDataFrame
    
    @staticmethod   
    def getBackground(minorFracs):
        if len(minorFracs) == 0: return np.nan,np.nan
        median = minorFracs.median()
        std = (np.quantile(minorFracs, 0.75) - median) / 0.675
        return median, std
        

    def scoreBarn(self):
        '''
        Labels each cell as belonging to mouse or human 
        '''
        cells = self.allCells
        cells['species'] = "None"
        cells['minorFrac'] = cells[['humanUmis', 'mouseUmis']].min(1) / cells[['humanUmis', 'mouseUmis']].sum(1)
        cells.loc[cells['pass'] & (cells.humanUmis > cells.mouseUmis), 'species'] = 'Human'
        cells.loc[cells['pass'] & (cells.mouseUmis > cells.humanUmis), 'species'] = 'Mouse'

        minHuman = minMouse = cells[cells['pass']].umis.min()
        if (cells.species=='Human').any() and (cells.species=='Mouse').any():
            minHuman = np.percentile(cells.humanUmis[cells.species=='Human'], 10)
            minMouse = np.percentile(cells.mouseUmis[cells.species=='Mouse'], 10)
        
        # Labelling low and doublet cells
        cells.loc[(cells.humanUmis < minHuman) & (cells.mouseUmis < minMouse), 'species'] = "None"
        cells.loc[(cells.humanUmis >= minHuman) & (cells.mouseUmis >= minMouse), 'species'] = "Mixed"
        # Labelling as Ambiguous
        humanBackMed, humanBackStd = self.getBackground(cells[(cells.species =='Mouse')].minorFrac)
        mouseBackMed, mouseBackStd = self.getBackground(cells[(cells.species =='Mouse')].minorFrac)
        cells.loc[(cells.species == 'Mixed') & (cells.humanUmis > cells.mouseUmis) & (cells.minorFrac < mouseBackMed + 2 * mouseBackStd), 'species'] = 'Ambiguous'
        cells.loc[(cells.species == 'Mixed') & (cells.mouseUmis > cells.humanUmis) & (cells.minorFrac < humanBackMed + 2 * humanBackStd), 'species'] = 'Ambiguous'
        cells.loc[(cells.species == 'Human') & (cells.minorFrac >= max(0.1, mouseBackMed + 3 * mouseBackStd)), 'species'] = 'Ambiguous'
        cells.loc[(cells.species == 'Mouse') & (cells.minorFrac >= max(0.1, humanBackMed + 3 * humanBackStd)), 'species'] = 'Ambiguous'

    def makeBarnyardStatsDf(self):
        passingBarnyardCells = self.allCells[(self.allCells['pass']) & (~self.allCells.species.isin(['Ambiguous', 'None']))]
        ncells = len(passingBarnyardCells.index)
        stats = []
        propHuman = (passingBarnyardCells.species=='Human').sum() / ncells
        stats.append(['Human Cells', f"{propHuman:.1%}"])
        propMouse = (passingBarnyardCells.species=='Mouse').sum() / ncells
        stats.append(['Mouse Cells', f"{propMouse:.1%}"])
        propMixed = (passingBarnyardCells.species=='Mixed').sum() / ncells
        stats.append(['Mixed Cells', f"{propMixed:.1%}"])
        stats.append(['Estimated Doublets', f"{min(1, (propMixed/(2*propHuman*propMouse))):.1%}"])
        stats.append(['Background', f"{passingBarnyardCells[passingBarnyardCells.species != 'Mixed'].minorFrac.median():.2%}"])
        statsAsDictList = [GeneralUtils.makeTableDict(['Metric','Value'], valuePair) for valuePair in stats]    
        result = pd.DataFrame(statsAsDictList)
        result['Category'] = 'Barnyard'
        return result

    def resolveSampleSpecificFilePaths(self, forFastqReport:bool) -> Dict[str, Union[Path, Dict[str, Path]]] :
        """
        Returns a dictionary of the needed file paths. 
        Differs based on whether -results argument was 
        provided during invocation  
        - isn't provided: (running from within ScaleAtac workflow)
        - provided: running independently 
        """
        outputDict = {}
        if (self.resultsDir is None):
            self.resultsDir = Path('.')
            soloOutName = f"{self.sampleName}.Solo.out" if forFastqReport else "Solo.out" 
            soloDirPrefix = Path(soloOutName)
        else:
            self.resultsDir = Path(self.resultsDir)
            soloDirPrefix = self.resultsDir / 'star' / f"{self.sampleName}.star" / 'Solo.out'
        exprMatrixPrefix = soloDirPrefix / self.featuresType 
        outputDict = dict(featuresFn=exprMatrixPrefix / "raw" / 'features.tsv',
                    barcodesFn=exprMatrixPrefix / "raw" / 'barcodes.tsv',
                    matrixFn=exprMatrixPrefix / "raw" / 'matrix.mtx', 
                    summaryFn=exprMatrixPrefix / "Summary.csv", 
                    cellReadsFn=exprMatrixPrefix / "CellReads.stats")
        return outputDict

