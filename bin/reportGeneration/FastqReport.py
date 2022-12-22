import functools
import pandas as pd 
import datapane as dp
import reportGeneration.myconstants as constants
import plotly.express as px
from reportGeneration.ReportUtil import GeneralUtils, CalculationUtils, DatapaneUtils
from reportGeneration.AbstractReport import AbstractReport
from reportGeneration.SampleReport import SampleReport
from pathlib import Path
from typing import Dict, List, Tuple


# Threshjson, filterDf (both in QC)

class FastqReport(AbstractReport):
    libName: str
    sampleNames: List[str]
    allCellsBetweenFiles: pd.DataFrame
    sampleFiles: List[Tuple[str,Dict[str,Path]]]
    sampleQCDf: pd.DataFrame
    sampleThresholds: Dict[str, Dict[str, float]]
    def __init__(self, libName: str, samplesheetFn: str, resultsDir: Path, featureType: str, makeAdvanced:bool, demuxJson: Dict[str,float], libDir:Path):
        super().__init__(resultsDir, makeAdvanced, demuxJson, libDir)
        self.libName = libName
        self.featureType = featureType
        self.sampleNames = sorted(FastqReport.getAssociatedSampleNames(samplesheetFn,libName))
        self.allCellsBetweenFiles= self.instantiateCellCountsDf()
        
    def instantiateCellCountsDf(self):
        '''
        Instantiate a multisample dataframe containing cell-level information for cells of all samples derived from same FASTQ
        '''
        sampleReportObjs = [SampleReport(sn,self.resultsDir, self.featureType, None, self.makeAdvanced, self.demuxJson, self.libDir, True) for sn in self.sampleNames]
        accumulatedDf = pd.DataFrame()
        thresholds = {}
        totalSamples = len(sampleReportObjs)
        for i,sampleObj in enumerate(sampleReportObjs):
            print(f"Processing Sample Data: {i + 1}/{totalSamples}")
            sampleName = self.sampleNames[i]
            sampleObj.instantiateGenesDf()
            sampleObj.instantiateBarcodesDf()
            sampleObj.instantiateCellCountsDf()
            sampleObj.allCells['sample'] = sampleName
            sampleObj.allCells['umis']  = sampleObj.allCells[["humanUmis", "mouseUmis"]].max(1)
            sampleObj.allCells.sort_values('umis', ascending=False, inplace=True)
            starStats = sampleObj.getStarStats()
            ncells = int(starStats['STAR Cells'])
            threshold = max(100, sampleObj.allCells.umis[ncells-1])
            thresholds[sampleName] = threshold
            sampleObj.allCells['Barcode'] = sampleObj.allCells.index
            index = range(len(sampleObj.allCells.index))
            # Assigning index prevents overwritting during concatenation 
            sampleObj.allCells.index = index
            accumulatedDf = pd.concat([accumulatedDf, sampleObj.allCells])
        self.sampleThresholds = thresholds
        return accumulatedDf

    def validateInput(self):
        GeneralUtils.makeDir(self.writeDir)
        for (_, sampleFilePaths) in self.sampleFiles.items():
            GeneralUtils.ensurePathsExist(sampleFilePaths)

    def build(self):
        (overallPassingStats, readsPage) = self.buildReadsPage()
        (barcodeStatsDf, cellsPage) = self.buildCellsPage()
        pages = [readsPage, cellsPage]
        report = dp.Report(
            blocks=pages
        )
        prefix = f"library_{self.libName}"
        barcodeStatsDf.to_csv(self.writeDir / f"{prefix}.typeLevelMatches.tsv", sep='\t', index=False)
        overallPassingStats.to_csv(self.writeDir / f"{prefix}.overallMatches.tsv", sep='\t', index=False)
        report.save(self.writeDir / f"{prefix}.report.html") 

    ###################### READS PAGE FUNCTIONS ##########################
    def buildMultisampleCellStatsDf(self):
        selectBlocks = []
        for sample in self.sampleNames:
            sampleSubsetDf = self.sampleQCDf[self.sampleQCDf['sample'] == sample]
            cellStatsDataFrame = SampleReport.makeCellStatsDf(sampleSubsetDf, self.sampleThresholds[sample]['UniqueReads'])
            cellStatsDataFrameStyled = cellStatsDataFrame.style.pipe(GeneralUtils.styleTable, title="Reads", hideColumnHeaders=True, boldColumn=['Metric'])
            cellStatsAsDpTable = DatapaneUtils.createTableIgnoreWarning(cellStatsDataFrameStyled, label=sample)
            selectBlocks.append(cellStatsAsDpTable)
        return dp.Select(type=dp.SelectType.TABS,blocks=selectBlocks)

    def buildDfFromDemuxSampleMetrics(self):
        totalCounts = []
        tgmtCounts = []
        for sampleName, sampleDict in self.demuxJson['samples'].items():
            readCount = sampleDict['reads'][0]
            totalCounts.append({'Sample':sampleName, 'TotalReads':readCount})
            tgmtBarcodeCounts = sampleDict['barcodes']
            if (len(tgmtBarcodeCounts) > 0):
                for tgmtId, metrics in tgmtBarcodeCounts.items():
                    tgmtCounts.append({'Sample': sampleName, 'tgmtBc': tgmtId, 'ReadCount': metrics['reads']})
        return (pd.DataFrame(totalCounts), pd.DataFrame(tgmtCounts))

    def getCellsAboveThreshold(self) -> Dict[str, int]:
        results = {}
        for sample, thresholds in self.sampleThresholds.items():
            cellsPassed = len(self.sampleQCDf[(self.sampleQCDf['UniqueReads'] >= thresholds['UniqueReads']) & (self.sampleQCDf['sample'] == sample)].index)
            results[sample] = cellsPassed
        return results

    def makeMultisampleKneePlot(self) -> dp.Plot:
        '''
        Makes a kneeplot using @field in @data; drawing a vertical line at @threshold
        '''  
        indiciesToInclude = set(CalculationUtils.getIndicesToInclude(max(self.allCellsBetweenFiles.index)))
        self.allCellsBetweenFiles['rankOrder'] = self.allCellsBetweenFiles.index
        plottingDf = self.allCellsBetweenFiles[self.allCellsBetweenFiles['rankOrder'].isin(indiciesToInclude)]
        fig = px.line(plottingDf, x=plottingDf.index, y='umis', color='sample', log_x=True, log_y=True, template=constants.DEFAULT_FIGURE_STYLE, labels={"index": "Cell Barcodes" })
        return dp.Plot(fig)

    def buildReadsPage(self):
        #sampleCellsAboveThreshold = self.getCellsAboveThreshold()
        multiSampleKneePlot = self.makeMultisampleKneePlot()
        barcodeReadsData = self.demuxJson['reads']
        barcodeReadsPerc = FastqReport.buildDfFromJSONDict(barcodeReadsData, "Type", "list")
        barcodeReadsTotal = FastqReport.buildDfFromJSONDict(barcodeReadsData, "Type", "list", 0)
        barcodeReadsTotal = barcodeReadsTotal[['Type', 'Reads']]
        barcodeReadsTotal.rename(columns={'Type':'Status'}, inplace=True)
        barcodeReadsTotal['Percent'] = barcodeReadsPerc['Reads']
        barcodeReadsTotalStyled = barcodeReadsTotal.style.pipe(GeneralUtils.styleTable, title="Barcode Read Status", numericCols=['Reads'])    
        barcodeReadStats = DatapaneUtils.createTableIgnoreWarning(barcodeReadsTotalStyled)
        (countsPerSampleDf, tgmtCountsPerSampleDf) = self.buildDfFromDemuxSampleMetrics()
        sortOrder = sorted(list(tgmtCountsPerSampleDf.tgmtBc.unique()), key=functools.cmp_to_key(CalculationUtils.wellStringComp))
        tgmtCountsPerSampleDf['tgmtBc'] = pd.Categorical(tgmtCountsPerSampleDf.tgmtBc, sortOrder)
        tgmtCountsPerSampleDf.sort_values(by=['tgmtBc'], inplace=True, ascending=False)
        sampleOrder = list(tgmtCountsPerSampleDf.Sample.unique())
        sampleOrder.reverse()
        countsPerSampleDf['Sample'] = pd.Categorical(countsPerSampleDf.Sample, sampleOrder + ["Unknown"])
        countsPerSampleDf.sort_values(by=['Sample'], inplace=True)
        colorMapToUse = FastqReport.matchColorsToNames(px.colors.qualitative.D3, list(countsPerSampleDf['Sample'].unique()))
        readsPerSample = px.bar(countsPerSampleDf, x='TotalReads', y='Sample', color='Sample',height=900, color_discrete_map=colorMapToUse,template=constants.DEFAULT_FIGURE_STYLE, title="Reads Per Sample", labels={"TotalReads":"Total Reads"})
        readsPerSample.update_layout(showlegend=False)
        readsPerWellBySample = px.bar(tgmtCountsPerSampleDf, x='ReadCount', y='tgmtBc', color='Sample',height=900,template=constants.DEFAULT_FIGURE_STYLE, color_discrete_map=colorMapToUse, labels={'tgmtBc': 'RT Barcodes', 'ReadCount':'Total Reads'}, title="Reads Per RT Well")
        return (barcodeReadsTotal,dp.Page(
                blocks=[
                    dp.Text(f"## libName: {self.libName}"),dp.Group(multiSampleKneePlot, barcodeReadStats,readsPerSample,readsPerWellBySample,columns=2)
                ],
                title='Reads')
        )    

    def matchColorsToNames(colorPalette, names) -> Dict[str,str]:
        '''
        Associated colors with categorical values @names
        for consistency between plots 
        '''
        colorMap = {"Unknown": 'rgb(179, 188, 201)'}
        if ("Unknown" in names):
            names.remove('Unknown')
        colors = px.colors.qualitative.D3
        for i,name in enumerate(sorted(names)):
            colorMap[name] = colors[i % len(colors)]
        return colorMap


    ###################### CELLS PAGE FUNCTIONS ##########################
    def buildMultisampleOverloadingFigure(self):
        barcodeNames = {}
        for sample in self.sampleNames:
            for name,bc in self.demuxJson["samples"][sample]["barcodes"].items():
                barcodeNames[bc['sequence']] =  name 
        tagmentBcPos = 0
        dropBcPos = 1
        tgmtBcs = []
        dropBcs = []
        for bc in self.sampleQCDf.Barcode:        
            bcs = bc.split('+')            
            tgmtBc = bcs[tagmentBcPos]
            tgmtBcs.append(barcodeNames.get(tgmtBc, tgmtBc))
            dropBcs.append(bcs[dropBcPos])
        self.sampleQCDf['tgmtBc'] = tgmtBcs
        self.sampleQCDf['dropBc'] = dropBcs
        cells = self.sampleQCDf[self.sampleQCDf['Filter'] == 'Pass']
        dropCnts = cells.groupby(['dropBc', 'sample']).size().reset_index()
        dropCnts.rename(columns={0: 'BarcodesInBead'}, inplace=True)
        cells['BarcodesInBead'] = cells.groupby(['dropBc', 'sample'])['dropBc'].transform('size')
        cellsPerDropBar = px.box(cells, y='UniqueReads', x='BarcodesInBead', log_y=True, color='sample', title='Unique reads per cell ',template="none", labels={"BarcodesInBead": "Cells Per Bead", "UniqueReads":"Unique Reads"}, points=False)
        DatapaneUtils.showAllTicks(cellsPerDropBar)
        barcodesPerBeadBySample = dropCnts.groupby(['BarcodesInBead', 'sample']).size().reset_index()
        barcodesPerBeadBySample.rename(columns={0: 'Counts'}, inplace=True)
        cellsPerDropLine = px.bar(barcodesPerBeadBySample, x='BarcodesInBead', y='Counts', color='sample',log_y=True,  template=constants.DEFAULT_FIGURE_STYLE, title="Cells per bead", labels={"BarcodesInBead": "Cells Per Bead"})
        DatapaneUtils.showAllTicks(cellsPerDropLine)
        return dp.Group(dp.Plot(cellsPerDropLine),dp.Plot(cellsPerDropBar), columns=2)
   
    def createBarcodeTypeMetricsTables(self):
        barcodes = self.demuxJson['barcodes']
        barcodesDf = FastqReport.buildDfFromJSONDict(barcodes, "Barcode", "dict")
        tableGroup = []
        allBarcodeTypes = list(barcodesDf['Barcode'].unique())
        for barcodeType in allBarcodeTypes:
            fullBarcodeTypeName = constants.BARCODE_SHORTHAND_TO_NAME[barcodeType]
            subset = barcodesDf[barcodesDf['Barcode'] == barcodeType][['Match', 'Reads']]
            styledDf = subset.style.pipe(GeneralUtils.styleTable, title=fullBarcodeTypeName)                       
            table = DatapaneUtils.createTableIgnoreWarning(styledDf)
            tableGroup.append(table)
        return (barcodesDf, dp.Group(blocks=tableGroup, columns=2))

    def buildCellsPage(self):    
        (barcodeTypeStatsDf, barcodeTypeStats) = self.createBarcodeTypeMetricsTables()
        blocks=[dp.Text(f"## libName: {self.libName}"),barcodeTypeStats]
        if self.makeAdvanced:
            overloadingFigures = self.buildMultisampleOverloadingFigure()
            blocks.append(overloadingFigures)
        return (barcodeTypeStatsDf,dp.Page(
                blocks=blocks,
                title='Cells')
    )

    ###################### STATIC HELPERS ##########################
    @staticmethod
    def getAssociatedSampleNames(sampleSheetFn: str, libName: str): 
        '''
        Returns a list of the values in the 'sample' column of @sampleSheetFn 
        for rows where 'libName' equals @libName
        '''
        sampleSheetDf = pd.read_csv(sampleSheetFn)
        associatedSampleMetaData = sampleSheetDf[sampleSheetDf['libName'] == libName]
        associatedSampleNames = list(associatedSampleMetaData['sample'])
        return associatedSampleNames

    @staticmethod
    def buildDfFromJSONDict(jsonDict: Dict, name: str, valueDataType: str, choiceIndex=1) -> pd.DataFrame:
        """
        Builds a dataframe from a json containing metrics  
        Used as helper to process demux metrics.json 
        """
        dictList = []
        for (keyName,valueObj) in jsonDict.items():
            if valueDataType == 'dict':
                for key,value in valueObj.items():
                    newDict={}
                    newDict[name]=keyName
                    newDict['Match']=key
                    newDict['Reads']=value[choiceIndex]
                    dictList.append(newDict)
            elif valueDataType == 'list':
                newDict = {"Reads": valueObj[choiceIndex]}
                newDict[name] = keyName
                dictList.append(newDict)
            elif valueDataType == 'int':
                newDict = {"Cells": valueObj}
                newDict[name] = keyName
                dictList.append(newDict)
        return pd.DataFrame(dictList)


