#!/usr/bin/env python
import pandas as pd
import argparse
from os.path import exists
from pathlib import Path
from typing import Union
from reportGeneration.ReportUtil import GeneralUtils
from reportGeneration.SampleReport import SampleReport
from reportGeneration.FastqReport import FastqReport

def validateArgs(sampleName:Union[str,None], libName:Union[str,None], sampleSheetFn:str, resultsDir:Union[str,None]) -> Union[str,None]:
    '''
    Validates:
    - @sampleSheetFn and @resultsDir exist
    - specified @sampleName, @libName exist in @sampleSheetFn if passed. 
    '''
    if (not exists(sampleSheetFn)):
        raise ValueError(f"A samplesheet does not exist at the specified path: {sampleSheetFn}")
    samplesheet = pd.read_csv(sampleSheetFn)
    if (sampleName is None and libName is None):
        raise ValueError(f"Specify a -sample or -libName for which to generate a report")
    # sample argument validation 
    availableSamples = list(samplesheet['sample'])
    if (sampleName is not None and not sampleName in availableSamples):
        raise ValueError(f"The sample provided:{sampleName} doesn't exist in specified samplesheet.\n Existing samples: {availableSamples}")
    # libName argument validation 
    availableFastqNames = list(samplesheet['libName'].unique())
    if (libName is not None and not libName in availableFastqNames):
        raise ValueError(f"The libName pvoided:{libName} doesn't exist in specified samplesheet\n Existing libNames: {availableFastqNames}")
    # results argument validation 
    if resultsDir is not None: 
        if (not exists(resultsDir)):
            raise FileExistsError(f"The results directory provided {resultsDir} doesn't exist")
        else: 
            resultsDir = resultsDir.rstrip("/")
    return resultsDir

def resolveFastqSpecificFilePath(resultsDir: Union[str,None],libName:Union[str,None]) -> Path:
    if (resultsDir is None):
        return Path("demuxMetrics.json")
    else:
        return Path(resultsDir , 'demux' , f'{libName}.demux' , 'metrics.json')

def getFastqName(sampleName:str,samplesheetFn:str):
    '''
    Retrieves the 'libName' field for the 'sample' equal to @sampleName 
    from the csv at @samplesheetFn

    Assumption: values in the 'sample' column of @samplesheetFn are unique 
    '''
    sampleSheet = pd.read_csv(samplesheetFn)
    sampleData = sampleSheet.loc[sampleSheet['sample'] == sampleName].iloc[0]
    return sampleData['libName']

def createReports(sampleName:Union[str,None], samplesheetFn:str, resultsDir:str, passedFastqName:Union[str,None], featureType:str, cells:Union[int,None], makeAdvanced:bool, libDir:Path) -> None:
    '''
    Writes Interactive HTML report(s):
    - If @sampleName is not None: for sample with @sampleName
    - If @passedFastqName is not Nont: for fastq with @libName
    '''
    makeSampleReport = sampleName is not None
    makeFastqReport = passedFastqName is not None
    libName = passedFastqName if makeFastqReport else getFastqName(sampleName, samplesheetFn)
    demuxMetricsPath = resolveFastqSpecificFilePath(resultsDir, libName)
    demuxJson = GeneralUtils.readJSON(demuxMetricsPath)
    if (makeSampleReport):
        reportObj = SampleReport(sampleName, resultsDir, featureType, cells, makeAdvanced, demuxJson, libDir)
        reportObj.build()
    if makeFastqReport:
        reportObj = FastqReport(libName, samplesheetFn, resultsDir, featureType, makeAdvanced, demuxJson, libDir)
        reportObj.build()
    
def main():
    parser = argparse.ArgumentParser(description="An interactive HTML report generator for the ScaleBio universal ATAC pipeline")
    parser.add_argument("--sample", metavar="sampleName", type=str, required=False, default=None,help="The name of the sample for which a report is being generated")
    parser.add_argument("--samplesheet", metavar="sampleSheet", required=True, default=None, help="Path to the samples.csv containing information about libName")
    parser.add_argument("--results", metavar="resultsFolder", type=str, help="Path to the output folder for which the report is to be generated", required=False, default=None)
    parser.add_argument("--libName", metavar="libName", default=None, type=str,  help="libName specified in passed samplesheet for which to generate a fastq report?")
    parser.add_argument("--libStruct", metavar="libraryStructureDirectory", required=False, default=None)
    parser.add_argument("--featureType", default="GeneFull_Ex50pAS", type=str,  help="STARSolo feature type used")
    parser.add_argument("--cells", default=None, type=int, help="Set a fixed number of cells to pass for UMI count threshold")
    parser.add_argument("--advanced", default=False, action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    resultsDir = validateArgs(args.sample, args.libName, args.samplesheet, args.results)
    createReports(args.sample, args.samplesheet, resultsDir, args.libName, args.featureType, args.cells, args.advanced, args.libStruct)

if __name__ == "__main__":
    main()
