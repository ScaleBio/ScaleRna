#!/usr/bin/env python

from reportGeneration.ReportUtil import GeneralUtils, CalculationUtils
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict

class AbstractReport(ABC):
    resultsDir: Path
    makeAdvanced: bool
    demuxJson: Dict[str, float]
    resultsDir: Path
    writeDir: Path
    libDir: Path

    @abstractmethod
    def __init__(self, resultsDir: Path, makeAdvanced:bool, demuxJson: Dict[str,float], libDir:Path):
        self.resultsDir = resultsDir
        self.makeAdvanced = makeAdvanced
        self.demuxJson = demuxJson
        self.resultsDir = resultsDir
        self.writeDir = GeneralUtils.resolveWriteDir(resultsDir)
        self.libDir = libDir

    @abstractmethod
    def build(self):
        '''
        Builds and saves an interactive HTML report
        '''
        pass

    @abstractmethod
    def buildReadsPage(self):
        '''
        Creates a dp.Page to be included in the HTML report labelled "Reads"
        '''
        pass

    