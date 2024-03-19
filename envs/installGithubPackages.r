#!/usr/bin/env Rscript

remotes::install_github('satijalab/azimuth', ref = 'master', upgrade = "never")
remotes::install_github("bnprks/BPCells", upgrade = "never")
SeuratData::InstallData("pbmcref")