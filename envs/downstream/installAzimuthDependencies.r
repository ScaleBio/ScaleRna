#!/usr/bin/env Rscript

remotes::install_github('satijalab/azimuth', ref = 'master', dependencies = TRUE, upgrade = "never")

SeuratData::InstallData("pbmcref")
