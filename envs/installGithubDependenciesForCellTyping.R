#!/usr/bin/env Rscript
packages <- c("remotes")
suppressMessages(invisible(lapply(packages, library, character.only = TRUE, quietly=TRUE)))
Sys.setenv(CONDA_BUILD_SYSROOT="/")
remotes::install_github("immunogenomics/presto",upgrade="always", build_opts=c("--no-build-vignettes"))
remotes::install_github("satijalab/seurat-data", build_opts=c("--no-build-vignettes"))
remotes::install_github("satijalab/azimuth", build_opts=c("--no-build-vignettes"))
library(SeuratData)
InstallData("pbmcref")
