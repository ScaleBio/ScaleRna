#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(BPCells))

#This function reads the UMI count matrix from the samples directory.
#Matrix must follow the naming convention as follows:
#features.tsv
#barcodes.tsv
readCountMatrix <- function(mtxDir, mtxName = "matrix.mtx"){

    message("Reading Count Matrix")

    mtx <- file.path(mtxDir, mtxName)
    genes <- file.path(mtxDir, "features.tsv")
    cells <- file.path(mtxDir, "barcodes.tsv")

    sampleId <- strsplit(x = mtx, split = "_")
    sampleId <- unlist(sampleId)[[1]]

    genes <- fread(file = genes, header = FALSE)[[2]]
    genes <- make.unique(names = genes, sep = "-")
    cells <- fread(file = cells, header = FALSE)[[1]]
    cells <- paste0(cells, "_", sampleId)

    mtx <- import_matrix_market(
      mtx_path = mtx,
      row_names = genes,
      col_names = cells
    )

    return(mtx)
}

main <- function(){
  
  argParser <- ArgumentParser(
    description = "This script creates an annData object."
  )

  argParser$add_argument("--matrixDir", required = TRUE, nargs = "+")
  argParser$add_argument("--sampleId", required = TRUE)

  argList <- argParser$parse_args()

  countMat <- lapply(argList$matrixDir, readCountMatrix)
  
  #If there are more than one matrix col bind into single matrix.
  if(length(countMat) > 1){
    countMat <- do.call(cbind, countMat)
  } else {
    countMat <- countMat[[1]]
  }

  message("Writing AnnData")
  outFn <- paste0(argList$sampleId, "_anndata.h5ad")
  write_matrix_anndata_hdf5(countMat, outFn)
}

main()
