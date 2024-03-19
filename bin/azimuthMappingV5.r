#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Azimuth))
suppressPackageStartupMessages(library(BPCells))
suppressPackageStartupMessages(library(Seurat))

#This function reads the UMI count matrix from the samples directory.
#Matrix must follow the naming convention as follows:
#features.tsv
#barcodes.tsv
readCountMatrix <- function(mtxDir, mtxName){

    message("Reading Count Matrix")

    sampleId <- unlist(strsplit(mtxDir, split = "_"))
    sampleId <- sampleId[[1]]

    mtx <- file.path(mtxDir, mtxName)
    genes <- file.path(mtxDir, "features.tsv")
    cells <- file.path(mtxDir, "barcodes.tsv")

    genes <- fread(file = genes, header = FALSE)[[2]]
    genes <- make.unique(names = genes, sep = "-")
    cells <- fread(file = cells, header = FALSE)
    if(nrow(cells) == 0){
      message("No Cells in Dataset")
      quit(save = "no")
    } else {
      message(paste0("Cells in Data: ", nrow(cells)))
      cells <- cells[[1]]
    }
    cells <- paste0(cells, "_", sampleId)

    #Directory name where BPCells will write mtx.
    bp_dir <- paste0(sampleId, "_bpcells")

    mtx <- import_matrix_market(
      mtx_path = mtx,
      row_names = genes,
      col_names = cells,
      tmpdir = file.path(bp_dir, "tmp"),
      outdir = file.path(bp_dir, "out")
    )
    #Modify dir value in mtx to relative path.
    #By default BPCells uses absolute paths but this causes issues when trying to move seurat objects.
    mtx@dir <- file.path(bp_dir, "out")

    return(mtx)
}

readAllCellsData <- function(filePath){
  
  message("Reading allCells.csv")
  
  dat <- fread(file = filePath)
  dat <- setnames(dat, old = 1, new = "cellBarcode")
  cellBarcodes <- paste0(dat[["cellBarcode"]], "_", dat[["sample"]])
  dat[["cellBarcode"]] <- cellBarcodes
  dat <- dat[dat$pass == TRUE, ]
  
  return(dat)
}

loadData <- function(argList){
  
  message("Loading Data")

  if(argList$workflow == "standard"){

    counts <- readCountMatrix(mtxDir = argList$matrixDir, mtxName = argList$star_matrix)
    meta <- readAllCellsData(filePath = argList$allCells)
    seurat <- CreateSeuratObject(
      counts = counts,
      meta.data = meta,
      project = argList$project
    )

  } else if(argList$workflow == "comparison"){
    
    matList <- lapply(argList$matrixDir, readCountMatrix, mtxName = "matrix.mtx")
    meta <- lapply(argList$allCells, readAllCellsData)
    meta <- rbindlist(meta)
    seurat <- CreateSeuratObject(
      counts = matList,
      meta.data = meta,
      project = argList$project
    )
    seurat <- JoinLayers(object = seurat)
    DefaultAssay(seurat) <- "RNA"
  }
  return(seurat)
}

parseArguments <- function(){

    message("Parsing Arguments")
    
    argParser <- ArgumentParser(
        description = "This script performs Azimuth reference based cell type annotation of a given sample."
    )
    argParser$add_argument("--reference", action = "store", required = TRUE,
        help = "Name of Azimuth Reference Dataset to use.")
    argParser$add_argument("--paramPath", action = "store", required = FALSE, default = FALSE,
        help = "Path to yaml file with Seurat function parameters.")
    argParser$add_argument("--star_matrix", action = "store", required = FALSE, default = "matrix.mtx",
        help = "Name of matrix file to use when loading UMI counts.")
    argParser$add_argument("--project", action = "store", required = TRUE)

    subParsers <- argParser$add_subparsers(
      help = "This subcommand determines how to load the data."
    )

    parserStandard <- subParsers$add_parser("standard", help = "standard help")
    parserStandard$set_defaults(workflow = "standard")
    parserStandard$add_argument("--matrixDir", required = TRUE, 
        help = "Path to directory which contains the UMI count matrix.")
    parserStandard$add_argument("--allCells", required = TRUE,
        help = "Path to the allCells.csv metric file.")

    parserComparison <- subParsers$add_parser("comparison", help = "comparison help")
    parserComparison$set_defaults(workflow = "comparison")
    parserComparison$add_argument("--matrixDir", required = TRUE, nargs = "+")
    parserComparison$add_argument("--allCells", required = TRUE, nargs = "+")

    argList <- argParser$parse_args()
    
    return(argList)
}


writeResults <- function(azimuthObject){
  
  azimuthObject <- AddMetaData(object = azimuthObject, metadata = azimuthObject@reductions$ref.umap@cell.embeddings[, "UMAP_1"], col.name = "Azimuth_1")
  azimuthObject <- AddMetaData(object = azimuthObject, metadata = azimuthObject@reductions$ref.umap@cell.embeddings[, "UMAP_2"], col.name = "Azimuth_2")
  
  dat <- azimuthObject[[]]

  sampleID <- azimuthObject@project.name
  
  cols2keep <- grep(pattern = "predicted.celltype.", x = colnames(dat), value = TRUE)
  cols2keep <- c("cellBarcode", "mapping.score", "Azimuth_1", "Azimuth_2", cols2keep)
  dat <- dat[, cols2keep]
  
  write.csv(x = dat, file = paste0(sampleID, "_azimuth_mapping_results.csv"), row.names = FALSE, quote = FALSE)
  saveRDS(object = azimuthObject, file = paste0(sampleID, "_AzimuthObject.rds"))
}

main <- function(){
  
  args <- parseArguments()

  seurat <- loadData(argList = args)

  message("Running Azimuth")
  seurat <- RunAzimuth(query = seurat, reference = args$reference)

  seurat@project.name <- args$project

  writeResults(azimuthObject = seurat)
}

main()