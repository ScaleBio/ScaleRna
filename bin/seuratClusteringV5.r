#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(data.table))
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
      message(paste0("Cells in Dataset: ", nrow(cells)))
      cells <- cells[[1]]
    }
    cells <- paste0(cells, "_", sampleId)

    #Directory name where BPCells will write mtx.
    bp_dir <- paste0(sampleId, "_bpcells")

    mtx <- import_matrix_market(
      mtx_path = mtx,
      row_names = genes,
      col_names = cells,
      tmpdir = file.path(bp_dir, "tmp"), #Set tmpdir used by bpcells to ./sampleId_bpcells/tmp
      outdir = file.path(bp_dir, "out") #Set outdir used by bpcells to ./sampleId_bpcells/out. This is the directory with the final matrix.
    )
    
    #Modify dir value in mtx to relative path.
    #By default BPCells uses absolute paths but this causes issues when trying to move seurat objects.
    mtx@dir <- file.path(bp_dir, "out")

    return(mtx)
}

#This function reads in the allCells csv file from the samples directory.
readAllCellsData <- function(filePath){
  
  message("Reading allCells.csv")
  
  dat <- fread(file = filePath)
  dat <- setnames(dat, old = 1, new = "cellBarcode")
  cellBarcodes <- paste0(dat[["cellBarcode"]], "_", dat[["sample"]])
  dat[["cellBarcode"]] <- cellBarcodes
  dat <- dat[dat$pass == TRUE, ]
  
  return(dat)
}

writeReportResults <- function(seuratObject, sketch){

  message("Writing Results")
  
  outDat <- seuratObject[[]]

  if(sketch){
    outDat[["Seurat_1"]] <- seuratObject@reductions$full.umap@cell.embeddings[, "fullumap_1"]
    outDat[["Seurat_2"]] <- seuratObject@reductions$full.umap@cell.embeddings[, "fullumap_2"]
  } else {
    outDat[["Seurat_1"]] <- seuratObject@reductions$umap@cell.embeddings[, "umap_1"]
    outDat[["Seurat_2"]] <- seuratObject@reductions$umap@cell.embeddings[, "umap_2"]
  }
  
  outPrefix <- seuratObject@project.name
  
  cols2keep <- c("Seurat_1", "Seurat_2", "cluster_full", "cluster_full.score", "leverage.score", "cellBarcode", "seurat_clusters")
  cols2keep <- cols2keep[which( cols2keep %in% colnames(outDat) )]

  outDat <- outDat[, cols2keep]

  fwrite(x = outDat, file = paste0(outPrefix, "_seurat_clustering_results.csv"))
  saveRDS(object = seuratObject, file = paste0(outPrefix, "_SeuratObject.rds"))
}

standardWorkflow <- function(argList){

  countMat <- readCountMatrix(mtxDir = argList$matrixDir, mtxName = argList$star_matrix)
  metaDat <- readAllCellsData(filePath = argList$allCells)

  sketch <- ifelse(nrow(metaDat) >= 50000, TRUE, FALSE)

  seurat <- CreateSeuratObject(
    counts = countMat,
    meta.data = metaDat,
    project = argList$project
  )

  seurat <- NormalizeData(object = seurat)
  seurat <- FindVariableFeatures(object = seurat)
  
  if(sketch){
    seurat <- SketchData(
      object = seurat,
      ncells = (.2 * ncol(seurat)),
      method = "LeverageScore",
      sketched.assay = "sketch",
    )
    DefaultAssay(seurat) <- "sketch"
  }

  seurat <- ScaleData(object = seurat)
  seurat <- RunPCA(object = seurat)
  seurat <- FindNeighbors(object = seurat)
  seurat <- FindClusters(object = seurat)
  #return.model = TRUE is necesarry for projecting umap results during sketch-based analysis.
  seurat <- RunUMAP(object = seurat, dims = 1:30, return.model = TRUE)

  #Projected PCA results are stored using value from full.reduction parameter.
  #Projected UMAP results are stored with the name full.umap by default.
  #return.model = TRUE must be set in RunUMAP command before this.
  if(sketch){
    seurat <- ProjectData(
      object = seurat,
      assay = "RNA",
      full.reduction = "full.pca",
      sketched.assay = "sketch",
      sketched.reduction = "pca",
      umap.model = "umap",
      dims = 1:50,
      refdata = list(cluster_full = "seurat_clusters")
    )
    DefaultAssay(seurat) <- "RNA"
  }

  clusterMarkers <- FindAllMarkers(object = seurat, logfc.threshold = 2, verbose = FALSE)
  fwrite(x = clusterMarkers, file = paste0(argList$project, "_seurat_cluster_markers.csv"))

  writeReportResults(seuratObject = seurat, sketch = sketch)
}

classicScaleWorkflow <- function(argList){

  countMat <- readCountMatrix(mtxDir = argList$matrixDir, mtxName = argList$star_matrix)
  metaDat <- readAllCellsData(filePath = argList$allCells)

  seurat <- CreateSeuratObject(
    counts = countMat,
    meta.data = metaDat,
    project = argList$project,
    min.features = 200
  )

  seurat <- NormalizeData(object = seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat <- FindVariableFeatures(object = seurat, selection.method = "vst", nfeatures = 2000)
  seurat <- ScaleData(object = seurat)
  seurat <- RunPCA(object = seurat, npcs = 50)
  seurat <- FindNeighbors(object = seurat)
  seurat <- FindClusters(object = seurat, resolution = 0.3)
  seurat <- RunUMAP(object = seurat, dims = 1:15, return.model = TRUE)

  writeReportResults(seuratObject = seurat, sketch = FALSE)
}

comparisonWorkflow <- function(argList){
  
  message("Comparison Workflow")

  matList <- lapply(argList$matrixDir, readCountMatrix, mtxName = "matrix.mtx")
  
  metaDat <- lapply(argList$allCells, readAllCellsData)
  metaDat <- rbindlist(metaDat)
  
  sketch <- ifelse(nrow(metaDat) >= 50000, TRUE, FALSE)
  
  seurat <- CreateSeuratObject(
    counts = matList,
    meta.data = metaDat,
    project = argList$project
  )
  
  seurat <- JoinLayers(object = seurat)
  seurat <- NormalizeData(object = seurat)
  seurat <- FindVariableFeatures(object = seurat)

  if(sketch){
    seurat <- SketchData(
      object = seurat,
      ncells = (.2 * ncol(seurat)),
      method = "LeverageScore",
      sketched.assay = "sketch",
    )
    DefaultAssay(seurat) <- "sketch"
  }
  seurat <- FindVariableFeatures(object = seurat)
  seurat <- ScaleData(seurat)
  seurat <- RunPCA(seurat)
  seurat <- FindNeighbors(seurat)
  seurat <- FindClusters(seurat)
  #return.model = TRUE is necesarry for projecting umap results during sketch-based analysis.
  seurat <- RunUMAP(seurat, dims = 1:30, return.model = TRUE)


  #Projected PCA results are stored using value from full.reduction parameter.
  #Projected UMAP results are stored with the name full.umap by default.
  #return.model = TRUE must be set in RunUMAP command before this.
  if(sketch){
    seurat <- ProjectData(
      object = seurat,
      assay = "RNA",
      full.reduction = "full.pca",
      sketched.assay = "sketch",
      sketched.reduction = "pca",
      umap.model = "umap",
      dims = 1:50,
      refdata = list(cluster_full = "seurat_clusters")
    )
    DefaultAssay(seurat) <- "RNA"
  }
  
  clusterMarkers <- FindAllMarkers(object = seurat, logfc.threshold = 2, verbose = FALSE)
  fwrite(x = clusterMarkers, file = paste0(argList$project, "_seurat_cluster_markers.csv"))
  
  writeReportResults(seuratObject = seurat, sketch = sketch)
}

main <- function(){
  
  argParser <- ArgumentParser(
    description = "This script performs Seurat based clustering and dimensionality reduction."
  )

  argParser$add_argument("--star_matrix", required = FALSE, default = "matrix.mtx")
  argParser$add_argument("--project", required = TRUE)

  subParsers <- argParser$add_subparsers(
    help = "This subcommand determines which version of the Seurat Workflow is performed."
  )

  parserStandard <- subParsers$add_parser("standard", help = "Performs standard Seurat analysis.")
  parserStandard$set_defaults(workflow = "standard")
  parserStandard$add_argument("--matrixDir", required = TRUE)
  parserStandard$add_argument("--allCells", required = TRUE)

  parserClassic <- subParsers$add_parser("classic", help = "Performs standard Seurat analysis.")
  parserClassic$set_defaults(workflow = "classic")
  parserClassic$add_argument("--matrixDir", required = TRUE)
  parserClassic$add_argument("--allCells", required = TRUE)

  parserComparison <- subParsers$add_parser("comparison", help = "Performs sketch based Seurat analysis comparing multiple samples.")
  parserComparison$set_defaults(workflow = "comparison")
  parserComparison$add_argument("--matrixDir", nargs = "+", required = TRUE)
  parserComparison$add_argument("--allCells", nargs = "+", required = TRUE)


  arguments <- argParser$parse_args()

  switch(
    arguments$workflow,
    standard = { standardWorkflow(argList = arguments) },
    classic = { classicScaleWorkflow(argList = arguments) },
    comparison = { comparisonWorkflow(argList = arguments)}
  )
}

main()
