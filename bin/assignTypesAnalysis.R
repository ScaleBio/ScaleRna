#!/usr/bin/env Rscript

#This is a script used to perform de novo Seurat analysis and/or Azimuth Cell Type classification of a given dataset.
#This script will output a Seurat object which only contains the given data if no options are given.
#The script will load, create, and filter a Seurat object with the given data. Default filtering parameters will be used.
#If a path to the output from Scrublet is given the script will also filter the data using this information.
#If given the script will also load meta data into the Seurat object. Such as information contained in the exampleSample.allCells.csv file.

#If the --seurat option is used the script will perform de novo Seurat Clustering and dimension reduction for the given data.
#The output from this de novo analysis will be included in the seurat object output by the script.

#If the --azimuth option is used the script will perform Azimuth Cell Type classification of the given data.
#The output from the azimuth cell type classification will be stored and output in the Seurat object produced by this script.
#The seurat object produced by this script can be used for further downstream analysis and visualization.
#The --seurat and --azimuth options can be used in combination.
#It is also possible to specify the normalization method (LogNormalize or SCT) used for the de novo seurat analysis.
#If the SCT option is used for normalization it is not possible to perform the azimuth analysis as well.

#Other options that can be used are --writeSeuratMarkers and --writeSeuratMeta.
#These options write marker genes for each seurat cluster and the meta data generated throughout the analysis to separate output files.
#A final option is --StarMtxFormat which helps specify the format of the data being loaded.

#The required arguments are --countDir and --projectName.
#The default parameters for the analysis performed in the script can be overwritten using the --paramYaml argument.

library(Seurat)
library(SeuratData)
library(Azimuth)
library(yaml)
library(patchwork)
library(ggplot2)
library(argparse)


parseArguments <- function(){
  
  argParser <- ArgumentParser()

  argParser$add_argument("--countDir", action = "store", required = TRUE)
  argParser$add_argument("--projectName", action = "store", required = TRUE)
  argParser$add_argument("--scrubletPath", action = "store", default = FALSE)
  argParser$add_argument("--paramYaml", action = "store", default = FALSE)
  argParser$add_argument("--metricsPath", action = "store", default = FALSE)
  argParser$add_argument("--clusterNormalization", choices = c("LogNormalize", "SCT"), action = "store")
  argParser$add_argument("--seurat", action = "store_true")
  argParser$add_argument("--azimuth", action = "store_true")
  argParser$add_argument("--writeSeuratMarkers", action = "store_true")
  argParser$add_argument("--writeSeuratMeta", action = "store_true")
  argParser$add_argument("--StarMtxFormat", action = "store_true")
  
  argList <- argParser$parse_args()
  
  if(argList$clusterNormalization == "SCT" & argList$azimuth == TRUE){
    stop("Can't perform SCT cluster normalization and Azimuth Classification at same time.\n  Please choose one or the other.")
  }
  
  return(argList)
}

#This function defines defaults parameters needed for running the various steps of the Seurat analysis.
#These parameters may need to be changed depending on the dataset being analyzed. In this case it is 
#possible to pass the path to a yaml file on the command line. Doing so will allow the user set these parameters
#to the desired values. If no yaml file is passed on the command line the defaults below are used for analysis.
parseParams <- function(paramYaml){
  
  if(isFALSE(paramYaml)){
    
    params <- list(
      
      maximumCounts = 1000000,
      maximumFeatures = 10000,
      minimumFeatures = 200,
      percentMito = 5,
      topNcells = 0,
      azimuthReference = "pbmcref",
      normalizationMethod = "LogNormalize",
      scaleFactor = 10000,
      selectionMethod = "vst",
      variableFeatures = 2000,
      numberPCs = 50,
      clusterResolution = 0.3,
      umapDims = 15
      
    )
    
  } else {
    
    params <- yaml.load_file(input = argList$paramYaml)
  }
  
  return(params)
}

#This function reads in the raw data output by the alignment pipeline.
seuratReadMatrix <- function(countDir, starFormat){

  stopifnot(dir.exists(countDir))
  
  if(starFormat){
    
    mtxPath <- paste0(countDir, "/matrix.mtx")
    featPath <- paste0(countDir, "/features.tsv")
    bcPath <- paste0(countDir, "/barcodes.tsv")
    
  } else{
    
    mtxPath <- paste0(countDir, "/matrix.mtx")
    featPath <- paste0(countDir, "/genes.tsv")
    bcPath <- paste0(countDir, "/barcodes.tsv")
    
  }
  
  mtx <- ReadMtx(
    mtx = mtxPath, 
    features = featPath, 
    cells = bcPath, 
    feature.column = 2
    )
  
  return(mtx)
}

#Function to load CellxGene matrix into a seurat object.
#If a path to the mapping metrics csv is given on the command line.
#These metrics will be added to the seurat object metadata.
seuratCreateObject <- function(countMatrix, projName, metPath){
  
  if(isFALSE(metPath)){
    
    seurat <- CreateSeuratObject(counts = countMatrix, project = projName)
    
  } else {

    stopifnot(file.exists(metPath))

    metDat <- read.csv(file = metPath)

    metDat <- metDat[metDat[["pass"]] == "True", ]

    rownames(metDat) <- metDat[[1]]
    
    seurat <- CreateSeuratObject(
      counts = countMatrix, 
      project = projName,
      meta.data = metDat
    )
    
  }
  
  seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^MT-")
  
  return(seurat)
  
}

#Function that adds the output from scrublet analysis to the metadata of a seurat object.
#This enables filtering the dataset based on the predicted doublets from scrublet.
seuratAddScrubletData <- function(seuratObject, scrubPath){

  stopifnot(file.exists(scrubPath))
  
  scrublet <- read.csv(file = scrubPath)
  seuratObject[["predicted_doublet"]] <- scrublet[["predicted_doublet"]]
  seuratObject[["doublet_score"]] <- scrublet[["doublet_score"]]
  
  return(seuratObject)
  
}

#Filters seurat object based on given scrublet output.
seuratFilterScrublet <- function(seuratObject){
  
  seuratObject <- subset(x = seuratObject, subset = predicted_doublet == "False")
  
  return(seuratObject)
  
}

#Filters seurat object based on various qc parameters. These parameters are defined in the parseParams function.
#The parameters can also be defined by passing a yaml file on the command line. This allows the user to change the
#parameters based on the dataset being analyzed.
seuratFilterMetrics <- function(seuratObject, maxCount, maxFeat, minFeat, ptMito){
  
  seuratObject <- subset(
    x = seuratObject,
    subset = nFeature_RNA > minFeat & nCount_RNA < maxCount & percent.mt < ptMito & nFeature_RNA < maxFeat
    )
  
  return(seuratObject)
  
}

#Filters seurat object. Chooses the top X cells based on nCount_RNA.
#X can be defined by user.
seuratChooseTopCells <- function(seuratObject, topNcells){
  
  topBCs <- names(sort(seuratObject$nCount_RNA, decreasing = TRUE))[seq_len(topNcells)]
  
  seuratObject <- seuratObject[, topBCs]
  
  return(seuratObject)
}

#This is a function that performs the normalization step of the Seurat analysis. 
#Log Normalization or SCTransform normalization can be performed.
#SCTransform normalization can not be performed at the same time as azimuth analysis.
seuratNormalize <- function(seuratObject, normMethod, scaleFactor, selectionMethod, variableFeats){
  
  if(normMethod == "SCT"){
    seuratObject <- SCTransform(object = seuratObject, conserve.memory = TRUE)
  } else if (normMethod == "LogNormalize") {
     seuratObject <- NormalizeData(
      object = seuratObject,
      normalization.method = normMethod,
      scale.factor = scaleFactor
     )

     seuratObject <- FindVariableFeatures(
      object = seuratObject,
      selection.method = selectionMethod,
      nfeatures = variableFeats
    )
    seuratObject <- ScaleData(object = seuratObject)
  }
  return(seuratObject)
}

#Function that performs dimension reduction of a given seurat object.
seuratPcaUmap <- function(seuratObject, umapDims){
  seuratObject <- RunPCA(object = seuratObject)
  seuratObject <- RunUMAP(object = seuratObject, dims = seq_len(umapDims))
  return(seuratObject)
}

#Function that clusters a given seurat object.
seuratNeighborsCluster <- function(seuratObject, clusterResolution){
  seuratObject <- FindNeighbors(seuratObject)
  seuratObject <- FindClusters(object = seuratObject, resolution = clusterResolution)
  return(seuratObject)
}

#This is a function that performs azimuth cell type classification analysis on a given seurat object.
#This code was mainly taken from the R script output by the Azimuth website when performing classification.
#This script is meant to be a means to reproduce the azimuth analysis without using the azimuth app.
azimuthMapping <- function(query, refName){

    reference <- LoadData(refName, type = "azimuth")
    
    query <- SCTransform(
        object = query,
        assay = "RNA",
        new.assay.name = "refAssay",
        residual.features = rownames(x = reference$map),
        reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
        method = 'glmGamPoi',
        ncells = 2000,
        n_genes = 2000,
        do.correct.umi = FALSE,
        do.scale = FALSE,
        do.center = TRUE,
        conserve.memory = TRUE
    )
    
    anchors <- FindTransferAnchors(
        reference = reference$map,
        query = query,
        k.filter = NA,
        reference.neighbors = "refdr.annoy.neighbors",
        reference.assay = "refAssay",
        query.assay = "refAssay",
        reference.reduction = "refDR",
        normalization.method = "SCT",
        features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
        dims = 1:50,
        n.trees = 20,
        mapping.score.k = 100
    )
    refdata <- lapply(
        X = c("celltype.l1", "celltype.l2", "celltype.l3"), 
        function(x) {reference$map[[x, drop = TRUE]]
    })   
    
    names(x = refdata) <- c("celltype.l1", "celltype.l2", "celltype.l3")
    
    query <- TransferData(
        reference = reference$map,
        query = query,
        dims = 1:50,
        anchorset = anchors,
        refdata = refdata,
        n.trees = 20,
        store.weights = TRUE
    )
    query <- IntegrateEmbeddings(
        anchorset = anchors,
        reference = reference$map,
        query = query,
        reductions = "pcaproject",
        reuse.weights.matrix = TRUE
    )
    query[["query_ref.nn"]] <- FindNeighbors(
        object = Embeddings(reference$map[["refDR"]]),
        query = Embeddings(query[["integrated_dr"]]),
        return.neighbor = TRUE,
        l2.norm = TRUE
    )
    query <- Azimuth:::NNTransform(
        object = query,
        meta.data = reference$map[[]]
    )
    query[["proj.umap"]] <- RunUMAP(
        object = query[["query_ref.nn"]],
        reduction.model = reference$map[["refUMAP"]],
        reduction.key = 'UMAP_'
    )
    query <- AddMetaData(
        object = query,
        metadata = MappingScore(anchors = anchors),
        col.name = "mapping.score"
    )
    rm(refdata, anchors, reference)
    return(query)
}

#Function that writes the Azimuth Cell Type classifications, prediction scores, and mapping scores to csv file.
writeAzimuthClusteringData <- function(azimuthComparisonObj, projName){
  
  outFile <- paste0(projName, "/Azimuth_ClusteringData.tsv")
  
  outFile <- file(outFile, "wb")
  
  cellTypeData <- azimuthComparisonObj@meta.data[,c("predicted.celltype.l1", "predicted.celltype.l1.score", "predicted.celltype.l2", "predicted.celltype.l2.score","nCount_RNA", 'mapping.score')]
  cellTypeData$CB <- rownames(cellTypeData)
  colnames(cellTypeData) <- c("L1_CellType","L1_CellType_Score", "L2_CellType", "L2_CellType_Score", "nCount_RNA", "Mapping_Score", "CB") # CB and mapping previously in the wrong order- check the report for issues fixing this
  cellTypeData <- format(x = cellTypeData, digits = 4)
  write.table(x = cellTypeData, file = outFile, row.names=FALSE, quote=FALSE,  sep='\t')
  
  close(outFile)
  
}

#Funtion that writes a table summarizing the Azimuth classifications.
#Writes number of cells classified as each cell type.
writeAzimuthCellTypeData <- function(azimuthComparisonObj, projName){
  
  outFile <- paste0(projName, "/Azimuth_CellTypes.csv")
  
  outFile <- file(outFile, "wb")
  
  dfForm <- as.data.frame(table(azimuthComparisonObj$predicted.celltype.l2))
  colnames(dfForm) <- c("CellType","Count")
  dfForm$Proportion = (dfForm$Count / sum(dfForm$Count)) * 100

  # L1 can be derived from L2
  write.csv(x = dfForm, file = outFile, row.names=FALSE,quote=FALSE)
  
  close(outFile)
  
}

# Calculates markers by seurat cluster and writes to CSV 
findAndWriteMarkers <- function(projName,seuratObj) {
  
  outFile <- paste0(projName, "/Seurat_ClusterMarkers.csv")
  
  all.markers <- FindAllMarkers(seuratObj)
  
  write.csv(x = all.markers, file = outFile, row.names=FALSE,quote=FALSE)
}

#Function to produce violin plots of qc metrics. These plots are used in a downstream dashboard.
#We generate these plots before and after filter the seurat object. Which requires us to plot in this script.
seuratPlotViolinQC <- function(seuratObj, outDir, outFile){
  
  featVec <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  
  plt <- VlnPlot(object = seuratObj, features = featVec, ncol = 3, pt.size = 0) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) & geom_boxplot(width = 0.1, fill = "white", outlier.colour = NA)
    
  outputFile = paste0("./", outDir, "/", outFile)
  
  ggsave(filename = outputFile, plot = plt, width = 7, height = 7, dpi = 300, bg = "white")
  
}

seuratPlotScatterQC <- function(seuratObj, outDir, outFile){
  
  outFile <- paste0("./", outDir, "/", outFile)
  
  plt1 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(title = NULL)
  
  plt2 <- FeatureScatter(object = seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + labs(title = NULL)
  
  ggsave(filename = outFile, plot = wrap_plots(plt1, plt2, ncol = 2), width = 10, height = 5, dpi = 300, bg = "white")
  
}

#Capture command line arguments passed by user.
argList <- parseArguments()

#Define default parameters used for seurat analysis.
params <- parseParams(argList$paramYaml)

if(! dir.exists(argList$projectName)){
  
  dir.create(argList$projectName)
}

seurat <- seuratCreateObject(
  countMatrix = seuratReadMatrix(countDir = argList$countDir, starFormat = argList$StarMtxFormat),
  projName = argList$projectName,
  metPath = argList$metricsPath
)

#Plot QC Metrics before filtering seurat object.
seuratPlotViolinQC(seuratObj = seurat, outDir = argList$projectName, outFile = "Seurat_QCboxplots.png" )

seuratPlotScatterQC(seuratObj = seurat, outDir = argList$projectName, outFile = "Seurat_QCscatter.png")

if(is.character(argList$scrubletPath)){
  
  stopifnot(file.exists(argList$scrubletPath))
  
  seurat <- seuratAddScrubletData(seuratObject = seurat, scrubPath = argList$scrubletPath)
  seurat <- seuratFilterScrublet(seuratObject = seurat)
}

#Filter Seurat object
seurat <- seuratFilterMetrics(
  seuratObject = seurat,
  maxCount = params$maximumCounts,
  maxFeat = params$maximumFeatures,
  minFeat = params$minimumFeatures,
  ptMito = params$percentMito
  )

#Plot QC Metrics after filtering
seuratPlotViolinQC(seuratObj = seurat, outDir = argList$projectName, outFile = "Seurat_QCboxplotsFilt.png")

seuratPlotScatterQC(seuratObj = seurat, outDir = argList$projectName, outFile = "Seurat_QCscatterFilt.png")

if(params$topNcells > 0){
  
  seurat <- seuratChooseTopCells(seuratObject = seurat, topNcells = params$topNcells)
  
}

#Perform Seurat analysis if --seurat option is used.
if(argList$seurat){
    
  seurat <- seuratNormalize(seuratObject = seurat, normMethod = argList$clusterNormalization, scaleFactor = params$scaleFactor, selectionMethod = params$selectionMethod, variableFeats = params$variableFeatures)
  seurat <- seuratPcaUmap(seuratObject = seurat, umapDims = params$umapDims)
  seurat <- seuratNeighborsCluster(seuratObject = seurat, clusterResolution = params$clusterResolution)
}

#Output marker genes for each cluster if --writeSeuratMarkers option is used.
if(argList$writeSeuratMarkers){
  findAndWriteMarkers(seuratObj = seurat, projName = argList$projectName)
}

#Perform azimuth reference mapping if --azimuth option is used.
if(argList$azimuth){

  featureMetaData <- seurat@assays$RNA[[]]

  seurat <- azimuthMapping(query = seurat, refName = params$azimuthReference)
  
  seurat@assays$RNA@meta.features <- featureMetaData

  seurat@assays$RNA@var.features <- rownames(featureMetaData)[featureMetaData[["vst.variable"]] == TRUE]
  
  rm(featureMetaData)

  writeAzimuthCellTypeData(azimuthComparisonObj = seurat,  projName = argList$projectName)

  writeAzimuthClusteringData(azimuthComparisonObj = seurat, projName = argList$projectName)

}

#Write meta data from seurat object if --writeSeuratMeta argument is used.
if(argList$writeSeuratMeta){
      write.csv(x = seurat@meta.data, file = paste0(argList$projectName, "/", "Seurat_object_meta.csv"))
}

#Save the final seurat object.
message("Saving Seurat Object")  
saveRDS(object = seurat, file = paste0(argList$projectName, "/Seurat_object.rds"))
