#!/usr/bin/env Rscript

#This is a script used to generate plots for Downstream dashboards and reports.
#The script takes in a Seurat object and Azimuth object and visualizes the data.
#The script only has one required argument the project name. Eg --projectName exampleSample
#This argument is used for naming output files and directories.
#It is also possible to pass the path of the Seurat Object directly to the script. Eg --seuratPath ~/path/to/seuratObject.rds
#You can also pass the path to the Azimuth Object directly as well. Eg --azimuthPath ~/path/to/azimuthObject.rds
#If the path to either of these objects is not given as a command line argument the script will assume the objects
#are in the current working directory. The assumed file names are Seurat_object.rds and Azimuth_object.rds

library(Seurat)
library(ggplot2)
library(patchwork)
library(glue)
library(magrittr)
library(argparse)


parseArguments <- function(){
  
  argParser <- ArgumentParser()
  
  argParser$add_argument("--projectName", action = "store", required = TRUE)
  argParser$add_argument("--seuratPath", action = "store", default = FALSE)
  argParser$add_argument("--sampleLevels", action = "store", default = 4)
  argParser$add_argument("--plotOrigSample", action = "store_true")
  argParser$add_argument("--plotRtBc", action = "store_true")
  argParser$add_argument("--qcPlots", action = "store_true")
  
  argList <- argParser$parse_args()
  
  return(argList)
  
}

#All of the functions below produce and save plots used for downstream dashboards/reports or general qc visualizations.
#The functions have the same general structure they have a projName argument which is used for naming outputs.
#The functions also take in the data object to be visualized. The functions which visualize data in a seurat object
#have an argument named seuratObj and functions which take in data from an azimuth object have an argument names
#azimuthComparisonObj. In general the function names should indicate the type of plot generated and which data is being visualized.

###RNA Sample Dashboard Functions###

azimuthCellTypeCountBarGraph <- function(azimuthComparisonObj, projName){
  
  dfForm <- as.data.frame(table(azimuthComparisonObj$predicted.celltype.l1))
  colnames(dfForm) <- c("CellType","Count")
  dfForm$Proportion = (dfForm$Count / sum(dfForm$Count)) * 100
  countsBarGraph <- ggplot(data=dfForm, aes(y=CellType, x=Count, fill=CellType)) +
    geom_bar(stat="identity", position=position_dodge())+ theme_minimal() + 
    theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylab("Cell Type") + xlab("Number of cells")
  ggsave(filename = glue("{projName}/Azimuth_CellTypeCountsL1.png"), plot = countsBarGraph, width=7, height=7, dpi=300, bg="white") 
    
  proportionsBarGraph <- ggplot(data=dfForm, aes(y=CellType, x=Proportion, fill=CellType)) +
    geom_bar(stat="identity", position=position_dodge()) + xlim(0,100) +
    theme_minimal() + theme(axis.text.x=element_text(angle = -45, hjust = 0)) + ylab("Cell Type") + xlab("Percent of cells")
  ggsave(filename = glue("{projName}/Azimuth_CellTypeProportionsL1.png"), plot = proportionsBarGraph, width=7, height=7, dpi=300, bg="white")
  
  dfForm <- as.data.frame(table(azimuthComparisonObj$predicted.celltype.l2))
  colnames(dfForm) <- c("CellType","Count")
  dfForm$Proportion = (dfForm$Count / sum(dfForm$Count)) * 100

  countsBarGraph <- ggplot(data=dfForm, aes(y=CellType, x=Count, fill=CellType)) +
    geom_bar(stat="identity", position=position_dodge())+ theme_minimal() + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylab("Cell Type") + xlab("Number of cells")
  ggsave(filename = glue("{projName}/Azimuth_CellTypeCountsL2.png"), plot = countsBarGraph, width=7, height=7, dpi=300, bg="white") 
  proportionsBarGraph <- ggplot(data=dfForm, aes(y=CellType, x=Proportion, fill=CellType)) +
    geom_bar(stat="identity", position=position_dodge())+ xlim(0,100) +
    theme_minimal() + theme(axis.text.x=element_text(angle = -45, hjust = 0)) + ylab("Cell Type") + xlab("Percent of cells")
  ggsave(filename = glue("{projName}/Azimuth_CellTypeProportionsL2.png"), plot = proportionsBarGraph, width=7, height=7, dpi=300, bg="white")  
}

azimuthCellTypesAzimuthUMAP <- function(azimuthComparisonObj, projName){
  
  p1 <- DimPlot(object = azimuthComparisonObj, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3, reduction = "proj.umap") + NoLegend()
  ggsave( filename = glue("{projName}/Azimuth_ClusteringUmapL1.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
  p1 <- DimPlot(object = azimuthComparisonObj, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, reduction = "proj.umap") + NoLegend()
  ggsave( filename = glue("{projName}/Azimuth_ClusteringUmapL2.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
  p1 <- DimPlot(object = azimuthComparisonObj, group.by = "predicted.celltype.l3", label = TRUE, label.size = 3, reduction = "proj.umap") + NoLegend()
  ggsave( filename = glue("{projName}/Azimuth_ClusteringUmapL3.png"), plot = p1, width=7, height=7, dpi=300, bg="white")  
}

azimuthScoresViolin <- function(azimuthComparisonObj, projName){
  
  p2 <- VlnPlot(object = azimuthComparisonObj, features = "predicted.celltype.l2.score", group.by = "predicted.celltype.l2", pt.size=0) & 
    theme(axis.text.x = element_text(angle = 90)) & ggtitle("Cell Type Prediction Scores") & NoLegend() & ylim(0,1)
  ggsave( filename = glue("{projName}/Azimuth_PredictionScoreL2.png"), plot = p2, width=7, height=5, dpi=300, bg="white") 
  p3 <- VlnPlot(object = azimuthComparisonObj, features = "mapping.score", group.by = "predicted.celltype.l2", pt.size=0) & 
    theme(axis.text.x = element_text(angle = 90)) & ggtitle("Cell Type Mapping Scores") & NoLegend() & ylim(0,1)
  ggsave( filename = glue("{projName}/Azimuth_MappingScoresL2.png"), plot = p3, width=7, height=5, dpi=300, bg="white")
  p2 <- VlnPlot(object = azimuthComparisonObj, features = "predicted.celltype.l1.score", group.by = "predicted.celltype.l1", pt.size=0) & 
    theme(axis.text.x = element_text(angle = 90)) & ggtitle("Cell Type Prediction Scores") & NoLegend() & ylim(0,1)
  ggsave( filename = glue("{projName}/Azimuth_PredictionScoreL1.png"), plot = p2, width=5, height=5, dpi=300, bg="white") 
  p3 <- VlnPlot(object = azimuthComparisonObj, features = "mapping.score", group.by = "predicted.celltype.l1", pt.size=0) & ylim(0,1)
  theme(axis.text.x = element_text(angle = 90)) & ggtitle("Cell Type Mapping Scores") & NoLegend() 
  ggsave( filename = glue("{projName}/Azimuth_MappingScoresL1.png"), plot = p3, width=5, height=5, dpi=300, bg="white")
}

azimuthScoreScatter <- function(azimuthComparisonObj, projName){
  
  #Can use Seurat::FeatureScatter  
  p4 <- ggplot(azimuthComparisonObj@meta.data, aes(x=predicted.celltype.l2.score, y=mapping.score, colour = predicted.celltype.l2)) + 
    geom_point(size=2, shape=23) + xlim(0,1) + ylim(0,1)
  ggsave( filename = glue("{projName}/Azimuth_ScoreScatterL2.png"), plot = p4, width=7, height=5, dpi=300, bg="white") 
  p5 <- ggplot(azimuthComparisonObj@meta.data, aes(x=predicted.celltype.l1.score, y=mapping.score, colour = predicted.celltype.l1)) + 
    geom_point(size=2, shape=23) + xlim(0,1) + ylim(0,1)
  ggsave( filename = glue("{projName}/Azimuth_ScoreScatterL1.png"), plot = p5, width=7, height=5, dpi=300, bg="white")
}

seuratElbowPlot <- function(seuratObj, projName){
  
  p1 <- ElbowPlot(object = seuratObj)
  ggsave( filename = glue("{projName}/Seurat_elbow.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
}

seuratClusterUmap <- function(seuratObj, projName){
  
  p1 <- DimPlot(seuratObj, reduction = "umap") 
  ggsave( filename = glue("{projName}/Seurat_clusterUmap.png"), plot = p1, width=7, height=5, dpi=300, bg="white")  
}

violinQCPlots <- function(seuratObj, projName){

  p1 <- VlnPlot(seuratObj, features = "nCount_RNA", ncol = 1, pt.size=0) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  ggsave( filename = glue("{projName}/Seurat_clsuterCountsVln.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
  p1 <- VlnPlot(seuratObj, features = "nFeature_RNA", ncol = 1, pt.size=0) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  ggsave( filename = glue("{projName}/Seurat_clsuterFeaturesVln.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
  p1 <- VlnPlot(seuratObj, features = "predicted.celltype.l2.score", ncol = 1, pt.size=0) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  ggsave( filename = glue("{projName}/Seurat_clsuterAzimuthPredictedScoreL2Vln.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
}

azimuthCellTypeonDeNovoUMAP <- function(seuratObj, projName){
  
  Idents(seuratObj) <- seuratObj$predicted.celltype.l2
  p1 <- DimPlot(seuratObj, reduction = "umap")
  ggsave( filename = glue("{projName}/Seurat_clusterAzimuthCellTypeL2.png"), plot = p1, width=10, height=5, dpi=300, bg="white")
  Idents(seuratObj) <- seuratObj$predicted.celltype.l1
  p1 <- DimPlot(seuratObj, reduction = "umap")
  ggsave( filename = glue("{projName}/Seurat_clusterAzimuthCellTypeL1.png"), plot = p1, width=9, height=5, dpi=300, bg="white")
  Idents(seuratObj) <- seuratObj$seurat_clusters
}


seuratPlotRTBarcode <- function(seuratObj, projName){
  
  seuratObj[["RT"]] <- sub(".*_", "", Cells(seuratObj))
  seuratObj[["RT"]] <- sapply(seuratObj$RT, function(x) {substr(x,10,19)})
  p1                <- DimPlot(seuratObj, group.by="RT", reduction = "umap") + ggtitle(paste(basename(projName), "RT"))
  ggsave(filename = glue("{projName}/Seurat_RT_Umap.png"), plot = p1, width=20, height=15, dpi=300, bg="white")
  # count RT per cluster
  dfRT            <- as.data.frame(table(seuratObj$RT, seuratObj$seurat_clusters))
  colnames(dfRT)  <- c("RT","seurat_clusters", "Count")
  dfRT$Proportion <- (dfRT$Count / sum(dfRT$Count)) * 100
  rtBarGraph <- ggplot(data=dfRT, aes(x=seurat_clusters, y=Proportion, fill=RT)) +
    geom_bar(stat="identity", position="fill") +
    theme_minimal() + theme(axis.text.x=element_text(angle = 45, hjust = 0),
                            axis.text = element_text(size=14, color="black"),
                            axis.title = element_text(size=14, color="black")
    ) + 
    xlab("Cluster") + ylab("Percent")
  ggsave(filename = glue("{projName}/Seurat_rtClusterProportions.png"), plot = rtBarGraph, width=20, height=15, dpi=300, bg="white")
}

###Mapping Metric qcPlots Functions###

heatmapQC <- function(seuratObj, projName){
  
  bcVec <- c("RT", "PCR", "Ligation")
  
  metricVec <- c("reads", "exonReads", "passingReads", "mitoReads", "genes", "mappedReads", "intronReads", "uniquePassingReads", "umis", "geneReads", "antisenseReads", "Saturation")
  
  for(bc in bcVec){
    
    for(met in metricVec){
      outFile <- paste0(projName, "/qcPlots/", bc, "/heatmap.", met, ".png")
      bcAlias <- paste0(bc, "_alias")
      plotDat <- seuratObj@meta.data[, c(met, bcAlias)]
      plotDat["Column"] <- stringr::str_extract(string = plotDat[[bcAlias]], pattern = "^[:digit:]{1,2}")
      plotDat[["Column"]] <- as.numeric(plotDat[["Column"]])
      plotDat[["Column"]] <- factor(x = plotDat[["Column"]])
      plotDat["Row"] <- stringr::str_extract(string = plotDat[[bcAlias]], pattern = "[:alpha:]$")
      plotDat[["Row"]] <- factor(plotDat[["Row"]])
      plotDat[["Row"]] <- factor(x = plotDat[["Row"]], levels = rev(levels(plotDat[["Row"]])))
      plt <- ggplot(data = plotDat, aes(x = Column, y = Row, fill = .data[[met]])) + geom_tile()
      ggsave(filename = outFile, plot = plt) 
    } 
  }
}

#Function used to visualize the proportion of each barcode by seurat clusters.
#Can group cells by variable other than cluster.
barcodeStackedBar <- function(seuratObj, projName){
  
  group <- "seurat_clusters"
  
  bcVec <- c("RT", "PCR", "Ligation")
  
  for(bc in bcVec){
    
    outFile <- paste0(projName, "/qcPlots/", bc, "/stackedBar.", bc, ".", group, ".png")
    
    if(length(unique(seuratObj@meta.data[[bc]])) > 10){
      plt <- ggplot(data = seuratObj@meta.data[, c(group, bc)], aes(x = .data[[group]], fill = .data[[bc]])) + geom_bar(position = "fill") + theme_minimal() + theme(legend.position = "none")
    } else {
      plt <- ggplot(data = seuratObj@meta.data[, c(group, bc)], aes(x = .data[[group]], fill = .data[[bc]])) + geom_bar(position = "fill") + theme_minimal() 
    }
    ggsave(filename = outFile, plot = plt, bg = "white") 
  }
}

seuratPlotMappingMetricsUMAP <- function(seuratObj, projName){
  
  metricVec <- c("reads", "exonReads", "passingReads", "mitoReads", "genes", "mappedReads", "intronReads", "uniquePassingReads", "umis", "geneReads", "antisenseReads", "Saturation")
  
  for(met in metricVec){
    outFile <- paste0(projName, "/qcPlots/umap.", met, ".png")
    umap <- FeaturePlot(object = seuratObj, features = met)
    ggsave(filename = outFile, plot = umap, device = "png") 
  }
}

seuratPlotMappingMetricsViolin <- function(seuratObj, projName){
  grp <- "seurat_clusters"
  
  metricVec <- c("reads", "exonReads", "passingReads", "mitoReads", "genes", "mappedReads", "intronReads", "uniquePassingReads", "umis", "geneReads", "antisenseReads", "Saturation")
  
  for(met in metricVec){
    outFile <- paste0(projName, "/qcPlots/violin.", grp, ".", met, ".png")
    vln <- VlnPlot(object = seurat, features = met, pt.size = 0, group.by = grp)
    ggsave(filename = outFile, plot = vln, device = "png") 
  }
}

seuratPlotBarcodeUMAP <- function(seuratObj, projName){
 
  bcVec <- c("RT", "PCR", "Ligation")
  
  for(bc in bcVec){
    
    outFile <- paste0(projName, "/qcPlots/", bc, "/umap.", bc, ".barcodes.png")
    
    if(length(unique(seuratObj@meta.data[[bc]])) > 10){
      umap <- DimPlot(object = seuratObj, group.by = bc) + theme(legend.position = "none")
    } else {
      umap <- DimPlot(object = seuratObj, group.by = bc) 
    }
    ggsave(filename = outFile, plot = umap, device = "png", width = 20, height = 15, dpi = 300, bg = "white") 
  }  
}

markerFeaturePlots <- function(seuratObj, projName){
  p1 <- FeaturePlot(seuratObj, features = c("CD4"))
  ggsave( filename = glue("{projName}/CD4_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
  p1 <- FeaturePlot(seuratObj, features = c("CD3E"))
  ggsave( filename = glue("{projName}/CD3E_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
  p1 <- FeaturePlot(seuratObj, features = c("CD8A"))
  ggsave( filename = glue("{projName}/CD8A_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
  p1 <- FeaturePlot(seuratObj, features = c("CD19"))
  ggsave( filename = glue("{projName}/CD19_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
  p1 <- FeaturePlot(seuratObj, features = c("NCAM1"))
  ggsave( filename = glue("{projName}/NCAM1_CD56_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
}

markerDotPlots <- function(seuratObj, projName){
  Idents(seuratObj) <- seuratObj$predicted.celltype.l1
  markers.to.plot <- c("CD4", "CD3E", "CD19", "CD8A", "NCAM1")
  p1 <- DotPlot(seuratObj, features = markers.to.plot, dot.scale = 8, assay = "RNA") +
    RotatedAxis()
  ggsave( filename = glue("{projName}/Azimuth_Level1PbmcMarkersDotPlot.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
  Idents(seuratObj) <- seuratObj$predicted.celltype.l2
  p1 <- DotPlot(seuratObj, features = markers.to.plot, dot.scale = 8, assay = "RNA") +
    RotatedAxis()
  ggsave( filename = glue("{projName}/Azimuth_Level2PbmcMarkersDotPlot.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
}

markerViolinPlots <- function(seuratObj, projName){
  markers.to.plot <- c("CD4", "CD3E", "CD19", "CD8A", "NCAM1")
  plots <- VlnPlot(seuratObj, features = markers.to.plot, combine = FALSE, pt.size=0)
  p1 <- wrap_plots(plots = plots, ncol = 1)
  ggsave( filename = glue("{projName}/Azimuth_Level2PbmcMarkersVlnPlot.png"), plot = p1, width=10, height=20, dpi=300, bg="white")
  Idents(seuratObj) <- seuratObj$predicted.celltype.l1
  plots <- VlnPlot(seuratObj, features = markers.to.plot, combine = FALSE, pt.size=0)
  p1 <- wrap_plots(plots = plots, ncol = 1)
  ggsave( filename = glue("{projName}/Azimuth_Level1PbmcMarkersVlnPlot.png"), plot = p1, width=10, height=20, dpi=300, bg="white")
}

featureQCPlots <- function(seuratObj, projName){  
  p1 <- FeaturePlot(seuratObj, features="nCount_RNA", reduction = "umap") + ggtitle("Cluster nCount_RNA")
  ggsave( filename = glue("{projName}/Seurat_nCount_Umap.png"), plot = p1, width=6, height=5, dpi=300, bg="white")
  p1 <- FeaturePlot(seuratObj, features="nFeature_RNA", reduction = "umap") + ggtitle("Cluster ")
  ggsave( filename = glue("{projName}/Seurat_nFeature_Umap.png"), plot = p1, width=6, height=5, dpi=300, bg="white")  
}

plotScrubletDat <- function(seuratObj, projName){
  outPath <- paste0(projName)  
  plt <- ggplot(seuratObj@meta.data, aes(x = predicted_doublet, y = doublet_score, color = predicted_doublet, fill = predicted_doublet)) + geom_boxplot(alpha = 0.2) + theme_bw()
  ggsave(path = outPath, filename = "scrublet_doublet_scores_boxplot.png", plot = plt, width=5, height=5, dpi=300, bg="white")
}

seuratVariableFeatures <- function(seuratObj, projName){
  if("SCT" %in% Assays(seuratObj)){
    top10 <- head(VariableFeatures(seuratObj, assay = "SCT"), 10)
    p1    <- VariableFeaturePlot(seuratObj, assay = "SCT")
    p2    <- LabelPoints(plot = p1, points = top10, repel = TRUE)
  } else{  
    top10 <- head(VariableFeatures(seuratObj, assay = "RNA", selection.method = "vst"), 10)
    p1    <- VariableFeaturePlot(seuratObj, selection.method = "vst", assay = "RNA")
    p2    <- LabelPoints(plot = p1, points = top10, repel = TRUE)
  }
  ggsave( filename = glue("{projName}/Seurat_varFeatures.png"), plot = p2, width=7, height=5, dpi=300, bg="white")
}

seuratPlotOriginalSample <- function(seuratObj, sampleLevels, projName){
  
  idDf <- matrix(unlist(strsplit( rownames(seuratObj[[]]) ,split = "_")) ,ncol=sampleLevels,byrow=T)[,1:sampleLevels] %>% as.data.frame
  cols                  <- colnames(idDf)[1:(ncol(idDf)-1)]
  if(length(cols) > 1){
    seuratObj[["sample"]] <- apply( idDf[ , cols ] , 1 , paste , collapse = "_" )
  }else{
    seuratObj[["sample"]] <- idDf[,cols]
  }
  p1 <- DimPlot(seuratObj, group.by = "sample")
  ggsave(filename = "Seurat_origSampleUmap.png", plot = p1, width=7, height=5, dpi=300, bg="white")
  p1 <- DimPlot(seuratObj, reduction = "pca", group.by = "sample")
  ggsave(filename = "Seurat_origSamplePCA.png", plot = p1, width=7, height=5, dpi=300, bg="white")
  ## TO DO add a plot with the number and proportion of sample Ids per cluster
  
}

azimuthQCFeaturePlot <- function(azimuthComparisonObj, projName){
  p1 <- FeaturePlot(azimuthComparisonObj, features = "nCount_RNA", reduction = "proj.umap") + ggtitle("nCount_RNA on Azimuth Umap") + theme(plot.title = element_text(hjust = 0.5))
  ggsave( filename = glue("{projName}/Azimuth_ClusteringUmap_nCounts.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
  
}

##########
#Script workflow starts here.
##########

params <- parseArguments()

#Loads data.
if(isFALSE(params$seuratPath)){
  seurat <- readRDS(file = paste0(params$projectName, "/Seurat_object.rds"))
} else {
  seurat <- readRDS(file = params$seuratPath)
}

if("refAssay" %in% Assays(object = seurat)){

  azimuthQCFeaturePlot(azimuthComparisonObj = seurat, projName = params$projectName)
  markerViolinPlots(seuratObj = seurat, projName = params$projectName)
  markerDotPlots(seuratObj = seurat, projName = params$projectName)
  azimuthCellTypeonDeNovoUMAP(seuratObj = seurat, projName = params$projectName)
  azimuthScoreScatter(azimuthComparisonObj = seurat, projName = params$projectName)
  azimuthScoresViolin(azimuthComparisonObj = seurat, projName = params$projectName)
  azimuthCellTypesAzimuthUMAP(azimuthComparisonObj = seurat, projName = params$projectName)
  azimuthCellTypeCountBarGraph(azimuthComparisonObj = seurat, projName = params$projectName)
  violinQCPlots(seuratObj = seurat, projName = params$projectName)

} else {
   message("No Azimuth Data to Plot.")
}

#Checks if RT barcode plots should be generated.
if(params$plotRtBc){
  seuratPlotRTBarcode(seuratObj = seurat, projName = params$projectName)  
}

#Checks if UMAP and PCA plot colored by original sample name should be generated.
if(params$plotOrigSample){
  seuratPlotOriginalSample(seuratObj = seurat, sampleLevels = params$sampleLevels, projName = params$projectName)
}

#Checks if mapping metric qc plots should be generated.
if(params$qcPlots){
  
  if(!dir.exists(paste0(params$projectName, "/qcPlots"))){
    dir.create(path = paste0(params$projectName, "/qcPlots"))
  }
  
  if(!dir.exists(paste0(params$projectName, "/qcPlots/Ligation"))){
    dir.create(path = paste0(params$projectName, "/qcPlots/Ligation")) 
  }

  if(!dir.exists(paste0(params$projectName, "/qcPlots/PCR"))){
    dir.create(path = paste0(params$projectName, "/qcPlots/PCR"))
  }

  if(!dir.exists(paste0(params$projectName, "/qcPlots/RT"))){
    dir.create(path = paste0(params$projectName, "/qcPlots/RT"))
  }
  
  barcodeStackedBar(seuratObj = seurat, projName = params$projectName)
  seuratPlotMappingMetricsUMAP(seuratObj = seurat, projName = params$projectName)
  seuratPlotMappingMetricsViolin(seuratObj = seurat, projName = params$projectName)
  seuratPlotBarcodeUMAP(seuratObj = seurat, projName = params$projectName)
  heatmapQC(seuratObj = seurat, projName = params$projectName)
}

###Misc. Plots###

markerFeaturePlots(seuratObj = seurat, projName = params$projectName)
featureQCPlots(seuratObj = seurat, projName = params$projectName)
plotScrubletDat(seuratObj = seurat, projName = params$projectName)
seuratVariableFeatures(seuratObj = seurat, projName = params$projectName)
seuratElbowPlot(seuratObj = seurat, projName = params$projectName)
seuratClusterUmap(seuratObj = seurat, projName = params$projectName)
