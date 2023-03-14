#!/usr/bin/env Rscript
args <- commandArgs()

help <- function(){
    cat("assignTypes.R : Run seurat and azimuth analysis\n")
    cat("Usage: \n")
    cat("--countDir     : Counts for which to perform analysis                              [required]
                                assumes file names are matrix.mtx, genes.tsv, barcodes.tsv\n")
    cat("--projName     : full path to save output dir                                      [required]\n")
    cat("--scrubletOut  : path to scrublet output for the matrix.mtx                        [required]\n")
    cat("--clustRes     : resolution to call clusters (larger value > 1 for more clusters)  [default = 0.4]\n")
    cat("--nCells       : top N cells to take by counts per cell                            [optional]\n")
    cat("--maxCounts    : max counts to allow                                               [default = max table]\n")
    cat("--minFeatures  : min features (genes) per cell to allow                            [default = 200]\n")
    cat("--MTpercent    : percent mito to filter cells on                                   [default = 5]\n")    
    cat("--plotOrigSample : If merged samples plot individual sample names (yes/no)         [default = no]\n")
    cat("--plotOrigSampleLevels : number of items in counts BCs that are seperated by _ from merge  [required with plotOrigSample]
                BC assumed to be last that will be trimmed\n")
    cat("--azimuthOnly  : Only run azimuth analysis (yes/no)                                [default = no]\n")
    cat("--StarMtxFormat: Mtx are either star (features.tsv) or droputils (genes.tsv) (yes/no) [default = no use genes.tsv]\n")
    cat("--plotRtBc     : Plot UMAP and stacked bar for cluster colored by RT (yes/no)         [default = no]\n")
    cat("--reference    : Azimuth ref string                                                   [default = pbmcref ]\n")
    cat("\n")
    q()
}

if(!is.na(charmatch("--help",args)) || !is.na(charmatch("-h",args))){
    help()
} else {
    countDir       <-sub('--countDir=', '',   args[grep('--countDir=', args)])
    projName       <-sub('--projName=',  '',  args[grep('--projName=', args)])
    scrubletOut    <-sub('--scrubletOut=','', args[grep('--scrubletOut=', args)])
    nCells         <-sub('--nCells=','',      args[grep('--nCells=', args)])
    maxCounts      <-sub('--maxCounts=','',   args[grep('--maxCounts=', args)])
    clustRes       <-sub('--clustRes=','',      args[grep('--clustRes=', args)])
    minFeatures    <-sub('--minFeatures=','', args[grep('--minFeatures=', args)])
    MTpercent      <-sub('--MTpercent=','',   args[grep('--MTpercent=', args)])
    plotOrigSample <-sub('--plotOrigSample=','', args[grep('--plotOrigSample=', args)])
    plotOrigSampleLevels <-sub('--plotOrigSampleLevels=','', args[grep('--plotOrigSampleLevels=', args)])
    azimuthOnly    <-sub('--azimuthOnly=','', args[grep('--azimuthOnly=', args)])
    StarMtxFormat  <-sub('--StarMtxFormat=','', args[grep('--StarMtxFormat=', args)])
    plotRtBc       <-sub('--plotRtBc=','', args[grep('--plotRtBc=', args)])
    reference       <-sub('--reference=','', args[grep('--reference=', args)])
}

packages <- c("Azimuth","dplyr","glue","ggplot2","HGNChelper","openxlsx","patchwork","scater","Seurat","SeuratData")
suppressMessages(invisible(lapply(packages, library, character.only = TRUE, quietly=TRUE)))

if( identical(MTpercent,character(0))){
    MTpercent <- 5
}else{
    MTpercent <- as.numeric(MTpercent)
}

if( identical(clustRes,character(0))){
    clustRes <- 0.3
    }else{
    clustRes <- as.numeric(clustRes)
}

if( identical(minFeatures,character(0))){
    minFeatures <- 200
}else{
    minFeatures <- as.numeric(minFeatures)
}

if( identical(maxCounts,character(0))){
    maxCounts <- 1000000
}else{
    maxCounts <- as.numeric(maxCounts)
}

if( identical(plotOrigSample,character(0))){
    plotOrigSample <- "no"
}

if(identical(plotOrigSampleLevels, character(0))){
    plotOrigSampleLevels <- 0
}else{
    plotOrigSampleLevels <- as.numeric(plotOrigSampleLevels)
}

if(identical(nCells, character(0))){
    nCells <- 0
}else{
    nCells <- as.numeric(nCells)
}

if(identical(azimuthOnly, character(0))){
    azimuthOnly <- "no"
}else{
    azimuthOnly <- azimuthOnly
}

if(identical(StarMtxFormat, character(0))){
    mtxVal      = "matrix.mtx"
    featuresVal = "genes.tsv"
    cellsVal    = "barcodes.tsv"
}else if(StarMtxFormat == "no"){
    mtxVal      = "matrix.mtx"
    featuresVal = "genes.tsv"
    cellsVal    = "barcodes.tsv"
}else if(StarMtxFormat == "yes"){
    mtxVal      = "matrix.mtx"
    featuresVal = "features.tsv"
    cellsVal    = "barcodes.tsv"
}

if( identical(plotRtBc,character(0))){
    plotRtBc <- "no"
}

if( identical(reference,character(0))){
    reference <- "pbmcref"
}

opts                      <- list()
opts$nCells               <- nCells
opts$clustRes             <- clustRes 
opts$MTpercent            <- MTpercent
opts$minFeatures          <- minFeatures
opts$maxCounts            <- maxCounts
opts$plotOrigSample       <- plotOrigSample
opts$plotOrigSampleLevels <- plotOrigSampleLevels # number strings seperated by underscores including BC assumed to be last that will be trimmed
opts$plotRtBc             <- plotRtBc
opts$reference            <- reference


main <- function(projName, countDir, scrubletOut, clustRes, nCells) {
    # Prevent writing of Rplots.pdf (which occurs by default)
    pdf(NULL)
    # Validate matrix dir existence
    stopifnot(file.exists(countDir))
    if(!(file.exists( file.path(projName) ))) {
        print(paste("mkdir", file.path(projName)))
        dir.create(file.path(projName),FALSE,TRUE)
    }
    # Using ReadMtx to be able to non gzipped tsv matrix representation
    mat_filt = ReadMtx(mtx            = glue("{countDir}/{mtxVal}"),
                       features       = glue("{countDir}/{featuresVal}"),
                       cells          = glue("{countDir}/{cellsVal}"),
                       feature.column = 2
    )
    seuratObj <- produceSeuratObj(mat_filt, projName)
    seuratObj <- filterSeuratObjDoublets(seuratObj, scrubletOut, projName)  
    if( nCells >0 ){
        print("filter nCells")
        seuratObj <- filterTopCellsByUmi(seuratObj, opts$nCells, projName)
    }
    # Compare to reference PBMC dataset using azimuth
    azimuthComparisonObj <- RunAzimuth(seuratObj, reference=opts$reference)
    produceAzimuthFigures(projName,azimuthComparisonObj)
    #saveRDS(azimuthComparisonObj, file = glue("{projName}/Azimuth_object.rds"))
    # Perform unsupervised clustering and cell type annotation
    # either replace seuratObj with azimuth meta
    if(azimuthOnly == "no"){
        stopifnot(rownames(seuratObj[[]]) == rownames(azimuthComparisonObj[[]]))
        seuratObj@meta.data <- azimuthComparisonObj@meta.data
        rm(azimuthComparisonObj)
        seuratObj <- normalizeAndCluster(seuratObj, projName, clustRes)
        makeSeuratClusterPlots(seuratObj, projName)
        makeMarkerPlotsPostClustering(projName=projName, seuratObj = seuratObj)
        saveRDS(seuratObj, file = glue("{projName}/Seurat_object.rds"))
        findAndWriteMarkers(projName,seuratObj) # add an option run time is long for this 
    }
}

produceAzimuthFigures <- function(projName,azimuthComparisonObj) {
    dfForm <- as.data.frame(table(azimuthComparisonObj$predicted.celltype.l1))
    colnames(dfForm) <- c("CellType","Count")
    dfForm$Proportion = (dfForm$Count / sum(dfForm$Count)) * 100
    countsBarGraph <- ggplot(data=dfForm, aes(y=CellType, x=Count, fill=CellType)) +
    geom_bar(stat="identity", position=position_dodge())+ theme_minimal() + 
    theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylab("Cell Type") + xlab("Number of cells")
    ggsave(glue("{projName}/Azimuth_CellTypeCountsL1.png"), plot = countsBarGraph, width=7, height=7, dpi=300, bg="white") 
    proportionsBarGraph <- ggplot(data=dfForm, aes(y=CellType, x=Proportion, fill=CellType)) +
    geom_bar(stat="identity", position=position_dodge()) + xlim(0,100) +
    theme_minimal() + theme(axis.text.x=element_text(angle = -45, hjust = 0)) + ylab("Cell Type") + xlab("Percent of cells")
    ggsave(glue("{projName}/Azimuth_CellTypeProportionsL1.png"), plot = proportionsBarGraph, width=7, height=7, dpi=300, bg="white")
    dfForm <- as.data.frame(table(azimuthComparisonObj$predicted.celltype.l2))
    colnames(dfForm) <- c("CellType","Count")
    dfForm$Proportion = (dfForm$Count / sum(dfForm$Count)) * 100
    write.csv(dfForm, glue("{projName}/Azimuth_CellTypes.csv"), row.names=FALSE,quote=FALSE) # L1 can be derived from L2
    countsBarGraph <- ggplot(data=dfForm, aes(y=CellType, x=Count, fill=CellType)) +
    geom_bar(stat="identity", position=position_dodge())+ theme_minimal() + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + ylab("Cell Type") + xlab("Number of cells")
    ggsave(glue("{projName}/Azimuth_CellTypeCountsL2.png"), plot = countsBarGraph, width=7, height=7, dpi=300, bg="white") 
    proportionsBarGraph <- ggplot(data=dfForm, aes(y=CellType, x=Proportion, fill=CellType)) +
    geom_bar(stat="identity", position=position_dodge())+ xlim(0,100) +
    theme_minimal() + theme(axis.text.x=element_text(angle = -45, hjust = 0)) + ylab("Cell Type") + xlab("Percent of cells")
    ggsave(glue("{projName}/Azimuth_CellTypeProportionsL2.png"), plot = proportionsBarGraph, width=7, height=7, dpi=300, bg="white")
    p1 <- DimPlot(azimuthComparisonObj, group.by = "predicted.celltype.l1", label = TRUE, label.size = 3) + NoLegend()
    ggsave(glue("{projName}/Azimuth_ClusteringUmapL1.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
    p1 <- DimPlot(azimuthComparisonObj, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) + NoLegend()
    ggsave(glue("{projName}/Azimuth_ClusteringUmapL2.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
    p1 <- DimPlot(azimuthComparisonObj, group.by = "predicted.celltype.l3", label = TRUE, label.size = 3) + NoLegend()
    ggsave(glue("{projName}/Azimuth_ClusteringUmapL3.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
    p2 <- VlnPlot(azimuthComparisonObj, features = "predicted.celltype.l2.score", group.by = "predicted.celltype.l2", pt.size=0) & 
    theme(axis.text.x = element_text(angle = 90)) & ggtitle("Cell Type Prediction Scores") & NoLegend() & ylim(0,1)
    ggsave(glue("{projName}/Azimuth_PredictionScoreL2.png"), plot = p2, width=7, height=5, dpi=300, bg="white") 
    p3 <- VlnPlot(azimuthComparisonObj, features = "mapping.score", group.by = "predicted.celltype.l2", pt.size=0) & 
    theme(axis.text.x = element_text(angle = 90)) & ggtitle("Cell Type Mapping Scores") & NoLegend() & ylim(0,1)
    ggsave(glue("{projName}/Azimuth_MappingScoresL2.png"), plot = p3, width=7, height=5, dpi=300, bg="white")
    p2 <- VlnPlot(azimuthComparisonObj, features = "predicted.celltype.l1.score", group.by = "predicted.celltype.l1", pt.size=0) & 
    theme(axis.text.x = element_text(angle = 90)) & ggtitle("Cell Type Prediction Scores") & NoLegend() & ylim(0,1)
    ggsave(glue("{projName}/Azimuth_PredictionScoreL1.png"), plot = p2, width=5, height=5, dpi=300, bg="white") 
    p3 <- VlnPlot(azimuthComparisonObj, features = "mapping.score", group.by = "predicted.celltype.l1", pt.size=0) & ylim(0,1)
    theme(axis.text.x = element_text(angle = 90)) & ggtitle("Cell Type Mapping Scores") & NoLegend() 
    ggsave(glue("{projName}/Azimuth_MappingScoresL1.png"), plot = p3, width=5, height=5, dpi=300, bg="white")
    p4 <- ggplot(azimuthComparisonObj@meta.data, aes(x=predicted.celltype.l2.score, y=mapping.score, colour = predicted.celltype.l2)) + 
    geom_point(size=2, shape=23) + xlim(0,1) + ylim(0,1)
    ggsave(glue("{projName}/Azimuth_ScoreScatterL2.png"), plot = p4, width=7, height=5, dpi=300, bg="white") 
    p5 <- ggplot(azimuthComparisonObj@meta.data, aes(x=predicted.celltype.l1.score, y=mapping.score, colour = predicted.celltype.l1)) + 
    geom_point(size=2, shape=23) + xlim(0,1) + ylim(0,1)
    ggsave(glue("{projName}/Azimuth_ScoreScatterL1.png"), plot = p5, width=7, height=5, dpi=300, bg="white")
    cellTypeData <- azimuthComparisonObj@meta.data[,c("predicted.celltype.l1", "predicted.celltype.l1.score", "predicted.celltype.l2", "predicted.celltype.l2.score","nCount_RNA", 'mapping.score')]
    cellTypeData$CB <- rownames(cellTypeData)
    colnames(cellTypeData) <- c("L1_CellType","L1_CellType_Score", "L2_CellType", "L2_CellType_Score", "nCount_RNA", "Mapping_Score", "CB") # CB and mapping previously in the wrong order- check the report for issues fixing this
    write.table(cellTypeData, glue("{projName}/Azimuth_ClusteringData.tsv"), row.names=FALSE, quote=FALSE,  sep='\t')
    p1 <- FeaturePlot(azimuthComparisonObj, features = "nCount_RNA") + ggtitle("nCount_RNA on Azimuth Umap") + theme(plot.title = element_text(hjust = 0.5))
    ggsave(glue("{projName}/Azimuth_ClusteringUmap_nCounts.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
}

# Calculates markers by cluster and writes to CSV 
findAndWriteMarkers <- function(projName,seuratObj) {
    all.markers <- FindAllMarkers(seuratObj)
    write.csv(all.markers, glue("{projName}/Seurat_ClusterMarkers.csv"), row.names=FALSE,quote=FALSE)
}

# 1. Creates seurat object
# 2. Normalizes and runs PCA on data
# 3. Clusters cells based on PCA output  
produceSeuratObj <- function(matrixObj, projName) {
    # Only include genes seen in min.cells, and that express atleast min.features unique features 
    seuratObj <- CreateSeuratObject(counts = matrixObj, project = projName) #, min.cells = 3, min.features = 200)
    Idents(seuratObj)    <- "all"
    seuratObj$orig.ident <- "all"
    seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
    #seuratObj[["percent.ribo"]] <- PercentageFeatureSet(seuratObj, "^RP[SL]")
    # save QC plots before filtering
    p1 <- VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) & geom_boxplot(width=0.1, fill="white", outlier.color = NA)
    ggsave(glue("{projName}/Seurat_QCboxplots.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
    p1 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme_bw() + theme(legend.position = "none")
    p2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme_bw() + theme(legend.position = "none")
    ggsave(glue("{projName}/Seurat_QCscatter.png"), plot = patchwork::wrap_plots(p1, p2, ncol = 2), width=10, height=5, dpi=300, bg="white")
    return(seuratObj)
}

filterSeuratObjDoublets <- function(seuratObj, scrubletOut, projName) {
    scrublet <- read.table(scrubletOut, sep=",", header=TRUE)
    # plot doublet score
    p5 <-ggplot(scrublet, aes(x=predicted_doublet, y=doublet_score, color=predicted_doublet, fill = predicted_doublet)) + geom_boxplot(alpha=0.2) + theme_bw()
    ggsave(glue("{projName}/scrublet_doublet_scores_boxplot.png"), plot = p5, width=5, height=5, dpi=300, bg="white")
    # remove doublets
    seuratObj[["predicted_doublet"]] <- scrublet$predicted_doublet
    seuratObj[["doublet_score"]]     <- scrublet$doublet_score
    rm(scrublet)
    seuratObj                        <- subset(x= seuratObj, subset= predicted_doublet == "False")
    # filter for low quality cells
    seuratObj <- subset(seuratObj, subset = nFeature_RNA > opts$minFeatures & nCount_RNA < opts$maxCounts & percent.mt < opts$MTpercent & nFeature_RNA < 10000) # add max nFeature option!
    p1 <- VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) & geom_boxplot(width=0.1, fill="white", outlier.color = NA)
    ggsave(glue("{projName}/Seurat_QCboxplotsFilt.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
    p1 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme_bw() + theme(legend.position = "none")
    p2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme_bw() + theme(legend.position = "none")
    ggsave(glue("{projName}/Seurat_QCscatterFilt.png"), plot = patchwork::wrap_plots(p1, p2, ncol = 2), width=10, height=5, dpi=300, bg="white")
    return(seuratObj)
}

filterTopCellsByUmi <- function(seuratObj, nCells, projName) {
    BCs <- head(seuratObj$nCount_RNA[order(seuratObj$nCount_RNA, decreasing=TRUE)], nCells) %>% names
    seuratObj[["topNcells"]] <- rownames(seuratObj[[]]) %in% BCs
    seuratObj <- subset(seuratObj, subset = topNcells == TRUE )
    seuratObj[["topNcells"]] <- NULL
    p1 <- VlnPlot(seuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) & geom_boxplot(width=0.1, fill="white", outlier.color = NA)
    ggsave(glue("{projName}/Seurat_QCboxplots_nCellsFilter.png"), plot = p1, width=7, height=7, dpi=300, bg="white")
    p1 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt") + theme_bw()
    p2 <- FeatureScatter(seuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + theme_bw()
    ggsave(glue("{projName}/Seurat_QCscatter_nCellsFilter.png"), plot = patchwork::wrap_plots(p1, p2, ncol = 2), width=9, height=5, dpi=300, bg="white")
    seuratObj
}

normalizeAndCluster <- function(seuratObj, projName, clustRes) {
    seuratObj <- normalizeScaleDimRed(seuratObj, projName)
    seuratObj <- clusterCells(seuratObj, projName, clustRes)
    return(seuratObj)
}

produceSCTypeFigures <- function(projName,seuratObj) {
    umap <- DimPlot(seuratObj, reduction = "umap")
    cellTypeUmap <- DimPlot(seuratObj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') & ggtitle('Cell Types UMAP')
    ggsave(glue("{projName}/ScType_Clustering.png"), plot= cellTypeUmap, width=7, height=7, dpi=300)
}

normalizeScaleDimRed <- function(seuratObj, projName) {
    seuratObj <- NormalizeData(seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
    seuratObj <- FindVariableFeatures(seuratObj, selection.method = "vst", nfeatures = 2000)
    # plot variable features 
    top10 <- head(VariableFeatures(seuratObj), 10)
    p1    <- VariableFeaturePlot(seuratObj)
    p2    <- LabelPoints(plot = p1, points = top10, repel = TRUE)
    ggsave(glue("{projName}/Seurat_varFeatures.png"), plot = p2, width=7, height=5, dpi=300, bg="white")
    seuratObj <- ScaleData(seuratObj)#, features = rownames(seuratObj))
    seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj))
    p1 <- ElbowPlot(object = seuratObj, ndims = 60)
    ggsave(glue("{projName}/Seurat_elbow.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
    return(seuratObj)
}

clusterCells <- function(seuratObj, projName, clustRes) {
    seuratObj <- FindNeighbors(seuratObj, dims = 1:15)
    seuratObj <- FindClusters(seuratObj, resolution = clustRes)
    seuratObj <- RunUMAP(seuratObj, dims = 1:15)
    p1 <- DimPlot(seuratObj, reduction = "umap") 
    ggsave(glue("{projName}/Seurat_clusterUmap.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
    return(seuratObj)
}

makeSeuratClusterPlots <- function(seuratObj, projName){
    # sample overlaid with orig sample name
    if(opts$plotOrigSample == "yes"){
        idDf <- matrix(unlist(strsplit( rownames(seuratObj[[]]) ,split = "_")) ,ncol=opts$plotOrigSampleLevels,byrow=T)[,1:opts$plotOrigSampleLevels] %>% as.data.frame
        cols                  <- colnames(idDf)[1:(ncol(idDf)-1)]
        if(length(cols) > 1){
            seuratObj[["sample"]] <- apply( idDf[ , cols ] , 1 , paste , collapse = "_" )
        }else{
            seuratObj[["sample"]] <- idDf[,cols]
        }
        p1 <- DimPlot(seuratObj, group.by = "sample")
        ggsave(glue("{projName}/Seurat_origSampleUmap.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
        p1 <- DimPlot(seuratObj, reduction = "pca", group.by = "sample")
        ggsave(glue("{projName}/Seurat_origSamplePCA.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
        ## TO DO add a plot with the number and proportion of sample Ids per cluster
    }
    if(opts$plotRtBc == "yes"){
        # add RT
        seuratObj[["RT"]] <- sub(".*_", "", Cells(seuratObj))
        seuratObj[["RT"]] <- sapply(seuratObj$RT, function(x) {substr(x,10,19)})
        p1                <- DimPlot(seuratObj, group.by="RT", reduction = "umap") + ggtitle(paste(basename(projName), "RT"))
        ggsave(glue("{projName}/Seurat_RT_Umap.png"), plot = p1, width=20, height=15, dpi=300, bg="white")
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
        ggsave(glue("{projName}/Seurat_rtClusterProportions.png"), plot = rtBarGraph, width=20, height=15, dpi=300, bg="white")
    }
    p1 <- FeaturePlot(seuratObj, features="nCount_RNA", reduction = "umap") + ggtitle("Cluster nCount_RNA")
    ggsave(glue("{projName}/Seurat_nCount_Umap.png"), plot = p1, width=6, height=5, dpi=300, bg="white")
    p1 <- FeaturePlot(seuratObj, features="nFeature_RNA", reduction = "umap") + ggtitle("Cluster ")
    ggsave(glue("{projName}/Seurat_nFeature_Umap.png"), plot = p1, width=6, height=5, dpi=300, bg="white")
    p1 <- VlnPlot(seuratObj, features = "nCount_RNA", ncol = 1, pt.size=0) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    ggsave(glue("{projName}/Seurat_clsuterCountsVln.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
    p1 <- VlnPlot(seuratObj, features = "nFeature_RNA", ncol = 1, pt.size=0) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    ggsave(glue("{projName}/Seurat_clsuterFeaturesVln.png"), plot = p1, width=7, height=5, dpi=300, bg="white")
    p1 <- VlnPlot(seuratObj, features = "predicted.celltype.l2.score", ncol = 1, pt.size=0) & theme_bw() & theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
    ggsave(glue("{projName}/Seurat_clsuterAzimuthPredictedScoreL2Vln.png"), plot = p1, width=7, height=5, dpi=300, bg="white")

    ## azimuth on UMAP
    Idents(seuratObj) <- seuratObj$predicted.celltype.l2
    p1 <- DimPlot(seuratObj, reduction = "umap")
    ggsave(glue("{projName}/Seurat_clusterAzimuthCellTypeL2.png"), plot = p1, width=10, height=5, dpi=300, bg="white")
    Idents(seuratObj) <- seuratObj$predicted.celltype.l1
    p1 <- DimPlot(seuratObj, reduction = "umap")
    ggsave(glue("{projName}/Seurat_clusterAzimuthCellTypeL1.png"), plot = p1, width=9, height=5, dpi=300, bg="white")
    Idents(seuratObj) <- seuratObj$seurat_clusters

}

makeMarkerPlotsPostClustering <- function(projName, seuratObj){
    Idents(seuratObj) <- seuratObj$predicted.celltype.l1
    markers.to.plot <- c("CD4", "CD3E", "CD19", "CD8A", "NCAM1")
    p1 <- DotPlot(seuratObj, features = markers.to.plot, dot.scale = 8, assay = "RNA") +
        RotatedAxis()
    ggsave(glue("{projName}/Azimuth_Level1PbmcMarkersDotPlot.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
    Idents(seuratObj) <- seuratObj$predicted.celltype.l2
    p1 <- DotPlot(seuratObj, features = markers.to.plot, dot.scale = 8, assay = "RNA") +
        RotatedAxis()
    ggsave(glue("{projName}/Azimuth_Level2PbmcMarkersDotPlot.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
    plots <- VlnPlot(seuratObj, features = markers.to.plot, combine = FALSE, pt.size=0)
    p1 <- wrap_plots(plots = plots, ncol = 1)
    ggsave(glue("{projName}/Azimuth_Level2PbmcMarkersVlnPlot.png"), plot = p1, width=10, height=20, dpi=300, bg="white")
    Idents(seuratObj) <- seuratObj$predicted.celltype.l1
    plots <- VlnPlot(seuratObj, features = markers.to.plot, combine = FALSE, pt.size=0)
    p1 <- wrap_plots(plots = plots, ncol = 1)
    ggsave(glue("{projName}/Azimuth_Level1PbmcMarkersVlnPlot.png"), plot = p1, width=10, height=20, dpi=300, bg="white")
    p1 <- FeaturePlot(seuratObj, features = c("CD4"))
    ggsave(glue("{projName}/CD4_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
    p1 <- FeaturePlot(seuratObj, features = c("CD3E"))
    ggsave(glue("{projName}/CD3E_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
    p1 <- FeaturePlot(seuratObj, features = c("CD8A"))
    ggsave(glue("{projName}/CD8A_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
    p1 <- FeaturePlot(seuratObj, features = c("CD19"))
    ggsave(glue("{projName}/CD19_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
    p1 <- FeaturePlot(seuratObj, features = c("NCAM1"))
    ggsave(glue("{projName}/NCAM1_CD56_clusterUmap.png"), plot = p1, width=6, height=6, dpi=300, bg="white")
}

main(projName, countDir, scrubletOut, clustRes, nCells)


