#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))

parseArguments <- function(){

    message("Parsing Arguments")
    
    argParser <- ArgumentParser(
        description = "This script takes in the two paths --resultsDir and --sampleID. Both arguments are required. The script will read in all of the result .csv files in the results directory for the given sample id. The script will then combine these data into a single data.frame and write out the combined results to a csv file"
    )

    argParser$add_argument("--resultsDir", action = "store", required = TRUE, 
        help = "Path to directory which contains the results from the Downstream Workflow.")
    argParser$add_argument("--sample", action = "store", required = TRUE)

    argList <- argParser$parse_args()
    return(argList)
}

readAllCellsData <- function(filePath){

    message("Reading allCells.csv")

    dat <- fread(filePath)
    dat <- setnames(dat, old = 1, new = "cellBarcode")
    cellBarcodes <- paste0(dat[["cellBarcode"]], "_", dat[["sample"]])
    dat[["cellBarcode"]] <- cellBarcodes
    dat <- dat[dat$pass == TRUE, ]

    return(dat)
}

combineMetaData <- function(dirPath, sampleID){

    message("Combining Data")

    resultPaths <- list.files(path = dirPath, pattern = "_results", full.names = TRUE)
    allCellsPath <- list.files(path = dirPath, pattern = "allCells.csv", full.names = TRUE)
    allCells <- lapply(allCellsPath, readAllCellsData)
    allCells <- rbindlist(allCells)
    fileList <- resultPaths
    datList <- lapply(fileList, read.csv)
    datList[["allCells"]] <- allCells
    dat <- Reduce(function(x,y) merge(x, y, all = TRUE), datList)

    write.csv(x = dat, file = paste0(sampleID, "_combinedResults.csv"), row.names = FALSE, quote = FALSE)
}

argList <- parseArguments()

combineMetaData(dirPath = argList$resultsDir, sampleID = argList$sample)