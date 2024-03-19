#!/usr/bin/env Rscript

library(rmarkdown)
library(argparse)

parseArguments <- function(){

    message("Parsing Arguments")
    
    argParser <- ArgumentParser(
        description = "This script takes in the two paths --resultsDir and --sampleID. Both arguments are required. The script will read in all of the result .csv files in the results directory for the given sample id. The script will then combine these data into a single data.frame and write out the combined results to a csv file"
    )

    argParser$add_argument("--rmdPath", action = "store", required = TRUE)
    argParser$add_argument("--results", action = "store", required = TRUE, 
        help = "Path to directory which contains the results from the Downstream Workflow.")
    argParser$add_argument("--sampleStats", action = "store", required = TRUE)
    argParser$add_argument("--comparison", action = "store", required = TRUE)

    argList <- argParser$parse_args()
    return(argList)
}

argList <- parseArguments()

paramList <- list(
    downstreamResults = argList$results,
    sample_stats = argList$sampleStats,
    is_comparison = as.logical(argList$comparison)
)

render(input = argList$rmdPath, params = paramList, output_file = "report.html", output_format = "html_document")
