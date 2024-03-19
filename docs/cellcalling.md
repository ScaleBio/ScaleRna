# Cell Calling Options

The ScaleBio Seq Suite: RNA Workflow has four different approaches to distinguish true cellular barcodes from ambient barcodes. 

## Default Thresholding

The default cell thresholding algorithm used by the ScaleRna workflow imposes a one-dimensional threshold on the UMI count per barcode ("unique transcript counts threshold"), determining that only barcodes with UMI counts above this threshold belong to real cells, while all other barcodes are associated with ambient RNA. This thresholding is performed according to the following logic: 

(1) First, barcodes with UMI counts less than the parameter *minUTC* are filtered from the data.

(2) Next, we calculate the UMI count associated with the *topCellPercent* percentile of this filtered dataset (using only top *expectedCells* number of cells if that parameter is > 0). The UMI count of the n-th percentile of these preliminary "cells" is then divided by the parameter *minCellRatio* to determine the unique transcript count threshold for the entire dataset. All barcodes with UMI counts greater than this unique transcript count threshold are called as cells.

### Parameters

parameter | description | default 
-- | -- | -- 
topCellPercent | Percentage of cells over minUTC to use as 'robust max.' | 99
minCellRatio | Ratio between transcript counts of top cells and the lower cell threshold | 10
minUTC | Minimum number of counts to consider a barcode as a potential cell | 100
expectedCells | Number of cells expected in the dataset | 0
UTC | Minimum number of unique transcript counts a barcode must be associated with to be considered a called cell. | 0

## CellFinder Thresholding

CellFinder is an [EmptyDrops](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)-like cell calling approach. This option is enabled when the parameter *cellFinder* is true. Cells that are above the *minUTC* but below the provided or calculated *UTC* can potentially be "rescued" based on expression differences from the ambient RNA. For a detailed description see the ScaleRna workflow handbook. 

### Parameters

parameter | description | default
-- | -- | --
cellFinder | Flag which enables CellFinder cell calling | TRUE
expectedCells | Number of cells expected in the dataset | 0
topCellPercent | Percentage of cells over minUTC to use as 'robust max.' | 99
minCellRatio | Ratio between transcript counts of top cells and the lower cell threshold | 10
minUTC | min. counts to consider a barcode as a potential cell | 100
UTC | The minimum number of unique transcript counts a barcode must be associated with to be considered a called cell. | 0
FDR | False discovery rate to use for rescuing cells based on deviation from ambient profile. | 0.01

## Unique Transcript Count Threshold

A unique transcript count threshold is used when the parameter *UTC* is greater than 0 and *cellFinder* is not enabled.
In this case every barcode with a UMI count greater than or equal to *UTC* is called as a cell.

### Parameters 

parameter | description | default
-- | -- | --
cellFinder | Flag which enables CellFinder cell calling | FALSE
UTC | If > 0 every barcode with UMIs >= *UTC* will be called as a cell. | > *minUTC* 

## Fixed Cell Number Threshold

A fixed number (*fixedCells*) of cell barcodes are returned when this parameter is true. In this case barcodes are sorted in descending order by UMI count and the top *expectedCells* barcodes are called as cells.

### Parameters 

parameter | description | default
-- | -- | --
fixedCells | Flag which enables fixed cell number thresholding | 0
expectedCells | Value which determines the number of cell barcodes to call | 0

# Outlier Identification

After identifying cell barcodes the workflow then identifies low quality (outlier) cells based on three metrics: number of umis, number of genes, and number of mitochondrial reads. These low quality cells are identified using Medain Absolute Deviation (MAD) thresholding for each of the three metrics. Different MAD thresholds can be set for each metric independently. By default cells classified as outliers are not filtered from the UMI count matrix. It is possible to remove cells classified as outliers using the parameter *filter_outliers*. It is important to note that this will filter all cells classified as outliers regardless of metric.

### Parameters

parameter | description | default
-- | -- | -- 
filter_outliers | Enables option to filter cells classified as outliers out of UMI count matrix | FALSE
num_mad_genes | Median absolute deviation threshold for number of genes | 8
num_mad_umis | Median absolute deviation threshold for number of umis | 8
num_mad_mito | Median absolute deviation threshold for mitochondrial read proportion | 5
