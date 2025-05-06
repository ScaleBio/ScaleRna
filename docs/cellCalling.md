# Cell Calling Options
For QuantumScale RNA, single cells are identified by a combination of RT (sample) barcode, bead barcode and library (PCR) index. However not all barcode combinations correspond to real cells, but instead to background, e.g. empty beads. The workflow implements multiple approaches to cells against background barcodes.

## CellFinder
_CellFinder_ is an [EmptyDrops](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y)-like cell calling approach used by default (`cellFinder` parameter). Cell-barcodes that are above a minimal transcript count (`minUTC`) but below the provided or calculated `UTC` threshold can be _rescued_ based on expression differences from the ambient RNA profile. For a detailed description see the ScaleRna workflow handbook. _Note_ that CellFinder is disabled for multi-species (_barnyard_) samples.

### Parameters
parameter | description | default
-- | -- | --
cellFinder | Enables CellFinder cell calling | True
UTC | Set a fixed unique transcript count thresholds for barcodes to be called a cell; Optional | 0
cellFinderFdr | False discovery rate to use for rescuing cells based on deviation from ambient profile. | 0.001

## TopCell Thresholding
_TopCell_ thresholding sets a hard threshold on the unique transcript count per cell-barcode; all barcodes with counts above this threshold are called as cells, while all other barcodes are background. This threshold is determined dynamically as follows: 

1. Barcodes with counts less than the parameter `minUTC` are filtered from the data.
2. Next, the *top cell*, the cell barcode  `topCellPercent` percentile count (if the `expectedCells` parameter is set, the percentile is applied to that number). 
3. The UTC *top cell* is then divided by the parameter `minCellRatio` to determine the unique transcript count threshold for the dataset.

### Parameters
parameter | description | default 
-- | -- | -- 
minUTC | Minimum number of counts to consider a barcode as a potential cell | 100
expectedCells | Approximate number of cells expected per sample; Optional; Can be set in `samples.csv` | 0
topCellPercent | Percentile to use for the *top cell* (_robust max._) | 99
minCellRatio | Ratio between transcript counts of top cell and the lower cell threshold | 10

## Fixed Threshold
A unique transcript count threshold is used when the parameter `UTC` is greater than 0 and *cellFinder* is not enabled.
In this case every barcode with a total count >= `UTC` is called as a cell.

## Fixed Cell Number Threshold
With this option a fixed number of cell barcodes with the highest UTCs are called. The number is given by `expectedCells` which can be set per sample in `samples.csv`

### Parameters 
parameter | description | default
-- | -- | --
fixedCells | Flag to enable fixed cell number calling | False
expectedCells | The number of cells to call | 0

# Outlier Identification
The workflow can flag potential outlier cells, based on multiple QC metrics. Barcodes are flagged in `allcells.csv` based on a Median Absolute Deviation (MAD) threshold for number of reads, fraction of usable reads (counted towards a gene), and fraction of mitochondrial reads. These parameters are defined in [create_mtx.config](../modules/create_mtx.config). By default cells classified as outliers are not filtered from the UMI count matrix. It is possible to remove cells classified as outliers using the parameter `filterOutliers`.
