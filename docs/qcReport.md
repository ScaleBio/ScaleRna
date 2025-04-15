# ScaleBio RNA QC reports

## Library report
The file called _library\_{LibraryName}.report.html_ contains the summary report at the library level, i.e. sequencing, barcode and demultiplexing information for all samples processed in one run of the ScaleRNA kit.

### Barcode Read Status
This table gives the barcode matching statistics for the full library. Reads that fail barcode matching are not assigned to any sample and are hence not included in any downstream analysis or metrics.
 
**Pass**: Reads for which all expected barcodes were found \
**Error**: Reads which were filtered because at least one barcode could not be found or matched against the expected sequences (whitelist). These reads are excluded from all further analysis.
  
### Reads per Sample
**Reads per Sample** shows the number of reads assigned to each sample based on the RT barcode match. These are total reads before alignment and duplicate removal.


### Barcodes tab
Plate-maps showing the total unique transcript counts (complexity) for each well (barcode) of the RT, ligation and PCR plates respectively.


## Sample Report
The files called _{SampleName}.{LibraryName}.report.html_ contain the summary report for a single sample, i.e. all or a subset of RT wells from a ScaleRNA library. It shows read, cell and barcode level summary metrics and plots for library and sample QC.

### Mapping Metrics
**Passing Reads**: Reads passing pre-alignment filters, specifically at least 16bp after Poly-A trimming. \
**Reads Mapped to Genome**: The fraction of _Passing Reads_ that are aligned anywhere to the genome \
**Reads Mapped to Transcriptome**: The fraction of _Reads Mapped to Genome_ that match one or more annotated genes (exon or intron, in sense direction) \
**Exonic Reads**: The fraction of _Reads Mapped to Genome_ overlapping an exon (by at least 50% of their length) in the sense direction \
**Antisense Reads**: The fraction of _Reads Mapped to Genome_ overlapping an exon in the antisense direction (opposite from gene annotation) \
**Mitochondrial Reads**: The fraction of reads mapping to the mitochondrial genome \
**Saturation**: The overall sequencing saturation level of this sample. _Saturation_ is defined as `1 - (UniqueReads / TotalReads)` on _Reads Mapped to Transcriptome_.


### Cell Metrics
**Note**: All numbers in this table depend on the number of cells called and hence on the _Unique Transcript Counts Threshold_. Check the threshold indicated on the _Barcode Rank plot_.

**Unique Transcript Counts Threshold**: The minimum number of unique reads mapped to transcripts required to separate cells from background barcodes \
**Cells above Threshold**: The number of cell barcodes passing the _Unique Transcript Counts Threshold_ \
**Mean Reads per cell**: The mean number of _Passing Reads_ for each cell \
**Median Unique Transcript Counts per cell**: The median number of unique reads for each cell matching a transcript; i.e. the number of transcripts detected \
**Median Genes per cell**: The median number of unique genes detected per cell \
**Reads in Cells**: The fraction of reads that come from a cell rather than a background barcode


### Plots
#### Barcode Rank plot
This shows the unique transcript counts for each barcode, sorted from high to low. When using Cell Finder cell calling, the points are colored by the proportion of called cells. This proportion is a rolling average of the proportion of called cells, calculated using a window size of 25. Otherwise when cell calling is performed using a _Unique Transcript Counts Threshold_ this threshold is indicated by a dashed line, separating (estimated) cells from background barcodes.

#### Complexity plot
This shows a statistical estimate for the unique transcript counts that would be observed at different shallower sequencing levels.

#### Genes Detected Per Cell
This shows the number of unique genes detected for each cell-barcode relative to the total reads for that barcode. _pass_ are cells passing the _Unique Transcript Counts Threshold_.

#### Saturation Per Cell
This shows the total number of reads vs. the sequencing saturation for each cell-barcode.

### Barcodes Tab
The plots on the left show the number of cells with each sample barcode. \
The plots of the right show the complexity (median unique transcript counts per cell) for each of these.

A breakdown of other barcodes is in the library report.
