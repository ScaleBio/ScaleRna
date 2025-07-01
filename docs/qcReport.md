# ScaleBio RNA QC reports

## Library report
The file called _library\_{libIndex2}.report.html_ contains the summary report at the library level, i.e. sequencing, barcode and demultiplexing information for all samples processed in one run of the Quantum Scale RNA kit.

### Read Status
This table gives the barcode matching statistics for the full library. Reads that fail barcode matching are not assigned to any sample.
 
**Barcode Pass**: Reads for which all expected barcodes were found. Includes TooShortError which are reads with less than 16 (`min_length`) bases of RNA sequence after adapter and Poly-A trimming. These are assigned to their sample based on RT barcode and included in "_Total Sample Reads_", but are filtered before genome alignment \
**Barcode Error**: Reads which were filtered because at least one barcode could not be found or matched against the expected sequences (whitelist). These reads are excluded from all further analysis
  
### Reads per Sample
**Reads per Sample** shows the number of reads assigned to each sample based on the RT (sample) barcode. These read counts are before adapter trimming, alignment and duplicate filtering.

### Barcodes tab
Plate-maps showing the total unique transcript counts for each RT barcode

## Sample Report
The files called _{SampleName}.{LibraryName}.report.html_ contain the summary report for a single sample, i.e. the set of RT wells (Sample Barcodes) from a Scale RNA dataset assigned to one sample in [samples.csv](samplesCsv.md). It shows read, cell and barcode level summary metrics and plots for sample QC. _Note_ that sample here could refer to a whole [ScalePlex pool](scalePlex.md).

### Read Metrics
**Total Sample Reads**: The number of reads assigned to the sample based on the RT (sample) barcode, after barcode error reads are removed. This is before adapter trimming, alignment and duplicate filters \
**Passing Sample Reads**: Reads passing pre-alignment filters, specifically RNA sequence length after Poly-A trimming ("_Too Short Error_") \
**Reads Mapped to Genome**: The fraction of _Passing Reads_ that are aligned anywhere to the genome. This includes multimapping reads \
**Passing Read Alignments**: The fraction of mapped reads retained after alignment filtering; specifically excluding multimappers to more than 6 (`starMaxLoci`) loci \
**Reads Mapped to Transcriptome**: The fraction of _Passing Read Alignments_ that match one or more annotated genes (exon or intron, in sense direction) \
**Exonic Reads**: The fraction of _Reads Mapped to Transcriptome_ overlapping an exon in the sense direction \
**Antisense Reads**: The fraction of _Passing Read Alignments_ overlapping an exon in the antisense direction (opposite from gene annotation) \
**Mitochondrial Reads**: The fraction of reads mapping to the mitochondrial genome (`chrM`) \
**Saturation**: The overall sequencing saturation level of this sample; defined as `1 - (UniqueReads / TotalReads)` on _Reads Mapped to Transcriptome_ \
**Average Trimmed Read Length**: The average length of RNA reads after adapter and Poly-A trimming \
**Average Mapped Length**: The average length of read alignments to the genome (after softclipping) \
**Mapped Mismatch Rate**: The mismatch rate between mapped reads and the genome (variants and sequencing errors)


### Cell Metrics
**Note**: All numbers in this table depend on the number of cells called for the sample. Before interpreting these numbers, check the _Barcode Rank plot_ to confirm that it matches expectations. 

**Cells called**: The number of cell barcodes passing filters to be [called as a cell](cellCalling.md) \
**Reads per cell**: The overall number of _Total Sample Reads_ divided by the number of cells called \
**Median Unique Transcript Counts per cell**: The median number of transcripts detected per cell; i.e. unique reads mapped to the transcriptome \
**Median Genes per cell**: The median number of unique genes detected per cell \
**Reads in Cells**: The fractions of transcriptome reads that belong to a cell rather than a background barcode


### Plots
#### Barcode Rank plot
This shows the unique transcript counts for each cell-barcode, sorted from high to low. Each point along the line indicates the proportion of barcodes called as cells (as opposed to background) with that transcript count.

#### Complexity plot
This shows a statistical estimate for the unique transcript counts per cell that would be observed at different shallower sequencing levels for the sample ("_Total Sample Reads_")

#### Genes Detected Per Cell
This scatterplot shows the number of unique genes detected for each cell-barcode relative to the reads for that barcode, separating cells from background barcodes. 

#### Saturation Per Cell
This scatterplot shows the number of reads vs. the sequencing saturation (`1 - (UniqueReads / TotalReads)`) for each cell-barcode.

### Barcodes Tab
The plot on the left show the number of cells called with each sample (RT) barcode. \
The plot of the right show the median unique transcript counts for cells with each barcode.

A breakdown of other barcodes is in the library report.
