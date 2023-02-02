# ScaleBio RNA QC reports

## Library report
The file called _library_<LibraryName>.report.html_ contains the summary report at the library level, i.e. sequencing, barcode and demultiplexing information for all samples processed in one run of the ScaleRNA kit.
  

### Barcode Read Status
This table gives the barcode matching statistics for the full library. Reads that fail barcode matching are not assigned to any sample and are hence not included in any downstream analysis or metrics.
 
**Pass**: Reads for which all expected barcodes were found \
**LinkerError**: Reads which were filtered because the fixed linker sequence between the barcodes could not be found \
**BarcodeError**: Reads for which at least one barcode could not be matched against the expected sequences (whitelist) \
**SequenceError**: Reads excluded from barcode matching, e.g. because they were too short
  
### Reads per Sample
**Reads per Sample** and **Reads Per RT Well** show the number of (Barcode passing) reads assinged to each sample based on the RT barcode. These are read counts before alignment and UMI collapsing.

### Cells tab
These tables show detailed metrics gor each barcode element (RT, Ligation, PCR barcodes and UMI).

**Exact**: The read contains an exact (no errors) match against one of the expected sequences (whitelist) for this barcode \
**Corrected**: The barcode sequence contains at least one mismatch, but could be corrected to an unique whitelist sequence \
**Ambiguous**: The barcode sequence has the same number of mismatches to two different whitelist sequences and is hence filtered \
**NoMatch**: The barcode sequence cannot be matched to any whitelist sequence \
**Error**: No barcode sequence can be extracted; typically because the linker sequence used to locate it in the read was not found \

The UMI is a random sequence with no whitelist. In that case all sequences containing 'N's or pure homopolymers are **Filtered**, other sequences are **Pass**.
  
## Sample Report

The files called _<SampleName>.report.html_ contain the summary report for a single sample, i.e. all or a subset of RT wells from a ScaleRNA library. It shows read, cell and barcode level summary metrics and plots for library and sample QC.

### STARSolo Metrics

**Total Reads**: The number of reads (pairs) input for this sample. This is after matching the barcodes (and possibly demultiplexing on the RT barcode), but before any alignent filters \
**Reads Mapped to Genome**: The fraction of _Total Reads_ that are aligned anywhere to the genome \
**Reads Mapped to Transcript**: The fraction of _Total Reads_ that are mapped to the genome overlapping an annotated transcript (exon or intron, in sense direction). \
**Reads Mapped to unique Transcript**: The fraction of _Total Reads_ that are mapped to the genome overlapping exactly one transcript uniquely \

### Cell Calling
All numbers in this table depend on the selected threshold to separate cells from background barcodes (automatic or manual)

**Cells above Threshold**: The number of cell barcodes passing the _Unique Reads Threshold_ \
**Unique Reads Threshold**: The UMI threshold to separate cells from background barcodes \
**Mean Reads per cell**: The mean number of reads for each cell barcode passing the _Unique Reads Threshold_ \
**Median UMIs per cell**: The median number of unique reads (UMIs) for each cell counted towards expression of a gene; i.e. the number of transcripts detected \
**Median Genes per cell**: The median number of unique genes detected per cell \
**Reads in Cells**: The fraction of reads that come from a cell passing threshold rather than a background barcode \
**Median Saturation**: The median sequencing saturation level in a cell. _Saturation_ is 1 - (_UMIs_ / _genic reads_) \

**Median Genic Reads**: The median per-cell fraction of reads mapping to a transcript (exon or intron) in sense direction \
**Median Exonic Reads**: The fraction of _Genic Reads_ overlapping an exon (by at least 50% of their length) \
**Median Antisense Reads**: The proportion of reads overlapping a transcript in the antisense orientation relative to *Genic Reads* (sense orientation) \
**Median Mitochondrial Reads**: The fraction of reads mapping to the mitochondrial genome

### Plots
#### Genes vs. UMIs detected
This shows the number of UMIs (transcripts detected) and unique genes detected for each cell-barcode. _pass_ are cells passing the threshold.

#### Saturation Per Cell
This shows the total number of reads vs. the sequencing saturation for each cell-barcode. _pass_ are cells passing the threshold.

### Barcodes Tab
The plots on the left show the number of cells (cell-barcodes passing the threshold) with each RT, ligation and PCR barcode respectively. \
The plots of the right show the number of UMIs per cell for each of these.
