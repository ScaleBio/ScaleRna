# Outputs

Analysis workflow results for all samples in the run will appear in the output directory defined by `--outDir` (`ScaleRna.out` by default). 
For detailed information about the library and sample level QC reports see [qcReport.md](qcReport.md)

## Key output files
| Directory | File | Description |
|-----------|------|-------------|
| `reports`| `multiqc_report.html` | [MultiQC](https://multiqc.info/) report for fastq generation, fastQC and trimming |
| | `<sample>.<libName>.report.html` | A standalone report including key QC metrics and figures for each sample; (`merged` for extended throughput runs)|
| | `<sample>_libraries` | For extended throughput runs, individual sample reports for each plate separately
| | `allSamples.reportStatistics.csv` | QC metrics from all samples in this analysis in one table
| | `csv/` | Summary and QC metrics for this sample in csv format |
| `reports/library` | `library_<libName>.report.html` | Barcode summary and demultiplexing statistics for the whole sequencing library |
| | `csv/` | Summary and QC metrics for this library in csv format | 
|  `samples` | `<sample>.<libName>.filtered.matrix/` | Pre-filtered gene expression matrix for passing cells; `merged` for results combined across multiple extended throughput plates |
| | `<sample>.<libName>.allCells.csv` | Metrics per called cell, including barcodes / well positions
| | `<sample>_libraries/` | For extended throughput runs, this contains output files per plate
| `fastq` | `fastqc/*_fastqc.html` | [fastqc](https://github.com/s-andrews/FastQC) report for each fastq file in the sequencing library |
| | `Reports/` | Fastq generation summary reports from [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) |
| `barcodes/<libName>.demux` | `<sample>.bam` | Sample unaligned bam files (Demultiplexed and barcode error-corrected); only included with `--bcParserBamOut true` |
| | `metrics.json` | Detailed barcode and demultiplexing information, including individual barcode error rates | 
| `alignment/<sample>.<libName>` | `<sample>.star.align/` | [STAR](https://github.com/alexdobin/STAR) alignment output, including BAM file, with single-cell barcode and UMI information in tags
|  | `<sample>.star.solo/` | [STARSolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) output for each sample, including unfiltered `.mtx` files


The gene expression matrix directories (.mtx) can be loaded into Seurat, ScanPy or similar tools for visualization or downstream analysis.

Columns in the `<sample>.<libName>.allCells.csv` file:

| Column name | Description |
|---------|-------------|
| cell_id | Unique id for cell composed of barcodes detected |
| counts | The number of unique transcript molecules detected |
| genes | The number of unique genes detected |
| totalReads | The number of reads demultiplexed to this cell barcode |
| countedReads | The number of reads contributing to counts in the expression matrix |
| mappedReads | The number of reads that aligned to the genome |
| geneReads | The number of reads that mapped to an annotated gene |
| exonReads | The number of reads that mapped to an exon |
| intronReads | The number of reads that mapped to an intron | 
| antisenseReads | The number of reads mapping antisense to annotated exons | 
| mitoReads | The number of reads mapping on mitochondrial genome
| countedMultiGeneReads | The number of multi-gene reads that contributed to counts in the expression matrix
| Saturation | `1 - (UniqueReads / TotalReads)` on _Reads Mapped to Transcriptome_ |
| mitoProp | Proportion of mapped reads that aligned to mitochondrial genome |
| PCR | The alias for the PCR (library) barcode |
| PBC | The RT well |
| bead_bc | The bead barcode (microwell) |
| bead_bcs_total | The number of cell-barcodes that are associated with this bead |
| bead_bcs_pass | The number of passing cells associated with this bead |
| sample | The sample name |
| flags | QC flags |