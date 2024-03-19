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
| | `<sample>.<libName>.allCells.csv` | Metrics per cell-barcode, including barcodes / well positions
| | `<sample>_libraries/` | For extended throughput runs, this contains output files per plate
| `fastq` | `fastqc/*_fastqc.html` | [fastqc](https://github.com/s-andrews/FastQC) report for each fastq file in the sequencing library |
| | `Reports/` | Fastq generation summary reports from [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) |
| `barcodes/<libName>.demux` | `<sample>.fastq.gz` | Sample fastq files (Demultiplexed and barcode error-corrected); only included with `--fastqOut true` |
| `alignment/<sample>.<libName>` | `<sample>.star.align/` | [STAR](https://github.com/alexdobin/STAR) alignment output, including BAM file, with single-cell barcode and UMI information in tags
|  | `<sample>.star.solo/` | [STARSolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) output for each sample, including unfiltered `.mtx` files


The gene expression matrix directories (.mtx) can be loaded into Seurat, ScanPy or similar tools for visualization or downstream analysis.