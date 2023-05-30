# Outputs

Analysis workflow results for all samples in the run will appear in the output directory defined by `--outDir` (`ScaleRna.out` by default). 
For detailed information about the library and sample level QC reports see [qcReport.md](qcReport.md)

## Key output files
| Directory | File | Description |
|-----------|------|-------------|
| `reports` | `<sample>.report.html` | A standalone report including key QC metrics and figures for each sample |
| | `<sample>.reportStatistics.csv` | Summary and QC metrics for this sample in csv format |
| | `library_<libName>.report.html` | Barcode summary and demultiplexing statistics for the whole sequencing library (potentially multiple samples) |
| | `multiqc_report.html` | [MultiQC](https://multiqc.info/) report for fastq generation, fastQC and trimming |
|  `samples` | `<sample>.filtered/` | Pre-filtered gene expression matrix for cells above the unique read threshold |
| | `<sample>.allCells.csv` | Metrics per cell-barcode, including barcodes / well positions
| `fastq` | `fastqc/*_fastqc.html` | [fastqc](https://github.com/s-andrews/FastQC) report for each fastq file in the sequencing library |
| | `Reports/` | Fastq generation summary reports from [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) |
| `barcodes/<libName>.demux` | `<sample>.fastq.gz` | Sample fastq files (Demultiplexed and barcode error-corrected); only included with `--fastqOut true` |
| `alignment` | `<sample>.star.align/` | [STAR](https://github.com/alexdobin/STAR) alignment output, including BAM file, with single-cell barcode and UMI information in tags
|  | `<sample>.star.solo/` | [STARSolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) output for each sample


The gene expression matrix directories (.mtx) can be loaded into Seurat, ScanPy or similar tools for visualization or downstream analysis.