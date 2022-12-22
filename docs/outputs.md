# Outputs

Analysis workflow results for all samples in the run will appear in the output directory defined by `--outDir` (`ScaleRna.out` by default). 

## Key output files
| Directory                   | File                                            | Description                                                                                                                                                                                                                                                                                  |
|-----------------------------|-------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `reports`                   | `<sample>.report.html`                          | An interactive standalone HTML report including key metrics/figures for each sample |
| | `<sample>.reportStatistics.tsv`                 | Summary and QC metrics for this sample in tsv (text) format                                                                                                                   |
| | `library_<libName>.report.html`                 | Barcode summary and demultiplexing statistics for the whole library (potentially multiple samples) |
| `fastq/fastqc/`             | `<libName>_fastqc.html`                     | [fastqc](https://github.com/s-andrews/FastQC) report for the sequencing library|
| `demux/<libName>.demux` | `<sample>.fastq.gz`                             | Sample fastq files (Demultiplexed) |
| `star/<sample>` | `Aligned.sortedByCoord.out.bam` | Read alignments to genome, with single-cell barcode and UMI information in tags
|  | `Solo.out` | STARSolo output for each sample. See https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md |
|  `star/<sample>/Solo.out/GeneFull_Ex50pAS` | `raw/matrix.mtx` | Cell-gene expression matrix (for Seurat/Scanpy etc.) |