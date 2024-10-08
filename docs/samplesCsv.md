# Sample table

A sample table (e.g. [samples.csv](examples/samples.csv)) file is used to list the samples included in an analysis run, their sample barcode (RT) sequences and optional sample-specific analysis parameters.

It is a comma separated file (csv), with a header line (column names), followed by one sample per line. 
The first column is required to be `sample` and contains the name of each sample. All other columns are optional. _Column names are case sensitive!_

 Column | Description | Example
:---- | ---- | :----:
sample | Sample name | Foobar-2
barcodes | RT-plate wells used for this sample (optional) | 1A-2H
libName | Name for the overall sequencing library / fastq files (optional) | ScaleRna
expectedCells | Approximate number of cells in this sample (optional) | 50000

* `sample` and `libName` should consist only of letters, numbers and dashes (`-`) and start with a letter.
* A single `libName` should be used for all samples in a ScaleRNA sequencing library (one final distribution plate); not a different one per sample.
    * When running from fastq file input, `libName` should match the first part of the fastq file name for this library, e.g.: `ScaleRNA` for `ScaleRNA_*.fastq.gz`.
    * See [samples.ext.csv](examples/extended-throughput/samples.ext.csv) for an example with two libraries for two extended-throughput plates
* `expectedCells` is optional. If it is left out or set to 0, the number will be estimated from the read count distribution.

## Demultiplexing samples
During analysis the sequencing data is first converted into library fastq files (`libName` column). If multiple samples were included in one sequencing library, these are then demultiplexed based on the sample (RT) barcodes. E.g.

sample | barcodes
-- | --
Foo | 1A-6H
Bar | 7A-12H

The RT wells used for each sample are given in `barcodes` as either
* An individual value (`1A`)
* A range of wells (`1A-2H`)
    * Wells are sorted first by number then by letter, i.e. `1A-1H`.
    * Note that all ranges are read in **column-wise** order; e.g. 1A-2C, refers to 1A-1H (all of column 1) plus 2A-2C.
* A list of values or ranges, separated by semicolon (`;`) (`1A;2A-2D`)
