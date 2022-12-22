# Sample table

A sample table (e.g. [samples.csv](examples/samples.csv)) file is used to list the samples included in an analysis run, their sample barcode (RT) sequences and optional sample-specific analysis parameters.

It is a comma separated file (csv), with a header line (column names), followed by one sample per line. 
The first column is required to be `sample` and contains the name of each sample. All other columns are optional. _Column names are case sensitive!_

 Column | Description | Example
:---- | ---- | :----:
sample | Sample name | Foobar-2
libName | Name for the sequencing library / fastq files | ScaleRna
barcodes | RT-plate wells used for this sample | 1A-2H
expectedCells | Approximate number of single cells in this sample | 50000

* `sample` and `libName` should consist only of letters, numbers, dash (-) and dot (.)
* When running from pre-existing fastq file input, `libName` should match the first part of the fastq file name for this sample, e.g.: `Foo1` for `Foo1_*.fastq.gz`.
* `expectedCells` is optional. If it is left out or set to 0, the number will be estimated from the read count distribution.

## Demultiplexing samples
During analysis the sequencing data is first converted into library fastq files (`libName` column). If multiple samples were included in one sequencing library, these are then be demultiplexed based on the sample (RT) barcodes. In that case the same library is repeated on multiple lines, once for each sample, with the specific tn5-wells for each listed in `barcodes`. E.g.

sample | barcodes
-- | --
Foo | 1A-6H
Bar | 7A-12H

The RT wells used for each sample are given in `barcodes` as either
* An individual value (`1A`)
* A range of wells (`1A-2H`)
    * Wells are sorted first by number then by letter, i.e. `1A-1H`.
    * Note that all ranges are read in column-wise order; e.g. 1A-2C, refers to 1A-1H (all of column 1) plus 2A-2C.
* A list of values or ranges, separated by semicolon (`;`) (`1A;2A-2D`)