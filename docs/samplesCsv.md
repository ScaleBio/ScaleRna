# Sample table

A sample table (e.g. [samples.csv](examples/samples.csv)) file is used to list the samples included in an analysis run, their sample barcode (RT) sequences and optional sample-specific analysis parameters.

It is a comma separated file (csv), with a header line (column names), followed by one sample per line. 
The first column is required to be `sample` and contains the name of each sample. All other columns are optional.

 Column | Description | Example
:---- | ---- | :----:
sample | Sample name | Foobar-2
barcodes | RT-plate wells used for this sample (optional) | 1A-2H
libIndex2 | PCR indices (semicolon separated) used for this experiment (optional) | QSR-1;QSR-2
expectedCells | Approximate number of cells in this sample (optional) | 50000

* `sample` should consist only of letters, numbers and dashes (`-`) and start with a letter.
* `libIndex2` is optional. If it is left out, all PCR indices will be defaulted in, and the workflow will deduce which ones are present.
* `expectedCells` is optional. If it is left out or set to 0, the number will be estimated from the read count distribution.

## Demultiplexing samples
During analysis the sequencing data is first converted into library fastq files. If multiple samples were included in one sequencing library, these are then demultiplexed based on the sample (RT) barcodes. E.g.

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
