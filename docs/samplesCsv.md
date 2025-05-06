# Sample table

A sample barcode table file is used to list the samples included in an analysis run and optional sample-specific analysis parameters. Each _sample_ in the workflow is defined by its set of well positions on the RT plate. For example, cells from Donor1 were fixed and loaded into the first column of the RT plate. That first column, 1A-1H, would then define the "Donor1" sample in the sample barcode table. 

A set of examples for the different Quantum RNA configurations is included in [docs/examples/quantum-sample-barcode-tables](examples/quantum-sample-barcode-tables). It is a comma separated file (csv), with a header line (column names), followed by one sample per line. 
The first column is required to be `sample` and contains the name of each sample. All other columns are optional.

 Column | Description | Example
:---- | ---- | :----:
sample | Sample name | PBMC-2
barcodes | RT-plate wells used for this sample | 1A-2H
libIndex2 | PCR indices (semicolon separated) used for this experiment (optional) | QSR-1;QSR-2
expectedCells | Approximate number of cells in this sample (optional) | 50000

*Note* 

* `sample` should consist only of letters, numbers and dashes (`-`) and start with a letter
* `barcodes` is optional if only a single sample was run. In that case all RT wells will be grouped together
* `libIndex2` is optional. If it is left out, the workflow will search for all possible PCR indicies in the input data
  * For data generated with the Scale RNA v1.1 kit, users should provide a `libIndex` column instead; for example "RNA-A-AP1"
* `expectedCells` is optional. If it is left out or set to 0, the number will be estimated from the read count distribution

## Sample Barcodes
If multiple samples were included in the experiment, they are effectively pooled in the sequencing library, and the analysis workflow will demultiplex them based on the sample (RT) barcodes. E.g

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
