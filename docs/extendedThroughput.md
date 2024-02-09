# Extended Throughput Kit v1.1

When using additional final distribution (PCR) plates from the Extended Throughput Kit v1.1, each plate has the same set of 96 i5-index sequences, but receives a different set of i7 index sequences from the P7 primer pool used (`RNA-A-AP1` to `RNA-A-AP4`). This data is analyzed by generating separate "fastq libraries" for each plate (i7 pool), analyzing each separately through the core workflow (_barcodeParser_, _STARSolo_, etc.) and then combining the outputs at the single-cell level (gene expression matrix and metrics).

Cells from different plates are distinguished by appending the library-name (P7 primer pool) to the cell barcode in the gene expression matrix outputs.

## Single Sequencing Run
If all plates are sequenced together in the same run, this can be done in one nextflow analysis run:
* Create a `samples.csv` for all plates and samples
    * List each sample on each plate (repeating the `sample` name for each plate used)
    * Set the `libIndex` column to the name of the barcode pool (e.g. `RNA-A-AP1`)
    * See [samples.ext.csv](examples/extended-throughput/samples.ext.csv) for an example with two samples on two plates
* Launch the workflow with the `--merge` option
    * This will produce outputs for each individual plate as well as merged outputs, combining all cells for each sample across all plates

A small test-case can be run with: \
`nextflow run /PATH/TO/ScaleRna -profile PROFILE -params-file /PATH/TO/ScaleRna/docs/examples/extended-throughput/runParams.yml --outDir output.ext`

## Merging Across Multiple Runs
If plates were sequenced on separate sequencing runs, each plate can first be analyzed individually and the results later merged in a separate nextflow run
* Create a `samples.csv` for each separate plate
    * Including the `libIndex` column with the P7 pool used for that plate
* Run the nextflow analysis workflow separately for each plate
* Create another `samples.csv` that lists all samples and plates to be merged
    * Add a new `resultDir` column that points to the workflow output directories for the individual plate runs
    * See [samples.ext-merge.csv](examples/extended-throughput/samples.ext-merge.csv)
* Run the workflow with the `--reporting` and `--merge` option
    * `--reporting` uses the previous alignment results instead of re-processing all reads
    * `--merge` creates merged outputs for each sample across all plates as before
    * Do not specify any input reads at this step (no `--runFolder` or `--fastqDir`)

## Fastq Generation
To run this workflow starting from fastq files, rather than a runFolder / BCLs, see [Fastq Generation](fastqGeneration.md)
