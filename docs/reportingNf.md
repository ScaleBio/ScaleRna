# Reporting-only Workflow

A separate sub-workflow, `reporting`, is included to run only the QC and reporting steps of the workflow, after alignment and quantification (i.e. after _STARsolo_). This includes
* Cell filtering (UMI threshold)
* Sample metric and report generation
* Library metric and report generation

This sub-workflow can be used to regenerate outputs after adjusting parameters, without re-running alignment. It is also used to combine (merge) results from multiple different extended throughput plates (libraries), sequenced on different sequencing runs; see [extendedThroughput](extendedThroughput.md).

## Usage
The reporting sub-workflow can be started by running

`nextflow run /PATH/TO/ScaleRna -profile ... --genome ... --samples samples.csv --reporting --resultDir <outDir from previous workflow run>`

where `samples.csv` and `genome` should be the same files as used in the previous (alignment) run of the workflow.


## Inputs
The `reporting` workflow will read the (raw) *STARSolo* output and a *STARSolo* log file from the previous pipeline run specified in `resultDir`. Specifically, the raw output is read from `resultDir/alignment/<sample>.<libName>.star.solo` and the log file is read from `resultDir/alignment/<sample>.<libName>.star.align`
It  also reads reference information from the `genome.json` and `library.json`

## Outputs
The `reporting` workflow produces
* A new filtered gene expression matrix (`samples/<sample>.<libName>.filtered/`)
* Cell metrics (`samples/<sample>.<libName>.allCells.csv`)
* New per-sample QC reports (`reports/<sample>/<sample>.<libName>.report.html`)
  - And `.csv` metric files

Currently the reporting sub-workflow does not aggregate read-trimming statistics from the original run, so the _Total Sample Reads_ (pre-trimming) and _Average Trimmed Read Length_ metrics will be missing from the re-generated reports.


## Combining Analysis Runs
The reporting sub-workflow can be used to combine results from multiple analysis workflow runs. To do this, multiple paths to previous workflow outputs can be specified in a `resultDir` column in `samples.csv`, instead of the `--resultDir` command-line option:

| sample | resultDir |
|--------|----------------------------|
| pbmc1  | /PATH/TO/RUN1/ScaleRna.out |
| pbmc2  | /PATH/TO/RUN2/ScaleRna.out |

See [samples.ext-merge.csv](examples/extended-throughput/samples.ext-merge.csv) for a complete example combining cells from multiple extended throughput plates. If `--merge` is set (default), the workflow produces combined outputs for each sample across all libraries (plates).

*Note* that this function only combines cells or samples from different runs, not reads for the same set of cells. Hence it cannot be used to combine two sequencing runs of the same library (multiple sets of reads for the same cells). Those need to be combined at the fastq level and go through alignment together in order to detect duplicate reads across fastq files.

