# Reporting-only workflow

A separate subworkflow, `reporting`, is included to run the downstream steps of the workflow, after alignment and quantification (after _STARsolo_). This includes
* Cell filtering (UMI threshold)
* Sample metric and report generation

This sub-workflow can be used to regenerate reports and filtered outputs after adjusting parameters without re-running alignment.

## Usage
The workflow can be started by running `nextflow run /PATH/TO/ScaleRna -profile ... --genome ... --samples samples.csv --reporting --resultDir <OutDir from previous pipeline run>

Where `samples.csv` includes a `resultDir` column, giving the base output directory from the preivous pipeline run for that sample. E.g.

| sample | resultDir |
| pbmc1 | /PATH/TO/RUN1/ScaleRna.out |
| pbmc2 | /PATH/TO/RUN2/ScaleRna.out |

## Inputs
The `reporting` workflow will read the (raw) *STARSolo* output from the previous pipeline run specified in `resultDir`. Specifically, `resultDir/alignment/<sample>.<libName>.star.solo`
It  also reads reference information from the `genome.json` and `library.json`

## Outputs
The `reporting` workflow produces
* A new filtered gene expression matrix (`samples/<sample>.<libName>.filtered/`)
* Cell metrics (`samples/<sample>.<libName>.allCells.csv`)
* New per-sample QC reports (`reports/<sample>/<sample>.<libName>.report.html`)
  - And `.csv` metric files

## Parameters
`--merge` Flag the determines wether or not to merge `<sample>.<libName>` libraries?
