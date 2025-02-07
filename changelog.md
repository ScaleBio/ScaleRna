# Version 1.6
## 1.6.3
* Fix remote URL in sample and library QC reports (HTML)
## 1.6.2
* Fix erroneous channel dump in inputReads.nf that was triggered by Nextflow 24.10
## 1.6.1
* Fix unnecessary parameter warning
## 1.6.0
* Add support for ScalePlex library analysis and assignment in tandem with RNA analysis, enabled by the `--scaleplex` parameter 

# Version 1.5
* New optional _EmptyDrops_-like cell calling method (`--cellFinder`)
* Correct up to one `N` basecall per barcode, like other sequencing errors 
* Rename resource limit parameters (_taskMaxMemory_, _taskMaxCpus_, _taskMaxTime_)
* Added option to output an anndata object of the cell-count matrix (`--annData`)
* Updated optional clustering workflow and report (seurat v5-based)
* Changed threshold settings (_expectedCells_, _fixedCells_, _minUTC_)
* Separated `pass` and `flags` columns for cell-calling in `allCells.csv`
* Make `--merge` the default for extended-throughput plate runs
* Fixed mouse background calculation for barnyard runs
* Compute mitochondrial fraction based on mapped reads
* Changed barcode rank plot to support _CellFinder_ calls
* Simplified parallelism strategy
* Improved task error retry strategy
* Added more input parameter validations
* All ScaleBio docker containers moved to AWS ECR

# Version 1.4
## 1.4.1
* Update `datapane` version to remove need to access users home dir.
* Always output sequencing-depth intrapolation results as (rounded) integers
* Update workflow test dataset to RNA kit v1.1

## 1.4.0
* Add support for RNA kit v1.1 (i5-indexed plates)
    - `libStructure` argument is now required; "libV1.json" or "libV1.1.json"
* Add support for Extended Throughput kit (multi-plate merging)
* Add 'reporting workflow' to re-generate outputs from previous alignments
* Change output file name and directory structure to support multiple per-library and merged outputs
* Add analysis metadata to sample QC report
* Remove read2 (RNA) read-length trimming default (previously was 48bp)
* Changed column label for PCR barcode in `allcells.csv` (`PCR` to `i7`/`i5`)
* Fix line-count in the header of the merged, unfiltered STARSolo `.mtx` files

# Version 1.3
## 1.3.3
* Updated bcl-convert dependency to version 3.9.3
* New tiny pipeline test dataset
## 1.3.2
* Output gene expression matrix (mtx) from STAR with resolved multimappers included
	- See `--starMulti` option in `nextflow.config`
* Added trimmed reads metric to sample QC report
* Exclude invalid i7 barcodes during fastq generation (`samplesheet.csv`)
* Switch metric file outputs to `.csv`
* Option for increased parallelization with `--splitFastq`
* New platemap color scheme

# Version 1.2
* Enable Poly-T trimming
	- Reads are now trimmed at Poly-T, in addition to Poly-A stretches, to remove a few antisense mismapping artifacts
* Metrics re-definition
	- Many QC metric definitions are updated. See [qcReport.md](docs/qcReport.md)
* Reorganized output (`--outDir`) directory. See [outputs.md](docs/outputs.md)
* Fix conda environment specification
* Remove Index reads from fastQC
* Performance improvements

# Version 1.1
* Enable Poly-A trimming
	- By default all cDNA reads are now trimmed at the first occurrence of a 8bp Poly-A stretch
	- This avoids mismapping of RT primer sequences

# Version 1.0
Initial Release
