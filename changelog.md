# Version 2.1
## 2.1.0
* Cell-calling is done for each sub-library separately, improving performance for large runs
* Changes to CellFinder to align more closely with EmptyDrops
  	- Cell-barcodes far below the median unique transcript count (`medianFraction`) are never called cells
* Beads exposed to a large ambient RNA signal are filtered by default
* Ultima .cram input supported without fixed filename pattern
* Short reads after adapter and Poly-A trimming are not longer reported in the library QC reports
 	- These reads are included in "Total Sample Reads" in the sample reports
* Saturation is always calculated on all transcriptome reads
  	- Unique and multimappers and both genomes for barnyard samples
* Barcode naming changes
	- The barcode from the RT plate is called `RT` and no longer `pbc` barcode
	- Bead barcode blocks are `01`-`96`, instead of `1A`-`12H`
* The workflow now checks minimal compute resource requirements after start-up
* Updated example `samplesheet.csv` [files](docs/examples/fastq-generation/)
	- The index orientation needed for NextSeq 2000 and NovaSeq X is now the default, rather than `_revComp`

# Version 2.0
## 2.0.1
* Remove trimming of cDNA (RNA) read to max 82bp length
* Make QuantumScale RNA the default library type
## 2.0.0
* Add support for QuantumScale RNA data
* Run cutadapt upstream of barcode demux to do adapter trimming
* Output unaligned BAM files post barcode demux, that are then sent to STAR as input
* Perform read length filtering during barcode demux
* Performance optimizations on ScalePlex assignment and reporting
* For barnyard experiments cell calling is now performed on a species-by-species basis
* Species specific metrics are now reported in the sample report for barnyard experiments
* Added sF, gn, and gx tags to BAM file output by STAR
* Updated read metric definitions, e.g. "Reads mapped to Transcriptome"; see [RNA QC reports](docs/qcReport.md)

# Version 1.6
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
