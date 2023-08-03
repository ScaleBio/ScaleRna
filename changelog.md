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
