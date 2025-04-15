# Analysis parameters

All analysis parameters can be set in a [runParams.yml](../docs/examples/runParams.yml) file, which is then passed to nextflow with `nextflow run -params-file runParams.yml`. 
Alternatively each parameter in this file can also be set on the nextflow command-line directly, overwriting the value in the parameter file. E.g.
`nextflow run --samples=samples.foo.csv`

*Note* that `nextflow` options such as `resume` or `params-file` are given with a single `-`, while workflow parameters (e.g. `samples`) are given with a double dash `--`.


## Inputs
The workflow can either start from an Illumina sequencer run folder (bcl files) or a directory with fastq files. Specify either
* runFolder : "path/to/runFolder" <br>
OR
* fastqDir : "path/to/fastqs" <br>
OR
* ultimaCramDir : "path/to/crams"

where `fastqDir` is a directory containing all input fastq files. See [Fastq Generation](fastqGeneration.md) for details on file names, etc.

`runFolder` is the top-level directory for a sequencer output (containg `RunInfo.xml`). The workflow uses Illumina [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) for automatic fastq generation.

`ultimaCramDir` is a directory containing all input unaligned cram files generated from an Ultima Genomics sequencer. The cram filenames need to follow the regular expression captured in [params.cramFilePattern](../modules/ultima.config)
### Sample Information
* samples : "samples.csv"

A [file](examples/samples.csv) listing all samples in the analysis with their names, sample barcodes (RT) and optional sample settings

### Reference Genome
* genome : "/genomes/grch38/genome.json"

Path to a [genome.json](genomes.md) file that contains the location of all sequence and index files as well as other parameters for the reference genome to use. 

### Kit Version
* libStructure: `libQuantumV1.0.json` or `libV1.json` or `libV1.1.json`

To analyze data generated with version 1.0 of the Quantum ScaleRNA kit use `libQuantumV1.0.json`. For version 1.0 of the RNA kit (i7-indexed PCR plates) use `libV1.json`. For the v1.1 kit and extended throughput kit (i5-indexed PCR plates), use `libV1.1.json`. 

## Optional and Advanced Parameters
See [nextflow.config](../nextflow.config) for a list of all available parameters. The file also includes nextflow system options (compute resource requirements, etc.).

### BAM Output
Setting `bamOut` to false will suppress alignment (.bam file) output from STAR. Gene expression results (.mtx), and all other workflow outputs, are still generated of course.
If alignments are not required for custom downstream analysis, disabling BAM output will save compute time and storage space.

### Parallel Execution for Large Datasets
The workflow is executed with extra parallelism by setting the `--splitFastq` parameter. This is on by default and is recommended for large datasets. Parallelization is implemented at different stages of the workflow's execution:
* If starting from a RunFolder (BCL): If using the v1.0 or v1.1 RNA kit, the workflow automatically splits the data into 96 fastq files (one per PCR barcode). For the v1.0 QuantumScale kit, the data is split into 384 fastq files (one file per combination of PCR and bead barcode)
* All sets of fastq files (`R1`, `R2`, `I1`, `I2`), either from --`fastqDir` or built-in bcl-convert, are processed in parallel through bcParser.
* After barcode parsing, the data is split by the sample barcode barcode through alignment.
* Split samples are merged after alignment (gene expression matrix, metric files, etc.).

*Note*: When starting from fastq files, the first workflow steps are parallelized per set of input fastq files, so it is important to have the data in multiple input fastq files for better performance. The easiest way is to use a `samplesheet.csv` to split the data by both indexes, `index2` (_I5_) and `index1` (_I7_), see [fastqGeneration](fastqGeneration.md)


### Library Structure Definition
The library structure JSON file defines 
* Where in the reads cell-barcodes are found
* What the list of allowed barcode sequences is
* Which parts of the reads represent cDNA for alignment

The workflow includes default files for version 1.0 of the Quantum ScaleRNA kit, and version 1.0 and 1.1 of the RNA kit [references](../references/).
