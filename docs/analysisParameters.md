# Analysis parameters

All analysis parameters can be set in a `runParams.yml` file, which is then passed to nextflow with `nextflow run -params-file runParams.yml`. 
Alternatively each option is this file can also be set on the nextflow command-line directly, overwriting the value in the parameter file. E.g.
`nextflow run --samples=samples.foo.csv`

*Note* that `nextflow` options such as `resume` or `params-file` are given with a single `-`, while workflow parameters (e.g. `samples`) are given with a double dash `--`.


## Inputs
The workflow can either start from an Illumina sequencer runFolder (bcl files) or a directory with fastq files. Specify either
* runFolder : "path/to/runFolder" <br>
OR
* fastqDir : "path/to/fastqs"

where fastqDir is a directory containing all input fastq files. See [Fastq Generation](fastqGeneration.md) for details on file names, etc.

When starting from a sequencer run folder the workflow uses Illumina [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) for automatic fastq generation.

### Sample Information
* samples : "samples.csv"

A [file](examples/samples.csv) listing all samples in the analysis with their names, barcode sequences and optional sample settings

### Reference Genome
* genome : "/genomes/grch38/genome.json"

Path to a [genome.json](genomes.md) file that contains the location of all sequence and index files as well as other parameters for the reference genome to use. 

### Kit Version
* libStructure: `libV1.json` or `libV1.1.json`

To analyze data generated with version 1.0 of the RNA kit (i7-indexed PCR plates) use `libV1.json`. For the v1.1 kit and extended throughput kit (i5-indexed PCR plates), use `libV1.1.json`

## Optional and Advanced Parameters
See [nextflow.config](../nextflow.config) for a list of all available parameters. The file also includes nextflow system options (compute resource requirements, etc.).

### BAM Output
Setting `bamOut` to false will suppress alignment (.bam file) output from STAR. Gene expression results (.mtx), and all other workflow outputs, are still generated of course.
If the user does not specifically need alignments for custom downstream analysis, disabling BAM output will save compute time and storage space.

### Parallel Execution
The workflow can be executed with extra parallelism by setting the `--splitFastq` parameter. Based on what type of input is provided, the workflow has two different modes of parallelism:
* *RunFolder*: Samples are split by PCR barcode in the `bcl-convert` step. These 96 fastq files (per library/lane) are processed separately through barcode parsing, trimming and alignment.
* *Fastq*: Each set of input fastq files (e.g. by lane) is processed separately through barcode parsing. During that step, samples are split further based on the RT barcode. Each resulting file subset goes through trimming and alignment in parallel.

Split samples are merged after alignment (gene expression matrix, metric files, etc.).

### Library Structure Definition
The library structure JSON file defines 
* Where in the reads cell-barcodes are found
* What the list of allowed barcode sequences is
* Which parts of the reads represent genomic DNA and which should be masked (e.g. RT and ligation barcodes)

The workflow includes default files for version 1.0 and 1.1 of the RNA kit in [references](../references/).
