# Analysis parameters

All analysis parameters can be set in a `runParams.yml` file, which is then passed to nextflow with `nextflow run -params-file runParams.yml`. 
Alternatively each option is this file can also be set on the nextflow command-line directly, overwriting the value in the parameter file. E.g.
`nextflow run --samples=samples.foo.csv`

*Note* that `nextflow` options such as `resume` or `dump-channels` are given with a single `-`, while workflow parameters (e.g. `samples`) are given with a double dash `--`.


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

### Reference genome
* genome : "/genomes/grch38/genome.json"

Path to a [genome.json](genomes.md) file that contains the location of all sequence and index files as well as other parameters for the reference genome to use. 

## Optional and Advanced parameters
Run `nextflow run path/to/ScaleRna --help` for a description of available options and see the example [runParams.yml](examples/runParams.yml) file.

System options (compute resource requirements, etc.) as well as all parameter defaults, are in the workflow [nextflow.config](../nextflow.config).

### Selected optional parameters
Setting `bamOut` to false will suppress alignment (.bam file) output from STAR. Gene expression results (.mtx), and all other workflow outputs, are still generated of course.
If the user does not specifically need alignments for custom downstream analysis, disabling BAM output will save compute time and storage space.

#### Distributed execution
The workflow can be executed with extra parallelism by setting the `splitFastq` parameter to true. Based on what type of input is provided, the workflow has two different modes of parallelism
* If starting from runFolder (bcl files) and `splitFastq` is set to true
  1. Samples are split by PCR barcode in the bcl convert step. So we get 96 fastq files for each sample that go into barcode parsing, sample demux and alignment. The split samples are merged downstream
* If starting from fastq files and `splitFastq` is set to true
  1. Samples are split by RT barcode in the barcode parsing step, which gives us 96 fastq files for each sample

All of this parallelism is abstracted from the user. All the user needs to do is set the `splitFastq` option in [nextflow.config](../nextflow.config) to true.

### Library Structure Definition
* libStructure : "libV1.json"

The library structure JSON file defines 
* Where in the reads cell-barcodes are found
* What the list of allowed barcode sequences is
* Which parts of the reads represent genomic DNA and which should be masked (e.g. RT and ligation barcodes)

The default file, for our standard product configuration, is included in [references/libV1.json](../references/libV1.json).
