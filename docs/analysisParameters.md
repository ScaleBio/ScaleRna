# Analysis parameters

All analysis parameters can be set in a [runParams.yml](../docs/examples/runParams.yml) file, which is then passed to nextflow with `nextflow run -params-file runParams.yml`. 
Alternatively each parameter in this file can also be set on the nextflow command-line directly, e.g.
`nextflow run --samples=samples.foo.csv`

*Note* that `nextflow` options such as `resume` or `params-file` are given with a single `-`, while workflow parameters (e.g. `samples`) are given with a double dash `--`. If both are used, command-line options overwrite values in the parameter file.


## Inputs
### Sequencing reads
The alignment workflow can either start from an Illumina sequencer run folder (_bcl files_) or a directory with sequencing read files (_fastq_). Specify one of these options:
* runFolder : "path/to/runFolder"
* fastqDir : "path/to/fastqs"
* ultimaCramDir : "path/to/crams"

`fastqDir` is a directory containing all input fastq files for this analysis. See [Fastq Generation](fastqGeneration.md) for details on file names, etc.

`runFolder` is the top-level directory for a sequencer output (containing `RunInfo.xml`). The workflow uses Illumina [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) for automatic fastq generation.

`ultimaCramDir` is a directory containing all input unaligned cram files pre-processed with Ultima `trimmer`. The cram filenames need to follow the regular expression captured in [params.cramFilePattern](../modules/ultima.config)

#### Reporting only runs
If the sequencing data has been analyzed previously, these outputs can be re-used to generate updated reports. In this case previous workflow outputs are re-used with `resultDir`; see [reportingNf](reportingNf.md).

### Sample Information
* samples : "samples.csv"

A [file](examples/samples.csv) listing all samples in the analysis with their names, sample barcodes (RT) and optional sample settings

### Reference Genome
* genome : "/genomes/grch38/genome.json"

Path to a [genome.json](genomes.md) file that contains the location of all sequence and index files as well as other parameters for the reference genome to use. 

### Scale Bio RNA Kit Version
* libStructure: "libQuantumV1.0.json"

This defines the version of the Scale Bio single-cell RNA kit used to generate the libraries. The default `libQuantumV1.0.json` matches version 1.0 of the Quantum ScaleRNA. For older Scale RNA (_3 level_) kits, set `quantum false` and `libStructure libV1.1.json` (Scale RNA v1.1 kit and extended throughput kit). 

## Outputs
* outDir: "ScaleRna.out"
The name of the output directory for all analysis results. See [outputs](outputs.md) for a list of output files. Pre-existing files will be overwritten

### BAM Output
By default `bamOut` is set to false, which will suppress alignment (.bam file) output from STAR. Gene expression results (.mtx), and all other workflow outputs, are still generated of course.
If alignments are not required for custom downstream analysis, disabling BAM output will save compute time and storage space.

## Optional and Advanced Parameters
See [nextflow.config](../nextflow.config) for a list of all available parameters and their default values. The file also includes nextflow system options (compute resource requirements, etc.).

### Multimapping Reads
The workflow uses STARsolo for alignments of reads to the genome and to transcripts. STARsolo handling of multimapping reads is controlled by these options
* starMaxLoci: 6
    * Reads that map to more genomic locations that filtered out, regardless of transcriptome overlap.
* starMulti: "PropUnique"
    * The algorithm used by handle reads that match multiple genes, either because of genomic multimapping or because of overlapping gene annotations; See [STARsolo documentation](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md#multi-gene-reads)
* roundCounts: "False"
    * If multi-gene reads are resolved, they lead to fractional unique transcript counts, e.g. "0.5" for a unique read that matches two genes equally well. If this setting is true, these counts are rounded to the closest full integer.

### Cell Calling
See [cellCalling.md](cellCalling.md) for a description of the method used to call cell barcodes against background and related options.

## Compute Resources
`splitFastq`, set to true by default, enables increased parallelization of the workflow by splitting the data based on input fastq files and RT barcodes. This can be set to false to reduce the number of compute jobs for small analysis.

`taskMaxMemory`, `taskMaxCpus` and `taskMaxTime` control the maximum amount of compute resources a single task within the workflow is allowed to reserve. These can be lowered to match the available resources or limits of a users compute cluster.

