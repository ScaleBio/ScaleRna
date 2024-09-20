# ScaleBio Seq Suite: RNA Workflow

This is a Nextflow workflow to run analysis of ScaleBio Single Cell RNA Sequencing libraries. It processes data from sequencing reads to alignments, single-cell outputs (gene-expression matrix, etc.), and QC reports.

## Getting started
* First install [Nextflow](http://www.nextflow.io) (version 23.10 or later)
* Download this workflow to your machine
* Setup [dependencies](docs/dependencies.md)
* Launch the small pipeline [test run](#workflow-test)
* Download / configure a reference [genome](docs/genomes.md) for your samples
* Create a [samples.csv](docs/samplesCsv.md) table for your samples
* Create [runParams.yml](docs/analysisParameters.md), specifying inputs and analysis options for your run
* Launch the workflow for your run

## Requirements
* Linux system with GLIBC >= 2.17 (such as CentOS 7 or later)
* Java 11 or later
* 64GB of RAM and 12 CPU cores
    * Smaller datasets can be run with 32GB of RAM and 6 CPUs

## Required Inputs
* Sequencing reads
    * Path to the Illumina sequencer RunFolder (`.bcl` files)
    * To instead start the workflow from `.fastq` files, generated outside (before) this workflow, see [Fastq generation](docs/fastqGeneration.md).
* Sample table
    * A `.csv` file listing all samples for this analysis run, optionally split by RT barcode. See [samples.csv](docs/samplesCsv.md).
* Reference genome
    * The workflow requires a reference genome, including a [STAR](https://github.com/alexdobin/STAR) index for alignment, and gene annotation. See [Reference Genomes](docs/genomes.md)
* Kit version / Library structure
    * Select the `libStructure` corresponding to the version of the ScaleBio RNA kit used. See [Analysis Parameters](docs/analysisParameters.md#kit-version)
    * Version 1.0 requires `Read1`, `Read2` and `Index1 or i7`
    * Version 1.1 requires `Read1`, `Read2`, `Index1 or i7` and `Index2 or i5`

## Outputs
The workflow produces per-sample and per-library ScaleBio demultiplexed reads (`fastqs`), alignments (`bam`), QC reports (`html`), a cell-by-gene count-matrix (`mtx`) and more; See [Outputs](docs/outputs.md) for a full list.

## Extended Throughput Kit
To analyze data from multiple final distribution plates in the extended throughput kit, see [Extended Throughput](docs/extendedThroughput.md)

## ScalePlex
To analyze your ScaleRNA sequencing data with ScalePlex labeling and assignment, see [ScalePlex](docs/scalePlex.md)

## Workflow Execution
### Workflow test
A small test run, with all input data stored online, can be run with the following command:

`nextflow run /PATH/TO/ScaleRna -profile PROFILE -params-file /PATH/TO/ScaleRna/docs/examples/runParams.yml --outDir output`

`-profile docker` is the preferred option if the system supports _Docker_ containers;  See [Dependency Management](#dependency-management) for alternatives.

With this command, nextflow will automatically download the input data from the internet (AWS S3), so please ensure that the compute nodes have internet access and storage space. Alternatively you can manually download the data first (using [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html)):
```
aws s3 sync s3://scale.pub/testData/rna/202312_rnaV1.1/fastqs/ fastqs --no-sign-request
aws s3 sync s3://scale.pub/testData/rna/GRCh38_chr1_genome GRCh38_chr1_genome --no-sign-request
```
and then run with
```
nextflow run /PATH/TO/ScaleRna/ -profile PROFILE -params-file /PATH/TO/ScaleRna/docs/examples/runParams.yml --genome GRCh38_chr1_genome/grch38.chr1.json --fastqDir fastqs --outDir /PATH/TO/OUTPUT_DIR
```

Note that this test run is merely a quick and easy way to verify that the pipeline executes properly and does not represent a complete dataset.

### Nextflow Command-line
`nextflow` options are given with a single `-` (e.g. `-profile`), while workflow parameters (e.g. `--outDir`) are given with a double dash `--`. See the [Nextflow command-line documentation](https://www.nextflow.io/docs/latest/cli.html) for further details.

A typical command to run the workflow on new data would then be:

```
nextflow run /PATH/TO/ScaleRna/ -profile docker --samples samples.csv --genome /PATH/TO/GRCh38/grch38.json --runFolder /PATH/TO/230830_A00525_1087_BHLFF3DSX7/ --splitFastq --outDir output
```

For large datasets (e.g. NovaSeq runs), setting `--splitFastq` increases the amount of parallelization, which can significantly reduce analysis time. See [Analysis Parameters](docs/analysisParameters.md#parallel-execution)


## Configuration
### Specifying Analysis Parameters
Analysis parameters (inputs, options, etc.) can be defined either in a [runParams.yml](docs/examples/runParams.yml) file or directly on the nextflow command-line. See [Analysis Parameters](docs/analysisParameters.md) for details on the options.

### Config File
In addition to the analysis parameters, a user-specific nextflow configuration file can be used for system settings (compute and storage resources, resource limits, storage paths, etc.):

`-c path/to/user.config`

See [Nextflow configuration](https://www.nextflow.io/docs/latest/config.html) for the way different configuration files, parameter files and the command-line interact.

## Dependency Management
The Nextflow workflow can automatically use pre-built docker containers with all dependencies included. Activating the included `-profile docker` enables the required Nextflow settings. For details and alternatives see [Dependencies](docs/dependencies.md).

## Running in the cloud
Nextflow itself supports running using [AWS](https://www.nextflow.io/docs/latest/aws.html), [Azure](https://www.nextflow.io/docs/latest/azure.html) and [Google Cloud](https://www.nextflow.io/docs/latest/google.html). 

In addition [Nextflow tower](https://tower.nf) offers another simple way to manage and execute nextflow workflows in Amazon AWS.

# Versions and Updates
See the [change log](changelog.md)

# License
By purchasing product(s) and downloading the software product(s) of ScaleBio, You accept all of the terms of the [License Agreement](LICENSE.md). If You do not agree to these terms and conditions, You may not use or download any of the software product(s) of ScaleBio.
