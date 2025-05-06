# ScaleBio Seq Suite: RNA Workflow

This is a Nextflow workflow to run analysis of ScaleBio Single Cell RNA Sequencing libraries. It processes data from sequencing reads to alignments, single-cell outputs (gene-expression matrix, etc.), and QC reports.

The QuantumScale RNA assay is supported with version 2.0 and later.

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

## Required Inputs
* Sequencing reads (one of the following)
    * An llumina sequencer RunFolder (`.bcl` files)
    * A directory with `.fastq` files including index reads; see [Fastq generation](docs/fastqGeneration.md).
    * A directory containing unaligned cram files from Ultima Genomics sequencing
* Sample table
    * A `.csv` file listing all samples for this analysis run, optionally split by RT barcode. See [samples.csv](docs/samplesCsv.md).
* Reference genome
    * The workflow requires a reference genome, including a [STAR](https://github.com/alexdobin/STAR) index for alignment, and gene annotation. See [Reference Genomes](docs/genomes.md)
* Kit version / Library structure definition
    * If using something other than Quantum RNA, see [Analysis Parameters](docs/analysisParameters.md#kit-version)

## Outputs
The workflow produces per-sample and per-library ScaleBio demultiplexed reads (`fastqs`), alignments (`bam`), QC reports (`html`), a cell-by-gene count-matrix (`mtx`) and more; See [Outputs](docs/outputs.md) for a full list.

## ScalePlex
To analyze your ScaleRNA sequencing data with ScalePlex labeling and assignment, see [ScalePlex](docs/scalePlex.md)

## Workflow Execution
### Workflow test
A small test run, with all input data stored online, can be run with the following command:

`nextflow run /PATH/TO/ScaleRna -profile PROFILE -params-file /PATH/TO/ScaleRna/docs/examples/runParams.yml --outDir output`

`-profile docker` is the preferred option if the system supports _Docker_ containers;  See [Dependency Management](#dependency-management) for alternatives.

With this command, nextflow will automatically download the input data from the internet (AWS S3), so please ensure that the compute nodes have internet access and storage space. Alternatively you can manually download the data first (using [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html)):
```
aws s3 sync s3://scale.pub/testData/rna/202504_quantumV1/fastq fastq --no-sign-request
aws s3 sync s3://scale.pub/testData/rna/GRCh38_chr20_genome GRCh38_chr20_genome --no-sign-request
```
and then run with
```
nextflow run /PATH/TO/ScaleRna/ -profile PROFILE -params-file /PATH/TO/ScaleRna/docs/examples/runParams.yml --genome GRCh38_chr1_genome/grch38.chr1.json --fastqDir fastqs --outDir /PATH/TO/OUTPUT_DIR
```

Note that this test run is merely a quick and easy way to verify that the pipeline executes properly and does not represent a complete or realistic dataset.

### Nextflow command line (CLI)
`nextflow` options are given with a single `-` (e.g. `-profile`), while workflow parameters (e.g. `--outDir`) are given with a double dash `--`. See the [Nextflow command line documentation](https://www.nextflow.io/docs/latest/cli.html) for further details.

A typical command to run the workflow would be:

```
nextflow run /PATH/TO/ScaleRna/ -profile docker --samples samples.csv --genome /PATH/TO/GRCh38/grch38.json --runFolder /PATH/TO/230830_A00525_1087_BHLFF3DSX7/ --outDir output
```

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
Nextflow itself supports execution on [AWS](https://www.nextflow.io/docs/latest/aws.html), [Azure](https://www.nextflow.io/docs/latest/azure.html), and [Google Cloud](https://www.nextflow.io/docs/latest/google.html). 

In addition [Seqera Platform](https://seqera.io/platform/) offers another simple way to manage and execute nextflow workflows in Amazon AWS.

# Versions and Updates
See the [change log](changelog.md)

# License
By purchasing product(s) and downloading the software product(s) of ScaleBio, You accept all of the terms of the [License Agreement](LICENSE.md). If You do not agree to these terms and conditions, You may not use or download any of the software product(s) of ScaleBio.
