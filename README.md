# ScaleBio Seq Suite: RNA Workflow

This is a Nextflow workflow to run analysis of ScaleBio Single Cell RNA Sequencing libraries. It processes data from sequencing reads to alignments, single-cell outputs (gene-count matrix, etc.), and QC reports.

## Getting started
* First install [Nextflow](http://www.nextflow.io) (22.04 or later)
* Download this workflow to your machine
* Install [dependencies](docs/dependencies.md)
* Launch the small pipeline [test run](#workflow-test)
* Download / configure a reference [genome](docs/genomes.md) for your samples
* Create a [samples.csv](docs/samplesCsv.md) table for your samples
* Create [runParams.yml](docs/analysisParameters.md), specifying inputs and analysis options for your run
* Launch the workflow for your run

## Inputs
* Sequencing reads
    * Path to the Illumina Sequencer RunFolder (bcl files)
    * If you prefer to start from fastq files, generated outside (before) this workflow, see [Fastq generation](docs/fastqGeneration.md).
* Sample Table
    * A .csv file listing all samples in the analysis, optionally split by RT barcode. See [samples.csv](docs/samplesCsv.md).
* Reference Genome
    * The workflow requires a reference genome, including a [STAR](https://github.com/alexdobin/STAR) index for alignment, and gene annotation. See [Reference Genomes](docs/genomes.md)

## Outputs
The workflow produces per-sample and per-library QC reports (`html`), alignments (`bam`), a cell-by-gene count-matrix (`mtx`) and more; See [Outputs](docs/outputs.md) for a full list.


## Workflow Execution
### Workflow test
A small test run, with all input data stored online or in this repository, can be run with the following command. Note that this test run is merely a quick and easy way to verify that the pipeline executes properly and does not represent a real assay: 

`nextflow run /PATH/TO/ScaleRna -profile PROFILE -params-file /PATH/TO/ScaleRna/docs/examples/runParams.yml --outDir output`

See [dependencies](docs/dependencies.md) for the best `PROFILE` to use on your system.

Note that this run will automatically download the example data from the internet (AWS S3), so please ensure that the compute nodes have internet access and storage space. Alternatively you can manually download the data first (using [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-welcome.html)):
```
aws s3 sync s3://scale.pub/testData/rna/202308_tinyPipelineTest/fastqs/ fastqs --no-sign-request
aws s3 sync s3://scale.pub/testData/rna/GRCh38_chr1_genome GRCh38_chr1_genome --no-sign-request
```
and then run with
```
nextflow run /PATH/TO/ScaleRna/ -profile docker,singularity --samples /PATH/TO/ScaleRna/docs/examples/samples.csv --genome GRCh38_chr1_genome/grch38.chr1.json --fastqDir fastqs --outDir /PATH/TO/OUTPUT_DIR
```

### Nextflow Command-line
**Note** that `nextflow` options are given with a single `-` (e.g. `-profile`), while workflow parameters (e.g. `--outDir`) are given with a double dash `--`.

See the [Nextflow command-line documentation](https://www.nextflow.io/docs/latest/cli.html) for the options to run `nextflow` on different systems (including HPC clusters and cloud compute).

## Configuration
### Specifying Analysis Parameters
Analysis parameters (inputs, options, etc.) can be defined either in a [runParams.yml](docs/examples/runParams.yml) file or directly on the nextflow command-line. See [analysisParameters](docs/analysisParameters.md) for details on the options.

### Config File
In addition to the analysis parameters, a user-specific nextflow configuration file can be used for system settings (compute and storage resources, resource limits, storage paths, etc.):

`-c path/to/user.config`

See [Nextflow configuration](https://www.nextflow.io/docs/latest/config.html) for the way different configuration files, parameter files and the command-line interact.

## Dependency Management
Different options to provide all required dependencies are described [here](docs/dependencies.md). Follow one approach there and then run nextflow with the corresponding `-profile`.

## Running in the cloud
Nextflow itself supports running using [AWS](https://www.nextflow.io/docs/latest/aws.html), [Azure](https://www.nextflow.io/docs/latest/azure.html) and [Google Cloud](https://www.nextflow.io/docs/latest/google.html). 

In addition [Nextflow tower](https://tower.nf) offers another simple way to manage and execute nextflow workflows in Amazon AWS.

# Versions and Updates
See the [Change log](changelog.md)

# License
For License information see [LICENSE.md](LICENSE.md).

By purchasing product(s) and downloading the software product(s) of ScaleBio, You accept all of the terms of the [License Agreement](LICENSE.md). If You do not agree to these terms and conditions, You may not use or download any of the software product(s) of ScaleBio. 


