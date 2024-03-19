# ScaleBio tools
In addition to third-party and open-source software the workflow also uses executable tools developed by ScaleBio:

- bc_parser
	- Extracts and error corrects cell-barcodes and UMIs from the original (input) fastq files
	- Splits (demultiplexes) the input fastq files into sample fastq files based on cell-barcodes (RT wells)
	- Barcode and read-level metrics

## Installation
These tools are included in the `scalerna` docker container image; when running the workflow with `-profile docker` (or another container engine, such as _singularity_, _podman_ etc.), they will be automatically available.

 These tools are however currently not available through Conda, so if running without a container system, they need to be installed first. A download script to get static pre-compiled binaries for linux (*x86_64*) is included. Simply run:

```/PATH/TO/ScaleRNA/envs/download-scale-tools.sh```
 
  This will install the binaries in `ScaleRNA/bin` (inside the Nextflow workflow directory), from where they will be available during workflow execution.

# bc_parser
The `bc_parser` output mode we are using for the RNA workflow puts all barcode information in the *_BC fastq file and the transcript sequence in *_R2. That makes the output easily compatible with STARSolo (and other tools). Specifically, the sequence in *_BC.fastq is 
`Ligation_BC (9bp) + RT_BC + i7(PCR) + UMI`
