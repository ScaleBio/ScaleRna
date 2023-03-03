# Fastq generation
The workflow can be started from a sequencer runFolder (bcl) [`--runFolder`]. In that case fastq files will be generated internally, using Illumina `bcl-convert`.

Alternatively it is also possible to generate fastq files upstream, for example when the ScaleRNA library is multiplexed together with other libraries during sequencing at a core facility. In that case we recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html), but older versions of the Illumina software are also possible.

An example [samplesheet.csv](examples/samplesheet.csv) with typical options is included. Here all 96 i7 barcode sequences from the PCR plate are merged into one set of fastq files.

## Index reads
For ScaleBio RNA libraries the RT and ligation barcodes are included in read 1, while the PCR barcode is in the index read. Hence we need to tell `bcl-convert` to generate index read fastqs using the `samplesheet.csv` setting: \
`CreateFastqForIndexReads,1`

# Using pre-generated fastq files as workflow input
Set `--fastqDir` to the directory containing the fastq files. 
The file names should follow the pattern `<LibraryName>_..._<Read>_...fastq.gz`, where
* `Name` is the library name (`ScaleRNA` by default, can be set in the `libName` column in `samples.csv`)
* `Read` is one of `R1`, `R2`, `I1`