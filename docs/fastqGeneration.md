# Fastq Generation
The workflow can be started from a sequencer runFolder (bcl) [`--runFolder`]. In that case fastq files will be generated internally, using Illumina `bcl-convert`.

Alternatively it is also possible to generate fastq files upstream using standard Illumina software, for example when the ScaleRNA library is multiplexed together with other libraries during sequencing at a core facility. We recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) version 3.9 or later.

## Samplesheet.csv
### RNA kit v1.0
An example [samplesheet.csv](examples/fastq-generation/samplesheet_v1.0.csv) with typical options is included. Here all 96 i7 barcode sequences from the PCR plate are merged into one set of fastq files. If an index2 (i5) read is used to demultiplex the ScaleBio RNA library with other libraries in the sequencing run, a `index2` column can be added with the constant i5 sequence of the ScaleRna library, which is `TGAACCTT[AC]` (8 or 10bp) in forward and `AAGGTTCA[GT]` in reverse orientation. See [samplesheet_v1.0_with_i5.csv](examples/fastq-generation/samplesheet_v1.0_with_i5.csv).

### RNA Kit v1.1 and Extended Throughput Kit
The v1.1 kit uses a different index read setup, one i5 index for each final distribution (PCR) plate well and P7-pools for each distinct plate. An example [samplesheet_v1.1.csv](examples/fastq-generation/samplesheet_v1.1.csv) is included. Some sequencers require the `index2` (i.e. `i5`) sequence in the opposite orientation; an example for that is included here: [samplesheet_v1.1_revComp.csv](examples/fastq-generation/samplesheet_v1.1_revComp.csv).

### Splitting Large Fastq Files
This creates a separate set of fastq files for each extended throughput plate (`RNA-A-AP1`, `RNA-A-AP2`, ...). See [Extended Throughput](extendedThroughput.md) for how to analyze these jointly.
To reduce workflow runtime through parallelization (`--splitFastq` mode), it is best to avoid combining the entire data in one fastq file. Generally the Illumina default of splitting per lane (i.e. not setting bcl-convert option `--no-lane-splitting`) is good. An alternative is to modify `samplesheet.csv` to assign a different `sample_ID` to each i5 barcode, following a pattern like `RNA-A-AP1_A1`, `RNA-A-AP1_A2`; as long as they all start with the same _libraryName_ followed by `_`, they will be grouped together for analysis in the ScaleRNA workflow.

## Index Reads
For ScaleBio RNA libraries the RT and ligation barcodes are included in read 1, while the PCR barcode is in the index read. Hence we need to tell `bcl-convert` to generate index read fastqs using the `samplesheet.csv` setting: \
`CreateFastqForIndexReads,1`

# Using Fastq Files as Workflow Input
Set `--fastqDir` to the directory containing the fastq files. 
The file names should follow the pattern `<LibraryName>_..._<Read>_...fastq.gz`, where
* `Name` is the library name (`ScaleRNA` by default, can be set in the `libName` or `libIndex` column in `samples.csv`)
* `Read` is one of `R1`, `R2`, `I1`

This can be achieved by matching the `sample_ID` in the Illumina _bcl-convert_ _samplesheet.csv_ to the `libName` in the ScaleBio `samples.csv`.

