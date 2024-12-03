# Fastq Generation
The workflow can be started from a sequencer runFolder (bcl) [`--runFolder`]. In that case fastq files will be generated internally, using Illumina `bcl-convert`.

Alternatively it is also possible to generate fastq files upstream using standard Illumina software, for example when the ScaleRNA library is multiplexed together with other libraries during sequencing at a core facility. We recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) version 3.9 or later.

## Samplesheet.csv
### RNA kit v1.0
An example [samplesheet.csv](examples/fastq-generation/samplesheet_v1.0.csv) with typical options is included. Here all 96 i7 barcode sequences from the PCR plate are merged into one set of fastq files. If an index2 (i5) read is used to demultiplex the ScaleBio RNA library with other libraries in the sequencing run, a `index2` column can be added with the constant i5 sequence of the ScaleRna library, which is `TGAACCTT[AC]` (8 or 10bp) in forward and `AAGGTTCA[GT]` in reverse orientation. See [samplesheet_v1.0_with_i5.csv](examples/fastq-generation/samplesheet_v1.0_with_i5.csv).

### RNA Kit v1.1
The v1.1 kit uses a different index read setup, one i5 sequence (`index2` read) for each final distribution (PCR) plate well and a pool of 4 i7 sequences (`index1` read) for the whole plate. An example Illumina[samplesheet.csv](examples/fastq-generation/ScaleRNA_3L_samplesheet_v1.1.csv) with all sequences is included. Some sequencers require the `index2` (i.e. `i5`) sequence in the opposite orientation; an example for that is included here: [ScaleRNA_3L_samplesheet_v1.1_revComp.csv](examples/fastq-generation/ScaleRNA_3L_samplesheet_v1.1_revComp.csv).

### Extended Throughput Kit
Each extended throughput plate uses the same set of i5 (`index2`) sequences, but a different pool of i7 (`index1`). A full samplesheet for all four plates can be found here: [ScaleRNA_3L_and_ET_samplesheet_v1.1.csv](examples/fastq-generation/ScaleRNA_3L_and_ET_samplesheet_v1.1.csv) (and with reverse index2: [ScaleRNA_3L_and_ET_samplesheet_v1.1_revComp.csv](examples/fastq-generation/ScaleRNA_3L_and_ET_samplesheet_v1.1_revComp.csv)).

### Splitting Large Fastq Files
This creates a separate set of fastq files for each extended throughput plate (`RNA-A-AP1`, `RNA-A-AP2`, ...). See [Extended Throughput](extendedThroughput.md) for how to analyze these jointly.
To reduce workflow runtime through parallelization (`--splitFastq` mode), it is best to avoid combining the entire data in one fastq file. Generally the Illumina default of splitting per lane (i.e. not setting bcl-convert option `--no-lane-splitting`) is good for smaller runs. An alternative is to modify `samplesheet.csv` to assign a different `sample_ID` to each i5 barcode, following a pattern like `RNA-A-AP1_A1`, `RNA-A-AP1_A2`; where all samples start with the same _library name_ followed by `_` and a unique tag. This way the fastq files will be grouped together for analysis in the ScaleRNA workflow.

## Index Reads
For ScaleBio RNA libraries the RT and ligation cell-barcodes are included in read 1, while the PCR cell-barcode is in the index read. Hence index read fastq files are a required input. When using `bcl-convert` include this `samplesheet.csv` setting: \
`CreateFastqForIndexReads,1`

# Using Fastq Files as Workflow Input
Set `--fastqDir` to the directory containing the fastq files. All the fastq files need to be present in this direcotyr or its subdirectories. 

FASTQ FILE NAMES: The names must follow this pattern `<libName>_..._<R1/R2/I1/I2>_...fastq.gz`

Read through the following checkpoints to ensure your fastq files names are correctly formatted:

* `libName` is the library name (`ScaleRNA` by default) and fastq files must contain the `libName` as a prefix up until the first underscore (`ScaleRNA_...<R1/R2/I1/I2>...fastq.gz`) of the `samples.csv`.

* Only the first part of `libName` in the fastq file-name has to match `libName` column from `samples.csv`
* e.g. fastq file name = `ScaleRNA_AXHGN18SWA_L001_R1.001.fastq.gz` and libName = `ScaleRNA`. This allows having multiple input fastq files for each library

* Generating fastqs with the `libName` prefix can be accomplished with the `sample_ID` column in the Illumina _bcl-convert_ _samplesheet.csv_. Copy the `sample_ID` name prefix up until the first underscore and add to the `libName` column in the ScaleBio `samples.csv`.

* `I1` and `I2` are index reads (i7 and i5, respectively) in fastq format.
