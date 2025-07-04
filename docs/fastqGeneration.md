# Fastq Generation
The workflow can be started from a sequencer runFolder (bcl) [`--runFolder`]. In that case fastq files will be generated internally, using Illumina `bcl-convert`.

Alternatively it is also possible to generate fastq files upstream using standard Illumina software, for example when the ScaleRNA library is multiplexed together with other libraries during sequencing at a core facility. We recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) version 3.9 or later.

## Samplesheet.csv
### QuantumScale RNA kit
Example [samplesheets](examples/fastq-generation/) for all kit sizes are provided. We use all 96 i7 bead1 barcode sequences along with all i5 PCR barcode sequences to enable a dual index demux strategy. Based on the QuantumScale kit used, a subset of PCR sequences will be used in the experiment. For the large and extra large kit we use PCR sequences QSR-1 to QSR-8; for the small and medium kit QSR-P; and for the modular kit QSR-1 to QSR-12. The example ScalePlex samplesheets can also be found in the same [directory](examples/fastq-generation/). If you are using a sequencer which produces index reads in the reverse orientation, example samplesheets for that are also provided.

### Splitting Large Fastq Files
To reduce workflow runtime, it is best to avoid combining the entire data in one fastq file. Generally the Illumina default of splitting per lane (i.e. not setting bcl-convert option `--no-lane-splitting`) is good for smaller runs. An alternative is to modify `samplesheet.csv` to assign a different `sample_ID` to each i5 * i7 barcode combination, following a pattern like `QSR-1_01`, `QSR-1_02`; where all samples start with the same _libIndex2_ followed by `_` and the well coordinate of the bead barcode (i7). This way the fastq files will be grouped together for analysis in the ScaleRNA workflow. This is the default way in which the workflow generates the samplesheet, if running the workflow from a runFolder.

### Index Reads
For ScaleBio QuantumScale RNA libraries the RT barcode is in read 2, the bead cell-barcodes are included in index 1, and the PCR cell-barcode is in the index 2 read. Hence, index read fastq files are a required input. When using `bcl-convert` include this `samplesheet.csv` setting: \
`CreateFastqForIndexReads,1`

### Using Fastq Files as Workflow Input
Set `--fastqDir` to the directory containing the fastq files. All the fastq files need to be present in this directory or its subdirectories. The only requirement when using fastq files as input to the workflow is to ensure that both index1 and index2 fastq files are present.

### ScaleRna v1 kits
#### RNA kit v1.0
An example [samplesheet.csv](examples/fastq-generation/samplesheet_v1.0.csv) with typical options is included. Here all 96 i7 barcode sequences from the PCR plate are merged into one set of fastq files. If an index2 (i5) read is used to demultiplex the ScaleBio RNA library with other libraries in the sequencing run, a `index2` column can be added with the constant i5 sequence of the ScaleRna library, which is `TGAACCTT[AC]` (8 or 10bp) in forward and `AAGGTTCA[GT]` in reverse orientation. See [samplesheet_v1.0_with_i5.csv](examples/fastq-generation/samplesheet_v1.0_with_i5.csv).

#### RNA Kit v1.1
The v1.1 kit uses a different index read setup, one i5 sequence (`index2` read) for each final distribution (PCR) plate well and a pool of 4 i7 sequences (`index1` read) for the whole plate. An example Illumina[samplesheet.csv](examples/fastq-generation/ScaleRNA_3L_samplesheet_v1.1.csv) with all sequences is included. Some sequencers require the `index2` (i.e. `i5`) sequence in the opposite orientation; an example for that is included here: [ScaleRNA_3L_samplesheet_v1.1_revComp.csv](examples/fastq-generation/ScaleRNA_3L_samplesheet_v1.1_revComp.csv).

#### Extended Throughput Kit
Each extended throughput plate uses the same set of i5 (`index2`) sequences, but a different pool of i7 (`index1`). A full samplesheet for all four plates can be found here: [ScaleRNA_3L_and_ET_samplesheet_v1.1.csv](examples/fastq-generation/ScaleRNA_3L_and_ET_samplesheet_v1.1.csv) (and with reverse index2: [ScaleRNA_3L_and_ET_samplesheet_v1.1_revComp.csv](examples/fastq-generation/ScaleRNA_3L_and_ET_samplesheet_v1.1_revComp.csv)).

#### Splitting Large Fastq Files for v1 kits
This creates a separate set of fastq files for each extended throughput plate (`RNA-A-AP1`, `RNA-A-AP2`, ...). See [Extended Throughput](extendedThroughput.md) for how to analyze these jointly.
To reduce workflow runtime through parallelization (`--splitFastq` mode), it is best to avoid combining the entire data in one fastq file. Generally the Illumina default of splitting per lane (i.e. not setting bcl-convert option `--no-lane-splitting`) is good for smaller runs. An alternative is to modify `samplesheet.csv` to assign a different `sample_ID` to each i5 barcode, following a pattern like `RNA-A-AP1_A1`, `RNA-A-AP1_A2`; where all samples start with the same _library name_ followed by `_` and a unique tag. This way the fastq files will be grouped together for analysis in the ScaleRNA workflow.
