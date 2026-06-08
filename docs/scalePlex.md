## ScalePlex Pool as a "Sample"

For the purposes of the pipeline inputs and outputs, the set of ScalePlex samples that were fixed on one ScalePlex fixation plate and pooled into a grouping of cells that were loaded into specific RT wells will be treated as a single "sample" in the pipeline. This includes defining them in the sample barcode table and for the generation of a "sample" report. When a pool is then defined, cell calling will first happen for that entire ScalePlex pool of cells independently of the associated ScalePlex data. Once the set of barcodes that refer to cells are determined, the ScalePlex portion of the workflow will attempt to assign those called cells to their ScalePlex fixation plate well of origin. When looking at the Sample Report that refers to the ScalePlex pool of cells, the Summary tab will show a single Barcode Rank Plot and a single value for Read Metrics and Cell Metrics for the combined pool. As described below, a ScalePlex tab will provide further details regarding that pool of cells performance and deconvolution into the individual ScalePlex-fixed samples, how many cells were successfully assigned, and their associated sensitivity metrics.

## ScalePlex Library

When using the ScalePlex Oligo Fixation Plate, a separate enriched library is generated from each index PCR reaction to allow further demultiplexing of cells from its paired RNA library. Additional analysis options to process the reads from this library are detailed here. 


## --scalePlex true

When this parameter is set to true, the workflow looks for the ScalePlex reads as specified in the Sample Barcode Table (samples.csv) portion either from BCL or FASTQ input. Example file here: [quantum sample barcode tables with ScalePlex](examples/quantum-sample-barcode-tables/)

- Note for Input Modes: When starting with fastq files or if the ScaleRNA and ScalePlex libraries were sequenced separately, fastq files will need to be generated for ScaleRNA and ScalePlex into separate files, then placed in one parent directory supplied in the `fastqDir` parameter (files can exist in subdirectories of this supplied path). See [Fastq Generation](fastqGeneration.md) and the example samplesheet.csv including ScalePlex libraries [Example Sample Sheets](examples/fastq-generation/) or for 3lvl RNA [ScaleRNA_3L_and_ET_with_ScalePlex_samplesheet_v1.1.csv](examples/fastq-generation/ScaleRNA_3L_v1.1/ScaleRNA_3L_and_ET_with_ScalePlex_samplesheet_v1.1.csv)

## Sample Barcode Table (samples.csv)

There are additional **optional** columns for the Sample Barcode Table when running the ScalePlex pipeline.
| Column              | Description         | Example |
|---------------------|---------------------|---------|
| scalePlexLibIndex2   | (Quantum Only) Index PCR sequences to associate with enriched library. Full sequence list here: [ScalePlex i5 Quantum](../references/quantum_scaleplex_pcr_pool.txt)   | QSR-1;QSR-2       |
| scalePlexLibIndex   | (Scale RNA v1.1 Only) i7 sequences to associate with enriched library. Full sequence list here: [ScalePlex i7](../references/scaleplex_p7.txt)   | ScalePlex-A-AP1        |
| scalePlexBarcodes   | Valid fixation plate wells for this sample, follows same specification scheme as the [RNA RT barcodes](samplesCsv.md), but is in reference to the ScalePlex fixation plate   | 1A-6H        |

## Analysis parameters

- `scalePlexAssignmentMethod` (Default: 'bg') Use background ('bg') or fold-change ('fc') algorithm for ScalePlex assignment
    - Note: The background based method will estimate a background profile per sample as defined in the Sample Barcode Table, and test counts of ScalePlex oligos vs that estimation. The counts that pass are then validated against the expected counts per ScalePlex fixation plate layout. The fold-changed based method is best suited to situations in which there are low overall ScalePlex oligo counts per  cell or if there is a low number of ScalePlex fixation plate wells used in a pool of samples and instead computes the fold-change of the expected enriched fraction of counts vs the next (third) highest and then assigns if that is above the specified threshold.
- `scalePlexPercentFromTopTwo` (Default: 0) Threshold percent of ScalePlex UMIs from top two unique to pass assignment, e.g. 50
- `scalePlexFCThreshold` (Default: 2) If using `fc` assignment method, set threshold for valid assignment based on fold change of second to third highest detected ScalePlex oligo per cell
## Outputs
- The `samples` directory will contain an amendment to the `<sample>.<libIndex2>.allCells.csv` file, where the `assigned_scaleplex` column (and `passing_scaleplex` when using the `bg` assignment method) will now be included.
- At the top-level of the output directory there is a `scaleplex` directory.
### ScalePlex directory contents
| Directory | File | Description |
|-----------|------|-------------|
| `scaleplex`| `demux/` | directory containing read level barcode level validation information|
| | `<sample>.<scalePlexLibIndex2>.raw.matrix/`| per-cell barcode, per-ScalePlex oligo UMI counts matrix for all potential cell barcodes detected in the ScalePlex library|
| | `<sample>.<scalePlexLibIndex2>.filtered.matrix/`| per-cell barcode, per-ScalePlex oligo UMI counts matrix for only cell barcodes detected in the ScalePlex library that were called as a cell in the RNA fraction of the analysis|
| | `<sample>.<scalePlexLibIndex2>.cellMetrics.parquet`| per-cell barcode ScalePlex library metadata for all barcodes (same cell barcodes as in the raw.matrix/ folder). Does not include assignment columns; see the amended `allCells.csv` in `samples/`|

### Samples directory additions with the inclusion of ScalePlex
- The `samples/<sample>.<libIndex2>.allCells.csv` file is replaced with an amended version that adds ScalePlex assignment results for passing RNA cells. The most important new column is `assigned_scaleplex`, which corresponds to the fixation plate well that the cell was fixed in prior to pooling. ScalePlex assignment is only performed for passing cells from the RNA workflow.

| Column | Description |
|--------|-------------|
| assigned_scaleplex | Final assignment of cell barcode, with successful assignment corresponding to the fixation plate well of sample origin |
| passing_scaleplex | (bg method only) Semicolon-delimited list of ScalePlex oligos that passed the background test for that cell, or `No_Pass` if none passed |

### ScalePlex library per-cell metrics
Per-cell ScalePlex library metrics for all detected barcodes are stored in `scaleplex/<sample>.<scalePlexLibIndex2>.cellMetrics.parquet`. This file covers all barcodes (not just RNA-called cells) and does not include assignment columns.

| Column | Description |
|--------|-------------|
| totalReads | The number of reads associated with that cell in the ScalePlex library |
| noScalePlex | The proportion of reads for that barcode that did not have a ScalePlex oligo sequence detected |
| counts | The total number of ScalePlex oligo UMIs detected |
| scaleplex | The number of distinct ScalePlex oligos detected (with at least one UMI) |
| max | The number of UMIs associated with the top / most highly detected ScalePlex oligo in that cell |
| second | The number of UMIs associated with the second highest detected ScalePlex oligo in that cell |
| third | The number of UMIs associated with the third highest detected ScalePlex oligo in that cell |
| | - Note: When using the FC based method for assignment, the workflow will check that the fold change of the second over the third value per cell is > `scalePlexFCThreshold` (default 2). If so, assignment proceeds |
| purity | Proportion of UMIs coming from the top ScalePlex oligo |
| topTwo | Proportion of UMIs coming from the top two ScalePlex oligos combined |
| minorFrac | Ratio of second to max |
| Saturation | `1 - (counts / usable reads)`, where usable reads are those with a ScalePlex oligo detected |
| topTwo_scaleplex | Which two ScalePlex oligos were the two highest detected |
| ALIAS_alias | columns that denote the alias of the well coordinate for the levels of cell barcoding |

### ScalePlex reports
- With the usage of ScalePlex in a workflow run, there are several amendments to the reporting structure that are worth noting. First and foremost is the generation of a library report for each ScalePlex library in the workflow. With ScalePlex, libraries are at the level of Index PCR reactions (or final distribution plates), such that each Index PCR reaction used for your analysis will have both an RNA library report as well as a Scaleplex library report. These capture the read attribution per sample within the library, barcode validation pass rates, and Scaleplex oligo detection pass rates per read of the library.

- Sample reporting, as defined by the individual rows of the Sample Barcode Table (samples.csv), will also have an updated "ScalePlex" tab that summarizes the performance of the ScalePlex data associated in the "ScalePlex Metrics" table. At a high level, these metrics are calculated much in the same way as the [RNA sample level metrics](qcReport.md), such as Reads Per Cell, Counts Per cell, Saturation, and Reads in Cells, but they now reference the ScalePlex library fraction rather than the RNA material. In addition, we also report the "Percent of Cells with Assigned ScalePlex". Critically, the cells referenced here are defined as barcodes that were "called" a cell in the RNA analysis, so "Percent of Cells with Assigned ScalePlex" says how many cells called by RNA also had a valid ScalePlex assignment.

- Other Reporting Figures:
    - Assigned ScalePlex Cell Counts: Bar chart of number of cells assigned each ScalePlex well of origin
    - Saturation Per Cell: Scatter plot of Saturation vs ScalePlex reads per barcode
    - Top ScalePlex Fraction: Histogram of the the percent of ScalePlex counts per cell that are originating from the top two ScalePlex oligos per cell. This ratio is important because in ideal scenarios, we are enriching for two ScalePlex oligos per cell, that in combination uniquely mark the fixation well of origin. Thus, a peak close to 1 is ideal. If the distribution shifts left, then using the "fc" method may improve your assignment.
    - Fixation Plate Assigned ScalePlex: Plate map of the ScalePlex fixation plate, with values corresponding to the number of cells assigned to each well
    - ScalePlex Assignment Errors:
        - Indeterminate
            - A cell is labeled “Indeterminate” in the “bg” method when no ScalePlex oligos in that cell pass the statistical background test (i.e. `passing_scaleplex` is `No_Pass`).
            - Solution: Verify fixation protocol.
        - Max_Fail
            - A cell is labeled “Max_Fail” in the “bg” method when at least one of the top two ScalePlex oligos detected in that cell were not among those that passed the statistical background test for that cell.
            - Solution: Check oligo quality and concentration.
        - Enrich_Fail
            - A cell is labeled “Enrich_Fail” in the “bg” method when `scalePlexPercentFromTopTwo` is enabled and the cell fails to pass that enrichment threshold.
            - A cell is labeled “Enrich_Fail” in the “fc” method when the fold change of the second highest ScalePlex oligo counts over the third highest ScalePlex oligo count is ≤ `scalePlexFCThreshold` (default: 2), or when no second oligo is detected.
            - Solution: Adjust quality parameters.
        - Unexpected
            - A cell is labeled “Unexpected” when it passed all other assignment criteria for either the “bg” or “fc” method, but the two top ScalePlex oligos in combination do not correspond to a valid fixation plate well. This can happen when:
                - The top two oligos are from the same row or same column rather than a valid row+column combination.
                - The user supplied a `scalePlexBarcodes` constraint for the sample and the assignment falls outside the specified wells.
             - Solution: Check plate layout and oligo design.
    - Per-well ScalePlex RNA Metrics: RNA library metrics and number of cells for each ScalePlex assignment group present within the sample. These are also summarized in the `<sample>.scaleplex_stats.csv` file in the `csv` folder of the reports directory. 
