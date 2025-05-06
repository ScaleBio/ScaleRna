## ScalePlex Pool as a "Sample"

For the purposes of the pipeline inputs and outputs, the ScalePlex pool of samples that were fixed on one ScalePlex fixation plate and pooled into a grouping of cells that were loaded into specific RT wells will be treated as a single "sample" in the pipeline. This includes defining them in the sample barcode table and for the generation of a "sample" report. When a pool is then defined, cell calling will first happen for that entire ScalePlex pool of cells as if there were no ScalePlex. Once the set of barcodes that refer to cells are determined, the ScalePlex portion of the workflow will attempt to assign those called cells to their ScalePlex fixation plate of origin. When looking at the Sample Report that refers to the ScalePlex pool of cells, the Summary tab will show a single Barcode Rank Plot and a single value for Read Metrics and Cell Metrics for the pool. As described below, a ScalePlex tab will go into greater detail regarding that pool of cells performance and deconvolution into the individual ScalePlex-fixed samples, how many cells were successfully assigned, and their associated sensitivity metrics.

## ScalePlex Library

When using the ScalePlex Oligo Fixation Plate a separate, paired enriched library is generated from each index PCR reaction to allow further demultiplexing of cells. Additional analysis options to process the reads from this library are detailed here. 


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
- The `samples` directory will contain an amendment to the `<sample>.allCells.csv` file, where the `passing_scaleplex` and `assigned_scaleplex` columns (described below) will now be included. In addition, there will be a `<sample>.ScalePlex.cellMetrics.csv` file for the additional ScalePlex per-cell metadata.
- At the top-level of the output directory there is a `scaleplex` directory. 
### ScalePlex directory contents
| Directory | File | Description |
|-----------|------|-------------|
| `scaleplex`| `demux/` | directory containing read level barcode level validation information|
| | `<sample>.<libIndex2>.ScalePlex.raw.matrix/`| per-cell barcode, per-ScalePlex oligo umi counts matrix for all potential cell barcodes detected in the ScalePlex library|
| | `<sample>.<libIndex2>.ScalePlex.filtered.matrix/`| per-cell barcode, per-scalePlex oligo umi counts matrix for only cell barocdes detected in the ScalePlex library that were called as a cell in the RNA fraction of the analysis|
| | `<sample>.<libIndex2>.ScalePlex.cellMetrics.csv`| using the same cell barcodes as in the raw.matrix/ folder, per-cell barcode ScalePlex library meta data information. Notably, this does not include assignment as this is largely used in our `--reporting` only pipeline run|

### Samples directory additions with the inclusion of ScalePlex
- In the `samples` directory there is a `<sample>.<libIndex2>.ScalePlex.allCells.csv` with per-cell ScalePlex-related metrics. The most important column is `assigned_scaleplex` that corresponds to the fixation well that the cell was fixed in prior to pooling. ScalePlex assignment is only performed for passing cells from the RNA workflow, and all per-cell metrics reported at a sample level are relative to those passing cells from the RNA workflow as the reference point. A full descriptor of the contents of the columns in this file can be found below:

| Measure | ​Description​ |
|---------|-------------|
| reads | The number of reads associated with that cell in the ScalePlex library |
| noScalePlex | The ratio of reads of the barcode in question that did not have a ScalePlex oligo sequence detected |
| umis | The number of ScalePlex oligo UMIs detected |
| scaleplex | The number of uniquely detected scalePlex oligos, defined as at least one read​ |
| max | The number of UMIs associated with the top / most highly detected ScalePlex oligo in that cell |
| second | The number of UMIs associated with the second highest detected ScalePlex oligo in that cell |
| third | The number of UMIs associated with the third highest detected ScalePlex oligo in that cell |
| | - Note: When using the FC based method for assignment, the workflow will check that the fold change enrichment of the second over the third value per cell is > `scalePlexFCThreshold` (default 2). If so, assignment proceeds |
| purity​ | Proportion of UMI coming from max ScalePlex oligo |​
| topTwo​ | Proportion of UMI coming from both max and second |​
| minorFrac​ | Ratio of second to max​ |
| Saturation​ | Ratio of UMI detected over Usable Reads |
| topTwo_scaleplex | Which two scalePlex oligos were the two highest detected |
| ALIAS_alias | columns that denote the alias of the well coordinate for the levels of cell barcoding |
| assigned_scaleplex | Final assignment of cell barocde, with successful assignment corresponding to the fixation plate well of sample origin |
| sample | Value for which sample as denoted by the samples.csv that the barcode originated from |

### reports
- With the usage of ScalePlex in a workflow run, there are several amendments to the reporting structure that are worth noting. First and foremost is the generation of a library report for each ScalePlex library in the workflow. With ScalePlex, libraries are at the level of Index PCR reactions (or final distribution plates), such that for each Index PCR reaction used for your analysis will have both an RNA library report as well as a Scaleplex library report. These capture the read attribution per sample within the library, barcode validation pass rates, and Scaleplex oligo detection pass rates per read of the library.

- Sample reporting, as defined by the individual rows of the Sample Barcode Table (samples.csv), will also have an updated "ScalePlex" tab that summarizes the performance of the ScalePlex data associated in the "ScalePlex Metrics" table. At a high level, these metrics are calculated much in the same way as the [RNA sample level metrics](qcReport.md), such as Reads Per Cell, Counts Per cell, Saturation, and Reads in Cells, but they now reference the ScalePlex library fraction rather than the RNA material. In addition, we also report the "Percent of Cells with Assigned ScalePlex". Critically, the cells in which we are referring to here are those barcodes that were "called" a cell in the RNA analysis, so "Percent of Cells with Assigned ScalePlex" says how many cells called by RNA also had a valid ScalePlex assignment.

- Other Reporting Figures:
    - Assigned ScalePlex Cell Counts: Bar chart of number of cells assigned each scaleplex well of origin
    - Saturation Per Cell: Scatter plot of Saturation vs ScalePlex reads per barcode
    - Top ScalePlex Fraction: Histogram of the the percent of ScalePlex counts per cell that are originating from the top two ScalePlex oligos per cell. This ratio is important because in ideal scenarios, we are enriching for two ScalePlex oligos per cell, that in combination uniquely mark the fixation well of origin. Thus, a peak close to 1 is ideal. If the distribution shifts left, then using the "fc" method may improve your assignment.
    - Fixation Plate Assigned ScalePlex: Plate map of the ScalePlex fixation plate, with values corresponding to the number of cells assigned to each well
    - ScalePlex Assignment Errors:
        - Max Fail
            - A cell fails assignment and is labeled "Max Fail” in the “bg” method when at least one of the top two ScalePlex oligos detected in that cell were not among those that passed the statistical background test for that cell
        - Enrich Fail
            - A cell fails assignment and is labeled “Enrich Fail” in the “bg” method in the case that "scalePlexPercentFromTopTwo" is enabled and fails to pass the enrichment threshold set by the user supplied value. 
            - A cell fails assignment and is labeled “Enrich Fail” in the “fc” method in the case that the fold change of the second highest ScalePlex oligo counts over the third highest ScalePlex oligo count was < 2.
        - Unexpected
            - A cell fails assignment and is labeled “Unexpected" when it passed all of the other assignment criteria for either the “bg” or “fc” method but the two passing ScalePlex oligos in combination do not match with a combination that either is:
                - not on the plate in the first place such as the top two being from two different columns rather than a combination of a row and column
                - If the user supplied an entry for “scalePlexBarcodes" for the sample in question and the assignment falls outside the constraint specified
        - Enrich Fail
            - A cell fails assignment and is labeled “Enrich_Fail” in the “bg” method in the case that " scalePlexPercentFromTopTwo" is enabled and fails to pass the enrichment threshold set by the user supplied value.
            - A cell fails assignment and is labeled “Enrich_Fail” in the “fc” method in the case that the fold change of the second highest ScalePlex oligo counts over the third highest ScalePlex oligo count was < 2.  
    - Per-well ScalePlex RNA Metrics: RNA library metrics and number of cells for each ScalePlex assignment group present within the sample. These are also summarized in the `<sample>.scaleplex_stats.csv` file in the `csv` folder of the reports directory. 
