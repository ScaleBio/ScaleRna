# P5 indexed PCR plates

To analyze data generated with the PCR-plate barcodes in the index2 read (P5 side), use the corresponding library structure file `libV1p5`; i.e.\
`--libStructure libV1p5.json`.

To generate fastq files for this data outside the workflow, see the example [samplesheet.csv](examples/p5-Plates/samplesheet.p5.csv).

## P7 index pools
If running multiple plates with different plate barcode pools (P7 pool), set the `libIndex` column to the name of the barcode pool (e.g. `RNA-A-AP1`); See [samples.p5.csv](examples/p5-Plates/samples.p5.csv)

Note that each plate is processed independently and samples will not be merged across plates. Hence the same sample name (ID) cannot appear twice in `samples.csv` for two different plates.