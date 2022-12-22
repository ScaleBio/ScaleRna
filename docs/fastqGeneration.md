# Fastq generation
We recommend using [bcl-convert](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) from Illumina for fastq generation.

An example [samplesheet.csv](examples/samplesheet.csv) with typical options is included.

## Index reads
For ScaleBio RNA libraries the RT and ligarion barcodes are included in read 1, while the PCR barcode is in the index read. Hence we need to tell `bcl-convert` to ignore the index read for sample demultiplexing. This can be achieved by declaring the _I1_ read as a `UMI` and generating index read fastqs; using `samplesheet.csv` settings:
```
CreateFastqForIndexReads,1
TrimUMI,0
OverrideCycles,Y34;U10;Y76
```

## Using pre-generated fastq files as workflow input
Set `--fastqDir` to the directory containing all fastq files for all samples in an analysis run. The filenames should follow the pattern `<Name>_..._<Read>_...fastq.gz`, where
* `Name` is the library name (`libName` column in `samples.csv`)
* `Read` is one of `R1`, `R2`, `I1`