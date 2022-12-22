# Reference genomes
The workflow requires a genome reference and annotation to run. All files and settings for a genome are defined in a [genome.json](examples/genomes.json) file. When launching the workflow, the reference genome is selected by passing the path to a specific `genome.json` in the `--genome` parameter.

The `genome.json` file includes

Field |  Description | Required? | Example
:-- | -- | -- | --
name | The name of the species / genome-version | Required | human 
star_index | Path to the STAR index directory | Required | `/PATH/TO/star.ann` 
gtf | Path to the gene annotation | Required | `filteredGTF/genes.gtf` 
speciesName | Full name of the species from [OrgDb](https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb) | optional | Homo sapiens

* All files (`star_index`, `gtf`, ...) can be specified either as
    - an absolute path (`/path/to/genome`)
    - a relative path starting from the location of the `genome.json` file (`genes/filtered.gtf`)
    - a AWS S3 url (s3://path/to/genome)

## STAR index
The provided STAR index needs to be compatible with STAR 2.7. It also has to include the gene annotation, e.g. be built with `--sjdbGTFfile /PATH/TO/genes.gtf` See the STAR [documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for additional options.

## Annotation
All annotated transcripts in the provided _GTF_ file are included in the analysis (i.e. in the output gene expression matrix). To exclude certain annotations (e.g. pseudo-genes), filter the GTF before generating the STAR index.

## Pre-built genomes
Pre-build reference genome for human is available for download:
* Human: http://scale.pub.s3.amazonaws.com/genomes/rna/grch38.tgz

Download these to your analysis server, unpack them and then use e.g.
`--genome /PATH/TO/grch38/grch38.json`

