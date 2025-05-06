# Reference genomes
The workflow requires a genome reference and gene annotation to run. All files and settings for a genome are defined in a [genome.json](examples/genome.json) file. When launching the workflow, the reference genome is selected by passing the path to a specific `genome.json` with the `--genome` parameter.

The `genome.json` file includes

Field |  Description | Required? | Example
:-- | -- | -- | --
name | The name of the species / genome-version | required | `GRCh38` 
star_index | Path to the STAR index directory | required | `/PATH/TO/starRef/` 
gtf | Path to the gene annotation | optional | `filteredGTF/genes.gtf` 
speciesName | Full name of the species from [OrgDb](https://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb) | optional | Homo sapiens
isBarnyard | Indicates that this is a combined multi-species genome | optional | false 

* All files (`star_index`, `gtf`, ...) can be specified either as
    - an absolute path (`/path/to/genome`)
    - a relative path starting from the location of the `genome.json` file (`genes/filtered.gtf`)
    - a AWS S3 url (s3://path/to/genome)

## STAR index
The provided STAR index needs to be built with STAR version `>= 2.7.4a`. It also has to include the gene annotation, e.g. be built with `--sjdbGTFfile /PATH/TO/genes.gtf` See the STAR [documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) for additional options. An example command would be
```
STAR --runMode genomeGenerate --runThreadN 16 --genomeDir star.ref --genomeFastaFiles Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile Homo_sapiens.GRCh38.103.biotypeFiltered.gtf
```

### Annotation
All transcripts in the _GTF_ file used to built the STAR index are included in the analysis (i.e. in the output gene expression matrix). To exclude certain annotations (e.g. pseudo-genes), filter the GTF before generating the STAR index. Only the `exon` entries in the GTF are relevant for STAR.

## Pre-built genomes
Pre-build reference genome for human and mouse are available for download:
* Human (GRCh38): http://scale.pub.s3.amazonaws.com/genomes/rna/grch38.tgz
* Mouse: http://scale.pub.s3.amazonaws.com/genomes/rna/mm39.tgz
* Human/Mouse Barnyard: http://scale.pub.s3.amazonaws.com/genomes/rna/grch38_mm39.tgz

Download these to your system, unpack them and then use e.g.
`--genome /PATH/TO/grch38/grch38.json`

