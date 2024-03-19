
# Cell Typing Workflow

This workflow enables the clustering and/or cell type annotation for each samples.

## Seurat

We use Seurat to perform initial normalization, dimension reduction, and clustering of the UMI count matrices for each sample.

A detailed description of the Seurat-based workflow can be found on the Seurat website [here](https://satijalab.org/seurat/articles/pbmc3k_tutorial).

If a dataset contains more than fifty-thousand cells the Seurat-based workflow will automatically perform a sketched-based analysis. This will select a "sketch" or subset of cells to analyze. The results from the analysis of the "sketched" dataset are then projected back to the full dataset. By default the size of the sketched dataset only contains 20% of the total cells.

A more detailed description of the sketch-based workflow can be found on the Seurat website [here](https://satijalab.org/seurat/articles/seurat5_sketch_analysis)

This clustering workflow can be enabled by setting the parameter `seurat = true`

## Azimuth

We also have the ability to perform cell type annotation using the Azimuth reference-based mapping approach for samples with a relevant Azimuth reference dataset. More details about the Azimuth workflow can be found on the Azimuth website [here](https://azimuth.hubmapconsortium.org/)

This cell type annotation workflow can be enabled by setting the parameter `azimuth = true`

## AnnData Output

It is possible to output UMI count matrices as an AnnData object. This can be enabled using the parameter `annData = true`. Doing so will create an AnnData object containing the UMI count matrix for each sample in the samplesheet. This parameter can also be used in combination with the `compSamples` parameter (described below). To create combined AnnData objects of multiple samples.

## Co-Analysis

We have added a feature that enables the co-analysis of multiple samples. Setting the parameter `compSamples = true` will concatenate the filtered UMI count matrices of each sample in the samplesheet into a single UMI count matrix. This combined UMI count matrix is then clustered using our seurat-based workflow. And optionally annotated with our azimuth-based cell typing workflow. Optionally, when `compSamples = true`, a `compSamples` column can be added to the samplesheet as well. The value of this column can be any string, and samples with the same `compSamples` value will be combined for co-analysis. Below are example samplesheets for passing a value in the `compSamples` column of the samplesheet. 

sample | barcodes | compSamples
-- | -- | --
Foo | 1A-6H | combine
Bar | 7A-12H | combine

This is an example samplesheet for combining samples from multiple runs using the reporting only workflow and `compSamples`.

sample | resultDir | compSamples
-- | -- | --
Foo1 | resultDir_1 | fooSamples
Bar1 | resultDir_2 | barSamples
Foo2 | resultDir_3 | fooSamples
Bar2 | resultDir_4 | barSamples