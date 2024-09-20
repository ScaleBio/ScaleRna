/**
* Create hash UMI matrices for each sample from the barcodes csv output by bcParser
* and use along with RNA passing cells to generate metrics and filtered UMI matrix
*
* Processes:
*     CountScalePlex
*     MergeCountRaw
*     MetricsGeneration
*/
include { SCALE_PLEX_REPORTING as REPORTING } from './scale_plex_metrics.nf'
include { SCALE_PLEX_REPORTING as MERGED_REPORTING } from './scale_plex_metrics.nf'

process CountScalePlex {
    tag "$meta.subsample"
    label 'report'

    input:
    path(libStructDir) // Directory containing the library structure definition
    val(libStructName) // Filename of the library structure definition .json
    tuple val(meta), path("barcodes_*.csv")

    output:
    tuple val(meta), path("${id}.cellMetrics.csv"), path("${id}.raw.umi.mtx"), emit: countRaw

    script:
	libStruct = "${libStructDir}/${libStructName}"
    id = meta.subsample + "." + meta.lib
    """
    # Concatenate output from multiple bc_parser jobs
    awk ' NR == 1 || FNR > 1 { print } ' barcodes_*.csv > barcodes.csv
    count_hashes.py barcodes.csv --references $libStructDir --lib_struct $libStruct --id $id --outDir .
    """
}

process MergeCountRaw {
    tag "$meta.sample"
    label 'report'

    input:
    tuple val(meta), path(cellMetrics), path(mtx)

    output:
    tuple val(meta), path("${id}.cellMetrics.csv"), path("${id}.raw.umi.mtx"), emit: countRaw

    script:
    id = meta.sample + "." + meta.lib
    """
    # Merge cellMetrics csv files only taking header from first one
    awk ' NR == 1 || FNR > 1 { print } ' $cellMetrics > ${id}.cellMetrics.csv
    merge_mtx.py $mtx --id $id
    """
}

process MergeSamples {
    tag "$meta.sample"
    label 'report'

    input:
    tuple val(meta), path(cellMetrics), path("raw_??????.mtx"), val(barcodes), val(rnaId) // pad filenames with 0s so they sort correctly

    output:
    tuple val(meta), path("${id}.cellMetrics.csv"), path("${id}.raw.umi.mtx"), emit: countRawMerged
    tuple val(meta), val(mergedBarcodes), emit: mergedBarcodes

    script:
    id = "${meta.sample}.merged"
    rnaId = rnaId.join(" ")
    // merge expected fixation plate wells for each sample in the group
    if (barcodes.every { it }) {
        mergedBarcodes = barcodes.join(';')
    } else {
        // if it's null or '' for any sample in the group, result is null, i.e. allow any valid combos
        mergedBarcodes = null
    }
    """
    merge_metrics.py $cellMetrics --id $id --rnaId $rnaId
    merge_mtx.py raw_*.mtx --id $id
    """
}

workflow SCALE_PLEX {
    take:
        libJson
        samples
        barcodesCsv
        allCells
        libCount
        perSampleCountRaw // provided if --reporting run

    main:
        if (!params.reporting) {
            // CountScalePlex on output of bcParser, dropping where sample was not identified
            CountScalePlex(libJson.getParent(), params.scalePlexLibStructure, barcodesCsv.filter { meta, _ -> meta.sample != 'Unknown' })
            CountScalePlex.out.countRaw
                .map { meta, cellMetrics, rawUmi ->
                    // drop subsample that was used for splitting
                    [[lib:meta.lib, sample:meta.sample], cellMetrics, rawUmi]
                }
                .set { perSampleCountRaw }
            if (params.splitFastq) {
                perSampleCountRaw
                    // collect all cellMetrics, barcodes, and umi.sparse.mtx for each lib and sample combination
                    .groupTuple()
                    | MergeCountRaw
                    | set { perSampleCountRaw }
            }
        }
        REPORTING(perSampleCountRaw, libJson, samples, allCells.id, libCount.id)
        if (params.merge) {
            perSampleCountRaw
                .join(
                    samples
                        .map { [[lib:it.libName, sample:it.sample], it.group, it.scalePlexBarcodes, it.rnaId] }
                )
                .groupTuple(by: 3) // grouping on samples.csv group column
                .filter { it[0].size() > 1 } // Only merge groups with more than one sample libname combination
                .map { 
                    meta, cellMetrics, rawUmi, group, scalePlexBarcodes, rnaId ->
                        [[lib:'merged', sample:group], cellMetrics, rawUmi, scalePlexBarcodes, rnaId]
                }
                .dump(tag:'perMergedSampleCountRaw')
                .set { perMergedSampleCountRaw }
            MergeSamples(perMergedSampleCountRaw)
            perMergedSampleCountRaw
                .join(MergeSamples.out.mergedBarcodes)
                .map { meta, cellMetrics, rawUmi, barcodes, rnaId, mergedBarcodes ->
                    // build samples channels for merging
                    ['libName':meta.lib, 'sample':meta.sample, 'rnaId':meta.sample, 'scalePlexBarcodes':mergedBarcodes]
                }
                .dump(tag:'mergedSamples')
                .set { mergedSamples }
            MERGED_REPORTING(MergeSamples.out.countRawMerged, libJson, mergedSamples, allCells.merged, libCount.merged)
            mergedMetrics = MERGED_REPORTING.out.metrics
        } else {
            mergedMetrics = Channel.empty() // initialize mergedMetrics
        }
    
    emit:
        metrics = REPORTING.out.metrics
        mergedMetrics = mergedMetrics
}
