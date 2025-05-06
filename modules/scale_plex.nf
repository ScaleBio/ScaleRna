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

    input:
    path(libStructDir) // Directory containing the library structure definition
    val(libStructName) // Filename of the library structure definition .json
    val(file_ext)
    path(mapping_file) // File containing the mapping of scale plex to RNA PCR barcodes
    tuple val(meta), path("barcodes_*.${file_ext}")

    output:
    tuple val(meta), path("${id}.cellMetrics.parquet"), path("${id}.raw.umi.mtx.gz"), emit: countRaw

    script:
	libStruct = "${libStructDir}/${libStructName}"
    id = meta.subsample + "." + meta.lib
    opts = ""
    if (mapping_file) {
        opts = "--mapping_file $mapping_file "
    }
    """
    count_hashes.py --bamDir . --references $libStructDir --lib_struct $libStruct --id $id --min_read_umi_thresh $params.scalePlexMinReadUmi --outDir . $opts
    """
}

process MergeCountRaw {
    tag "$meta.sample"
	label 'large'

    input:
    tuple val(meta), path(cellMetrics), path(mtx)

    output:
    tuple val(meta), path("${id}.cellMetrics.parquet"), path("${id}.raw.umi.mtx.gz"), emit: countRaw

    script:
    id = meta.sample + "." + meta.lib
    """
    merge_metrics.py $cellMetrics --id $id
    merge_mtx.py $mtx --id $id
    """
}

process MergeSamples {
    tag "$meta.sample"
	label 'large'

    input:
    tuple val(meta), path(cellMetrics), path("raw_??????.mtx.gz"), val(barcodes), val(rnaId) // pad filenames with 0s so they sort correctly

    output:
    tuple val(meta), path("${id}.cellMetrics.parquet"), path("${id}.raw.umi.mtx.gz"), emit: countRawMerged
    tuple val(meta), val(mergedBarcodes), emit: mergedBarcodes

    script:
    id = "${meta.sample}.merged"
    // merge expected fixation plate wells for each sample in the group
    if (barcodes.every { it }) {
        mergedBarcodes = barcodes.join(';')
    } else {
        // if it's null or '' for any sample in the group, result is null, i.e. allow any valid combos
        mergedBarcodes = null
    }
    """
    merge_metrics.py $cellMetrics --id $id --rnaId ${rnaId.join(" ")}
    merge_mtx.py raw_*.mtx.gz --id $id
    """
}

workflow SCALE_PLEX {
    take:
        libJson
        scalePlexToRnaMapping
        samples
        ubam
        allCells
        libCount
        perSampleCountRaw // provided if --reporting run

    main:
        if (!params.reporting) {
            if (params.ultimaCramDir) {
                file_ext = 'cram'
            } else {
                file_ext = 'bam'
            }
            // CountScalePlex on output of bcParser, dropping where sample was not identified
            CountScalePlex(libJson.getParent(), libJson.getName(), file_ext, scalePlexToRnaMapping, ubam.filter { meta, _ubam_file -> meta.sample != 'Unknown' })
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
                    _meta, cellMetrics, rawUmi, group, scalePlexBarcodes, rnaId ->
                        [[lib:'merged', sample:group], cellMetrics, rawUmi, scalePlexBarcodes, rnaId]
                }
                .dump(tag:'perMergedSampleCountRaw')
                .set { perMergedSampleCountRaw }
            MergeSamples(perMergedSampleCountRaw)
            perMergedSampleCountRaw
                .join(MergeSamples.out.mergedBarcodes)
                .map { meta, _cellMetrics, _rawUmi, _barcodes, _rnaId, mergedBarcodes ->
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
