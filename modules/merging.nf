include { SampleStats } from './create_mtx.nf'
include { SampleStats as ScalePlexStats } from './scale_plex_metrics.nf'
include { SAMPLE_REPORTING } from './sample_reporting.nf'
include { DOWNSTREAM } from './downstream.nf'

// Reorder list based on what would sort refList
// In practice, refList would be sample names and the other lists would be workdir paths
def reorderList(list, refList) {
    def argsort = refList.sort(false).collect { name -> refList.indexOf(name) } // false means don't sort in-place
    argsort.collect { index -> list[index] }
}

// Merge cell metrics and filtered mtx for individual libraries
// Assumes different cells in each subsample
process MergeMetrics {
	tag "$group"
	label 'large'

	publishDir { file(params.outputDir) / "samples" }, pattern:"${group}_filtered_star_output", saveAs: { "${group}.merged.filtered.matrix" }, mode: 'copy'
	publishDir { file(params.outputDir) / "samples" }, pattern: "${group}.merged.allCells.csv", mode: 'copy'

	input:
	tuple val(group), path(allCells), path(allBarcodes), path(filteredMtx)
	
	output:
    tuple val(group), val("merged"), path("${group}.merged.allCells.csv"), emit: allCells
    tuple val(group), val("merged"), path("${group}.merged.allBarcodes.parquet"), emit: allBarcodes
    tuple val(group), val(libCount), emit: libCount
	tuple val(group), path("${group}_filtered_star_output"), emit: filtered_star_mtx

	script:
    libCount = allCells.size()
	"""
    concat_columnar.py $allCells --outputFile ${group}.merged.allCells.csv
    concat_columnar.py $allBarcodes --columns pass "COLUMNS('counts\$|species|bead_bc')" totalReads genes Saturation --outputFile ${group}.merged.allBarcodes.parquet
    concat_mtx.py $filteredMtx --id $group --threads $task.cpus --memory '$task.memory'
	"""
}

// Merge ScalePlex metrics for reporting
// Assumes different cells in each subsample
process MergeScalePlex {
	tag "$group"
	label 'large'
        
	input:
	tuple val(group), path(allBarcodes)
	
	output:
    tuple val(group), path("${group}.merged.ScalePlex.allBarcodes.parquet"), emit: allBarcodes

	script:
	"""
    concat_columnar.py $allBarcodes --columns pass counts totalReads Saturation topTwo assigned_scaleplex --outputFile ${group}.merged.ScalePlex.allBarcodes.parquet
	"""
}

// Merge STAR logs for individual libraries
process MergeReadStats {
	tag "$group"
	label 'large'

	input:
	tuple val(group), path("Log_??????.final.out")
	
	output:
    tuple val(group), path("Log.final.out"), emit: mergedStarLog

	script:
	"""
	merge_raw_star_output.py --log_out_dir . --star_log Log_*.final.out
	"""
}

workflow MERGING {
take:
    samples // each row of samples.csv parsed into a channel 
    libJson // library structure json file
    isBarnyard // Different matrix filename for barnyard
    allCells // Cell metrics for passing cells
    allBarcodes // Metrics for all cell-barcode combinations
    filteredMtx // Gene expression matrix for passing cells
    soloLog // STAR Solo log files
    scalePlexResults // ScalePlex barcode metrics and umi matrix if --scalePlex
    totalSampleReadsBySampleID // Total reads per sample from bcParser
main:
    samples
        .map { tuple(it.id, it.group) }
        .join(allCells)
        .join(allBarcodes)
        .join(filteredMtx)
        .groupTuple(by:1) // All samples grouped for merging by library/group
        .filter{ it[0].size() > 1 } // Only merge samples with multiple libraries
        .map {
            ids, group, _libs, cellMetrics, _libs2, barcodeMetrics, mtxDirs ->
            // Compute indices that would sort the samples list for deterministic output
            [
                group,
                reorderList(cellMetrics, ids),
                reorderList(barcodeMetrics, ids),
                reorderList(mtxDirs, ids),
            ]
        }
        .dump(tag:'mergeMetricsInput')
        .set{ mergeMetricsInput }
    MergeMetrics(mergeMetricsInput)

    samples
        .map { tuple(it.id, it.group) }
        .join(soloLog)
        .groupTuple(by:1) // All samples grouped for merging by library/group
        .filter{ it[0].size() > 1 } // Only merge samples with multiple libraries
        .map {
            ids, group, starLog ->
            [
                group,
                reorderList(starLog, ids)
            ]
        }
        .dump(tag:'mergeReadStatsInput')
        .set{ mergeReadStatsInput }
    MergeReadStats(mergeReadStatsInput)

    samples
        .map { tuple(it.id, it.group) }
        .join(allBarcodes)
        .join(totalSampleReadsBySampleID)
        .groupTuple(by:1) // All samples grouped for merging by library/group
        .filter{ it[0].size() > 1 } // Only merge samples with multiple libraries
        .map {
            it ->
                // unpack tuple
                def (ids, group, _libs, cellMetrics, totalLibReads) = it
            def totalReads = 0
            // Sum up total reads of samples belonging to same group
            totalLibReads.each { reads ->
                totalReads += reads.toLong()
            }
            [
                group,
                'merged',
                reorderList(cellMetrics, ids),
                totalReads,
            ]
        }
        .join(MergeReadStats.out.mergedStarLog)
        .dump(tag:'sampleStatsInput')
        .set{ sampleStatsInput }
        SampleStats(sampleStatsInput, isBarnyard)
        
    if (params.scalePlex) {
        samples
            .map { tuple(it.id, it.group) }
            .join(scalePlexResults)
            .join(allCells)
            .groupTuple(by:1) // All samples grouped for merging by library/group
            .filter{ it[0].size() > 1 } // Only merge samples with multiple libraries
            .map {
                ids, group, scalePlexBarcodes, _stats, _metrics, _lib, _lib2, allCellsWithAssignment ->
                [
                    [sample: group, lib: "merged"],
                    reorderList(allCellsWithAssignment, ids),
                    reorderList(scalePlexBarcodes, ids),
                    ids.size(), // number of libraries
                ]
            }
            .dump(tag:'scalePlexInput')
            .set{ scalePlexInput }
        ScalePlexStats(scalePlexInput)
        scalePlexInput.map {
            meta, _allCells, allScalePlexBarcodes, _libCount ->
                [meta.sample, allScalePlexBarcodes]
        }
        .dump(tag:'mergeScalePlex')
        .set{ mergeScalePlex }
        MergeScalePlex(mergeScalePlex)
        MergeScalePlex.out.allBarcodes
            .join(ScalePlexStats.out.metrics
                .map { meta, stats, metrics ->
                    [meta.sample, stats, metrics, meta.lib]
                })
        .dump(tag:'scalePlexReportingResults')
        .set { scalePlexReportingResults }
    }
    SAMPLE_REPORTING(
        samples,
        libJson,
        isBarnyard,
        MergeMetrics.out.allCells,
        MergeMetrics.out.allBarcodes,
        MergeMetrics.out.libCount,
        SampleStats.out.sample_metrics,
        params.scalePlex ? scalePlexReportingResults : [],
        true,
    )
    if (params.seurat || params.azimuth || params.annData) {
        DOWNSTREAM(samples, MergeMetrics.out.allCells, MergeMetrics.out.filtered_star_mtx, SampleStats.out.sample_metrics, libJson, false)
    }
    emit:
        sampleStats = SAMPLE_REPORTING.out.sampleStats
}