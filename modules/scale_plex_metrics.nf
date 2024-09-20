/**
* Use RNA passing cells and hash UMI matrix to generate metrics and filtered UMI matrix
*
* Processes:
*     MetricsGeneration
*/

process MetricsGeneration {
    tag "$meta.sample"
    label 'report'
    publishDir { outDir }, pattern: "${id}.filtered.matrix", mode: 'copy'
    publishDir { outDir }, pattern: "${id}.raw.matrix", mode: 'copy'
    publishDir { outDir }, pattern: "${id}.cellMetrics.csv", mode: 'copy'
    publishDir { samplesDir }, pattern: "${id}.ScalePlex.allCells.csv", mode: 'copy'
	publishDir { reportDir / "csv" }, pattern: "${id}.scaleplex_stats.csv", mode: 'copy'

    input:
    path(libStructDir) // Directory containing the library structure definition
    val(libStructName) // Filename of the library structure definition .json
    tuple val(meta), path(cellMetrics), path(rawMtx), path(allCells), val(libCount)

    output:
    tuple val(meta), path("${id}.ScalePlex.allCells.csv"), path("${id}.filtered.matrix"), path("${id}.scaleplex_stats.csv"), path("metrics.csv"), emit: metrics
    tuple val(meta), path("${id}.raw.matrix")
    tuple val(meta), path("${id}.cellMetrics.csv")

    script:
	libStruct = "${libStructDir}/${libStructName}"
    id = meta.sample + "." + meta.lib
    outDir = file(params.outDir) / "scaleplex"
	samplesDir = file(params.outDir) / "samples"
    reportDir = file(params.outDir) / "reports"
	if (libCount > 1 && meta.lib != "merged") {
		samplesDir = samplesDir / "${meta.sample}_libraries"
		reportDir = reportDir / "${meta.sample}_libraries"
	}
	opts = ""
	if (meta.barcodes) {
		opts += "--expected_combos '${meta.barcodes}'"
    }
    """
    metrics_generation.py --umi_matrix $rawMtx --all_cells $allCells --cell_stats $cellMetrics --references $libStructDir --lib_struct $libStruct --id $id --outDir . --assignment_method $params.scalePlexAssignmentMethod --toptwo_frac ${params.scalePlexPercentFromTopTwo / 100} $opts --fc_threshold $params.scalePlexFCThreshold
    # move raw matrix to matrix directory
    mv $rawMtx ${id}.raw.matrix/matrix.mtx
    """
}

workflow SCALE_PLEX_REPORTING {
    take:
        perSampleCountRaw
        libJson
        samples
        allCells
        libCount

    main:
        perSampleCountRaw
            .join(
                // recover matching rnaId
                samples
                    .filter { it.rnaId != null }
                    .map {
                        meta = [lib:it.libName, sample:it.sample]
                        [meta, it.rnaId, it.scalePlexBarcodes]
                    }
            )
            .map {
                // move rnaId to first item in tuple to join with allcells
                meta, cellMetrics, rawUmi, rnaId, barcodes ->
                [rnaId, meta, cellMetrics, rawUmi, barcodes]
            }
            .join(
                // attach allCells based on rnaId
                allCells
                    // match CountScalePlex out for join
                    .map { id, lib, allCells ->
                        [id, allCells]
                    }
            )
            .join(libCount)
            .map {
                // can add rnaId and barcodes to meta now that allCells is joined
                rnaId, meta, cellMetrics, rawUmi, barcodes, allCells, libCount ->
                [meta + [rnaId: rnaId, barcodes: barcodes], cellMetrics, rawUmi, allCells, libCount]
            }
            .dump(tag:'perSampleCountRaw')
            .set { perSampleCountRaw }
        MetricsGeneration(libJson.getParent(), params.scalePlexLibStructure, perSampleCountRaw)
    
    emit:
        metrics = MetricsGeneration.out.metrics
}
