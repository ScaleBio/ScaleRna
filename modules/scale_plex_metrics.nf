/**
* Use RNA passing cells and hash UMI matrix to generate metrics and filtered UMI matrix
*
* Processes:
*     MetricsGeneration
*/

process MetricsGeneration {
    tag "$meta.sample"
	label 'large'
    publishDir { file(params.outputDir) / "scaleplex" }, pattern: "${id}.filtered.matrix", mode: 'copy'
    publishDir { file(params.outputDir) / "scaleplex" / "${id}.raw.matrix" }, pattern: "${id}.raw.matrix/features.tsv.gz", mode: 'copy', saveAs: { "features.tsv.gz"}
    publishDir { file(params.outputDir) / "scaleplex" / "${id}.raw.matrix" }, pattern: "${id}.raw.matrix/barcodes.tsv.gz", mode: 'copy', saveAs: { "barcodes.tsv.gz"}
    publishDir { file(params.outputDir) / "scaleplex" / "${id}.raw.matrix" }, pattern: "matrix.mtx.gz", mode: 'copy'
    publishDir { file(params.outputDir) / "scaleplex" }, pattern: "${id}.cellMetrics.parquet", mode: 'copy'
	publishDir { file(params.outputDir) / samplesDir }, pattern: "${id}_allCellsWithAssignment.csv", mode: 'copy', saveAs: { meta.lib == "merged" ? "${rnaId}.merged.allCells.csv" : "${rnaId}.allCells.csv" }
    publishDir { file(params.outputDir) / reportDir / "csv" }, pattern: "${id}.scaleplex_stats.csv", mode: 'copy'

    input:
    path(libStructDir) // Directory containing the library structure definition
    val(libStructName) // Filename of the library structure definition .json
    tuple val(meta), path(cellMetrics), path('matrix.mtx.gz'), path(allCells), val(libCount)

    output:
    tuple val(meta), path("${id}.ScalePlex.allBarcodes.parquet"), path("${id}.filtered.matrix"), path("${id}.scaleplex_stats.csv"), path("metrics.csv"), emit: metrics
    tuple val(meta), path("${id}.raw.matrix/features.tsv.gz")
    tuple val(meta), path("${id}.raw.matrix/barcodes.tsv.gz")
    tuple val(meta), path("matrix.mtx.gz", includeInputs: true)
    tuple val(meta), path("${id}.cellMetrics.parquet")
    tuple val(meta), path("${id}_allCellsWithAssignment.csv")

    script:
    libStruct = "${libStructDir}/${libStructName}"
    id = meta.sample + "." + meta.lib
    rnaId = meta.rnaId
    samplesDir = "samples"
    reportDir = "reports"
    if (libCount > 1 && meta.lib != "merged") {
        samplesDir += "/${meta.sample}_libraries"
        reportDir += "/${meta.sample}_libraries"
    }
    opts = ""
    if (meta.barcodes) {
        opts += "--expected_combos '${meta.barcodes}'"
    }
    """
    scaleplex_assignment.py --umi_matrix matrix.mtx.gz --all_cells $allCells --cell_stats $cellMetrics --references $libStructDir --lib_struct $libStruct --id $id --outDir . --assignment_method $params.scalePlexAssignmentMethod --toptwo_frac ${params.scalePlexPercentFromTopTwo / 100} $opts --fc_threshold $params.scalePlexFCThreshold
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
                        def meta = [lib:it.libName, sample:it.sample]
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
                    .map { id, _lib, allCells_file ->
                        [id, allCells_file]
                    }
            )
            .join(libCount)
            .map {
                // can add rnaId and barcodes to meta now that allCells is joined
                rnaId, meta, cellMetrics, rawUmi, barcodes, allCells_file, libraryCount ->
                [meta + [rnaId: rnaId, barcodes: barcodes], cellMetrics, rawUmi, allCells_file, libraryCount]
            }
            .dump(tag:'perSampleCountRaw')
            .set { perSampleCountRaw }
        MetricsGeneration(libJson.getParent(), libJson.getName(), perSampleCountRaw)
    
    emit:
        metrics = MetricsGeneration.out.metrics
}
