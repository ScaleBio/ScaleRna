/*
* Compute per cell metrics from STARsolo output
*
*/

// Generate cell-level metrics per sample from STARsolo output
process CellMetricsGeneration {
	tag "$id"
	label 'large'

	input:
	tuple(val(id), val(sampleName), val(libName), val(expectedCells), path("Solo.out"))
	val(libStructName)
	val(isBarnyard)
	path(libStructDir)
	
	output:
	tuple val(id), val(libName), path("${id}_metrics/${id}_allBarcodes.parquet"), emit: raw_cell_metrics

	script:
	opts = ""
	starMultiParam = params.starMulti
	if (isBarnyard) {
		opts = opts + "--isBarnyard "
		starMultiParam = params.starMultiBarnyard
	}

	libStruct = "${libStructDir}/${libStructName}"
	matrixFn = Utils.starMatrixFn(starMultiParam)
"""
	export POLARS_MAX_THREADS=$task.cpus
	build_all_cells.py --STARsolo_out Solo.out --feature_type $params.starFeature --matrix_type $matrixFn \
	--sample $id --libStruct $libStruct $opts --threads $task.cpus --memory '$task.memory'
"""
}

// Analyze detected barcodes from library for bead filtering
process BeadFiltering {
	tag "$libName"
	label 'large'

	input:
	tuple val(libName), path(allBarcodes)
	val(sampleBarcode)
	
	output:
	tuple val(libName), path("bead_scores.parquet"), emit: beadScores

	script:
"""
	filter_beads.py --barcodesDir . --minUTC $params.minUTC --threads $task.cpus --memory '$task.memory' --sampleBarcode $sampleBarcode
"""
}

// Run cell finding algorithm and generate filtered expression matrix,
// also generate sample-level statistics for use in reporting module.
// libCount is the number of libraries associated with a specific sampleName.
process FilteredMatrixGeneration {
	tag "$id"
	label 'large'

	// If libCount > 1, we organize outputs by library
	publishDir { file(params.outputDir) / samplesDir }, pattern:"${id}_filtered_star_output", saveAs: { "${id}.filtered.matrix" }, mode: 'copy'
	// allCells.csv will be published after scaleplex assignment if --scalePlex is set
	publishDir { file(params.outputDir).resolve(samplesDir) }, pattern: "*/${id}_allCells.csv", saveAs: { "${id}.allCells.csv" }, mode: 'copy'
	publishDir { file(params.outputDir).resolve(samplesDir) }, pattern: "*/${id}_allBarcodes.parquet", saveAs: { "${id}.allBarcodes.parquet" }, mode: 'copy'
	publishDir { file(params.outputDir).resolve(samplesDir) }, pattern: "*/${id}_cellcalling_stats.csv", saveAs: { "internal/${id}_cellcalling_stats.csv" }, mode: 'copy', enabled: params.internalReport

	input:
	tuple val(id), val(sampleName), val(libName), val(expectedCells), path("Solo.out"), val(libCount), path("allBarcodes.parquet"), path("bead_scores.parquet")
	val(isBarnyard)
	val(quantum)
	
	output:
	tuple(val(id), val(libName), path("${id}_metrics/${id}_allCells.csv"), emit: cell_metrics)
	tuple(val(id), val(libName), path("${id}_metrics/${id}_allBarcodes.parquet"), emit: all_barcodes_metrics)
	tuple(val(id), path("${id}_filtered_star_output"), emit: filtered_star_mtx)
	tuple(val(id), path("${id}_metrics/${id}_cellcalling_stats.csv"), emit: stats)

	script:
	samplesDir = "samples"
	if (libCount > 1) {
		samplesDir += "/${sampleName}_libraries"
	}
	opts = ""
	starMultiParam = params.starMulti
	if (params.cellFinder) {
		if (!isBarnyard) {
			opts = opts + "--cellFinder --FDR $params.cellFinderFdr --medianFraction $params.medianFraction "
		} else {
			log.warn("CellFinder is not supported for barnyard samples. Default cell thresholding will be used.")
		}
		if(params.fixedCells){
			log.warn("Both of the parameters fixedCells and cellFinder are set to true. CellFinder will be used for cell calling.")
		}
	}
	if (params.filterOutliers){
		opts = opts + "--filter_outliers "
	}
	if (params.internalReport) {
		opts = opts + "--internalReport "
	}
	if (isBarnyard) {
		opts = opts + "--isBarnyard "
		starMultiParam = params.starMultiBarnyard
	}
	if (params.fixedCells) {
		opts = opts + "--fixedCells "
	}
	if (quantum){
		opts = opts + "--isQuantum --beadScores bead_scores.parquet --minDivergence ${params.minBeadDivergence} "
	}
	if (params.filterAmbientBeads){
		opts = opts + "--filterBeads "
	}
	if (params.roundCounts){
		opts = opts + "--roundCounts "
	}	
	matrixFn = Utils.starMatrixFn(starMultiParam)
"""
	export POLARS_MAX_THREADS=$task.cpus
	call_cells.py --sampleMetrics allBarcodes.parquet --STARsolo_out Solo.out --feature_type $params.starFeature --matrix_type $matrixFn \
	--sample $id --expectedCells $expectedCells --topCellPercent $params.topCellPercent --minCellRatio $params.minCellRatio \
	--minUTC $params.minUTC --UTC $params.UTC \
	--madsReads $params.madsReads --madsPassingReads $params.madsPassingReads --madsMito $params.madsMito $opts \
	--threads $task.cpus --memory '$task.memory'
	pigz -p $task.cpus ${id}_filtered_star_output/matrix.mtx
"""
}

// Generate sample-level statistics
process SampleStats {
	tag "$id"
	label 'large'

	input:
	tuple val(id), val(libName), path(allBarcodes), val(totalSampleReads), path("Log.final.out")
	val(isBarnyard)
	
	output:
	tuple val(id), path("${id}_metrics/${id}_sample_stats.csv"), emit: sample_metrics

	script:
	opts = ""
	if (params.cellFinder && !isBarnyard) {
		opts = opts + "--cellFinder "
	}
	if (params.filterOutliers){
		opts += "--filterOutliers "
	}
	if (params.internalReport) {
		opts += "--internalReport "
	}
	if (isBarnyard) {
		opts += "--isBarnyard "
	}
"""
	sample_stats.py --sample $id --cellMetrics $allBarcodes --starLog Log.final.out --totalSampleReads $totalSampleReads \
		--threads $task.cpus --memory '$task.memory' $opts
"""
}


workflow CELL_CALLING {
take:
	samplesDict // each row from samples.csv parsed into a channel
	libraryInfo // map of library information (lib jsons, adapters, etc.)
	soloResults // STARSolo outputs for each sample
	isBarnyard // Flag for reporting
	soloLog // STAR log files
	totalSampleReadsBySampleID // Total sample reads computed from bcParser metrics
main:

	libStructName = libraryInfo.rnaLibraryStructureFile.getName()
	libStructDir = libraryInfo.rnaLibraryStructureFile.getParent()
	// Convert 'expectedCells' column to int
	samplesDict = samplesDict.map{ it.plus(['expectedCells': Utils.toIntOr0(it.expectedCells)])}
	soloLogResults = soloLog

	samplesDict.map{ [it.id, it.sample, it.libName, it.expectedCells] }
	.join(soloResults)
	.set{ cellMetricsGenerationInput }

	// libCount determines the number of libNames associated with each sampleName in the samplesheet.
	// This allows us to determine whether or not multiple libraries are being merged.
	// libCount -> (id, sample, libName)
	samplesDict
		.filter{ it.rnaId == null } // Only consider RNA samples
		.map{ tuple(it.id, it.sample, it.libName) }
		.groupTuple(by: 1) // Gather all libNames for each sample
		.map{ tuple(it[0], it[2].size()) }.transpose()
		.set { libCount }

	// cellMetricsGenerationInput -> ([id, sample, group]_name, library_name, expectedCells, star_output, libCount)
	CellMetricsGeneration(cellMetricsGenerationInput, libStructName, isBarnyard, libStructDir)
	if (libraryInfo.quantum) {
		BeadFiltering(
			CellMetricsGeneration.out.raw_cell_metrics
				.groupTuple(by: 1) // collect by library
				.map{ _id, lib, metrics -> [lib, metrics] },
			libraryInfo.sampleBarcode
		)
		CellMetricsGeneration.out.raw_cell_metrics
			.combine(
				BeadFiltering.out.beadScores
					// add null id to match cardinality of raw_cell_metrics
					.map{ libName, beadScores -> [null, libName, beadScores] },
				by: 1,
			)
			.map{
				_lib, id, metrics, _null, beadScores ->
					[id, metrics, beadScores]

			}
			.dump(tag:'cellMetricsGenerationOutput')
			.set{ cellMetricsGenerationOutput }
	} else {
		CellMetricsGeneration.out.raw_cell_metrics
			.map{ id, _lib, metrics -> [id, metrics, []] } // drop lib, add empty beadScores
			.dump(tag:'cellMetricsGenerationOutput')
			.set{ cellMetricsGenerationOutput }
	}
	FilteredMatrixGeneration(
		cellMetricsGenerationInput
			.join(libCount)
			.join(cellMetricsGenerationOutput),
		isBarnyard,
		libraryInfo.quantum
	)
	SampleStats(
		FilteredMatrixGeneration.out.all_barcodes_metrics
			.join(totalSampleReadsBySampleID)
			.join(soloLogResults),
		isBarnyard
	)
		
	
emit:
	allCells = FilteredMatrixGeneration.out.cell_metrics
	allBarcodes = FilteredMatrixGeneration.out.all_barcodes_metrics
	sampleStats = SampleStats.out.sample_metrics
	filteredMtx = FilteredMatrixGeneration.out.filtered_star_mtx
	libCount = libCount
	beadScores = libraryInfo.quantum ? BeadFiltering.out.beadScores : Channel.empty() // Only emit bead scores if quantum is true
	cellCallingStats = FilteredMatrixGeneration.out.stats
}
