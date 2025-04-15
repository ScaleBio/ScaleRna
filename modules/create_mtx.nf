/*
* Compute per cell metrics from STARsolo output
*
* Processes:
*     CellMetricsGeneration
*     FilteredMatrixGeneration
*     MergeSamples
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
	val(mergedSamples)
	
	output:
	tuple(val(id), path("${id}_metrics/${id}_allBarcodes.parquet"), emit: raw_cell_metrics)

	script:
	opts = ""
	starMultiParam = params.starMulti
	if (isBarnyard) {
		opts = opts + "--isBarnyard "
		starMultiParam = params.starMultiBarnyard
	}
	if (mergedSamples){
		opts = opts + "--isMerge "
	}
	libStruct = "${libStructDir}/${libStructName}"
	matrixFn = Utils.starMatrixFn(starMultiParam)
"""
	export POLARS_MAX_THREADS=$task.cpus
	build_all_cells.py --STARsolo_out Solo.out --feature_type $params.starFeature --matrix_type $matrixFn \
	--sample $id --libStruct $libStruct $opts --threads $task.cpus --memory '$task.memory'
"""
}

// Run cell finding algorithm and generate filtered expression matrix,
// also generate sample-level statistics for use in reporting module.
// libCount is the number of libraries associated with a specific sampleName.
process FilteredMatrixGeneration {
	tag "$id"
	label 'large'

	// If libCount > 1, we organize outputs by library
	// merged samples are published as  <sampleName>.merged.*
	publishDir { file(params.outputDir) / samplesDir }, pattern:"${id}_filtered_star_output", saveAs: { mergedSamples ? "${id}.merged.filtered.matrix" : "${id}.filtered.matrix" }, mode: 'copy'
	// allCells.csv will be published after scaleplex assignment if --scalePlex is set
	publishDir { file(params.outputDir).resolve(samplesDir) }, pattern: "*/${id}_allCells.csv", saveAs: { mergedSamples ? "${id}.merged.allCells.csv" : "${id}.allCells.csv" }, mode: 'copy', enabled: !params.scalePlex
	publishDir { file(params.outputDir).resolve(samplesDir) }, pattern: "*/${id}_allBarcodes.parquet", saveAs: { mergedSamples ? (params.internalReport ? "${id}.merged.allBarcodes.parquet" : null) : "${id}.allBarcodes.parquet" }, mode: 'copy'

	input:
	tuple(val(id), val(sampleName), val(libName), val(expectedCells), path("Solo.out"), val(libCount), path("allBarcodes.parquet"), val(totalSampleReads))
	val(isBarnyard)
	val(mergedSamples)
	
	output:
	tuple(val(id), val(libName), path("${id}_metrics/${id}_allCells.csv"), emit: cell_metrics)
	tuple(val(id), val(libName), path("${id}_metrics/${id}_allBarcodes.parquet"), emit: all_barcodes_metrics)
	tuple(val(id), path("${id}_metrics/${id}_sample_stats.csv"), emit: sample_metrics)
	tuple(val(id), path("${id}_filtered_star_output"), emit: filtered_star_mtx)
	tuple(val(id), val(libCount), emit: libCount)

	script:
	samplesDir = "samples"
	if (libCount > 1) {
		samplesDir += "/${sampleName}_libraries"
	}
	opts = ""
	starMultiParam = params.starMulti
	if (params.cellFinder) {
		if (!isBarnyard) {
			opts = opts + "--cellFinder "
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
	if (params.quantum){
		opts = opts + "--isQuantum "
	}
	if (params.filterBeads){
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
	--minUTC $params.minUTC --UTC $params.UTC --FDR $params.cellFinderFdr --madsReads $params.madsReads --totalSampleReads $totalSampleReads \
	--madsPassingReads $params.madsPassingReads --madsMito $params.madsMito --maxBeadBcs $params.maxBeadBcs $opts --threads $task.cpus --memory '$task.memory'
	pigz -p $task.cpus ${id}_filtered_star_output/matrix.mtx
	"""
}


// Merge STARSolo outputs (.mtx and CellReads Metrics)
// Used with --merge
// Assumes different cells in each subsample
process MergeSamples {
	tag "$group"
	label 'large'

	publishDir file(params.outputDir) / "alignment", mode: 'copy'

	input:
	tuple( val(samples), val(group), path("Solo.out*"), path("log.final.out*") )
	val(isBarnyard)
	
	output:
	tuple(val(group), path(outDir), emit: mergeOut)
	tuple(val(group), path("$logOutDir/Log.final.out"), emit: merge_log)

	script:
	starMultiParam = params.starMulti
	if(isBarnyard){
		starMultiParam = params.starMultiBarnyard
	}
	sampleIDs = String.join(" ", samples)
	matrixFn = Utils.starMatrixFn(starMultiParam)
	outDir = "${group}.merged/${group}.merged.star.solo/"
	logOutDir = "${group}.merged/${group}.merged.star.align"
	"""
	merge_raw_star_output.py --star_dirs Solo.out* --star_feature $params.starFeature --star_matrix $matrixFn --sample_ids $sampleIDs --out_dir $outDir --log_out_dir $logOutDir --star_log log.final.out*
	pigz -p $task.cpus $outDir/${params.starFeature}/raw/$matrixFn
	"""
}

workflow FILTER_MTX {
take:
	samplesDict // each row from samples.csv parsed into a channel
	libJson // library structure json file
	soloResults // STARSolo outputs for each sample (or sample group if --merge)
	isBarnyard // Flag for reporting
	soloLog // STAR log files
	totalSampleReadsBySampleID // Total sample reads computed from bcParser metrics
	mergedSamples // If merged run
main:

	libStructName = libJson.getName()
	libStructDir = libJson.getParent()
	// Convert 'expectedCells' column to int
	samplesDict = samplesDict.map{ it.plus(['expectedCells': Utils.toIntOr0(it.expectedCells)])}
	soloLogResults = soloLog

	if (mergedSamples) {
		samplesDict
			.map{ tuple(it.id, it.group) }
			.join(totalSampleReadsBySampleID)
			// Group by it.group
			.groupTuple(by:1)
			.map { it -> // [[sample1.group, sample2.group], group, [sample1 total_reads, sample2 total_reads]]
				def total_reads = 0
				// Sum up total reads of samples belonging to same group
				it[2].each { reads ->
					total_reads += reads.toLong()
				}
				tuple(it[1], total_reads)
			}
			.dump(tag:'mergedTotalSampleReadsBySampleID')
			.set{ totalSampleReadsBySampleID }
		// soloResults is already merged by group. We just need to create a conbined samples.csv entry per group
		samplesDict
			.map{ tuple(it.id, it.group) }
			.join(soloResults)
			.join(soloLog)
			.groupTuple(by:1) // All samples grouped for merging, with each column as a list of values
			.filter{ it[0].size() > 1 } // Only merge samples with multiple libraries
			.set{ mergeSamplesInput }
		MergeSamples(mergeSamplesInput, isBarnyard)

		sampleGroupsDict = samplesDict.map{[it.group, it]}.groupTuple()
		// Collapse sample.csv entries for one merged sample
		samplesDict = sampleGroupsDict.map{ ['id':it[0], 'sample':it[1].sample[1], 'libName':'merged', 'expectedCells':it[1].expectedCells.sum()] }

		// Set soloResults as merged STARSolo output
		soloResults = MergeSamples.out.mergeOut
		soloLogResults = MergeSamples.out.merge_log
	}

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
	CellMetricsGeneration(cellMetricsGenerationInput, libStructName, isBarnyard, libStructDir, mergedSamples)
	FilteredMatrixGeneration(cellMetricsGenerationInput.join(libCount).join(CellMetricsGeneration.out.raw_cell_metrics).join(totalSampleReadsBySampleID), isBarnyard, mergedSamples)
	
emit:
	allCells = FilteredMatrixGeneration.out.cell_metrics
	allBarcodes = FilteredMatrixGeneration.out.all_barcodes_metrics
	sampleStats = FilteredMatrixGeneration.out.sample_metrics
	filteredMtx = FilteredMatrixGeneration.out.filtered_star_mtx
	libCount = FilteredMatrixGeneration.out.libCount
	rawMtx = soloResults
	soloLog = soloLogResults
}
