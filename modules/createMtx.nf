// Return 0 for null or non-integer strings
def toIntOr0(str) {
	if ((str != null) && str.isInteger())
		return str as int;
	else
		return 0
}

// Generate cell-level metrics per sample from STARsolo output
process cellMetricsGeneration {
input:
	tuple(val(id), val(sampleName), val(libName), val(expectedCells), path("Solo.out"))
	val(libStructName)
	val(isBarnyard)
	path(libStructDir)
output:
	tuple(val(id), path("${id}_metrics/${id}_allCells.csv"), emit: raw_cell_metrics)
label "report"
tag "$id"
script:
	opts = ""
	if (isBarnyard) {
		opts = opts + "--isBarnyard "
	}
	if (params.mergedSamples){
		opts = opts + "--isMerge "
	}
	libStruct = "${libStructDir}/${libStructName}"
	matrixFn = Utils.starMatrixFn(params.starMulti)
"""
	buildAllCells.py --STARsolo_out Solo.out --feature_type $params.starFeature --matrix_type $matrixFn \
	--sample $id --libStruct $libStruct $opts
"""
}

// Run cell finding algorithm and generate filtered expression matrix,
// also generate sample-level statistics for use in reporting module.
// libCount is the number of libraries associated with a specific sampleName.
// This helps us determine if samples are being merged. libCount > 1 == merging. libCount = 1 == no merging.
process filteredMatrixGeneration {
tag "$id"
// If libCount is equal to 1 publishDir = $params.OutDir/samples
// If libCount is not equal to 1 publishDir = $params.outDir/samples/<sampleName>
// If params.mergedSamples is true files are renamed <sampleName>.merged.filtered.matrix or <sampleName>.merged.allCells.csv
// If params.mergedSamples if false files are name <sampleName>.<libName>.filtered.matrix or <sampleName>.<libName>.allCells.csv
publishDir path: {libCount == 1 ? "${params.outDir}/samples/" : "${params.outDir}/samples/${sampleName}_libraries/"}, pattern:"${id}_filtered_star_output", saveAs: {params.mergedSamples ? "${id}.merged.filtered.matrix" : "${id}.filtered.matrix"}, mode: 'copy'
publishDir path: {libCount == 1 ? "${params.outDir}/samples/" : "${params.outDir}/samples/${sampleName}_libraries/"}, pattern: "*/${id}_allCells.csv", saveAs: {params.mergedSamples ? "${id}.merged.allCells.csv" : "${id}.allCells.csv"}, mode: 'copy'

input:
	tuple(val(id), val(sampleName), val(libName), val(expectedCells), path("Solo.out"), val(libCount), path("allCells.csv"))
	
output:
	tuple(val(id), val(libName), path("${id}_metrics/${id}_allCells.csv"), emit: cell_metrics)
	tuple(val(id), path("${id}_metrics/${id}_sample_stats.csv"), emit: sample_metrics)
	tuple(val(id), path("${id}_filtered_star_output"), emit: filtered_star_mtx)
	tuple(val(id), val(libCount), emit: libCount)

script:
	opts = ""
	if (params.cellFinder) {
		opts = opts + "--cellFinder "
		if(params.fixedCells){
			log.warn("Both of the parameters fixedCells and cellFinder are set to true. CellFinder will be used for cell calling.")
		}
	}
	if (params.filter_outliers){
		opts = opts + "--filter_outliers "
	}
	if (params.internalReport) {
		opts = opts + "--internalReport "
	}
	if (params.fixedCells){
		opts = opts + "--fixedCells "
	}
	matrixFn = Utils.starMatrixFn(params.starMulti)
"""
	callCells.py --sampleMetrics allCells.csv --STARsolo_out Solo.out --feature_type $params.starFeature --matrix_type $matrixFn \
	--sample $id --expectedCells $expectedCells --topCellPercent $params.topCellPercent --minCellRatio $params.minCellRatio \
	--minUTC $params.minUTC --UTC $params.UTC --FDR $params.cellFinderFdr --num_mad_genes $params.num_mad_genes \
	--num_mad_umis $params.num_mad_umis --num_mad_mito $params.num_mad_mito $opts
"""
}


// Merge STARSolo outputs (.mtx and CellReads Metrics)
// Used with --merge
// Assumes different cells in each subsample
process mergeSamples {
input:
	tuple( val(samples), val(group), path(starDirs) )
output:
	tuple(val(group), path(outDir), emit: mergeOut)
publishDir "$params.outDir/alignment", mode: 'copy'
label "report"
tag "$group"
script:
	sampleIDs = String.join(" ", samples)
	matrixFn = Utils.starMatrixFn(params.starMulti)
	outDir = "${group}.merged/${group}.merged.star.solo/"
"""
	mergeRawSTARoutput.py --star_dirs $starDirs --star_feature $params.starFeature --star_matrix $matrixFn --sample_ids $sampleIDs --out_dir $outDir
"""
}

workflow filterMtx {
take:
	samplesDict // each row from samples.csv parsed into a channel
	libJson // library structure json file
	soloResults // STARSolo outputs for each sample (or sample group if --merge)
	filterBarnyard // Flag for reporting
main:

	libStructName = libJson.getName()
	libStructDir = libJson.getParent()
	// Convert 'expectedCells' column to int
	samplesDict = samplesDict.map{ it.plus(['expectedCells': toIntOr0(it.expectedCells)])}


	if(params.mergedSamples) {
		// soloResults is already merged by group. We just need to create a conbined samples.csv entry per group
		samplesDict.map{ tuple(it.id, it.group) }
		.join(soloResults)
		.groupTuple(by:1) // All samples grouped for merging, with each column as a list of values
		.filter{ it[0].size() > 1 } // Only merge samples with multiple libraries
		.set{ mergeSamplesInput }
		mergeSamples(mergeSamplesInput)

		sampleGroupsDict = samplesDict.map{[it.group, it]}.groupTuple()
		// Collapse sample.csv entries for one merged sample
		samplesDict = sampleGroupsDict.map{ ['id':it[0], 'sample':it[1].sample[1], 'libName':'merged', 'expectedCells':it[1].expectedCells.sum()] }

		// Set soloResults as merged STARSolo output
		soloResults = mergeSamples.out.mergeOut
	}
	
	samplesDict.map{ [it.id, it.sample, it.libName, it.expectedCells] }
	.join(soloResults)
	.set{ cellMetricsGenerationInput }

	// libCount determines the number of libNames associated with each sampleName in the samplesheet.
	// This allows us to determine whether or not multiple libraries are being merged.
	// libCount -> (id, sample, libName)
	libCount = samplesDict.map{tuple(it.id, it.sample, it.libName)}
	// libCount -> (id, [sample], [libName])
	libCount = libCount.groupTuple(by : 1)
	// libCount -> (id, numberOfLibraries)
	libCount = libCount.map{tuple(it[0], it[2].size())}.transpose()

	// cellMetricsGenerationInput -> ([id, sample, group]_name, library_name, expectedCells, star_output, libCount)
	cellMetricsGeneration(cellMetricsGenerationInput, libStructName, filterBarnyard, libStructDir)
	filteredMatrixGeneration(cellMetricsGenerationInput.join(libCount).join(cellMetricsGeneration.out.raw_cell_metrics))
	
emit:
	allCells = filteredMatrixGeneration.out.cell_metrics
	sampleStats = filteredMatrixGeneration.out.sample_metrics
	filteredMtx = filteredMatrixGeneration.out.filtered_star_mtx
	libCount = filteredMatrixGeneration.out.libCount
	rawMtx = soloResults
}
