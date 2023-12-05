// Return 0 for null or non-integer strings
def toIntOr0(str) {
	if ((str != null) && str.isInteger())
		return str as int;
	else
		return 0
}

// Generate metrics per sample from STARsolo output
// libCount is the number of libraries associated with a specific sampleName.
// This helps us determine if samples are being merged. libCount > 1 == merging. libCount = 1 == no merging.
process sampleMetricsGeneration {
input:
	tuple(val(id), val(sampleName), val(libName), val(expectedCells), path("Solo.out"), val(libCount))
	val(libStructName)
	val(isBarnyard)
	path(libStructDir)
output:
	tuple(val(id), val(libName), path("${id}_metrics/allCells.csv"), emit: cell_metrics)
	tuple(val(id), path("${id}_filtered_star_output"), emit: filtered_star_mtx)

// If libCount is equal to 1 publishDir = $params.OutDir/samples
// If libCount is not equal to 1 publishDir = $params.outDir/samples/<sampleName>
// If params.mergedSamples is true files are renamed <sampleName>.merged.filtered.matrix or <sampleName>.merged.allCells.csv
// If params.mergedSamples if false files are name <sampleName>.<libName>.filtered.matrix or <sampleName>.<libName>.allCells.csv
publishDir path: {libCount == 1 ? "${params.outDir}/samples/" : "${params.outDir}/samples/${sampleName}_libraries/"}, pattern:"${id}_filtered_star_output", saveAs: {params.mergedSamples ? "${id}.merged.filtered.matrix" : "${id}.filtered.matrix"}, mode: 'copy'
publishDir path: {libCount == 1 ? "${params.outDir}/samples/" : "${params.outDir}/samples/${sampleName}_libraries/"}, pattern: "*/allCells.csv", saveAs: {params.mergedSamples ? "${id}.merged.allCells.csv" : "${id}.allCells.csv"}, mode: 'copy'
label 'report'
tag "$id"
script:
	opts = ""
	if (isBarnyard) {
		opts = opts + "--isBarnyard "
	}
	if (!params.useSTARthreshold) {
		opts = opts + "--calcCellThreshold "
	}
	libStruct = "${libStructDir}/${libStructName}"
	matrixFn = Utils.starMatrixFn(params.starMulti)
"""
	getSampleMetrics.py --sample $id --libStruct $libStruct \
	--topCellPercent $params.topCellPercentage --minCellRatio $params.minCellRatio --minReads $params.minReads \
	--star_out Solo.out --starFeature $params.starFeature --starMatrix $matrixFn --expectedCells $expectedCells \
	$opts
"""
}

// Generate report for each sample from metrics generated in sampleMetricsGeneration process
process sampleReportGeneration {
input:
	tuple(val(id), val(libName), path("${id}_metrics/allCells.csv"), val(sampleInfo), path("trim_stats*.json"), val(libCount))
	val(libStructName)
	path(libStructDir)
	val(isBarnyard)
output:
	path("$outDir/${sampleId}.report.html")
	path("$outDir/csv/${sampleId}.reportStatistics.csv"), emit: stats
	path("$outDir/${sampleId}_figures"), optional: true
	path("$outDir/csv/${sampleId}_unique_transcript_counts_by_RT_Index_well.csv")
	path("$outDir/csv/${sampleId}_num_cells_by_RT_Index_well.csv")
publishDir path: "$params.outDir", mode: 'copy'
errorStrategy 'ignore'
label 'report'
tag "$id"
script:
	opts = ""
	sampleId = id
	outDir = "reports"
	if(params.mergedSamples){
		sampleId = id + ".merged"
		sampleName = id
	} else {
		// This regex splits a string at (.). In this process sampleId == PBMC1.ScaleRNA.
		// The regex splits PBMC1.ScaleRNA into PBMC1 and ScaleRNA
		sampleName = sampleId =~ /(.+)\.(.+)/
		sampleName = sampleName[0][1]
	}
	if(libCount > 1){
		outDir = outDir + "/" + sampleName + "_libraries"
	}
	if (isBarnyard) {
		opts = opts + "--isBarnyard "
    }
	if (params.internalReport) {
		opts = opts + "--internalReport "
	}
	if (sampleInfo.barcodes) {
		opts = opts + "--barcodes '$sampleInfo.barcodes' "
	}
	if (sampleInfo.libName) {
		opts = opts + "--libName $sampleInfo.libName "
	}
	libStruct = "${libStructDir}/${libStructName}"
	// trimStats is optional. If empty (no files), skip the --tripStats option
"""
	export TMPDIR=\$(mktemp -p `pwd` -d)
	if test -n "\$(shopt -s nullglob; echo trim_stats*)"
	then
		trimOpt="--trim_stats trim_stats*" 
	else
		trimOpt=""
	fi
	generateSampleReport.py --outDir $outDir --sampleId $sampleId --sampleName $sampleName --libraryStruct $libStruct --sampleMetrics ${id}_metrics \$trimOpt \
	--workflowVersion $workflow.manifest.version $opts
"""
}

process multiSampleReport {
// Generate multi sample metric .csv from per-sample csv's
// Called once on all samples in a run (from outside `sampleReport` workflow)
// Takes output from sampleReportGeneration process as input.
input:
	path(sampleStats)
output:
	path(outFn)
publishDir "$params.outDir/reports", mode: 'copy'
errorStrategy 'ignore'
label 'report'
script:
	outFn = "allSamples.reportStatistics.csv"
"""
	mergeReportStats.py $sampleStats > $outFn
"""
}

// Run sample metrics generation and per-sample outputs (.mtx, HTML report, ...)
// This is run separately for original samples and merged (grouped) samples when using --merge
workflow sampleReport {
take:
	samples		// each row from samples.csv parsed into a channel
	samplesCsv  // samples.csv file
	libJson  // library structure json file
    starSoloOuts // STARsolo outputs for each sample (or sample group if --merge)
	trimFqStats // list of cutadapt stats files for each sample (can be empty list)
	isBarnyard // Flag for reporting

main:
	starSoloOuts.dump(tag:'sampleReport/soloOuts')
	trimFqStats.dump(tag:'sampleReport/trimFqStats')

	// Convert 'expectedCells' column to int
	samples = samples.map{ it.plus(['expectedCells': toIntOr0(it.expectedCells)]) }

	if (params.mergedSamples){ 
		// starSoloOuts is already merged by group. We just need to create a conbined samples.csv entry per group
		sampleGroups = samples.map{ [it.group, it] }.groupTuple() // All samples grouped for merging, with each column as a list of values
		samples = sampleGroups.map{ // Collapse sample.csv entries for one merged sample 
			['id':it[0], 'sample':it[1].sample[1], 'libName':'merged', 'expectedCells':it[1].expectedCells.sum()] }
		samples.dump(tag:'sampleReport/mergedSamples')
	}
	sampleInfo = samples.map{ [it.id, it] } // SampleId -> Dict with all columns from samples.csv
	// sampleStatsInput -> (id, sample, libName, expectedCells, starSoloOuts)
	sampleStatsInput = samples.map{ [it.id, it.sample, it.libName, it.expectedCells] }.join(starSoloOuts)
	
	// libCount determines the number of libNames associated with each sampleName in the samplesheet.
	// This allows us to determine whether or not multiple libraries are being merged.
	// libCount -> (id, sample, libName)
	libCount = sampleStatsInput.map{tuple(it[0], it[1], it[2])}
	// libCount -> (id, [sample], [libName])
	libCount = libCount.groupTuple(by : 1)
	// libCount -> (id, numberOfLibraries)
	libCount = libCount.map{tuple(it[0], it[2].size())}.transpose()
	// sampleStatsInput -> (id, sample, libName, expectedCells, starSoloOuts, libCount)
	sampleStatsInput = sampleStatsInput.join(libCount)

	sampleInfo.dump(tag:'sampleReport/samples')
	sampleStatsInput.dump(tag:'sampleReport/sampleStatsInput')

	// sampleStatsInput -> ([id, sample, group]_name, library_name, expectedCells, star_output)
	sampleStatsInput.dump(tag:'sampleStatsInput')
	sampleMetricsGeneration(sampleStatsInput, libJson.getName(), isBarnyard, libJson.getParent())

	sampleMetrics = sampleMetricsGeneration.out.cell_metrics.join(sampleInfo)
		.join(trimFqStats, remainder:true).map { [*it[0..-2], it[-1] ?: []] } // If no trimStats (null), pass an empty list
	sampleMetrics = sampleMetrics.join(libCount)
	sampleMetrics.dump(tag:'sampleMetrics')
	sampleReportGeneration(sampleMetrics, libJson.getName(), libJson.getParent(), isBarnyard)

emit:
	cellMetrics = sampleMetricsGeneration.out.cell_metrics // allCells.csv per sample
	filtered_star_mtx = sampleMetricsGeneration.out.filtered_star_mtx // mtx directory per sample
	sampleStats = sampleReportGeneration.out.stats // metrics csv per sample
}
