// Return 0 for null or non-integer strings
def toIntOr0(str) {
	if ((str != null) && str.isInteger())
		return str as int;
	else
		return 0
}

// Generate report for each sample from metrics generated in cellMetricsGeneration process
process sampleReportGeneration {
input:
	tuple(val(id), val(libName), path("${id}_metrics/allCells.csv"), val(sampleInfo), path("trim_stats*.json"), val(libCount), path("${id}_metrics/sample_stats.csv"))
	val(libStructName)
	path(libStructDir)
	val(isBarnyard)
output:
	path("$outDir/${sampleId}.report.html")
	path("$outDir/csv/${sampleId}.reportStatistics.csv"), emit: stats
	path("$outDir/figures_internal/*"), optional: true
	path("$outDir/csv/${sampleId}_unique_transcript_counts_by_RT_Index_well.csv")
	path("$outDir/csv/${sampleId}_num_cells_by_RT_Index_well.csv")
publishDir params.outDir, mode: 'copy'
label 'optional'
label 'report'
tag "$id"
script:
	opts = ""
	sampleId = id
	outDir = "reports"
	if (params.mergedSamples) {
		sampleId = id + ".merged"
		sampleName = id
	} else {
		// This regex splits a string at (.). In this process sampleId == PBMC1.ScaleRNA.
		// The regex splits PBMC1.ScaleRNA into PBMC1 and ScaleRNA
		sampleName = sampleId =~ /(.+)\.(.+)/
		sampleName = sampleName[0][1]
	}
	if (libCount > 1) {
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
publishDir file(params.outDir) / "reports", mode: 'copy'
label 'optional'
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
	libJson  // library structure json file
	trimFqStats // list of cutadapt stats files for each sample (can be empty)
	isBarnyard // Flag for reporting
	allCells // allCells.csv generated by callCells
	libCount // value equal to number of sequencing libraries for sample
	metricsSampleLevel // sample-level metrics generated by callCells
main:
	trimFqStats.dump(tag:'sampleReport/trimFqStats')
	allCells.dump(tag:'sampleReport/allCells')

	// Convert 'expectedCells' column to int
	samples = samples.map{ it.plus(['expectedCells': toIntOr0(it.expectedCells)]) }

	if (params.mergedSamples) { 
		// starSoloOuts is already merged by group. We just need to create a combined samples.csv entry per group
		samples.map{ [it.group, it] }
			.groupTuple() // All samples grouped for merging, with each column as a list of values
			.map{ // Collapse sample.csv entries for each merged sample (group)
				['id':it[0], 'sample':it[1].sample[1], 'libName':'merged', 
				'expectedCells':it[1].expectedCells.sum(),
				'barcodes': Utils.combineSampleBarcodes(it[1].barcodes)] }
			.dump(tag:'sampleReport/mergedSamples')
			.set{ samples }
	}
	sampleInfo = samples.map{ [it.id, it] } // SampleId -> Dict with all columns from samples.csv
	sampleInfo.dump(tag:'sampleReport/samples')
	
	// sampleMetrics -> [id, libName, allCells.csv, sampleInfo, trimStatsJsonPath, libCount, sampleLevelMetrics]
	allCells
		.join(sampleInfo)
		.join(trimFqStats)
		.join(libCount)
		.join(metricsSampleLevel)
		.dump(tag:'sampleMetrics')
		.set{ sampleMetrics }
	sampleReportGeneration(sampleMetrics, libJson.getName(), libJson.getParent(), isBarnyard)

emit:
	sampleStats = sampleReportGeneration.out.stats // metrics csv per sample
}
