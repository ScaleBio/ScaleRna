nextflow.enable.dsl=2

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include { inputReads } from './inputReads'

// Load .json file into object
// (reference genomes, library design, etc.)
def loadJson(json) {
	def jsonSlurper  = new JsonSlurper()
	return jsonSlurper.parse( json )
}

def starMatrixFn(multimappers) {
	if (multimappers == "Unique") {
		return "matrix.mtx"
	} else {
		return "UniqueAndMult-${multimappers}.mtx"
	}
}

// Create a 'file()' from string 'path'
// 'path' could be null, a s3 URL, an absolute file-path or relative to 'baseDir'
def expandPath(path, baseDir) {
	if (path == null) { return null}
	if (path =~ /^s3/) { return file(path)}
	return baseDir.resolve(path)
}

// Return 0 for null or non-integer strings
def toIntOr0(str) {
	if ((str != null) && str.isInteger())
		return str as int;
	else
		return 0
}

// Reference genome files and parameters
// Paths can be absolute or relative to location of json
def loadGenome(json) {
	def baseDir = json.getParent()
	genome = loadJson(json)
	genome.star_index = expandPath(genome.star_index, baseDir)
	genome.gtf = expandPath(genome.gtf, baseDir)
	return genome
}

// Load per-sample read-counts from bcParser output
def loadDemuxReadCounts(demuxMetrics) {
	def jsonSlurper  = new JsonSlurper()
	def counts = []
		json = jsonSlurper.parse(demuxMetrics)
		for (sample in json["samples"]) {
			counts.add(tuple(sample.value["name"], sample.value["reads"][0]))
		}
	return counts
}

// Prepare samples.csv with defaults, rename legacy columns, etc.
process regularizeSamplesCsv {
input: 
	path("samples.in.csv")
output: 
	path("samples.csv")
publishDir "${params.outDir}", mode: 'copy'
label 'small'
script:
	opts=""
	if (params.splitFastq) {
		if (params.runFolder == null) {
			opts += "--splitSample"
		}
	}
"""
	regularizeSamplesCsv.py samples.in.csv $opts > samples.csv
"""
}

// Run the STAR executable on demuxed fastq files
process starsolo {
input: 
	path(indexDir) // STAR genome index
	val(library) // library.json (contains STAR barcode parameters)
	tuple(val(sample), val(count), path("transcript*.fastq.gz"), path("barcode*.fastq.gz")) // transcript and BC fastq file
output: 
	tuple(val(sample), path("$starDir/Aligned.sortedByCoord.out.bam*"), emit: bam, optional: true)
	tuple(val(sample), path("$starDir/Log*"), emit: log)
	tuple(val(sample), path("${sample}.${count}.star.solo"), emit: solo)
publishDir "$pubDir", pattern: "*.solo", mode: 'copy'
publishDir "$pubDir", pattern: "$starDir/*bam*"
publishDir "$pubDir", pattern: "$starDir/Log*", mode: 'copy'
tag "$sample"
script:
	if (params.splitFastq) {
		pubDir = "${params.outDir}/alignment/split"
		starDir = "${sample}.${count}.star.align"
	} else {
		pubDir = "${params.outDir}/alignment"
		starDir = "${sample}.star.align"
	}
	barcodeParam = library["star_barcode_param"]
	if (params.bamOut) {
		bamOpts = "--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM --outSAMunmapped Within"
    } else {
        bamOpts = "--outSAMtype None" 
    }
"""
	STAR --runThreadN $task.cpus --genomeDir $indexDir \
	$barcodeParam ${params.starTrimming} \
    $bamOpts --outSJtype None --soloCellReadStats Standard \
	--soloStrand ${params.starStrand} --soloFeatures ${params.starFeature} --soloMultiMappers ${params.starMulti} \
	--readFilesIn <(cat transcript*.fastq.gz) <(cat barcode*.fastq.gz) --readFilesCommand zcat \
	--outFileNamePrefix $starDir/
	mv $starDir/Solo.out ${sample}.${count}.star.solo
"""
}

// Merge smaller raw STAR outputs for a single sample into a single merged raw STAR output
// Run only when splitFastq is set to true
process merge_star {
input:
    tuple(val(sample), path("Solo.out*"))
output:
    tuple(val(sample), path("${sample}.star.solo/"), emit: merge)
publishDir "$params.outDir/alignment", mode:'copy'
script:
"""
    mergeRawSTARoutput.py --star_dirs Solo.out* --star_feature ${params.starFeature} --star_matrix ${starMatrixFn(params.starMulti)} --sample_name ${sample}
"""
}

// Run cell typing on STARsolo output
process cellTyping {
input:
	tuple(val(sample), val(libName), path("cellReadMetrics"), path("filteredMatrix"))
output:
	path("${sample}")
tag "$sample"
publishDir "$params.outDir/cellTyping", mode: 'copy'
script:
"""
	runScrublet.py --counts filteredMatrix --outDir .

	assignTypesAnalysis.R --projectName ${sample} --countDir filteredMatrix --scrubletPath scrublet_output_table.csv --metricsPath cellReadMetrics --clusterNormalization LogNormalize --StarMtxFormat --seurat --azimuth --writeSeuratMeta --writeSeuratMarkers
	assignTypesPlots.R --projectName ${sample} --plotRtBc --qcPlots
"""
}

// Generate metrics per sample from STARsolo output
process sampleMetricsGeneration {
input:
	tuple(val(sample), val(libName), path("demuxMetrics.json"), path("Solo.out"), val(expectedCells))
	path(samplesCsv)
	val(libStructName)
	val(isBarnyard)
	path(libStructDir)
output:
	tuple(val(sample), path("${sample}_metrics/sample_metrics"), emit: sample_metrics_json_for_report)
	tuple(val(sample), val(libName), path("${sample}_metrics/allCells.csv"), emit: cell_metrics)
	tuple(val(sample), path("${sample}_filtered_star_output"), emit: filtered_star_mtx)
	
publishDir "$params.outDir/samples", pattern:"${sample}_filtered_star_output", saveAs:{"${sample}.filtered"}, mode: 'copy'
publishDir "$params.outDir/samples", pattern:"*/allCells.csv", saveAs:{"${sample}.allCells.csv"}, mode: 'copy'

label 'report'
tag "$sample"
script:
	if (isBarnyard) {
		opts = "--isBarnyard"
	} else {
		opts = ""
	}
	if (params.useSTARthreshold) {
		opts = opts + " --useSTARthreshold"
	}

"""
	getSampleMetrics.py --sample ${sample} --samplesheet ${samplesCsv} --libJsonName ${libStructName} \
	--topCellPercent ${params.topCellPercentage} --libraryStructPath ${libStructDir}\
	--minCellRatio ${params.minCellRatio} --minReads ${params.minReads} --star_out Solo.out \
	--starFeature ${params.starFeature} --starMatrix ${starMatrixFn(params.starMulti)} \
	$opts
"""
}

// Generate report for each sample from metrics generated in sampleMetricsGeneration process
process sampleReportGeneration {
input:
	tuple(val(sample), val(libName), path("${sample}_metrics/allCells.csv"), path("${sample}_metrics/sample_metrics"), path("trim_stats*.json"))
	val(libJson)
	val(isBarnyard)
	path(libStructDir)
output:
	path("reports/${sample}.report.html")
	path("reports/csv/${sample}.reportStatistics.csv"), emit: stats
	path("reports/${sample}_figures"), optional: true
	path("reports/csv/${sample}_unique_transcript_counts_by_RT_Index_well.csv")
	path("reports/csv/${sample}_num_cells_by_RT_Index_well.csv")
publishDir "$params.outDir", mode: 'copy'
errorStrategy 'ignore'
label 'report'
tag "$sample"
script:
	if (isBarnyard) {
                opts = "--isBarnyard "
    } else {
                opts = ""
    }
	if (params.internalReport) {
		opts = opts + "--internalReport "
	}
"""
	export TMPDIR=\$(mktemp -p `pwd` -d)
	generateSampleReport.py --sampleName ${sample} --libJsonName ${libJson} --libraryStructPath ${libStructDir} --sample_metrics ${sample}_metrics --trim_stats trim_stats* $opts
"""
}

// Generate metrics for each library from the metrics of all samples in the library
process fastqMetricsGeneration {
input:
	tuple(val(sample), val(libName), path("allCells*.csv"), path("sample_metrics"))
output:
	tuple(val(libName), path("library_${libName}_metrics"))
errorStrategy 'ignore'
label 'report'
tag "$libName"
script:
"""
	getFastqMetricsFromSampleMetrics.py --sample_metrics allCells* --libName ${libName}
"""
}

// Generate report for each library from metrics generated in fastqMetricsGeneration
process fastqReportGeneration {
input:
	tuple(val(libName), path("metrics"))
	val(libJson)
	path(libStructDir)
output:
	path("reports/library_${libName}.report.html")
	path("reports/csv/library_${libName}*.csv")
publishDir "$params.outDir", mode: 'copy'
label 'report'
tag "$libName"
script:
	if (params.internalReport) {
		opts = "--internalReport"
	} else {
		opts = ""
	}
"""
	export TMPDIR=\$(mktemp -p `pwd` -d)
	generateFastqReport.py --libName ${libName} --libJsonName ${libJson} --libMetrics metrics --libraryStructPath ${libStructDir} $opts
"""
}

// Generate multi sample report from output of sampleReportGeneration
process multiSampleReport {
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

//// Main entry point
// Run the workflow for one or multiple samples
// either from one runFolder or one / multiple sets of fastq files
workflow {
	// Call initialise function in Helper.groovy to pretty print important
	// parameters along with some additional information about the pipeline
	Helper.initialise(workflow, params, log)
	
	regularizeSamplesCsv(file(params.samples))
	samplesCsv = regularizeSamplesCsv.out
	samples = samplesCsv.splitCsv(header:true, strip:true)
	
	libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
	libStructure = loadJson(libJson)
	genome = loadGenome(file(params.genome))

	isBarnyard = genome.get('isBarnyard', false)
	
	// Call the inputReads workflow
	inputReads(samples, samplesCsv, libJson, params.runFolder, params.fastqDir, params.fastqSamplesheet, params.trimFastq, params.fastqc, params.outDir, params.bclConvertParams, params.splitFastq, params.trimAdapt)
	
	// sampleDemuxMetrics -> (sample_name, library_name, demux_json for the library containing the sample)
	sampleDemuxMetrics = inputReads.out.metrics.cross(samples.map{[it.libName,it.sample]}).map({[it[1][1], it[0][0], it[0][1]]})
	sampleDemuxMetrics.dump(tag:'sampleDemuxMetrics')

	// Sort the input fq files from inputReads to not break caching
	// fqsSorted -> (sample_name.well_coordinate, transcript, barcode)
	fqsSorted = inputReads.out.fqs.toSortedList().flatMap()
	fqsSorted.dump(tag:'fqsSorted')
	
	if (params.splitFastq) {
		// Group by sample_name.well_coordinate(first entry of fqSorted)
		// alignFqGroupedBySampleAndWell -> (sample_name.well_coordinate, [transcript], [barcode])
		alignFqGroupedBySampleAndWell = fqsSorted.groupTuple()

		// alignFqJobGroups -> (sample_name, [[transcript]], [[barcode]])
		alignFqJobGroups = alignFqGroupedBySampleAndWell.map {
			// Get sample name from {sample_name.well_coordinate}
			// Not including last element of tokenized because it will be the well coordinate
			def tokenized = it[0].tokenize(".")
			tokenized.removeLast()
			// Get the actual sample name
			def sample = tokenized.join(".")
			def transcript = it[1]
			def barcode = it[2]
			tuple(sample, transcript, barcode)
		}.groupTuple(size:params.starGroupSize, remainder:true)
		
		// Flatten list of list of transcripts and barcodes
		// alignFqGroups -> (sample_name, [transcript], [barcode])
		alignFqJobGroups = alignFqJobGroups.map {
			def sample = it[0]
			def transcript = it[1].flatten()
			def barcode = it[2].flatten()
			tuple(sample, transcript, barcode)
		}
	} else {
		// Group by sample name
		alignFqJobGroups = fqsSorted.groupTuple()
	}

	// Sort transcript and barcode files to ensure input to star is deterministic
	alignFqJobGroups = alignFqJobGroups.map {
		tuple(it[0], it[1].sort(), it[2].sort())
	}
	alignFqJobGroups = alignFqJobGroups.toSortedList().flatMap()	
	
	// Add counter to alignFqJobGroups to aid output naming when splitFastq is true
	def count = 1
	alignFqJobGroups = alignFqJobGroups.map{
		tuple(it[0], count++, it[1], it[2])
	}
	alignFqJobGroups.dump(tag:'alignFqJobGroups')
	
	starsolo(genome.star_index, libStructure, alignFqJobGroups)
	starsolo.out.solo.dump(tag:'solo')
	
	if (params.splitFastq) {
		// Group on sample name to ensure star outputs of the same sample are grouped together
		star_solo_by_sample = starsolo.out.solo.groupTuple()
		star_solo_by_sample.dump(tag: 'star_solo_by_sample')
		
		// Merge star outputs
		merge_star(star_solo_by_sample)
		solo_out = merge_star.out.merge
	} else {
		solo_out = starsolo.out.solo
	}

	expectedCells = samples.map{[it.sample, toIntOr0(it.expectedCells)]}

	if (params.splitFastq) {
		trimFqStatsBySample = inputReads.out.trimFqStats.map {
			// Get sample name from {sample_name.well_coordinate}
			// Not including last element of tokenized because it will be the well coordinate
			def tokenized = it[0].tokenize(".")
			tokenized.removeLast()
			// Get the actual sample name
			def sample = tokenized.join(".")
			def file = it[1]
			tuple(sample, file)
		}.groupTuple()
	}
	else {
		trimFqStatsBySample = inputReads.out.trimFqStats
	}
	
	//trimFqStatsBySample -> (sample_name, cutadapt output files)
	trimFqStatsBySample.dump(tag:'trimFqStatsBySample')
	
	// stats -> (sample_name, library_name, demux_json, star_output, expectedCells)
	stats =  sampleDemuxMetrics.join(solo_out).join(expectedCells)
	stats.dump(tag:'stats')

	sampleMetricsGeneration(stats, samplesCsv, libJson.getName(), isBarnyard, libJson.getParent())
	sampleMetricsBySample = sampleMetricsGeneration.out.cell_metrics.join(sampleMetricsGeneration.out.sample_metrics_json_for_report).join(trimFqStatsBySample)
	sampleReportGeneration(sampleMetricsBySample, libJson.getName(), isBarnyard, libJson.getParent())

	sampleMetricsBySample = sampleMetricsBySample.map{
		def sample = it[0]
		def libName = it[1]
		def allCells = it[2]
		def metrics = it[3]
		tuple(sample, libName, allCells, metrics)
	}
	
	// Group by libName
	metricsByLib = sampleMetricsBySample.groupTuple(by: 1)
	fastqMetricsGeneration(metricsByLib)
	multiSampleReport(sampleReportGeneration.out.stats.collect())
	
	fastqReportGeneration(fastqMetricsGeneration.out, libJson.getName(), libJson.getParent())
	cellTypingInput = sampleMetricsGeneration.out.cell_metrics.join(sampleMetricsGeneration.out.filtered_star_mtx)

	if (params.cellTyping) {
		cellTyping(cellTypingInput)
	}

}

// Function that gets invoked when workflow completes
// Publishes a file named workflow_info.json that contains information about the pipeline execution
workflow.onComplete {
	testing = false
	def data = ["Workflow Information": ["Execution status": "${ workflow.success ? 'OK' : 'failed' }",
								         "Pipeline completion timestamp": "$workflow.complete",
								         "Git repo URL": "$workflow.repository",
								         "Configuration files": "$workflow.configFiles",
         		 						 "Container": "$workflow.container",
		         						 "Command line executed": "$workflow.commandLine",
				         				 "Configuration profile": "$workflow.profile",
						         		 "Start timestamp": "$workflow.start",
         								 "Stop timestamp": "$workflow.complete",
										 "Run name": "$workflow.runName",
										 "Revision": "$workflow.revision",
										 "Commit ID": "$workflow.commitId",
										 "Duration": "$workflow.duration"]]
	def params_data = ["Parameters": [:]]
	for (p in params) {
		if (!p.key.contains('-')) {
			params_data."Parameters".put("$p.key", "$p.value")
		}
		if (p.key.equals('testing_run')) {
			testing = true
			data."Workflow Information".put("Exit status", "$workflow.exitStatus")
			data."Workflow Information".put("Error message", "$workflow.errorMessage")

		}
        
    }
	def reference_data = ["Reference Genome": [:]]
	for (p in genome) {
		reference_data."Reference Genome".put("$p.key", "$p.value")
	}
	def manifest_data = ["Manifest": [:]]
	for (p in workflow.manifest.getProperties()) {
		p = p.toString()
		def split_str = p.split("=")
		if (split_str[0].equals("name") || split_str[0].equals("version")) {
			manifest_data."Manifest".put(split_str[0], split_str[1])
		}
	}

	def json_str = JsonOutput.toJson(data+params_data+reference_data+manifest_data)
	def json_beauty = JsonOutput.prettyPrint(json_str)
	workflow_info = file("$params.outDir/workflow_info.json")
	workflow_info.write(json_beauty)
}
