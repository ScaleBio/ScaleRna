nextflow.enable.dsl=2

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include { inputReads } from './modules/inputReads.nf'
include { alignment } from './modules/alignment.nf'
include { multiSampleReport} from './modules/sampleReport.nf'
include { sampleReport } from './modules/sampleReport.nf' addParams(merge:false)
include { sampleReport as mergeSampleReport } from './modules/sampleReport.nf' addParams(merge:true)

// Load .json file into object
// (reference genomes, library design, etc.)
def loadJson(json) {
	def jsonSlurper  = new JsonSlurper()
	return jsonSlurper.parse( json )
}

// Create a 'file()' from string 'path'
// 'path' could be null, a s3 URL, an absolute file-path or relative to 'baseDir'
def expandPath(path, baseDir) {
	if (path == null) { return null}
	if (path =~ /^s3/) { return file(path)}
	return file(baseDir.resolve(path), checkIfExists: true)
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

def validateParams() {
	if (params.samples == null) {
		ParamLogger.throwError(log, "Must specify --samples (e.g. samples.csv)")
	}
	if (params.libStructure == null) {
		ParamLogger.throwError(log, "Must specify --libStructure (e.g. libV1.1.json)")
	}
	if (params.genome == null) {
		ParamLogger.throwError(log, "Must specify --genome")
	}
	if (params.reporting) {
		if (params.runFolder || params.fastqDir) {
			ParamLogger.throwError(log, "Cannot specify --runFolder or --fastqDir when running reporting-only (--reporting)")
		}
	}
}

// Prepare samples.csv with defaults, rename legacy columns, etc.
process regularizeSamplesCsv {
input: 
	path("samples.in.csv")
output: 
	path("samples.csv")
publishDir "${params.outDir}", mode: 'copy'
label 'small'
cache 'deep'
script:
	opts=""
	if (params.splitFastq && (params.runFolder == null)) {
		opts += "--splitSample "
	}
	if (params.reporting) {
		opts += "--reporting "
	}
	if (params.resultDir) {
		opts += "--resultDir ${params.resultDir} "
	}
"""
	regularizeSamplesCsv.py samples.in.csv $opts > samples.csv
"""
}

// Generate metrics for each library from the metrics of all samples in the library
process fastqMetricsGeneration {
input:
	tuple(val(libName), val(samples), path("allCells*.csv"))
output:
	tuple(val(libName), path("library_${libName}_metrics"))
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
	tuple(val(libName), path(demuxJson), path("metrics"))
	val(libJsonFn)
	path(libStructDir)
output:
	path("reports/library_${libName}.report.html")
	path("reports/csv/library_${libName}*.csv")
publishDir "$params.outDir", mode: 'copy'
label 'report'
label 'optional'
tag "$libName"
script:
	if (params.internalReport) {
		opts = "--internalReport"
	} else {
		opts = ""
	}
	libJson = "${libStructDir}/${libJsonFn}"
"""
	export TMPDIR=\$(mktemp -p `pwd` -d)
	generateFastqReport.py --libName $libName --demuxMetrics $demuxJson --libStruct $libJson --libMetrics metrics $opts
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
label "mergeRawStarOutput"
tag "$group"
script:
	sampleIDs = String.join(" ", samples)
	matrixFn = Utils.starMatrixFn(params.starMulti)
	outDir = "merged/${group}.star.solo/"
"""
	mergeRawSTARoutput.py --star_dirs $starDirs --star_feature $params.starFeature --star_matrix $matrixFn --sample_ids $sampleIDs --out_dir $outDir
"""
}

//// Sub-workflow to run reporting and output steps downstream of STARSolo
// Runs cell-filtering, sample metrics, sample report and cell-typing/clustering
// This does not run the library-level (fastq) report, since that depends on the 
// bcParser metrics and the samples might come from more than 1 library (pipeline / bcParser run).
workflow sampleReporting {
take:
	samples		// each row from samples.csv parsed into a channel
	samplesCsv  // samples.csv file
	libStructure // parsed library structure information
	libJson  // library structure json file
	genome // Reference genome
	soloOut // STARsolo output from previous pipeline run
	trimStats // Read trimming statistics from previous pipeline run
main:
	soloOut.dump(tag:'soloOut')
	trimStats.dump(tag:'trimStats')
	isBarnyard = genome.get('isBarnyard', false)

	// Group cutadapt statistics by sample (all subsamples)
	if (params.splitFastq) {
		trimStatsBySample = trimStats.map {
			// Get sample name from "sampleName_wellCoordinate"
			def sample = it[0].tokenize("_")[0]
			def file = it[1]
			tuple(sample, file)
		}.groupTuple()
	} else {
		trimStatsBySample = trimStats
	}

	// Run reporting on original (un-merged) samples
	sampleReport(samples, samplesCsv, libJson, soloOut, trimStatsBySample, isBarnyard)
	sampleStats = sampleReport.out.sampleStats

	if (params.merge) {
		// Merge samples (alignment results) by group
		sampleGroups = samples.map{ tuple(it.id, it.group) }.join(soloOut).groupTuple(by:1)
		sampleGroups = sampleGroups.filter { it[0].size() > 1} // Don't create merged samples from just one single sample
		sampleGroups.dump(tag:'sampleGroups')
		mergeSamples(sampleGroups)
		mergeTrimStats = sampleGroups.map{ tuple(it[1], []) }
		// Run reporting on merged samples
		mergeSampleReport(samples, samplesCsv, libJson, mergeSamples.out.mergeOut, mergeTrimStats, isBarnyard)
		sampleStats = sampleStats.concat(mergeSampleReport.out.sampleStats)
	}	
	multiSampleReport(sampleStats.collect()) // Combined original and merged samples

emit:
	cellMetrics = sampleReport.out.cellMetrics //allcells.csv outputs per sample
}

//// Full workflow
// Get input reads (fastq or bcl), run barcode parsing, alignment and reporting
workflow endToEnd {
take:
	samples		// each row from samples.csv parsed into a channel
	samplesCsv  // samples.csv file
	libStructure // parsed library structure information
	libJson  // library structure json file
	genome // Reference genome
main:
	// inputReads gets/generates fastq files and runs bcParser to demux per sample
	inputReads(samples, samplesCsv, libJson, params.runFolder, params.fastqDir)
	// STARSolo
	alignment(inputReads.out.fqs, libStructure, genome)
	//// REPORTING
	// Per sample QC report
	sampleReporting(samples, samplesCsv, libStructure, libJson, genome, alignment.out.soloOut, inputReads.out.trimFqStats)

	// Per library QC report
	metricsByLib = sampleReporting.out.cellMetrics.map{
		tuple(it[1], it[0], it[2]) // [libName, sample, allCells]
	}.groupTuple()
	metricsByLib.dump(tag:'metricsByLib')
	fastqMetricsGeneration(metricsByLib)
	fastqReportMetrics = inputReads.out.metrics.join(fastqMetricsGeneration.out)
	fastqReportGeneration(fastqReportMetrics, libJson.getName(), libJson.getParent())
}

//// Main entry point
// Run the workflow for one or multiple samples
// either from reads (--runFolder/ --fastqDir) or pre-existing alignment results (--resultDir)
workflow {
	// Pretty print important parameters along with some additional information about the pipeline
	ParamLogger.initialise(workflow, params, log)
	validateParams()
	// Prepare and load samples.csv
	regularizeSamplesCsv(file(params.samples, checkIfExists:true))
	samplesCsv = regularizeSamplesCsv.out
	samples = samplesCsv.splitCsv(header:true, strip:true)
	samples.dump(tag:'samples')
	// Load library structure json and reference genome
	libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
	libStructure = loadJson(libJson)
	genome = loadGenome(file(params.genome))

	if (params.reporting) {
		// Load STARsolo output from previous pipeline output directory
		soloOut = samples.map{tuple(it.id, file("${it.resultDir}/alignment/${it.id}.star.solo", checkIfExists:true))}
		// We are skipping read trimming statistics when re-running reporting
		// Doesn't matter for anything other than the field in the HTML report
		trimStats = Channel.empty()
		sampleReporting(samples, samplesCsv, libStructure, libJson, genome, soloOut, trimStats)
	} else {
		endToEnd(samples, samplesCsv, libStructure, libJson, genome)
	}
}


//// When the workflow completes, publish 'workflow_info.json' containing information pipeline metadata
workflow.onComplete {
	testing = false
	def workflow_data = ["Nextflow Information": 
		["Execution status": "${ workflow.success ? 'OK' : 'failed' }",
		"Run name": "$workflow.runName",
		"Pipeline completion timestamp": "$workflow.complete",
		"Git repo URL": "$workflow.repository",
		"Configuration files": "$workflow.configFiles",
		"Container": "$workflow.container",
		"Command line executed": "$workflow.commandLine",
		"Configuration profile": "$workflow.profile",
		"Start timestamp": "$workflow.start",
		"Stop timestamp": "$workflow.complete",
		"Exit status": "$workflow.exitStatus",
		"Error message": "$workflow.errorMessage",
		"Revision": "$workflow.revision",
		"Commit ID": "$workflow.commitId",
		"Duration": "$workflow.duration"]
	]
	def params_data = ["Parameters": [:]]
	for (p in params) {
		if (!p.key.contains('-')) {
			params_data."Parameters".put("$p.key", "$p.value")
		}
    }
	def reference_data = ["Reference Genome": [:]]
	for (p in genome) {
		reference_data."Reference Genome".put("$p.key", "$p.value")
	}
	def manifest_data = ["Workflow Manifest": [:]]
	for (p in workflow.manifest.getProperties()) {
		p = p.toString()
		def split_str = p.split("=")
		if (split_str[0].equals("name") || split_str[0].equals("version")) {
			manifest_data["Workflow Manifest"].put(split_str[0], split_str[1])
		}
	}

	def json_str = JsonOutput.toJson(manifest_data + workflow_data +params_data+reference_data)
	def json_beauty = JsonOutput.prettyPrint(json_str)
	workflow_info = file("$params.outDir/workflow_info.json")
	workflow_info.write(json_beauty)
}
