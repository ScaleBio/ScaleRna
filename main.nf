nextflow.enable.dsl=2

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

include { inputReads } from './modules/inputReads.nf'
include { alignment } from './modules/alignment.nf'
include { multiSampleReport} from './modules/sampleReport.nf'
include { sampleReport } from './modules/sampleReport.nf' addParams(mergedSamples:false)
include { sampleReport as mergeSampleReport } from './modules/sampleReport.nf' addParams(mergedSamples:true)
include { filterMtx } from './modules/createMtx.nf' addParams(mergedSamples:false)
include { filterMtx as mergedFilterMtx } from './modules/createMtx.nf' addParams(mergedSamples:true)
include { downstream } from './modules/downstream.nf' addParams(comparison:false)
include { downstream as mergedDownstream } from './modules/downstream.nf' addParams(comparison:false)
include { downstream as compareDownstream } from './modules/downstream.nf' addParams(comparison:true)

// Load .json file into object
// (reference genomes, library design, etc.)
def loadJson(json) {
	def jsonSlurper  = new JsonSlurper()
	return jsonSlurper.parse( json )
}

// Create a 'file()' from string 'path'
// 'path' could be null, a URL, an absolute file-path or relative to 'baseDir'
def expandPath(path, baseDir) {
	if (path == null) { return null}
	uri = new URI(path)
	if (uri.isAbsolute()) {
		return file(path)
	} else {
		return baseDir.resolve(path)
	}
}

// Reference genome files and parameters
// Paths can be absolute or relative to location of json
def loadGenome(json) {
	if (json.isDirectory()) {
		ParamLogger.throwError("--genome should be json file, not a directory: ${json}")
	}
	def baseDir = json.getParent()
	genome = loadJson(json)
	genome.star_index = expandPath(genome.star_index, baseDir)
	if (genome.star_index == null || !genome.star_index.exists()) {
		ParamLogger.throwError("Genome star_index not found: ${genome.star_index?.toUriString() ?: ''}")
	}
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
	val(libStructName)
	path(libStructDir)
output: 
	path("samples.csv")
publishDir params.outDir, mode: 'copy'
label 'small'
cache 'deep'
script:
	opts=""
	if (params.splitFastq) {
		opts += "--splitSample "
	}
	if (params.reporting) {
		opts += "--reporting "
	}
	if (params.resultDir) {
		opts += "--resultDir ${params.resultDir} "
	}
	libStruct = "${libStructDir}/${libStructName}"
"""
	regularizeSamplesCsv.py samples.in.csv --libraryStruct $libStruct $opts > samples.csv
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
	path("${outDir}/library_${libName}.report.html")
	path("${outDir}/csv/library_${libName}*.csv")
publishDir params.outDir, mode: 'copy'
label 'report'
label 'optional'
tag "$libName"
script:
	outDir = "reports/library"
	if (params.internalReport) {
		opts = "--internalReport"
	} else {
		opts = ""
	}
	libJson = "${libStructDir}/${libJsonFn}"
"""
	export TMPDIR=\$(mktemp -p `pwd` -d)
	generateFastqReport.py --libName $libName --outDir $outDir --demuxMetrics $demuxJson --libStruct $libJson --libMetrics metrics $opts
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
	libJson  // library structure json file
	genome // Reference genome
	soloOut // STARsolo output from previous pipeline run
	trimStats // Read trimming statistics from previous pipeline run
main:
	soloOut.dump(tag:'soloOut')
	trimStats.dump(tag:'trimStats')
	isBarnyard = genome.get('isBarnyard', false)

	// Group cutadapt statistics by sample
	// All split-fastq subsamples of the form "sampleName_wellCoordinate" are grouped together
	// This is a no-op if splitFastq false
	trimStats.map {
		def sample = it[0].tokenize("_")[0]
		def file = it[1]
		tuple(sample, file)
	}.groupTuple()
	.map { tuple(it[0], it[1].flatten())} // Flatten list of empty lists in case of --reporting and --trimStats
	.dump(tag:'trimStatsBySample')
	.set{ trimStatsBySample }

	filterMtx(samples, libJson, soloOut, isBarnyard)

	if(params.seurat || params.azimuth || params.annData){
		downstream(samples, filterMtx.out.allCells, filterMtx.out.filteredMtx, filterMtx.out.rawMtx, filterMtx.out.sampleStats)
		if(params.compSamples){
			compareDownstream(samples, filterMtx.out.allCells, filterMtx.out.filteredMtx, filterMtx.out.rawMtx, filterMtx.out.sampleStats)
		}
	}
	// Run reporting on original (un-merged) samples
	sampleReport(samples, libJson, trimStatsBySample, isBarnyard, filterMtx.out.allCells, filterMtx.out.libCount, filterMtx.out.sampleStats)
	sampleStats = sampleReport.out.sampleStats

	if (params.merge) {
		// Collect all trimStats for one sample (group) into a single list.
		// The reporting step will sum them up.
		samples.map{ tuple(it.id, it.group) }
		.join(trimStatsBySample)
		.groupTuple(by:1) 
		.map { tuple(it[1], it[2].flatten())} // Flatten here for splitFastq where each sub-sample already is a list
		.dump(tag:'mergeTrimStats')
		.set { mergeTrimStats }
		// This returns merged samples only for samples with multiple libraries
		// (Same sample name, different libName)
		mergedFilterMtx(samples, libJson, soloOut, isBarnyard)

		if(params.seurat || params.azimuth || params.annData){
			mergedDownstream(samples, mergedFilterMtx.out.allCells, mergedFilterMtx.out.filteredMtx, mergedFilterMtx.out.rawMtx, mergedFilterMtx.out.sampleStats)
		}
		// Run reporting on merged samples (if any)
		mergeSampleReport(samples, libJson, mergeTrimStats, isBarnyard, mergedFilterMtx.out.allCells, mergedFilterMtx.out.libCount, mergedFilterMtx.out.sampleStats)
		sampleStats = sampleStats.concat(mergeSampleReport.out.sampleStats)
	}	
	multiSampleReport(sampleStats.collect()) // Combined original and merged samples

emit:
	cellMetrics = filterMtx.out.allCells //allcells.csv outputs per sample
}

//// Full workflow
// Get input reads (fastq or bcl), run barcode parsing, alignment and reporting
workflow endToEnd {
take:
	samples		// each row from samples.csv parsed into a channel
	samplesCsv  // samples.csv file
	libJson  // library structure json file
	genome // Reference genome
main:
	// inputReads gets/generates fastq files and runs bcParser to demux per sample
	inputReads(samples, samplesCsv, libJson, params.runFolder, params.fastqDir)
	libJsonContents = loadJson(libJson)
	// STARSolo
	alignment(inputReads.out.fqs, libJsonContents, genome)
	//// REPORTING
	// Per sample QC report
	sampleReporting(samples, samplesCsv, libJson, genome, alignment.out.soloOut, inputReads.out.trimFqStats)

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
	// Validate and print key pipeline parameters
	ParamLogger.initialise(workflow, params, log)
	// Compute library structure json path to enable staging in of json file to downstream processes
	libJson = expandPath(params.libStructure, file(projectDir) / "references")
	// Fill in samples csv with defaults and workflow specific parameters
	regularizeSamplesCsv(file(params.samples, checkIfExists:true), libJson.getName(), libJson.getParent())
	samplesCsv = regularizeSamplesCsv.out
	samples = samplesCsv.splitCsv(header:true, strip:true)
	samples.dump(tag:'samples')
	genome = loadGenome(file(params.genome, checkIfExists:true))

	if (params.reporting) {
		// Load STARsolo output from previous pipeline output directory
		soloOut = samples.map{tuple(it.id, file("${it.resultDir}/alignment/${it.id}/${it.id}.star.solo", checkIfExists:true))}
		// We are skipping read trimming statistics when re-running reporting
		// Doesn't matter for anything other than the field in the HTML report
		trimStats = samples.map{ tuple(it.id, []) }
		sampleReporting(samples, samplesCsv, libJson, genome, soloOut, trimStats)
	} else {
		endToEnd(samples, samplesCsv, libJson, genome)
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
	workflow_info = file(params.outDir) / "workflow_info.json"
	workflow_info.write(json_beauty)
}
