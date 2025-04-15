nextflow.enable.dsl=2

// Validate and print key pipeline parameters
ParamLogger.initialise(workflow, params, log)
include { INPUT_READS } from './modules/input_reads.nf'
include { SCALE_PLEX } from './modules/scale_plex.nf'
include { SCALE_PLEX as SCALE_PLEX_REPORTING} from './modules/scale_plex.nf'
include { ALIGNMENT } from './modules/alignment.nf'
include { MultiSampleReport } from './modules/sample_report.nf'
include { MultiLibraryReport } from './modules/internal.nf'
include { ULTIMA } from './modules/ultima.nf'
include { SAMPLE_REPORT } from './modules/sample_report.nf'
include { SAMPLE_REPORT as MERGE_SAMPLE_REPORT } from './modules/sample_report.nf'
include { FILTER_MTX } from './modules/create_mtx.nf'
include { FILTER_MTX as MERGED_FILTER_MTX } from './modules/create_mtx.nf'
include { DOWNSTREAM } from './modules/downstream.nf'
include { DOWNSTREAM as MERGED_DOWNSTREAM } from './modules/downstream.nf'
include { DOWNSTREAM as COMPARE_DOWNSTREAM } from './modules/downstream.nf'

def shouldRunLibDetection(samples) {
	def header = samples.withReader { reader ->
		reader.readLine()
	}
	if (!header.contains("libName") && params.quantum && params.fastqDir) {
		return true
	}
	return false
}

// Create a 'file()' from string 'path'
// 'path' could be null, a URL, an absolute file-path or relative to 'baseDir'
def expandPath(path, baseDir) {
	if (path == null) { return null}
	def uri = new URI(path)
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
	def genome = Utils.loadJson(json)
	genome.star_index = expandPath(genome.star_index, baseDir)
	if (genome.star_index == null || !genome.star_index.exists()) {
		ParamLogger.throwError("Genome star_index not found: ${genome.star_index?.toUriString() ?: ''}")
	}
	genome.gtf = expandPath(genome.gtf, baseDir)
	return genome
}

def resolveScalePlexToRnaMapping(scalePlexLibJson) {
	if (params.scalePlexToRnaMapping) {
		scalePlexToRnaMapping = file(params.scalePlexToRnaMapping, checkIfExists:true)
	} else {
		scalePlexLibStructure = Utils.loadJson(scalePlexLibJson)
		if (scalePlexLibStructure.get("scaleplex_to_rna_mapping")) {
			scalePlexToRnaMapping = file(scalePlexLibJson.getParent()) / scalePlexLibStructure["scaleplex_to_rna_mapping"]
		} else {
			// path input cannot be null, no mapping file used in 3lvl pipeline
			scalePlexToRnaMapping = []
		}
	}
	return scalePlexToRnaMapping
}

// Prepare samples.csv with defaults, rename legacy columns, etc.
process RegularizeSamplesCsv {
	label 'small'
	cache 'deep'

	// Publish user provided samples csv to output directory
	publishDir params.outputDir, pattern: "samples.in.csv", mode: 'copy', saveAs: {'samples.csv'}

	input: 
	path("samples.in.csv")
	// Directory containing the library structure definition (barcode sequence lists are loaded from here)
	path("references")
	// Filename of the library structure definition .json
	val(libStructName) 
	// File containing the mapping of scale plex to RNA PCR barcodes
    path(mapping_file) 
	
	output: 
	path("samples.csv"), emit: samples_csv
	path("barcode_range_unwrapped_samples.csv"), emit: barcode_unwrapped_samples_csv, optional: true
	path("samples.in.csv")

	script:
	libStruct = "references/$libStructName"
	opts=""
	if (params.scalePlex) {
		opts += "--scalePlex "
		if (mapping_file) {
			opts += "--scalePlexToRnaMapping $mapping_file "
		}
	}
	if (params.splitFastq) {
		opts += "--splitSample "
	}
	if (params.reporting) {
		opts += "--reporting "
	}
	if (params.resultDir) {
		opts += "--resultDir ${params.resultDir} "
	}
	if (params.quantum) {
		opts += "--quantum "
	}
	if (params.ultimaCramDir) {
		opts += "--ultima "
	}
	if (params.fastqDir || params.ultimaCramDir) {
		opts += "--fastq "
	}
	"""
	regularize_samples_csv.py samples.in.csv --libraryStruct $libStruct $opts
	"""
}

// Generate metrics for each library from the metrics of all samples in the library
process LibraryMetricsGeneration {
	tag "$libName"
	label 'large'
	
	input:
	tuple(val(libName), val(samples), path("allCells*.parquet"))
	
	output:
	tuple(val(libName), path("library_${libName}_metrics"))

	script:
	"""
	get_library_metrics_from_sample_metrics.py --sample_metrics allCells* --libName ${libName}
	"""
}

// Generate report for each library from metrics generated in LibraryMetricsGeneration
process LibraryReportGeneration {
	tag "$libName"
	label 'report'
	label 'optional'

	publishDir params.outputDir, mode: 'copy'

	input:
	tuple(val(libName), path(demuxJson), path("metrics"))
	val(libJsonFn)
	path(libStructDir)
	val(scalePlexLibJsonFn)
	path("scalePlexLibStructDir") // Cannot use a variable here since it will cause input name collision if libStructDir and scalePlexLibStructDir are the same
	val(hashLibNames)
	
	output:
	path("${outDir}/library_${libName}.report.html")
	path("${outDir}/csv/library_${libName}.overallMatches.csv"), emit: overallMatches, optional: true
	path("${outDir}/csv/library_${libName}.typeLevelMatches.csv"), emit: typeLevelMatches, optional: true
	path("${outDir}/csv/library_${libName}_unique_reads*.csv")

	script:
	outDir = "reports/library"
	opts= "--minPassingSampleReads ${params.minPassingSampleReads} "
	if (params.internalReport) {
		opts = opts + "--internalReport "
	}
	if (params.ultimaCramDir) {
		opts = opts + "--ultima "
	} else {
		opts = opts + "--demuxMetrics $demuxJson "
	} 
	if (params.scalePlex && libName in hashLibNames) {
		libJson = "scalePlexLibStructDir/${scalePlexLibJsonFn}"
		opts += " --scalePlexLib"
	} else {
		libJson = "${libStructDir}/${libJsonFn}"
	}
	"""
	export DATAPANE_CDN_BASE="https://d3j2ibc7el7s1k.cloudfront.net/v0.17.0"
	export TMPDIR=\$(mktemp -p `pwd` -d)
	generate_library_report.py --libName $libName --outDir $outDir --libStruct $libJson --libMetrics metrics $opts
"""
}

//// Sub-workflow to run reporting and output steps downstream of STARSolo
// Runs cell-filtering, sample metrics, sample report and cell-typing/clustering
// This does not run the library-level (fastq) report, since that depends on the 
// bcParser metrics and the samples might come from more than 1 library (pipeline / bcParser run).
workflow SAMPLE_REPORTING {
take:
	samples	     // each row from samples.csv parsed into a channel
	libJson      // library structure json file
	isBarnyard   // Bool indicating if the genome is barnyard
	filterMtx    // Cell metrics for individual samples ('id') and ('merged') samples
	hashResults  // Hash metrics if --scalePlex
main:
	if(params.seurat || params.azimuth || params.annData){
		DOWNSTREAM(samples, filterMtx.id.allCells, filterMtx.id.filteredMtx, filterMtx.id.sampleStats, libJson, false)
		if(params.compSamples){
			COMPARE_DOWNSTREAM(samples, filterMtx.id.allCells, filterMtx.id.filteredMtx, filterMtx.id.sampleStats, libJson, true)
		}
	}
	// Run reporting on original (un-merged) samples
	SAMPLE_REPORT(samples, libJson, isBarnyard, filterMtx.id.allCells, filterMtx.id.allBarcodes, filterMtx.id.libCount, filterMtx.id.sampleStats, filterMtx.id.soloLog, hashResults.id, false)
	sampleStats = SAMPLE_REPORT.out.sampleStats

	if (params.merge) {
		if(params.seurat || params.azimuth || params.annData){
			MERGED_DOWNSTREAM(samples, filterMtx.merged.allCells, filterMtx.merged.filteredMtx, filterMtx.merged.sampleStats, libJson, false)
		}
		// Run reporting on merged samples
		MERGE_SAMPLE_REPORT(samples, libJson, isBarnyard, filterMtx.merged.allCells, filterMtx.merged.allBarcodes, filterMtx.merged.libCount, filterMtx.merged.sampleStats, filterMtx.merged.soloLog, hashResults.merged, true)
		sampleStats = sampleStats.concat(MERGE_SAMPLE_REPORT.out.sampleStats)
	}	
	MultiSampleReport(sampleStats.collect()) // Combined original and merged samples
}

//// Full workflow
// Get input reads (fastq or bcl), run barcode parsing, alignment and reporting
workflow END_TO_END {
take:
	samples		 // each row from samples.csv parsed into a channel
	samplesCsv   // samples.csv file
	libStructure // parsed library structure information
	libJson      // library structure json file
	scalePlexLibJson // ScalePlex library structure json file
	scalePlexToRnaMapping // File with mapping of ScalePlex to RNA PCR barcodes
	genome       // Reference genome
	barcodeUnwrappedSamplesCsv // samples csv with barcode ranges unwrapped
	runLibDetection // Run library detection
main:
	if (params.ultimaCramDir) {
		ULTIMA(samples, barcodeUnwrappedSamplesCsv)
		alignmentInput = ULTIMA.out.alignmentInput
		ubamHash = ULTIMA.out.ubamHash
		// Construct totalSampleReadsBySampleID channel and set it to 0 for all samples
		samples
			.filter { it.rnaId == null } // Only RNA samples are passed on to FILTER_MTX
			.map {
				tuple(it.id, 0)
			}
			.set { totalSampleReadsBySampleID }
	}
	else {
		// INPUT_READS gets/generates fastq files and runs bcParser to demux per sample
		INPUT_READS(samples, samplesCsv, libStructure, libJson, scalePlexLibJson, params.runFolder, params.fastqDir, runLibDetection)
		alignmentInput = INPUT_READS.out.ubam
		ubamHash = INPUT_READS.out.ubamHash
		totalSampleReadsBySampleID = INPUT_READS.out.totalSampleReadsBySampleID
		samples = INPUT_READS.out.samples
	}
	// STARSolo
	ALIGNMENT(alignmentInput, genome)
	// Cell calling and metrics generation
	isBarnyard = genome.get('isBarnyard', false)
	// Cell calling and metrics generation
	FILTER_MTX(samples, libJson, ALIGNMENT.out.soloOut, isBarnyard, ALIGNMENT.out.soloLog, totalSampleReadsBySampleID, false)
	if (params.merge) {
		// This returns merged samples only for samples with multiple libraries
		// (Same sample name, different libName)
		MERGED_FILTER_MTX(samples, libJson, ALIGNMENT.out.soloOut, isBarnyard, ALIGNMENT.out.soloLog, totalSampleReadsBySampleID, true)
	}
	//// ScalePlex
	if (params.scalePlex) {
		SCALE_PLEX(scalePlexLibJson, scalePlexToRnaMapping, samples, ubamHash,
			[id: FILTER_MTX.out.allCells, merged: params.merge ? MERGED_FILTER_MTX.out.allCells : null],
			[id: FILTER_MTX.out.libCount, merged: params.merge ? MERGED_FILTER_MTX.out.libCount : null],
			[]) // perSampleCountRaw is empty if not a --reporting run
		SCALE_PLEX
			.out.metrics
			// [meta, cellMetrics, filteredMtx, scaleplex_stats.csv, metrics.csv]
			.map { [it[0].rnaId, it[1], it[3], it[4], it[0].lib] }
			// [rnaId, cellMetrics, scaleplex_stats.csv, metrics.csv, hashLibName]
			.set { hashResults }
			if (params.merge) {
				SCALE_PLEX
					.out.mergedMetrics
					// [meta, cellMetrics, filteredMtx, scaleplex_stats.csv, metrics.csv]
					.map { [it[0].rnaId, it[1], it[3], it[4], it[0].lib] }
					// [rnaId, cellMetrics, scaleplex_stats.csv, metrics.csv, hashLibName]
					.set { mergedHashResults }
			}
	} else {
		hashResults = []
		mergedHashResults = []
	}
	//// REPORTING
	// Per sample QC report
	SAMPLE_REPORTING(samples, libJson, isBarnyard,
		[id: FILTER_MTX.out, merged: params.merge ? MERGED_FILTER_MTX.out : null],
		[id: hashResults, merged: params.merge ? mergedHashResults : null])

	// Per library QC report
	FILTER_MTX
		.out.allBarcodes
		.map {
			tuple(it[1], it[0], it[2]) // [libName, sample, allCells]
		}.groupTuple()
		.set { metricsByLib }

	if (params.scalePlex) {
		metricsByLib.concat(
			SCALE_PLEX
				.out.metrics
				.map { 
					// [meta, cellMetrics, filteredMtx, scaleplex_stats.csv, metrics.csv]
					tuple(it[0].lib, it[0].sample, it[1])
				}.groupTuple()
		).set { metricsByLib }
	}

	LibraryMetricsGeneration(metricsByLib)

	if (params.ultimaCramDir) {
		// Empty list indicates bcParser demux metrics, which do not exist when starting from cram files
		LibraryMetricsGeneration.out
			.map { tuple(it[0], [], it[1]) }
			.set { libraryReportMetrics }
	}
	else {
		INPUT_READS
			.out.metrics.join(LibraryMetricsGeneration.out)
			.set { libraryReportMetrics }
	}
	libraryReportMetrics.dump(tag:'libraryReportMetrics')
	// Get hash lib names
	samples
		.filter { it.rnaId != null }
		.collect { it.libName }
		.ifEmpty { [] }
		.dump(tag:'hashLibNames')
		.set { hashLibNames }
	LibraryReportGeneration(libraryReportMetrics, libJson.getName(), libJson.getParent(), scalePlexLibJson.getName(), scalePlexLibJson.getParent(), hashLibNames)
	
	// Barcode stats do not exist when starting from ultima cram files 
	if (params.internalReport && !params.ultimaCramDir) {
		MultiLibraryReport(LibraryReportGeneration.out.typeLevelMatches.collect(), LibraryReportGeneration.out.overallMatches.collect())
	}
}

//// Main entry point
// Run the workflow for one or multiple samples
// either from reads (--runFolder/ --fastqDir) or pre-existing alignment results (--resultDir)
workflow {
	runLibDetection = shouldRunLibDetection(file(params.samples, checkIfExists:true))
	// Load library structure json and reference genome
	libJson = expandPath(params.libStructure, file(projectDir) / "references")
	libStructure = Utils.loadJson(libJson)
	// Load scalePlex library structure json
	scalePlexLibJson = expandPath(params.scalePlexLibStructure, file(projectDir) / "references")
	scalePlexToRnaMapping = resolveScalePlexToRnaMapping(scalePlexLibJson)
	// Prepare and load samples.csv
	RegularizeSamplesCsv(file(params.samples, checkIfExists:true), libJson.getParent(), libJson.getName(), scalePlexToRnaMapping)
	samplesCsv = RegularizeSamplesCsv.out.samples_csv
	samples = samplesCsv.splitCsv(header:true, strip:true)
	samples.dump(tag:'samples')

	if (params.ultimaCramDir) {
		RegularizeSamplesCsv.out.barcode_unwrapped_samples_csv.splitCsv(header:true, strip:true)
			.map { [it.id, it.libName, it.barcode] }
			.dump(tag:'barcodeUnwrappedSamplesCsv')
			.set {barcodeUnwrappedSamplesCsv}
	} else {
		barcodeUnwrappedSamplesCsv = []
	}
	
	genome = loadGenome(file(params.genome, checkIfExists:true))

	if (params.reporting) {
		samples
			.map { it.sample }
			.unique()
			.collect()
			.dump(tag:'sampleNames')
			.set { sampleNames }
		
		samples
			.filter { it.resultDir == null }
			.map { 
				ParamLogger.throwError("Sample ${it.id} does not have a resultDir")
			}

		// Filter out sample-library combinations not present in the resultDir
		samples
			.filter { it.rnaId == null &&  !(file(it.resultDir) / "alignment" / "${it.id}" / "${it.id}.star.solo").exists() }
			.map { it.id }
			.collect()
			.dump(tag:'missingLibs')
			.set { missingLibs }
			
		samples
			.filter {
				// Keep sample-library combinations present in the resultDir
				def id = it.rnaId ?: it.id
				!(id in missingLibs.val)
			}
			.dump(tag:'filteredSamples')
			.set { samples }
			
		// Check that no samples were completely filtered out
		samples
			.map { it.sample }
			.unique()
			.collect()
			.ifEmpty { [] }
			.dump(tag:'filteredSampleNames')
			.set { filteredSampleNames }
			
		sampleNames
			.flatMap()
			.filter { !(it in filteredSampleNames.val) }
			.map {
				ParamLogger.throwError("Sample ${it} was not found in the resultDir")
			}

		// Load STARsolo output from previous pipeline output directory
		samples
			.filter { it.rnaId == null } // Only RNA samples have alignment results
			.map {
				def alignmentDir = file(it.resultDir) / "alignment"
				tuple(it.id, file(alignmentDir / "${it.id}" / "${it.id}.star.solo", checkIfExists:true))
			}
			.set { soloOut }
		// Load STARsolo log file from previous pipeline output directory
		samples
			.filter { it.rnaId == null } // Only RNA samples have alignment results
			.map {
				def alignmentDir = file(it.resultDir) / "alignment"
				tuple(it.id, file(alignmentDir / "${it.id}" / "${it.id}.star.align" / "Log.final.out", checkIfExists:true))
			}
			.set { soloLog }
		// Get totalSampleReads per sample from allSamples.reportStatistics.csv
		samples
			.filter { it.rnaId == null }
			.map { file(file(it.resultDir) / "reports" / "allSamples.reportStatistics.csv", checkIfExists:true) }
			.unique() // Multiple samples can have same resultDir
			.map { reportStatistics ->
				def lines = reportStatistics.text.readLines()
				def header = lines[0].split(',')
				// Get row with "Total Sample Reads"
				def totalSampleReadsRow = lines.find { it.split(',')[1] == 'Total Sample Reads' }.split(',')
				// Ignore column names with merged, Category and Metric in the header
				// Remaining column names are sample names (sampleName.libName)
				def totalSampleReads = header.findAll { !it.contains('merged') && it != 'Category' && it != 'Metric' }
					.collect { sampleName -> 
						[sampleName, totalSampleReadsRow[header.toList().indexOf(sampleName)]]
					}
				totalSampleReads
			}
			// flatMap because totalSampleReads is a list of lists
			.flatMap { totalSampleReads -> 
				totalSampleReads.collect { sampleReadPair -> 
					tuple(sampleReadPair[0], sampleReadPair[1])
				}
			}
			.dump(tag:'totalSampleReadsBySampleID')
			.set { totalSampleReadsBySampleID }
		// Cell calling and metrics generation
		isBarnyard = genome.get('isBarnyard', false)
		FILTER_MTX(samples, libJson, soloOut, isBarnyard, soloLog, totalSampleReadsBySampleID, false)
		if (params.merge) {
			// This returns merged samples only for samples with multiple libraries
			// (Same sample name, different libName)
			MERGED_FILTER_MTX(samples, libJson, soloOut, isBarnyard, soloLog, totalSampleReadsBySampleID, true)
		}
		//// ScalePlex
		if (params.scalePlex) {
			samples
				.filter { it.rnaId != null } // Hash samples
				.map {
					def scalePlexDir = file(it.resultDir) / "scaleplex"
					[[lib: it.libName, sample: it.sample], file(scalePlexDir / "${it.id}.cellMetrics.parquet", checkIfExists: true), file(scalePlexDir / "${it.id}.raw.matrix/matrix.mtx.gz", checkIfExists: true)]
				}
				.set { perSampleCountRaw}
			SCALE_PLEX_REPORTING(scalePlexLibJson, scalePlexToRnaMapping, samples, [], // ubam is not used for reporting
				[id: FILTER_MTX.out.allCells, merged: params.merge ? MERGED_FILTER_MTX.out.allCells : null],
				[id: FILTER_MTX.out.libCount, merged: params.merge ? MERGED_FILTER_MTX.out.libCount : null],
				perSampleCountRaw)
			SCALE_PLEX_REPORTING
				.out.metrics
				// [meta, cellMetrics, filteredMtx, scaleplex_stats.csv, metrics.csv]
				.map { [it[0].rnaId, it[1], it[3], it[4], it[0].lib] }
				// [rnaId, cellMetrics, scaleplex_stats.csv, metrics.csv, hashLibName]
				.set { hashResults }
			if (params.merge) {
				SCALE_PLEX_REPORTING
					.out.mergedMetrics
					// [meta, cellMetrics, filteredMtx, scaleplex_stats.csv, metrics.csv]
					.map { [it[0].rnaId, it[1], it[3], it[4], it[0].lib] }
					// [rnaId, cellMetrics, scaleplex_stats.csv, metrics.csv, hashLibName]
					.set { mergedHashResults }
			}
		} else {
			hashResults = []
			mergedHashResults = []
		}
		SAMPLE_REPORTING(samples, libJson, isBarnyard,
			[id: FILTER_MTX.out, merged: params.merge ? MERGED_FILTER_MTX.out : null],
			[id: hashResults, merged: params.merge ? mergedHashResults : null])
	} else {
		END_TO_END(samples, samplesCsv, libStructure, libJson, scalePlexLibJson, scalePlexToRnaMapping, genome, barcodeUnwrappedSamplesCsv, runLibDetection)
	}
}


//// When the workflow completes, publish 'workflow_info.json' containing information on pipeline metadata
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

	def json_str = groovy.json.JsonOutput.toJson(manifest_data + workflow_data +params_data+reference_data)
	def json_beauty = groovy.json.JsonOutput.prettyPrint(json_str)
	def workflow_info = file(params.outputDir) / "workflow_info.json"
	workflow_info.write(json_beauty)
}

