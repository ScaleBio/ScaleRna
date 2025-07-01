nextflow.enable.dsl=2

// Validate and print key pipeline parameters
ParamLogger.initialise(workflow, params, log)
include { INPUT_READS } from './modules/input_reads.nf'
include { SCALE_PLEX } from './modules/scale_plex.nf'
include { SCALE_PLEX as SCALE_PLEX_REPORTING} from './modules/scale_plex.nf'
include { ALIGNMENT } from './modules/alignment.nf'
include { MultiSampleReport } from './modules/sample_reporting.nf'
include { MultiLibraryReport } from './modules/internal.nf'
include { ULTIMA } from './modules/ultima.nf'
include { SAMPLE_REPORTING } from './modules/sample_reporting.nf'
include { CELL_CALLING } from './modules/create_mtx.nf'
include { DOWNSTREAM } from './modules/downstream.nf'
include { DOWNSTREAM as COMPARE_DOWNSTREAM } from './modules/downstream.nf'
include { MERGING } from './modules/merging.nf'
include { getTotalSampleReads } from './modules/utils.nf'
include { shouldRunLibDetection } from './modules/utils.nf'
include { expandPath } from './modules/utils.nf'
include { loadGenome } from './modules/utils.nf'
include { resolveScalePlexToRnaMapping } from './modules/utils.nf'

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
	// Whether data was generated using QuantumScale assay
	val(quantum)
	
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
	if (quantum) {
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
	val(hashLibNames)
	val(quantum)
	
	output:
	tuple(val(libName), path("library_${libName}_metrics"))

	script:
	opts = ""
	if (quantum && params.internalReport && !(params.scalePlex && libName in hashLibNames)) {
		opts += " --beadMetrics "
	}
	"""
	get_library_metrics_from_sample_metrics.py --sample_metrics allCells* --libName ${libName} $opts
	"""
}

// Generate report for each library from metrics generated in LibraryMetricsGeneration
process LibraryReportGeneration {
	tag "$libName"
	label 'report'
	label 'optional'

	publishDir params.outputDir, mode: 'copy'

	input:
	tuple(val(libName), path(demuxJson), path("metrics"), path(beadScores))
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
	path("${outDir}/figures_internal"), optional: true

	script:
	outDir = "reports/library"
	opts= "--minPassingSampleReads ${params.minPassingSampleReads} "
	if (params.internalReport) {
		opts = opts + "--internalReport --minDivergence ${params.minBeadDivergence} "
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

//// Main entry point
// Run the workflow for one or multiple samples
// either from reads (--runFolder/ --fastqDir) or pre-existing alignment results (--resultDir)
workflow {
	// Load rna library structure json
	libJson = expandPath(params.libStructure, file(projectDir) / "references")
	libJsonContents = Utils.loadJson(libJson)
	// Load scalePlex library structure json
	scalePlexLibJson = expandPath(params.scalePlexLibStructure, file(projectDir) / "references")
	scalePlexLibJsonContents = Utils.loadJson(scalePlexLibJson)
	scalePlexToRnaMapping = resolveScalePlexToRnaMapping(scalePlexLibJson)
	libraryInfo = [
		rnaLibraryStructureFile: libJson,
		rnaTrimAdapter: libJsonContents["trimAdapter"],
		rnaGenomicR1: libJsonContents.get('genomic_r1', false),
		rnaGenomicR2: libJsonContents.get('genomic_r2', false),
		scalePlexLibraryStructureFile: scalePlexLibJson,
		scalePlexTrimAdapter: scalePlexLibJsonContents["trimAdapter"],
		scalePlexToRnaMappingFile: scalePlexToRnaMapping,
		quantum: libJsonContents.get('quantum', false),
		sampleBarcode: libJsonContents["sample_barcode"],
	]
	// Prepare and load samples.csv
	RegularizeSamplesCsv(
		file(params.samples, checkIfExists:true),
		libraryInfo.rnaLibraryStructureFile.getParent(),
		libraryInfo.rnaLibraryStructureFile.getName(),
		libraryInfo.scalePlexToRnaMappingFile,
		libraryInfo.quantum
	)
	samplesCsv = RegularizeSamplesCsv.out.samples_csv
	samples = samplesCsv.splitCsv(header:true, strip:true)
	samples.dump(tag:'samples')
	runLibDetection = shouldRunLibDetection(file(params.samples, checkIfExists:true), libraryInfo.quantum)
	if (params.ultimaCramDir) {
		RegularizeSamplesCsv.out.barcode_unwrapped_samples_csv.splitCsv(header:true, strip:true)
			.dump(tag:'barcodeUnwrappedSamplesCsv')
			.set {barcodeUnwrappedSamplesCsv}
	} else {
		barcodeUnwrappedSamplesCsv = []
	}
	
	genome = loadGenome(file(params.genome, checkIfExists:true))
	isBarnyard = genome.get('isBarnyard', false)
	
	if (params.reporting) {
		// For reporting run, validate that all samples are present in the resultDir
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
		samples
			.filter { it.rnaId == null }
			.map { sample -> getTotalSampleReads(sample).flatten() }
			.dump(tag:'totalSampleReadsBySampleID')
			.set { totalSampleReadsBySampleID }
		// Load STARsolo output from previous pipeline output directory
		samples
			.filter { it.rnaId == null } // Only RNA samples have alignment results
			.map {
				def alignmentDir = file(it.resultDir) / "alignment"
				tuple(it.id, file(alignmentDir / "${it.id}" / "${it.id}.star.solo", checkIfExists:true))
			}
			.dump(tag:'soloOut')
			.set { soloOut }
		// Load STARsolo log file from previous pipeline output directory
		samples
			.filter { it.rnaId == null } // Only RNA samples have alignment results
			.map {
				def alignmentDir = file(it.resultDir) / "alignment"
				tuple(it.id, file(alignmentDir / "${it.id}" / "${it.id}.star.align" / "Log.final.out", checkIfExists:true))
			}
			.dump(tag:'soloLog')
			.set { soloLog }
		if (params.scalePlex) {
			samples
				.filter { it.rnaId != null } // Hash samples
				.filter {
					def scalePlexDir = file(it.resultDir) / "scaleplex"
					def cellMetricsFile = scalePlexDir.resolve("${it.id}.cellMetrics.parquet")
					def rawMatrixFile = scalePlexDir.resolve("${it.id}.raw.matrix/matrix.mtx.gz")
					cellMetricsFile.exists() && rawMatrixFile.exists()
				}
				.map {
					def scalePlexDir = file(it.resultDir) / "scaleplex"
					[
						[lib: it.libName, sample: it.sample],
						scalePlexDir.resolve("${it.id}.cellMetrics.parquet"),
						scalePlexDir.resolve("${it.id}.raw.matrix/matrix.mtx.gz"),
					]
				}
				.set { perSampleCountRaw }
		}
	} else {
		if (params.ultimaCramDir) {
			ULTIMA(
				samples,
				barcodeUnwrappedSamplesCsv,
				libraryInfo,
			)
			alignmentInput = ULTIMA.out.alignmentInput
			ubamHash = ULTIMA.out.ubamHash
			samples = ULTIMA.out.samples
			// Construct totalSampleReadsBySampleID channel and set it to 0 for all samples
			samples
				.filter { it.rnaId == null } // Only RNA samples are passed on to CELL_CALLING
				.map {
					tuple(it.id, 0)
				}
				.set { totalSampleReadsBySampleID }
		}
		else {
			// INPUT_READS gets/generates fastq files and runs bcParser to demux per sample
			INPUT_READS(
				samples,
				samplesCsv,
				libraryInfo,
				params.runFolder,
				params.fastqDir,
				runLibDetection
			)
			alignmentInput = INPUT_READS.out.ubam
			ubamHash = INPUT_READS.out.ubamHash
			totalSampleReadsBySampleID = INPUT_READS.out.totalSampleReadsBySampleID
			samples = INPUT_READS.out.samples
		}
		// STARSolo
		ALIGNMENT(alignmentInput, genome)
		ALIGNMENT.out.soloOut
			.dump(tag:'soloOut')
			.set { soloOut }
		ALIGNMENT.out.soloLog
			.dump(tag:'soloLog')
			.set { soloLog }
	}
	// Cell calling and metrics generation
	CELL_CALLING(
		samples,
		libraryInfo,
		soloOut,
		isBarnyard,
		soloLog,
		totalSampleReadsBySampleID,
	)
	//// ScalePlex
	allCellsWithAssignment = Channel.empty()
	if (params.scalePlex) {
		SCALE_PLEX(
			libraryInfo,
			samples,
			params.reporting ? [] : ubamHash, // No ubamHash for reporting run
			CELL_CALLING.out.allCells,
			CELL_CALLING.out.libCount,
			params.reporting ? perSampleCountRaw : [] // perSampleCountRaw is empty if not a --reporting run
		)
		SCALE_PLEX
			.out.metrics
			// [meta, cellMetrics, filteredMtx, scaleplex_stats.csv, metrics.csv]
			.map { [it[0].rnaId, it[1], it[3], it[4], it[0].lib] }
			// [rnaId, cellMetrics, scaleplex_stats.csv, metrics.csv, hashLibName]
			.set { scalePlexResults }
		SCALE_PLEX
			.out.allCells
			// [meta, cellMetrics]
			.map { [it[0].rnaId, it[1]] }
			// [rnaId, library, cellMetrics]
			.dump(tag:'allCellsWithAssignment')
			.set { allCellsWithAssignment }
	}
	//// REPORTING
	// Downstream analysis and per sample QC report
	if (params.seurat || params.azimuth || params.annData) {
		DOWNSTREAM(
			samples,
			CELL_CALLING.out.allCells,
			CELL_CALLING.out.filteredMtx,
			CELL_CALLING.out.cellCallingStats,
			libJson,
			false,
		)
		if (params.compSamples) {
			COMPARE_DOWNSTREAM(
				samples,
				CELL_CALLING.out.allCells,
				CELL_CALLING.out.filteredMtx,
				CELL_CALLING.out.cellCallingStats,
				libJson,
				true,
			)
		}
	}
	SAMPLE_REPORTING(
		samples,
		libJson,
		isBarnyard,
		CELL_CALLING.out.allCells,
		CELL_CALLING.out.allBarcodes,
		CELL_CALLING.out.libCount,
		CELL_CALLING.out.sampleStats,
		params.scalePlex ? scalePlexResults : [],
		false,
	)
	sampleStats = SAMPLE_REPORTING.out.sampleStats
	CELL_CALLING.out.allCells
		.join(allCellsWithAssignment, remainder: true)
		.map {
			id, lib, allCells, allCellsWithAssignment ->
			// If allCellsWithAssignment is null, use allCells
			[id, lib, allCellsWithAssignment ?: allCells]
		}
		.dump(tag:'allCells')
		.set { allCells }
	
	if (params.merge) {
		// Concatenate cell metrics, filtered mtx, and STAR log for individual libraries
		// Create merged sample reports and do downstream analysis
		MERGING(
			samples,
			libJson,
			isBarnyard,
			allCells,
			CELL_CALLING.out.allBarcodes,
			CELL_CALLING.out.filteredMtx,
			soloLog,
			params.scalePlex ? scalePlexResults : [],
			totalSampleReadsBySampleID,
		)
		sampleStats = sampleStats.concat(MERGING.out.sampleStats)
	}
	sampleStats.dump(tag:'sampleStats')
	MultiSampleReport(sampleStats.collect()) // Combined original and merged samples

	// Per library QC report
	// Usually not re-run for reporting because there are no changes but if --internalReport,
	// that tab has beads per passing cell which can change in reporting run
	if (!params.reporting || params.internalReport) {
		CELL_CALLING
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
		metricsByLib.dump(tag:'metricsByLib')

		// Get hash lib names
		samples
			.filter { it.rnaId != null }
			.collect { it.libName }
			.ifEmpty { [] }
			.dump(tag:'hashLibNames')
			.set { hashLibNames }
		LibraryMetricsGeneration(
			metricsByLib,
			hashLibNames,
			libraryInfo.quantum
		)

		// attach read metrics to library metrics for reporting
		if (params.ultimaCramDir) {
			// Empty list indicates bcParser demux metrics, which do not exist when starting from cram files
			LibraryMetricsGeneration.out
				.map { tuple(it[0], [], it[1]) }
				.set { libraryReportMetrics }
		}
		else if (params.reporting) {
			// Load metrics.json from previous pipeline output directory
			samples
				.map { it.subMap("libName", "resultDir") }
				.unique()
				.map {
					def barcodesDir = file(it.resultDir) / "barcodes"
					if (it.libName in hashLibNames.val) {
						barcodesDir = file(it.resultDir).resolve("scaleplex/demux")
					}
					// libName, metrics.json
					[it.libName, file(barcodesDir / "${it.libName}.metrics.json")]
				}
				// It's possible for the ScalePlex library to be dropped in the original run
				.filter { lib, metricsFile ->
					metricsFile.exists()
				}
				.join(LibraryMetricsGeneration.out)
				.set { libraryReportMetrics }
		} else {
			INPUT_READS
				.out.metrics.join(LibraryMetricsGeneration.out)
				.set { libraryReportMetrics }
		}
		libraryReportMetrics.dump(tag:'libraryReportMetrics')
		LibraryReportGeneration(
			libraryReportMetrics.join(CELL_CALLING.out.beadScores, remainder: true).map {
				libName, demuxJson, metrics, beadScores ->
				// [libName, demuxJson, metrics, beadScores]
				[libName, demuxJson, metrics, beadScores ?: []] // path value cannot be null
			},
			libraryInfo.rnaLibraryStructureFile.getName(),
			libraryInfo.rnaLibraryStructureFile.getParent(),
			libraryInfo.scalePlexLibraryStructureFile.getName(),
			libraryInfo.scalePlexLibraryStructureFile.getParent(),
			hashLibNames
		)
		
		// Barcode stats do not exist when starting from ultima cram files 
		if (params.internalReport && !params.ultimaCramDir) {
			MultiLibraryReport(
				LibraryReportGeneration.out.typeLevelMatches.collect(),
				LibraryReportGeneration.out.overallMatches.collect()
			)
		}
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

