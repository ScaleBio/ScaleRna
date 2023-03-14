nextflow.enable.dsl=2

import groovy.json.JsonOutput
import groovy.json.JsonSlurper

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

def getReadFromFqname(fqName) {
    m = ( fqName =~ /_([RI][12])[_.]/ )
    if (!m) { return null }
    return m[0][1]
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

def throwError(errMessage) {
	log.error errMessage
	sleep(200)
	System.exit(1)
}

// Prepare samples.csv with defaults, rename legacy columns, etc.
process regularizeSamplesCsv {
input: 
	path("samples.in.csv")
output: 
	path("samples.csv")
publishDir "${params.outDir}", mode: 'copy'
label 'small'

"""
	regularizeSamplesCsv.py samples.in.csv > samples.csv
"""
}

// Create a bcl-convert samplesheet for libraries in samples.json
process makeBclConvertSamplesheet {
input: 
	path(samplesCsv)
	path(lib)
	path(runinfo)
output: 
	path("samplesheet.csv")
publishDir "${params.outDir}/fastq", mode: 'copy'
label 'small'

"""
	bclConvertSheet.py $samplesCsv $lib $runinfo > samplesheet.csv
"""
}

/*Run bcl-convert, used when starting from a sequencing run-folder
  Requires a separate bcl-convert samplesheet*/
process bclconvert {
input: 
	path(run)
	path(samplesheet)
output: 
	path("fastq/*fastq.gz"), emit: fastq
	path("fastq/Reports/*"), emit: stats
publishDir "${params.outDir}/", pattern: 'fastq/Reports/*', mode: 'copy'
publishDir "${params.outDir}/", pattern: 'fastq/*.fastq.gz', enabled: params.fastqOut


script:
	pthreads = ((task.cpus)/3 - 0.5).round()
"""
	bcl-convert --sample-sheet $samplesheet --bcl-input-directory $run \
    --bcl-num-conversion-threads $pthreads --bcl-num-compression-threads $pthreads --bcl-num-decompression-threads $pthreads \
    $params.bclConvertParams --output-directory fastq
"""
}

process trimFq {
input:
	tuple(val(name), path(transFastq), path(bcFastq))
output: 
	tuple(val(name),path("trimmed/${transName}.fastq.gz"), path("trimmed/${bcName}.fastq.gz"), emit: fastq)
	path("${name}.trim_stats"), emit: stats
tag "$name"

script:
	transName = transFastq.getSimpleName()
	bcName = bcFastq.getSimpleName()
"""
	mkdir trimmed
	cutadapt -j${task.cpus} -e0.15 -O2 -m16 ${params.trimAdapt} -o trimmed/${transName}.fastq.gz -p trimmed/${bcName}.fastq.gz $transFastq $bcFastq > ${name}.trim_stats
"""
}

// Run fastQC on (a pair of) fastq files from inputs or bcl-convert
// Note that these are the input fastq files, not outputs of bcParser
process fastqc {
input:
	path(fq)
output:
	path("fastqc/*.html"), emit: html
	path("fastqc/*.zip"), emit: zip
publishDir "${params.outDir}/fastq/", mode: 'copy'
label 'small'
tag "${getReadFromFqname(fq)}"
"""
	mkdir fastqc
	fastqc -o fastqc $fq
"""
}

// Use multiQC to report a single QC report for all fastQC and bcl-convert reports
process multiqc {
input:
	path(reports)
output:
	path("multiqc_report.html")
publishDir "${params.outDir}/reports", mode: 'copy'
label "small"

"""
	multiqc .
"""
}

// Run bcParser to extract and correct cell-barcodes
// Optionally demuxes fastq files based on some barcodes
process barcodeParser {
input:
	path(sheet) // Samplesheet json
	path(libStructDir) // Directory containing the library structure definition (barcode sequence lists are loaded from here)
	val(libStructName) // Filename of the library structure definition .json
	tuple(val(libName), path(fqFiles)) // Input fastq file
output:
	tuple(val(libName), path("$outDir/*_S[1-9]*_R?[._]*fastq.gz"), emit: readFq)
	tuple(val(libName), path("$outDir/*_S[1-9]*_BC[._]*fastq.gz"), emit: bcFq)
	path("$outDir/*_S0_*.fastq.gz"), emit: unknown, optional: true
	path("$outDir/*.tsv")
	tuple(val(libName), path("$outDir/metrics.json"), emit: metrics)
publishDir "${params.outDir}/barcodes/", pattern: "$outDir/*gz", enabled: params.fastqOut
publishDir "${params.outDir}/barcodes/", mode: 'copy', pattern: "$outDir/*{txt,tsv,json}"
tag "$libName"

script:
	outDir = "${libName}.demux"
	libStruct = "$libStructDir/$libStructName"
"""
	bc_parser --lib-struct $libStruct --demux $sheet --lib-name $libName -v --reads ${fqFiles.join(" ")} --write-fastq --write-barcode-fq --out $outDir
"""
}

// Fastq generation, trimming, QC, barcode extraction and sample demux
workflow inputReads {
take:
	samples
	samplesCsv
	libJson
main:
	runDir = params.runFolder
	fqDir = params.fastqDir
	fqs = null
	if (runDir != null) {
		if (params.fastqSamplesheet == null) {
			makeBclConvertSamplesheet(samplesCsv, libJson, 
				file("$runDir/RunInfo.xml", checkIfExists:true))
			fqSheet = makeBclConvertSamplesheet.out
		} else {
			fqSheet = file(params.fastqSamplesheet)
		}
		bclconvert(file(runDir), fqSheet)
		fqs = bclconvert.out.fastq.flatten()
	} else if (fqDir != null) {
		fqs = Channel.fromPath("$fqDir/*fastq.gz", checkIfExists: true)
	} else {
		throwError("Must specify either 'runFolder' or 'fastqDir'")
	}

	// Organize fastq files by sample
	fqFiles = fqs.map { file ->
		def ns = file.getName().toString().tokenize('_')
		return tuple(ns.get(0), file)
	}.groupTuple()
	fqSamples = samples.map({it.libName}).unique().join(fqFiles)
	fqSamples.dump(tag:'fqSamples')

	// Process cell-barcodes and (optionally) split fastqs into samples based on tagmentation barcode
	barcodeParser(samplesCsv, libJson.getParent(), libJson.getName(), fqSamples)
	readFqs = barcodeParser.out.readFq.flatMap({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('_')
		return tuple(ns.get(0), file)
	}.groupTuple(size:2) // We still get R1 and R2, even if only R2 is used
	bcFqs = barcodeParser.out.bcFq.flatMap({it[1]}).map { file ->
		def ns = file.getName().toString().tokenize('_')
		return tuple(ns.get(0), file)
	}
	demuxFqs = readFqs.join(bcFqs).map({tuple(it[0], it[1][1], it[2])})
	demuxFqs.dump(tag:'demuxFqs')

	if (params.trimFastq) {
		trimFq(demuxFqs)
		demuxFqs = trimFq.out.fastq
		demuxFqs.dump(tag:'trimmedFqs')
	}
	if (params.fastqc) {
		// Exclude I1 and I2 from fastqc
		mainReads = fqSamples.flatMap{it.get(1)}.filter{getReadFromFqname(it.getName())[0] == 'R'}
		fastqc(mainReads)
		reports = fastqc.out.zip
		if (runDir != null) {
			reports = reports.mix(bclconvert.out.stats)
		}
		if (params.trimFastq) {
			reports = reports.mix(trimFq.out.stats)
		}
		multiqc(reports.collect())
	}

emit:
	fqs = demuxFqs
	metrics = barcodeParser.out.metrics
}

// Run the STAR executable on demuxed fastq files
process starsolo {
input: 
	path(indexDir) // STAR genome index
	val(library) // library.json (contains STAR barcode parameters)
	tuple(val(sample), path(transcriptFq), path(barcodeFq)) // transcript and BC fastq file
output: 
	tuple(val(sample), path("$starDir/Aligned.sortedByCoord.out.bam*"), emit: bam, optional: true)
	tuple(val(sample), path("$starDir/Log*"), emit: log)
	tuple(val(sample), path("${sample}.star.solo"), emit: solo)
publishDir "${params.outDir}/alignment", pattern: "*.solo", mode: 'copy'
publishDir "${params.outDir}/alignment", pattern: "$starDir/*"
//memory { ( indexDir.directorySize() < 32000000000 ? 32.GB : 64.GB ) * task.attempt }
tag "$sample"
script:
	starDir = "${sample}.star.align"
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
	--soloStrand ${params.starStrand} --soloFeatures ${params.starFeature} --soloMultiMappers PropUnique \
	--readFilesIn $transcriptFq $barcodeFq --readFilesCommand zcat \
	--outFileNamePrefix $starDir/
	mv $starDir/Solo.out ${sample}.star.solo
"""
}

// Run cell typing on STARsolo output
process cellTyping {
input:
	tuple(val(sample), val(libName), path("demuxMetrics.json"), path("Solo.out"), val(expectedCells), path("filtered_matrix"))
output:
	path("${sample}")
tag "$sample"
publishDir "$params.outDir/cellTyping", mode: 'copy'
script:
"""
	runScrublet.py --counts filtered_matrix --outDir .
	assignTypes.R --projName=${sample} --countDir=filtered_matrix --scrubletOut=scrublet_output_table.csv --StarMtxFormat=yes
"""
}

// Generate metrics per sample from STARsolo output
process sampleMetricsGeneration {
input:
	tuple(val(sample), val(libName), path("demuxMetrics.json"), path("Solo.out"), val(expectedCells))
	path(samplesCsv)
	path(libJson)
	path("references")
	val(isBarnyard)
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
	}
	else {
		opts = ""
	}
	if (params.useSTARthreshold) {
		opts = opts + "--useSTARthreshold"
	}
	else {
		opts = opts + ""
	}

"""
	getSampleMetrics.py --sample ${sample} --samplesheet ${samplesCsv} --libStruct ${libJson} \
	--topCellPercent ${params.topCellPercent} \
	--minCellRatio ${params.minCellRatio} --minReads ${params.minReads} --star_out Solo.out \
	$opts
"""
}

// Generate report for each sample from metrics generated in sampleMetricsGeneration process
process sampleReportGeneration {
input:
	tuple(val(sample), val(libName), path("${sample}_metrics/allCells.csv"), path("${sample}_metrics/sample_metrics"))
	path(libJson)
	val(isBarnyard)
output:
	path("reports/${sample}.report.html")
	path("reports/${sample}.reportStatistics.tsv"), emit: stats
	path("reports/${sample}_figures"), optional: true
	path("reports/${sample}_unique_transcript_counts_by_RT_Index_well.csv")
	path("reports/${sample}_num_cells_by_RT_Index_well.csv")
publishDir "$params.outDir", mode: 'copy'
errorStrategy 'ignore'
label 'report'
tag "$sample"
script:
	if (isBarnyard) {
                opts = "--isBarnyard"
        }
    else {
                opts = ""
    }
	if (params.internalReport) {
		opts = opts + "--internalReport"
	}
"""
	export TMPDIR=\$(mktemp -p `pwd` -d)
	generateSampleReport.py --sampleName ${sample} --libDir ${libJson} --sample_metrics ${sample}_metrics $opts
"""
}

// Generate metrics for each library from the metrics of all samples in the library
process fastqMetricsGeneration {
input:
	tuple(val(sample), val(libName), path("allCells*.csv"), path("sample_metrics"))
	path(samplesCsv)
output:
	tuple(val(libName), path("library_metrics"))
errorStrategy 'ignore'
label 'report'
tag "$libName"
script:
"""
	getFastqMetricsFromSampleMetrics.py --sample_metrics allCells* --samplesheet ${samplesCsv} --libName ${libName}
"""
}

// Generate report for each library from metrics generated in fastqMetricsGeneration
process fastqReportGeneration {
input:
	tuple(val(libName), path("metrics"))
	path(libJson)
	path("references")
output:
	path("reports/library_${libName}.report.html")
	path("reports/library_${libName}*.tsv")
	path("reports/unique_reads_ligation_well.csv")
	path("reports/unique_reads_pcr_well.csv")
	path("reports/unique_reads_rt_well.csv")
publishDir "$params.outDir", mode: 'copy'
label 'report'
tag "$libName"
script:
	if (params.internalReport) {
		opts = "--internalReport"
	}
	else {
		opts = ""
	}
"""
	export TMPDIR=\$(mktemp -p `pwd` -d)
	generateFastqReport.py --libName ${libName} --libDir ${libJson} --libMetrics metrics $opts
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
	outFn = "allSamples.reportStatistics.tsv"
"""
	mergeReportStats.py $sampleStats > $outFn
"""
}

//// Main entry point
// Run the workflow for one or multiple samples
// either from one runFolder or one / multiple sets of fastq files
workflow {
	// Inputs
	regularizeSamplesCsv(file(params.samples))
	samplesCsv = regularizeSamplesCsv.out
	samples = samplesCsv.splitCsv(header:true, strip:true)
	
	libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
	libStructure = loadJson(libJson)
	genome = loadGenome(file(params.genome))

	isBarnyard = genome.get('isBarnyard', false)
	
	inputReads(samples, samplesCsv, libJson)
	
	sampleDemuxMetrics = inputReads.out.metrics.cross(samples.map{[it.libName,it.sample]}).map({[it[1][1], it[0][0], it[0][1]]})

	starsolo(genome.star_index, libStructure,  inputReads.out.fqs)

	expectedCells = samples.map{[it.sample, toIntOr0(it.expectedCells)]}
	stats =  sampleDemuxMetrics.join(starsolo.out.solo).join(expectedCells)
	
	sampleMetricsGeneration(stats, samplesCsv, libJson, "${projectDir}/references", isBarnyard)
	sampleMetricsBySample = sampleMetricsGeneration.out.cell_metrics.join(sampleMetricsGeneration.out.sample_metrics_json_for_report)
	sampleReportGeneration(sampleMetricsBySample, libJson, isBarnyard)
	
	// Group by libName
	metricsByLib = sampleMetricsBySample.groupTuple(by: 1)
	fastqMetricsGeneration(metricsByLib, samplesCsv)
	multiSampleReport(sampleReportGeneration.out.stats.collect())
	
	fastqReportGeneration(fastqMetricsGeneration.out, libJson, "${projectDir}/references")
	cellTypingInput = stats.join(sampleMetricsGeneration.out.filtered_star_mtx)
	if (params.cellTyping) {
		cellTyping(cellTypingInput)
	}
}
workflow.onComplete {
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
		         						 "Duration": "$workflow.duration"]]
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
