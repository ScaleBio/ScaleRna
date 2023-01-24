nextflow.enable.dsl=2

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
	path("fastq/Reports"), emit: stats
publishDir "${params.outDir}/", pattern: 'fastq/Reports/*', mode: 'copy'
publishDir "${params.outDir}/", pattern: 'fastq/*.fastq.gz'

script:
	pthreads = ((task.cpus-4)/3).round()
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
	tuple(val(name), path("trimmed/${transName}.fastq.gz"), path("trimmed/${bcName}.fastq.gz"), emit: fastq)
	path("${name}.trim_stats"), emit: stats
tag "$name"

script:
	transName = transFastq.getSimpleName()
	bcName = bcFastq.getSimpleName()
"""
	mkdir trimmed
	cutadapt -j${task.cpus} -a A{8}N{100} -e0.15 -O2 -m16 -o trimmed/${transName}.fastq.gz -p trimmed/${bcName}.fastq.gz $transFastq $bcFastq > ${name}.trim_stats
"""
}

// Run fastQC on (a pair of) fastq files from inputs or bcl-convert
// Note that these are the input fastq files, not outputs of bcParser
process fastqc {
input:
	path(fqs)
output:
	path("fastqc/*.html"), emit: html
	path("fastqc/*.zip"), emit: zip
publishDir "${params.outDir}/fastq/fastqc", mode: 'copy'
label 'small'

"""
	mkdir fastqc
	fastqc -o fastqc $fqs
"""
}

// Use multiQC to report a single QC report for all fastQC and bcl-convert reports
process multiqc {
input:
	path(reports)
output:
	path("multiqc_report.html")
publishDir "${params.outDir}/reports", mode: 'copy'
tag "small"

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
	path("$outDir/*_S0_*.fastq.gz"), emit: unknown optional true
	path("$outDir/*.tsv")
	tuple(val(libName), path("$outDir/metrics.json"), emit: metrics)
publishDir "${params.outDir}/demux/", pattern: "$outDir/*gz"
publishDir "${params.outDir}/demux/", mode: 'copy', pattern: "$outDir/*{txt,tsv,json}" //, saveAs: {"${libName}.${it.getName()}"}
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
		fastqc(fqSamples.flatMap{it.get(1)})
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

process starsolo {
input: 
	path(indexDir) // STAR genome index
	val(library) // library.json (contains STAR barcode parameters)
	tuple(val(sample), path(transcriptFq), path(barcodeFq)) // transcript and BC fastq file
output: 
	tuple(val(sample), path("$starDir/Aligned.sortedByCoord.out.bam*"), emit: bam)
	tuple(val(sample), path("$starDir/Log*"), emit: log)
	tuple(val(sample), path("${sample}.Solo.out"), emit: solo)
publishDir "${params.outDir}/star", pattern: "${sample}.Solo.out", mode: 'copy', saveAs: {'Solo.out'}
publishDir "$params.outDir/star"
//memory { ( indexDir.directorySize() < 32000000000 ? 32.GB : 64.GB ) * task.attempt }
tag "$sample"
script:
	starDir = "${sample}.align"
	barcodeParam = library["star_barcode_param"]
"""
	STAR --runThreadN $task.cpus --genomeDir $indexDir  --outSAMtype BAM SortedByCoordinate --soloCellReadStats Standard \
	$barcodeParam ${params.starTrimming} \
	--soloStrand ${params.starStrand} --soloFeatures ${params.starFeature} --soloMultiMappers PropUnique \
	--readFilesIn $transcriptFq $barcodeFq --readFilesCommand zcat \
    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
	--outFileNamePrefix $starDir/
	mv $starDir/Solo.out ${sample}.Solo.out
"""
}

process sampleReport {
input:
	tuple(val(sample), val(libName), path("demuxMetrics.json"), path("Solo.out"), val(expectedCells))
	path(samplesCsv)
	path(libJson)
output:
	path("reports/${sample}.report.html")
	path("reports/${sample}.reportStatistics.tsv"), emit: stats
	path("reports/${sample}_figures")
publishDir "$params.outDir", mode: 'copy'
errorStrategy 'ignore'
label 'report'
tag "$sample"
script:
"""
	generateReport.py --sample ${sample} --samplesheet ${samplesCsv} --libStruct ${libJson}
"""
}

process libraryReport {
input:
	tuple(val(libName), path("demuxMetrics.json"), path(soloFiles))
	path(samplesCsv)
	path(libJson)
output:
	path("reports/library_${libName}.report.html")
	path("reports/library_${libName}*.tsv")
publishDir "$params.outDir", mode: 'copy'
errorStrategy 'ignore'
label 'report'
tag "$libName"
script:
"""
	generateReport.py --libName ${libName} --samplesheet ${samplesCsv} --libStruct ${libJson}
"""
}

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
	samples.dump(tag:'samples')
	libJson = expandPath(params.libStructure, file("${projectDir}/references/"))
	libStructure = loadJson(libJson)
	genome = loadGenome(file(params.genome))

	inputReads(samples, samplesCsv, libJson)
	inputReads.out.metrics.dump(tag:'demuxMetrics')
	sampleDemuxMetrics = inputReads.out.metrics.cross(samples.map{[it.libName,it.sample]}).map({[it[1][1], it[0][0], it[0][1]]})

	starsolo(genome.star_index, libStructure,  inputReads.out.fqs)
	starsolo.out.solo.dump(tag:'solo')

	reportTemplate = file("${projectDir}/report/reportRna.ipynb")
	expectedCells = samples.map{[it.sample, toIntOr0(it.expectedCells)]}
	stats =  sampleDemuxMetrics.join(starsolo.out.solo).join(expectedCells)
	stats.dump(tag:'stats')
	sampleReport(stats, samplesCsv, libJson)
	multiSampleReport(sampleReport.out.stats.collect())

	// Library level report (all samples with same libName in one report)
	statsByLib = stats.map{ row -> 
		return tuple(row[1], row[2], row[3]) // (libName, libMetrics, soloOut)
	}.groupTuple(by:[0,1]) // groupBy (libName, libMetrics)
	statsByLib.dump(tag:'statsGroupedFastq')
	libraryReport(statsByLib, samplesCsv, libJson)
}
