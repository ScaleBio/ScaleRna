/*
inputReads module that deals with preprocessing, demuxing and qc of input data
*/
def getReadFromFqname(fqName) {
    m = ( fqName =~ /_([RI][12])[_.]/ )
    if (!m) { return null }
    return m[0][1]
}

// @channel_contents -> library name, matching key, fastq file
// Function to add well coordinate to sample name to facilitate grouping downstream
// Sample name is computed from filename and the well coordinate we're getting from the matching key
// Filename is of the form {sampleName}_{wellCoordinate}*.fastq.gz
// Function only called when splitFastq set to True
def constructSampleNameWithWellCoordinate (channel_contents) {
	def libName = channel_contents[0]
	def matching_key = channel_contents[1]
	def file = channel_contents[2]
	def sample = file.getName().toString().tokenize('_')[0] + "." + channel_contents[1].tokenize("_")[1]
	return tuple(sample, matching_key, file)
}

// @channel_contents -> library name, matching key, fastq file
// Function to compute sample name from filename
// Filename is of the form {sampleName}_*.fastq.gz
// Function only called when splitFastq set to False or when splitFastq is True but we're starting from fastq
// We do not allow underscores in the sample name so we can take the
// first element after tokenizing on "_"
def constructSampleName (channel_contents) {
	def libName = channel_contents[0]
	def matching_key = channel_contents[1]
	def file = channel_contents[2]
	def sample = file.getName().toString().tokenize('_')[0]
	return tuple(sample, matching_key, file)
}

def throwError(errMessage) {
	log.error errMessage
	sleep(200)
	System.exit(1)
}

def constructMatchingKey(fname) {
	def identifying_char_removed = fname.replaceAll(/_[RI][12][_\.]/, "")
	return identifying_char_removed
}

// Create a bcl-convert samplesheet for libraries in samples.json
process makeBclConvertSamplesheet {
input: 
	path(samplesCsv)
	path(libStructDir) // Directory containing the library structure definition (barcode sequence lists are loaded from here)
	val(libStructName) // Filename of the library structure definition .json
	path(runinfo) // RunInfo.xml from Sequencer RunFolder (Read-lengths etc.)
	val(outDir)
	val(splitFastq)
output: 
	path("samplesheet.csv")
publishDir "${outDir}/fastq", mode: 'copy'
label 'small'
script:
	libJson = "$libStructDir/$libStructName"
	opts = ""
	if (splitFastq) {
		opts += "--splitFastq"
	}
"""
	bclConvertSheet.py $samplesCsv $libJson $runinfo $opts > samplesheet.csv
"""
}

/*Run bcl-convert, used when starting from a sequencing run-folder
  Requires a separate bcl-convert samplesheet*/
process bclconvert {
input: 
	path(run)
	path(samplesheet)
	val(outDir)
	val(bclConvertParams)
output: 
	path("fastq/*fastq.gz"), emit: fastq
	path("fastq/Reports/*"), emit: stats
publishDir "${outDir}/", pattern: 'fastq/Reports/*', mode: 'copy'
publishDir "${outDir}/", pattern: 'fastq/*.fastq.gz', enabled: params.fastqOut


script:
	wthreads = (task.cpus/2 - 0.5).round()
	dthreads = (task.cpus/3 - 0.5).round()
"""
	bcl-convert --sample-sheet $samplesheet --bcl-input-directory $run \
    --bcl-num-conversion-threads $wthreads --bcl-num-compression-threads $wthreads --bcl-num-decompression-threads $dthreads \
    $bclConvertParams --output-directory fastq
"""
}

// Run cutadapt on demuxed fastq files
process trimFq {
input:
	tuple(val(name), path(transFastq), path(bcFastq), val(count))
	val(trimAdapt)
output: 
	tuple(val(name), path("trimmed/${transName}"), path("trimmed/${bcName}"), emit: fastq)
	tuple(val(name), path("*.trim_stats.json"), emit: stats)
	path("*.trim_stats"), emit: stats_for_multiqc
tag "$name $count"

script:
	transName = transFastq.getName()
	bcName = bcFastq.getName()
"""
	mkdir trimmed
	cutadapt -j${task.cpus} -e0.15 -O2 -m16 ${trimAdapt} -o trimmed/${transName} -p trimmed/${bcName} $transFastq $bcFastq --json ${name}.trim_stats.json | tee ${name}.${count}.trim_stats
"""
}

// Run fastQC on (a pair of) fastq files from inputs or bcl-convert
// Note that these are the input fastq files, not outputs of bcParser
process fastqc {
input:
	path(fq)
	val(outDir)
output:
	path("fastqc/*.html"), emit: html
	path("fastqc/*.zip"), emit: zip
publishDir "${outDir}/fastq/", mode: 'copy'
label 'small'
tag "${fq.getSimpleName()}"
"""
	mkdir fastqc
	fastqc -o fastqc $fq
"""
}

// Use multiQC to report a single QC report for all fastQC and bcl-convert reports
process multiqc {
input:
	path(reports)
	val(outDir)
output:
	path("multiqc_report.html")
publishDir "${outDir}/reports", mode: 'copy'
label 'optional'
"""
	multiqc .
"""
}

// Run bcParser to extract and correct cell-barcodes
// Optionally demuxes fastq files based on some barcodes
process barcodeParser {
input:
	// samples.csv
	path(sheet)
	// Directory containing the library structure definition (barcode sequence lists are loaded from here)
	path(libStructDir)
	// Filename of the library structure definition .json
	val(libStructName) 
	// Input fastq file
	tuple(val(libName), val(matching_key), path(fqFiles), val(count))
	val(outputDir)
output:
	//Output fastqs might not exist if the input has no passing reads
	tuple(val(libName), val(matching_key), path("$outDir/*_S[1-9]*_R?[._]*fastq.gz"), emit: readFq, optional: true)
	tuple(val(libName), val(matching_key), path("$outDir/*_S[1-9]*_BC[._]*fastq.gz"), emit: bcFq, optional: true)
	path("$outDir/*_S0_*.fastq.gz"), emit: unknown, optional: true
	path("$outDir/*.tsv")
	tuple(val(libName), path("$outDir/metrics.json"), emit: metrics)
publishDir "${outputDir}/${dir_name}", pattern: "$outDir/*gz", enabled: params.fastqOut
publishDir "${outputDir}/${dir_name}", mode: 'copy', pattern: "$outDir/*{txt,tsv,json}"
tag "$libName $count"

script:
	outDir = "${libName}.demux"
	libStruct = "$libStructDir/$libStructName"
	if (params.splitFastq) {
		dir_name = "barcodes/bcparser.${count}/"
	} else {
		dir_name = "barcodes/"
	}
"""
	bc_parser --lib-struct $libStruct --demux $sheet --lib-name $libName -v --reads ${fqFiles.join(" ")} --write-fastq --write-barcode-fq --out $outDir
"""
}

process merge_demux {
input:
    tuple(val(libName), path("demux_metrics*.json"))
    path(libJson)
	val(outDir)
output:
    tuple(val(libName), path("metrics.json"))
    publishDir "${outDir}/barcodes", mode:'copy'
tag "$libName"
script:
    """
    mergeBCparserOutput.py --bc_jsons demux_metrics* --lib_json ${libJson}
    """
}

// Fastq generation, trimming, QC, barcode extraction and sample demux
workflow inputReads {
take:
	samples
	samplesCsv
	libJson
    runFolder
    fastqDir
    fastqSamplesheet
    trimFastq
    fastqc
    outDir
    bclConvertParams
	splitFastq
	trimAdapt
main:
	runDir = runFolder
	fqDir = fastqDir
	fqs = null
	if (runDir != null) {
		if (fastqSamplesheet == null) {
			makeBclConvertSamplesheet(samplesCsv, libJson.getParent(), libJson.getName(), 
				file("$runDir/RunInfo.xml", checkIfExists:true), outDir, splitFastq)
			fqSheet = makeBclConvertSamplesheet.out
		} else {
			fqSheet = file(fastqSamplesheet)
		}
		bclconvert(file(runDir), fqSheet, outDir, bclConvertParams)
		fqs = bclconvert.out.fastq.flatten()
	} else if (fqDir != null) {
		fqs = Channel.fromPath("$fqDir/*fastq.gz", checkIfExists: true)
	} else {
		throwError("Must specify either 'runFolder' or 'fastqDir'")
	}

	// fqs -> (fastq file)
	fqs.dump(tag:'fqs')
	
	// Organize fastq files by sample
	// fqFiles_with_fastq_matching_key -> (library_name, fastq_matching_key, file)
	fqFiles = fqs.map { file ->
		def fname = file.getName().toString()
		def libName = fname.tokenize('_')[0]
		// Construct matching key to join on later on in the workflow
		if (splitFastq) {
			def matching_key = constructMatchingKey(fname)
			return tuple(libName, constructMatchingKey(fname), file)
		}
		// Matching key is just the library name if splitFastq is false
		else {
			return tuple(libName, libName, file)
		}
	}
	pairedFqFiles = fqFiles.groupTuple(by:[0,1])
	
	// This keeps only those sets of fastq files that correspond to a sample from our samples.csv
	// it.libName takes only the library name from samples csv
	// unique because samples csv has multiple entries for multiple samples, but only one library name
	// cross ensures that the right library names are paired with the right fastq files
	// The last map is to remove the extra library name element(it[0]) in the tuple that we get after the cross
	fqSamples = samples.map({it.libName}).unique().cross(pairedFqFiles).map{it[1]}
	
	// The counter aids in output directory naming when splitFastq is set to true
	def count = 1
	fqSamples = fqSamples.map{
		def libName = it[0]
		def fastq_matching_key = it[1]
		def fqFiles = it[2]
		tuple(libName, fastq_matching_key, fqFiles, count++)
	}
	// fqSamples -> (library name, matching key, [I1, R1, R2], counter) 
	fqSamples.dump(tag:'fqSamples')
	
	// Check that for each library (unique libName in samples.csv), we have a complete set of fastq files
	// checkFastq -> (libName, sampleName, matchingKey, Fastqs, counter)
	// This slightly strange channel construction (including it.sample) works around the lack of 'left-join'
        // in nextflow (this way we can distinguish 'left-only' and 'right-only' join results)
	checkFastq = samples.map{[it.libName, it.sample]}.unique({it[0]}).join(fqSamples, remainder: true)
	checkFastq.dump(tag:'checkFastq')
	checkFastq.map {
		// If samples.csv has a library name which does not have corresponding fastq files
		// then 4th element of checkFastq will be null
		if (it[3] == null) {
			throwError("Library ${it[0]} does not have matching fastq files. None of the provided fastq filenames start with ${it[0]}")
		}
		// Check that each library has index1, read1 and read2 fastq files
		if (it[3].any {fq -> fq.getName().contains("_R1_") or fq.getName().contains("_R1.")} == false) {
			throwError("Library ${it[0]} does not have read1 fastq file")
		}

		if (it[3].any {fq -> fq.getName().contains("_R2_") or fq.getName().contains("_R2.")} == false) {
			throwError("Library ${it[0]} does not have read2 fastq file")
		}

		if (it[3].any {fq -> fq.getName().contains("_I1_") or fq.getName().contains("_I1.")} == false) {
			throwError("Library ${it[0]} does not have index fastq file")
		}
	}

	// Process cell-barcodes and (optionally) split fastqs into samples based on tagmentation barcode
	barcodeParser(samplesCsv, libJson.getParent(), libJson.getName(), fqSamples, outDir)

	// If starting from runfolder and splitFastq is true, we do not split in bcparser,
	// so the sample name does not have the well coordinate in it
	// We add the well coordinate to the sample name by calling the constructSampleNameWithWellCoordinate function
	if (splitFastq && runDir != null)  {
		readFqs = barcodeParser.out.readFq.transpose().map{
			constructSampleNameWithWellCoordinate(it)
		}.groupTuple(by:[0,1], size:2)

		bcFqs = barcodeParser.out.bcFq.transpose().map{
			constructSampleNameWithWellCoordinate(it)
		}
	}
	// Indicates splitFastq is False, or splitFastq is True and we're starting the workflow from
	// fastq files. In that case we split in bcParser and the well coordinate is already attached
	else {
		// Each barcodeParser run outputs fastqs for each sample in one list
		// readFqs -> ([[sample_name]], fastq_matching_key, file)
		readFqs = barcodeParser.out.readFq.transpose().map{
			constructSampleName(it)
		}.groupTuple(by:[0,1], size:2) // We still get R1 and R2, even if only R2 will be used downstream
		readFqs.dump(tag:'readFqs')
		
		// Same as above but for the barcode reads (_BC_)
		// No group tuple because only one barcode read file is produced per sample
		bcFqs = barcodeParser.out.bcFq.transpose().map{
			constructSampleName(it)
		}
	}
	bcFqs.dump(tag:'bcFqs')
	// Combine Read2 and barcode read for each sample / sample x RT combination
	// demuxFqs -> (sample_name, read2, barcode read)
	// matching key is dropped here
	// if splitFastq true, then sample_name="{sample_name}.{rt_well_coordinate}"
	demuxFqs = readFqs.join(bcFqs, by:[0,1]).map({tuple(it[0], it[2][1], it[3])})
	demuxFqs.dump(tag:'demuxFqs')

	if (trimFastq) {
		demuxFqs = demuxFqs.toSortedList().flatMap()
		def counter=1
		demuxFqs = demuxFqs.map {
			def libName = it[0]
			def transFastq = it[1]
			def bcFastq = it[2]
			tuple(libName, transFastq, bcFastq, counter++)
		}
		trimFq(demuxFqs, trimAdapt)
		demuxFqs = trimFq.out.fastq
		
		//demuxFqs -> (sample name, R2, BC)
		demuxFqs.dump(tag:'trimmedFqs')
	}
	if (fastqc) {
		// Exclude I1 and I2 from fastqc
		mainReads = fqSamples.flatMap{it.get(2)}.filter{getReadFromFqname(it.getName())[0] == 'R'}
		channel.dump(tag:'fastqc')
		fastqc(mainReads, outDir)
		reports = fastqc.out.zip
		if (runDir != null) {
			reports = reports.mix(bclconvert.out.stats)
		}
		if (trimFastq) {
			reports = reports.mix(trimFq.out.stats_for_multiqc)
		}
		multiqc(reports.collect(), outDir)
	}

	if (splitFastq) {
		// Need to merge the bcParser outputs for a single library together
        organized_demux = barcodeParser.out.metrics.groupTuple()
        merge_demux(organized_demux, libJson, outDir)
		metrics = merge_demux.out
    }
	else {
		metrics =  barcodeParser.out.metrics
	}

emit:
	fqs = demuxFqs
	metrics = metrics
	trimFqStats = trimFq.out.stats
}
