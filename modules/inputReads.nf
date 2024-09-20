/*
Subworkflow for SaleRNA that deals with fastq generation, barcode processing, sample demux and qc of input data
bcl-convert, fastqc, bcParser, trimFq
*/

def getReadFromFqname(fqName) {
    m = ( fqName =~ /_([RI][12])[_.]/ )
    if (!m) { return null }
    return m[0][1]
}

// Extract sampleName (and optionally subsample, split by RT well) for one bcParser output
// @bcParserOut -> library name, fastq file
def constructSampleName (bcParserOut, splitFastq) {
	def libName = bcParserOut[0]
	def file = bcParserOut[1]
	def tok = file.getName().toString().tokenize('_') // bcParser outputs are "sampleName_well_...fastq.gz"
	def sampleName = tok[0]
	def subsampleName = sampleName
	if (splitFastq) {
		subsampleName = "${sampleName}_${tok[1]}"
	}
	return tuple(libName, sampleName, subsampleName, file)
}

def constructMatchingKey(fname) {
	def identifying_char_removed = fname.replaceAll(/_[RI][12][_\.]/, "")
	return identifying_char_removed
}

// Create a bcl-convert samplesheet for libraries in samples.json
process makeBclConvertSamplesheet {
input: 
	path(samplesCsv)
	val(hashLibNames) // Hash libraries
	path(libStructDir) // Directory containing the library structure definition (barcode sequence lists are loaded from here)
	val(libStructName) // Filename of the library structure definition .json
	path(runinfo) // RunInfo.xml from Sequencer RunFolder (Read-lengths etc.)
output: 
	path("samplesheet.csv")
publishDir file(params.outDir) / "fastq", mode: 'copy'
label 'small'
script:
	numHashRows = hashLibNames.size()
	libStruct = "$libStructDir/$libStructName"
	hashlibStruct = "$libStructDir/${params.scalePlexLibStructure}"
	opts = ""
	if (params.splitFastq) {
		opts += "--splitFastq "
	}
	if (params.scalePlex)
		"""
			# create RNA samplesheet
			awk -F',' 'NR <= $numHashRows + 1' $samplesCsv > rna_samples.csv
			bclConvertSheet.py rna_samples.csv $libStruct $runinfo $opts > samplesheet.csv
			# create Hash samplesheet
			awk -F',' 'NR==1 || NR > $numHashRows + 1' $samplesCsv > hash_samples.csv
			# run bclConvertSheet.py and append the output without header info to samplesheet.csv
			bclConvertSheet.py hash_samples.csv $hashlibStruct $runinfo $opts | awk '/\\[Data\\]/,EOF { print \$0 }' | sed '1,2d' >> samplesheet.csv
		"""
	else
		"""
			bclConvertSheet.py $samplesCsv $libStruct $runinfo $opts > samplesheet.csv
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
publishDir params.outDir, pattern: 'fastq/Reports/*', mode: 'copy'
publishDir params.outDir, pattern: 'fastq/*.fastq.gz', enabled: params.fastqOut


script:
	wthreads = (task.cpus/2 - 0.5).round()
	dthreads = (task.cpus/3 - 0.5).round()
"""
	bcl-convert --sample-sheet $samplesheet --bcl-input-directory $run \
    --bcl-num-conversion-threads $wthreads --bcl-num-compression-threads $wthreads --bcl-num-decompression-threads $dthreads \
    $params.bclConvertParams --output-directory fastq
"""
}

// Run cutadapt on demuxed fastq files
process trimFq {
input:
	tuple(val(id), val(subsample), path(transFastq), path(bcFastq), val(count))
output: 
	tuple(val(id), val(subsample), path("trimmed/${transName}"), path("trimmed/${bcName}"), emit: fastq)
	tuple(val(id), path("*.trim_stats.json"), emit: stats)
	path("*.trim_stats"), emit: stats_for_multiqc
tag "$id $count"

script:
	transName = transFastq.getName()
	bcName = bcFastq.getName()
	trimAdaptParam = params.scalePlex ? "${params.trimAdapt} ${params.scalePlexTrimAdapt}" : params.trimAdapt
"""
	mkdir trimmed
	cutadapt -j$task.cpus -e0.15 -O2 -m16 $trimAdaptParam -o trimmed/$transName -p trimmed/$bcName $transFastq $bcFastq --json ${id}.trim_stats.json | tee ${id}.${count}.trim_stats
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
publishDir file(params.outDir) / "fastq", mode: 'copy'
label 'small'
tag "${fq.getSimpleName()}"
"""
	mkdir fastqc
	fastqc -o fastqc $fq
"""
}

// Use multiQC to report a single QC report for all fastQC and bcl-convert reports
process multiqc {
memory { (reports.size() > 100 ? 12.GB : 6.GB) * task.attempt}
input:
	path(reports)
output:
	path("multiqc_report.html")
publishDir file(params.outDir) / "reports", mode: 'copy'
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
	// List of hash libraries, [] if not --scalePlex
	val(hashLibNames)
	// Directory containing the library structure definition (barcode sequence lists are loaded from here)
	path(libStructDir)
	// Filename of the library structure definition .json
	val(libStructName) 
	// Input fastq file
	tuple(val(libName), path(fqFiles), val(count))
output:
	//Output fastqs might not exist if the input has no passing reads
	tuple(val(libName), path("$demuxDir/*_S[1-9]*_R?[._]*fastq.gz"), emit: readFq, optional: true)
	tuple(val(libName), path("$demuxDir/*_S[1-9]*_BC[._]*fastq.gz"), emit: bcFq, optional: true)
	path("$demuxDir/*_S0_*.fastq.gz"), emit: unknown, optional: true
	path("$demuxDir/*.tsv")
	tuple(val(libName), path("$demuxDir/metrics.json"), emit: metrics)
	tuple(val(libName), path("$demuxDir/*.barcodes.csv"), emit: barcodesCsv, optional: true)

	
publishDir { file(params.outDir) / dir_name }, pattern: "$demuxDir/*gz", enabled: params.fastqOut
publishDir { file(params.outDir) / dir_name }, mode: 'copy', pattern: "$demuxDir/*{txt,tsv,json}"
tag "$libName $count"

script:
	demuxDir = "${libName}.demux"
	libStruct = "$libStructDir/$libStructName"
	opts = ""
	if (params.scalePlex && libName in hashLibNames) {
		opts += "--barcode-info"
		libStruct = "$libStructDir/${params.scalePlexLibStructure}"
		dir_name = "scaleplex/demux"
	} else {
		opts += "--write-fastq --write-barcode-fq"
		dir_name = "barcodes"
	}
	if (params.splitFastq) {
		dir_name = "$dir_name/bcparser.${libName}.${count}"
	}
"""
	bc_parser --lib-struct $libStruct --demux $sheet --lib-name $libName -v --reads ${fqFiles.join(" ")} --out $demuxDir $opts
"""
}

process merge_demux {
input:
    tuple(val(libName), path("demux_metrics*.json"))
	// List of hash libraries, [] if not --scalePlex
	val(hashLibNames)
	path(libStructDir)
    val(libStructName)
output:
    tuple(val(libName), path("*metrics.json"))
publishDir { file(params.outDir) / dir_name }, mode:'copy'
tag "$libName"
script:
	libStruct = "$libStructDir/$libStructName"
	if (params.scalePlex && libName in hashLibNames) {
		libStruct = "$libStructDir/${params.scalePlexLibStructure}"
		dir_name = "scaleplex/demux"
	} else {
		dir_name = "barcodes"
	}
    """
    mergeBCparserOutput.py --bc_jsons demux_metrics* --lib_json $libStruct --libName $libName
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
main:
	runDir = runFolder
	fqDir = fastqDir
	fqs = null

	// Get hash lib names
	samples
		.filter { it.rnaId != null }
		.collect { it.libName }
		.ifEmpty { [] }
		.dump(tag:'hashLibNames')
		.set { hashLibNames }
	
	if (runDir != null) {
		if (params.fastqSamplesheet == null) {
			makeBclConvertSamplesheet(samplesCsv, hashLibNames, libJson.getParent(), libJson.getName(), 
					          file("$runDir/RunInfo.xml", checkIfExists:true))
			fqSheet = makeBclConvertSamplesheet.out
		} else {
			fqSheet = file(params.fastqSamplesheet)
		}
		bclconvert(file(runDir), fqSheet)
		fqs = bclconvert.out.fastq.flatten()
	} else if (fqDir != null) {
		fqs = Channel.fromPath("$fqDir/**fastq.gz", checkIfExists: true)
	} else {
		ParamLogger.throwError("Must specify either '--runFolder' or '--fastqDir' when running alignment")
	}

	// fqs -> (fastq file)
	fqs.dump(tag:'fqs')
	
	// Organize fastq files by sample
	// fqFiles_with_fastq_matching_key -> (library_name, fastq_matching_key, file)
	fqFiles = fqs.map { file ->
		def fname = file.getName().toString()
		def libName = fname.tokenize('_')[0]
		// Construct matching key to join on later on in the workflow
		if (params.splitFastq) {
			return tuple(libName, constructMatchingKey(fname), file)
		}
		// Matching key is just the library name if params.splitFastq is false
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
	
	// Counter that we use for grouping bcParser jobs with the same libName together
	// Intention is to create params.bcParserJobs number of groups of fastq files belonging to the same libName
	// count%params.bcParserJobs will give us a number between 0 and params.bcParserJobs-1 which we can then group on
	def count = 0
	fqSamples = fqSamples.map{
		def libName = it[0]
		def fqFiles = it[2]
		count++
		tuple(libName, fqFiles, count % params.bcParserJobs)
	}
	// fqSamples -> (library name, [I1, R1, R2], counter)
	fqSamples.dump(tag:'fqSamples')

	// Group on libName and counter to create params.bcParserJobs number of groups of fastq files
	// map statement is to join all lists of fastq files together after grouping
	fqGroups = fqSamples.groupTuple(by:[0,2]).map{
		def fqFiles = []
		for (sublist in it[1]) {
			fqFiles += sublist
		}
		tuple(it[0], fqFiles, it[2])
	}
	fqGroups.dump(tag:'fqGroups')
	 
	// Check that for each library (unique libName in samples.csv), we have a complete set of fastq files
	// checkFastq -> (libName, sampleName, Fastqs, counter)
	// This slightly strange channel construction (including it.id) works around the lack of 'left-join'
        // in nextflow (this way we can distinguish 'left-only' and 'right-only' join results)
	checkFastq = samples.map{ [it.libName, it.id] }.unique{ it[0] }.join(fqSamples, remainder: true)
	checkFastq.dump(tag:'checkFastq')
	checkFastq.map {
		// If samples.csv has a library name which does not have corresponding fastq files
		// then 3rd element of checkFastq will be null
		if (it[2] == null) {
			ParamLogger.throwError("Library ${it[0]} does not have matching fastq files. None of the provided fastq filenames start with ${it[0]}")
		}
		// Check that each library has index1, read1 and read2 fastq files
		if (it[2].any {fq -> fq.getName().contains("_R1_") or fq.getName().contains("_R1.")} == false) {
			ParamLogger.throwError("Library ${it[0]} does not have read1 fastq file")
		}

		if (it[2].any {fq -> fq.getName().contains("_R2_") or fq.getName().contains("_R2.")} == false) {
			ParamLogger.throwError("Library ${it[0]} does not have read2 fastq file")
		}

		if (it[2].any {fq -> fq.getName().contains("_I1_") or fq.getName().contains("_I1.")} == false) {
			ParamLogger.throwError("Library ${it[0]} does not have index fastq file")
		}
	}

	// Process cell-barcodes and (optionally) split fastqs into samples based on tagmentation barcode
	barcodeParser(samplesCsv, hashLibNames, libJson.getParent(), libJson.getName(), fqGroups)

	// Check that there exists fastq files(R1, R2) post bcParser, for every sample that corresponds to a libName in samples.csv
	samples
		// check post bcParser fastq generation only for RNA libraries
		.filter { it.rnaId == null }
		.map{ [it.libName, it.id] }
		// get all samples corresponding to a libName
		.groupTuple()
		.join(barcodeParser.out.readFq.groupTuple(), remainder: true)
		// bcParserOutCheck -> (libName, [sample names corresponding to that libName], [[R1, R2 fastq files]])
		.dump(tag:'bcParserOutCheck')
		.map {
			if (it[2] == null) {
				ParamLogger.throwError("Library ${it[0]} does not have any passing R1 and R2 fastq files after bcParser")
			}
			// iterate over every sample name corresponding to a libName
			for (sample in it[1]) {
				// extract the sample name from the id
				sampleName = sample.tokenize('.')[0]
				exists = false
				// iterate over every fastq file for that libName and check that the sample name is present in the fastq file name
				for (file in it[2].flatten()) {
					if (file.getName().toString().contains(sampleName)) {
						exists = true
					}
				}
				if (!exists) {
					log.warn("Sample ${sampleName} does not have passing R1 and R2 files from bcParser.")
				}
			}
		}
	
	// Each barcodeParser run outputs fastqs for one libName
	// There can be multiple barcodeParser jobs for one libName if splitFastq is true
	// readFqs -> ([[sample_name]], fastq_matching_key, file)
	readFqs = barcodeParser.out.readFq.transpose().map{
		constructSampleName(it, params.splitFastq)
	}.groupTuple(by:[0,1,2], size:2) // We still get R1 and R2, even if only R2 will be used downstream
	readFqs.dump(tag:'readFqs')
	
	// Same as above but for the barcode reads (_BC_)
	// No group tuple because only one barcode read file is produced per sample
	bcFqs = barcodeParser.out.bcFq.transpose().map{
		constructSampleName(it, params.splitFastq)
	}
	bcFqs.dump(tag:'bcFqs')

	// Combine Read2 and barcode read for each sample / sample x RT combination
	// demuxFqs -> (sample_id, read2, barcode read)
	// matching key is dropped here
	// if params.splitFastq true, then sample_name="{sample_name}.{rt_well_coordinate}"
	readFqs.
		join(bcFqs, by:[0,1,2])
		.map{ [it[0], it[1], it[2], it[3][1], it[4]] }
		.dump(tag:'demuxFqs1')
		// bcParser outputs (demuxFqs) are by libName and sampleName.
		// We need to add the sample id (can occur multiple times for splitFastq subsamples)
		.combine(samples.map{ [it.libName, it.sample, it.id] }, by:[0,1])
		.dump(tag:'demuxFqs2')
		.map{ [it[5], it[2], it[3], it[4]] }
		//demuxFqs -> (sampleId, subsample, R2, BC)
		.dump(tag:'demuxFqs3')
		.set { demuxFqs }

	if (params.trimFastq) {
		demuxFqs = demuxFqs.toSortedList().flatMap()
		// We need a unique count on trimJobs to make output file names distinct
		// across split fastqs (when running from --fastqDir)
		def counter=1
		trimJobs = demuxFqs.map { it + (counter++) } 
		trimFq(trimJobs)
		demuxFqs = trimFq.out.fastq
		demuxFqs.dump(tag:'trimmedFqs')
	}
	if (params.fastqc) {
		// Exclude I1 and I2 from fastqc
		mainReads = fqSamples.flatMap{it.get(1)}.filter{getReadFromFqname(it.getName())[0] == 'R'}
		channel.dump(tag:'fastqc')
		fastqc(mainReads)
		reports = fastqc.out.zip
		if (runDir != null) {
			reports = reports.mix(bclconvert.out.stats)
		}
		if (params.trimFastq) {
			reports = reports.mix(trimFq.out.stats_for_multiqc)
		}
		multiqc(reports.collect())
	}

	if (params.splitFastq) {
		// Need to merge the bcParser outputs for a single library together
		organized_demux = barcodeParser.out.metrics.groupTuple()
		merge_demux(organized_demux, hashLibNames, libJson.getParent(), libJson.getName())
		metrics = merge_demux.out
    	}
	else {
		metrics =  barcodeParser.out.metrics
	}
	// create matching key for barcodesCsv output from bcParser
	barcodeParser
		// barcodesCsv -> [libName, [sample_1A_S1.barcodes.csv, sample_1B_S1.barcodes.csv, ...]]
		.out.barcodesCsv
		.transpose()
		.map { constructSampleName(it, params.splitFastq) }
		// collect barcodesCsv file for each lib [0], sample [1], and subsample [2] (or well)
		// subsample will equal sample if splitFastq is false 
		.groupTuple(by:[0, 1, 2])
		.map { lib, sample, subsample, file ->
			meta = [lib:lib, sample:sample, subsample:subsample]
			[meta, file]
		}
		// barcodesCsv -> [[lib, sample, sample_1A], sample_1A_S1.barcodes.csv, sample_1A_S1.barcodes.csv, ...]]
		.dump(tag:'barcodesCsv')
		.set { barcodesCsv }

emit:
	fqs = demuxFqs
	metrics = metrics
	trimFqStats = trimFq.out.stats
	barcodesCsv = barcodesCsv 
	hashLibNames = hashLibNames
}
