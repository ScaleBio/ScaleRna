// Run the STAR executable on demuxed fastq files
process starsolo {
input: 
	path(indexDir) // STAR genome index
	val(library) // library.json (contains STAR barcode parameters)
	tuple(val(sample), val(count), path("transcript*.fastq.gz"), path("barcode*.fastq.gz")) // transcript and BC fastq file
output: 
	tuple(val(sample), path("Aligned.sortedByCoord.out.bam"), emit: bam, optional: true)
	tuple(val(sample), path("Log.*"), emit: log)
	tuple(val(sample), path("${sample}"), emit: solo)

// If params.splitFastq is true we publish individual split outputs into params.outDir/alignment/<sample>/split/
// with output filenames <sample>.<count>.*
// Otherwise just params.outDir/alignment/<sample>/<sample>.*
// BAM file is only published with --bamOut
publishDir { outDir }, pattern: "${sample}", mode: 'copy', saveAs: { "${outName}.star.solo" }
publishDir { outDir / "${outName}.star.align" }, pattern: "*bam", mode: 'copy', saveAs: { "${outName}.bam" }, enabled: params.bamOut
publishDir { outDir / "${outName}.star.align" }, pattern: "Log*", mode: 'copy'
tag "$sample $count"
script:
	outDir = file(params.outDir) / "alignment" / sample
	outName = sample
	if (params.splitFastq) {
		outDir = outDir / "split"
		outName = "${sample}.${count}"
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
	--readFilesIn <(cat transcript*.fastq.gz) <(cat barcode*.fastq.gz) --readFilesCommand zcat
	mv Solo.out/ $sample
"""
}

// Concatinate raw STAR outputs for multiple sub-sample into a single merged raw STAR output
// Run only when splitFastq is set to true
process merge_star {
input:
    tuple(val(sample), path("Solo.out*"))
output:
    tuple(val(sample), path(outDir), emit: merge)
label "report"
tag "$sample"
publishDir file(params.outDir) / "alignment", mode:'copy'
script:
	outDir = "${sample}/${sample}.star.solo"
	matrixFn = Utils.starMatrixFn(params.starMulti)
"""
    mergeRawSTARoutput.py --star_dirs Solo.out* --star_feature $params.starFeature --star_matrix $matrixFn --out_dir $outDir
"""
}


workflow alignment {
take:
	fastqs		// (sampleId, subsample, transcript, barcode)
	libStructure // parsed library structure information
	genome // Reference genome
main:
	// Sort the input fq files from inputReads to not break caching
	fqsSorted = fastqs.toSortedList().flatMap()
	fqsSorted.dump(tag:'fqsSorted')

	if (params.splitFastq) {
		// Group all fastq files for one subsample
		// alignFqGroupedBySampleAndWell -> (sampleID, subsample, [transcript], [barcode])
		alignFqGroupedBySampleAndWell = fqsSorted.groupTuple(by:[0,1])

		// Batch subsamples into groups for alignment
		// alignFqJobGroups -> (sampleID, [[transcript]], [[barcode]])
		alignFqJobGroups = alignFqGroupedBySampleAndWell.map {
			def sample = it[0]
			def transcript = it[2]
			def barcode = it[3]
			tuple(sample, transcript, barcode)
		}.groupTuple(size:params.starGroupSize, remainder:true)
		
		// Flatten list of list of transcripts and barcodes
		// alignFqGroups -> (sampleID, [transcript], [barcode])
		alignFqJobGroups = alignFqJobGroups.map {
			def sample = it[0]
			def transcript = it[1].flatten()
			def barcode = it[2].flatten()
			tuple(sample, transcript, barcode)
		}
	} else {
		// Group by sample name (sampleID, [transcript], [barcode])
		alignFqJobGroups = fqsSorted.groupTuple().map { [it[0], it[2], it[3]]}
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
	
	// If we split work, merge STAR results for subsamples together again
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

emit:
    soloOut = solo_out
}