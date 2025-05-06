/*
* Perform alignment of demultiplexed reads using STAR
*
* Processes:
*     StarSolo
*     MergeStar
*/

// Run the STAR executable on demuxed bam files
process StarSolo {
	tag "$sample ${task.index}"

	// Outputs published under params.outDir/alignment/<sample>/<sample>.*
	// BAM file is only published with --bamOut under params.outDir/alignment/aligned_bams if splitFastq true
	// Otherwise under params.outDir/alignment/<sample>/<sample>.star.align
	publishDir { file(params.outputDir) / starDir }, pattern: "Solo.out", mode: 'copy', saveAs: { "${outName}.star.solo" }, enabled: !params.splitFastq
	publishDir { file(params.outputDir) / bamDir }, pattern: "*bam", mode: 'copy', saveAs: { "${outName}.bam" }, enabled: params.bamOut
	publishDir { file(params.outputDir) / starDir / "${outName}.star.align" }, pattern: "Log*", mode: 'copy', enabled: !params.splitFastq

	input: 
	path(indexDir) // STAR genome index
	tuple(val(sample), path(demuxed_ubam, name:"demuxed_ubam*.bam")) // bcParser output bam file
	val(isBarnyard)

	output: 
	tuple(val(sample), path("Aligned.sortedByCoord.out.bam"), emit: bam, optional: true)
	tuple(val(sample), path("Log.final.out"), emit: log)
	tuple(val(sample), path("Solo.out"), emit: solo)

	script:
	starDir = "alignment/$sample"
	outName = sample
	bamDir = "alignment/${sample}/${sample}.star.align"
	if (params.splitFastq) {
		starDir += "/split"
		outName = "${sample}.${task.index}"
		bamDir += "/split_bams"
	}
	if (params.ultimaCramDir) {
		strand = "Reverse"
	} else {
		strand = params.starStrand
	}
	if (params.bamOut) {
		bamOpts = "--outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI nM AS GX GN gx gn sF --outSAMunmapped Within"
    } else {
        bamOpts = "--outSAMtype None" 
    }
	if(isBarnyard){
		starMultiParam = params.starMultiBarnyard
	} else {
		starMultiParam = params.starMulti
	}
	matrixFn = Utils.starMatrixFn(starMultiParam)
	ubam_files = demuxed_ubam.join(',')
	"""
	STAR --runThreadN $task.cpus --genomeDir $indexDir --soloType CB_UMI_Simple --soloBarcodeReadLength 0 --soloCBwhitelist None --soloCBtype String $bamOpts --outSJtype None\
	--soloCellReadStats Standard --soloStrand $strand ${params.starTrimming} --soloFeatures ${params.starFeature} --soloMultiMappers $starMultiParam --readFilesIn $ubam_files\
	--readFilesType SAM SE --readFilesCommand samtools view --soloInputSAMattrBarcodeSeq CB UM --soloCellFilter None --outFilterMultimapNmax ${params.starMaxLoci}
	pigz -p $task.cpus Solo.out/${params.starFeature}/raw/$matrixFn
	pigz -p $task.cpus Solo.out/${params.starFeature}/raw/barcodes.tsv
	pigz -p $task.cpus Solo.out/${params.starFeature}/raw/features.tsv
	"""
}

// Concatinate raw STAR outputs for multiple sub-sample into a single merged raw STAR output
// Run only when splitFastq is set to true
process MergeStar {
	tag "$sample"
	label 'large'

	publishDir file(params.outputDir) / "alignment", mode:'copy'
	
	input:
    tuple(val(sample), path("Solo.out*"), path("log.final.out*"))
	val(isBarnyard)

	output:
    tuple(val(sample), path(outDir), emit: merge)
	tuple(val(sample), path("$logOutDir/Log.final.out"), emit: merge_log)

	script:
	outDir = "${sample}/${sample}.star.solo"
	// To preserve original STAR output directory structure, we write the log file in the star.align directory
	logOutDir = "${sample}/${sample}.star.align"
	if(isBarnyard){
		starMultiParam = params.starMultiBarnyard
	} else {
		starMultiParam = params.starMulti
	}
	matrixFn = Utils.starMatrixFn(starMultiParam)
	"""
    merge_raw_star_output.py --star_dirs Solo.out* --star_log log.final.out* --star_feature $params.starFeature --star_matrix $matrixFn --out_dir $outDir --log_out_dir $logOutDir
	pigz -p $task.cpus $outDir/${params.starFeature}/raw/$matrixFn
	"""
}


workflow ALIGNMENT {
take:
	bams   // (sampleId, subsample, bam)
	genome // Reference genome
main:
	
	isBarnyard = genome.get('isBarnyard', false)

	// Ultima files will always be split by sample barcode
	if (params.splitFastq || params.ultimaCramDir) {
		// Group all bam files for one subsample
		// alignFqGroupedBySampleAndWell -> (sampleID, subsample, [bam])
		alignFqGroupedBySampleAndWell = bams.groupTuple(by:[0,1], sort: true)
		// Batch subsamples into groups for alignment
		// alignFqJobGroups -> (sampleID, [bam])
		alignFqGroupedBySampleAndWell
			// Sort channel so results are cached
			.toSortedList { a, b ->
				// Sort by id, subsample and then bam file names
				def result = a[0] <=> b[0]
				if (result == 0) {
					result = a[1] <=> b[1]
					if (result == 0) {
						result = a[2].join(",") <=> b[2].join(",")
					}
				} 
				return result
			}
			.flatMap()
			// Drop subsample so that we can sort on bam files in subsequent groupTuple
			.map{ id, _subsample, bam ->
				[id, bam]
			}
			.groupTuple(size: params.rtBarcodesPerStarJob, sort: { a, b -> a.join(",") <=> b.join(",") }, remainder: true)
			.map { id, bam ->
				[id, bam.flatten()]
			}
			.set{ alignFqJobGroups }
	} else {
		// Group by sample name (sampleID, [bam])
		bams
			.groupTuple()
			.map { [it[0], it[2].flatten().sort()] }
			.set{ alignFqJobGroups }
	}
	alignFqJobGroups.dump(tag:'alignFqJobGroups')
	
	StarSolo(genome.star_index, alignFqJobGroups, isBarnyard)
	StarSolo.out.solo.dump(tag:'solo')
	
	// If we split work, merge STAR results for subsamples together again
	// Ultima files will always be split by sample barcode
	if (params.splitFastq || params.ultimaCramDir) {
		// Group on sample name to ensure star outputs of the same sample are grouped together
		star_solo_by_sample = StarSolo.out.solo.join(StarSolo.out.log).groupTuple()
		star_solo_by_sample.dump(tag: 'star_solo_by_sample')
		
		// Merge star outputs
		MergeStar(star_solo_by_sample, isBarnyard)
		solo_out = MergeStar.out.merge
		solo_log = MergeStar.out.merge_log
	} else {
		solo_out = StarSolo.out.solo
		solo_log = StarSolo.out.log
	}
	solo_out.dump(tag:'soloOut')

emit:
    soloOut = solo_out
	soloLog = solo_log
}
