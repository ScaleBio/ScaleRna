/*
* Perform downstream analysis of ScaleBio scRNA-seq data
* It performs seurat-based clustering, azimuth-based cell type annotation, and optionally makes an anndata object
*
* Processes:
*     SeuratClustering
*     AzimuthMapping
*     RnaDashboard
*     CombineResults
*     MakeAnnData
*/

process SeuratClustering {
	tag "$sample"
	label "seuratv5"
	// If comparison is true "_allSamples" is added as a suffix to the comparison group Id. This value is then used as the output directory.
	// If comparison is false the sample Id is used as the output directory.
	publishDir { file(params.outputDir) / cellTypingDir }, pattern: "*_SeuratObject.rds", mode: 'copy'
	publishDir { file(params.outputDir) / cellTypingDir }, pattern: "*_seurat_clustering_results.csv", mode: 'copy'
	publishDir { file(params.outputDir) / cellTypingDir }, pattern: "*_seurat_cluster_markers.csv", mode: 'copy'
	publishDir { file(params.outputDir) / cellTypingDir }, pattern: "*_bpcells", mode: 'copy', type: 'dir'

	input:
	tuple(val(sample), path(soloOut), path(cellReadMetrics))
	val(comparison)
	
	output:
	tuple(val(sample), path("*_seurat_clustering_results.csv"), emit: results)
	tuple(val(sample), path("*_SeuratObject.rds"), emit: object)
	tuple(val(sample), path("*_seurat_cluster_markers.csv"))
	tuple(val(sample), path("*_bpcells", type: 'dir'))

	script:
	if(comparison){
		wf = "comparison"
		cellTypingDir = "cellTyping/${sample}_allSamples"
	} else {
		wf = "standard"
		cellTypingDir = "cellTyping/$sample"
	}
	"""
	seuratClusteringV5.r --project $sample $wf \
	--matrixDir $soloOut --allCells $cellReadMetrics
	"""
}

process AzimuthMapping {
	tag "$sample"
	label "seuratv5"
	// If comparison is true "_allSamples" is added as a suffix to the comparison group Id. This value is then used as the output directory.
	// If comparison is false the sample Id is used as the output directory.
	publishDir { file(params.outputDir) / azimuthDir }, pattern: "*_AzimuthObject.rds", mode: 'copy'
	publishDir { file(params.outputDir) / azimuthDir }, pattern: "*_azimuth_mapping_results.csv", mode: 'copy'
	publishDir { file(params.outputDir) / azimuthDir }, pattern: "*_bpcells", mode: 'copy', type: 'dir'

	input:
	tuple(val(sample), path(soloOut), path(cellReadMetrics))
	val(comparison)
	
	output:
	tuple(val(sample), path("*_azimuth_mapping_results.csv"), emit: results)
	tuple(val(sample), path("*_AzimuthObject.rds"), emit: object)
	tuple(val(sample), path("*_bpcells", type: 'dir'))

	script:
	if(comparison){
		wf = "comparison"
		azimuthDir = "cellTyping/${sample}_allSamples"
	} else {
		wf = "standard"
		azimuthDir = "cellTyping/$sample"
	}
	"""
	azimuthMappingV5.r --reference ${params.azimuthRef} --project $sample \
	$wf --matrixDir $soloOut --allCells $cellReadMetrics
	"""
}

process RnaDashboard {
	tag "$sample"
	label "reporting"
	// If comparison is true "_allSamples" is added as a suffix to the comparison group Id. This value is then used as the output directory.
	// If comparison is false the sample Id is used as the output directory.
	publishDir { file(params.outputDir) / cellTypingDir }, pattern: "report.html", mode: 'copy', saveAs: {_filename -> "${sample}.html"}
	publishDir { file(params.outputDir) / cellTypingDir }, pattern: "*_summary.csv", mode: 'copy'

	input:
	tuple(val(sample), path("results.csv"), path("sampleStats.csv"), path("report.rmd"))
	path("lib.json")
	val(comparison)
	
	output:
	path("report.html")
	path("*_summary.csv")

	script:
	cellTypingDir = "cellTyping"
	if (comparison) {
		cellTypingDir += "/${sample}_allSamples"
	} else {
		cellTypingDir += "/${sample}"
	}
	"""
	renderRNADash.r --rmdPath report.rmd --results results.csv --sampleStats sampleStats.csv --comparison $comparison
	"""
}

process CombineResults {
	tag "$sample"
	label "reporting"
	// If comparison is true "_allSamples" is added as a suffix to the comparison group Id. This value is then used as the output directory.
	// If comparison is false the sample Id is used as the output directory.
	publishDir { file(params.outputDir) / "cellTyping" / (comparison ? "${sample}_allSamples" : sample) }, pattern: "*_combinedResults.csv", mode: 'copy', saveAs: {_filename -> "${sample}_cellTypingResults.csv"}

	input:
	tuple(val(sample), path(results))
	val(comparison)

	output:
	tuple(val(sample), path("*_combinedResults.csv"), emit: results)

	script:
	"""
	combineDownstreamResults.r --resultsDir ./ --sample $sample
	"""
}

process MakeAnnData {
	tag "$id"
	label "anndata"
	// If comparison is true search for file named "merged_anndata.h5ad", otherwise search for any file with the suffix _anndata.h5ad
	publishDir file(params.outputDir) / "samples", pattern: params.comparison ? "merged_anndata.h5ad" : "*_anndata.h5ad", mode: 'copy', saveAs: {filename -> "${id}_anndata.h5ad" }

	input:
	tuple(val(id), path("filtered_mtx_dir"), path("all_cells"), val(name))
	
	output:
	tuple(val(id), path("*_anndata.h5ad"), emit: annData)
	
	script:
	mtx_dirs = filtered_mtx_dir.join(" ")
	all_cells_files = all_cells.join(" ")
	sample_names = name.join(" ")
"""
	build_anndata.py --matrix_dir $mtx_dirs --all_cells $all_cells_files --sample_ids $sample_names
"""
}

workflow DOWNSTREAM {
take:
	samples
	allCells
	filteredMtx
	sampleStats
	libJson
	comparison
main:
	// allCells -> [sampleId, path(allCells.csv)]
	allCells = allCells.map{ tuple(it[0], it[2]) }
	resultList = allCells
	filteredMtxMeta = filteredMtx.join(allCells)

	// If samples are being "compared" we need to group samples by their compSamples value.
	if (comparison){
		// If compSamples is given in the samplesheet comparisonDict -> [sampleId, compSamples]
		// Else comparisonDict -> [sampleId, ScaleRNA]. This gives all samples a common value to group on.
		comparisonDict = samples.map{tuple(it.id, it.compSamples ? it.compSamples : "ScaleRNA")}
		// comparisonDict -> [sampleId, compSamples, path(filteredMtx), path(allCells.csv)]
		comparisonDict = comparisonDict.join(filteredMtxMeta)
		// comparisonDict -> [[sampleIds], compSamples, [path(filteredMtx)], [path(allCells.csv)]]
		comparisonDict = comparisonDict.groupTuple(by:1)
		comparisonDict.dump(tag:'comparisonDict')
		// filteredMtxMeta -> [ compSamples, [path(filteredMtx)], [path(allCells.csv)]]
		filteredMtxMeta = comparisonDict.map{tuple(it[1], it[2], it[3])}
		// resultList -> [compSamples, [path(allCells.csv)]]
		resultList = comparisonDict.map{tuple(it[1], it[3])}
		resultList.dump(tag:'resultList')
	}

	if (params.seurat){
		SeuratClustering(filteredMtxMeta, comparison)
		resultList = resultList.join(SeuratClustering.out.results)
		resultList.dump(tag:'resultListSeurat')
	}

	if (params.azimuth){
		AzimuthMapping(filteredMtxMeta, comparison)
		resultList = resultList.join(AzimuthMapping.out.results)
		resultList.dump(tag:'resultListAzimuth')
	}

	if (params.seurat || params.azimuth){
		combineResultsInput = resultList.map{tuple(it[0], it.tail().flatten())}
		combineResultsInput.dump(tag:'combineResultsInput')
		CombineResults(combineResultsInput, comparison)
		if (comparison){
			// When comparing samples there are multiple mad_threshold.csv files.
			// Currently the downstream report does not handle this well.
			// So we take the first file as "dummy" input and do not plot mad info in comparison reports.
			// mad_thresh = path(mad_threshold.csv)
			sample_stats = sampleStats.first().map{it[1]}
			// dashInput = [sampleId, path(sampleId_DownstreamResults.csv), path(sampleId_mad_threshold.csv)]
			dashInput = CombineResults.out.results.combine(sample_stats)
		} else {
			dashInput = CombineResults.out.results.join(sampleStats)
		}
		// dashInput -> [sampleId/compSamples, path(DownstreamResults.csv), path(mad_threshold_info.csv), path(report.rmd) ]
		dashInput = dashInput.combine(Channel.fromPath("$projectDir/bin/downstream_clustering_report.rmd", checkIfExists: true))
		dashInput.dump(tag:'dashInput')
		RnaDashboard(dashInput, libJson, comparison)
	}

	if(params.annData){

		if(params.comparison){
			comparisonDict
				.map{ tuple it[1], it[2], it[3], it[0]}
				.dump(tag:"makeAnnDataInput")
				.set{ makeAnnDataInput } // makeAnnDataInput -> [ compSamples, [path(filteredMtx)], [path(allCells.csv)], [sampleId] ]
			
		} else {
			filteredMtxMeta
				.map{ tuple it[0], it[1], it[2], [it[0]]}
				.dump(tag:"makeAnnDataInput")
				.set{ makeAnnDataInput } // makeAnnDataInput -> [ sampleId, path(filteredMtx), path(allCells.csv), sampleId ]
		}

		MakeAnnData(makeAnnDataInput)
	}
}
