/// This module perform downstream analysis of scRNA-seq data. 
/// It performs seurat-based clustering, azimuth-based cell type annotation, and optionally scrublet doublet detection.
/// The reults from these tools are visualized with an html report.

nextflow.enable.dsl=2

process seuratClustering {
input:
	tuple(val(sample), path(soloOut), path(cellReadMetrics))
output:
	tuple(val(sample), path("*_seurat_clustering_results.csv"), emit: results)
	tuple(val(sample), path("*_SeuratObject.rds"), emit: object)
	tuple(val(sample), path("*_seurat_cluster_markers.csv"))
	tuple(val(sample), path("*_bpcells", type: 'dir'))
tag "$sample"
label "seuratv5"
// If comparison is true "_allSamples" is added as a suffix to the comparison group Id. This value is then used as the output directory.
// If comparison is false the sample Id is used as the output directory.
publishDir { out_dir }, pattern: "*_SeuratObject.rds", mode: 'copy'
publishDir { out_dir }, pattern: "*_seurat_clustering_results.csv", mode: 'copy'
publishDir { out_dir }, pattern: "*_seurat_cluster_markers.csv", mode: 'copy'
publishDir { out_dir }, pattern: "*_bpcells", mode: 'copy', type: 'dir'
script:
	out_dir = file(params.outDir) / "cellTyping"
	if (params.comparison) {
		workflow = "comparison"
		out_dir = out_dir / "${sample}_allSamples"
	} else {
		workflow = "standard"
		out_dir = out_dir / sample
	}
"""
	seuratClusteringV5.r --project $sample $workflow \
	--matrixDir $soloOut --allCells $cellReadMetrics
"""
}

process azimuthMapping {
input:
	tuple(val(sample), path(soloOut), path(cellReadMetrics))
output:
	tuple(val(sample), path("*_azimuth_mapping_results.csv"), emit: results)
	tuple(val(sample), path("*_AzimuthObject.rds"), emit: object)
	tuple(val(sample), path("*_bpcells", type: 'dir'))
tag "$sample"
label "seuratv5"
// If comparison is true "_allSamples" is added as a suffix to the comparison group Id. This value is then used as the output directory.
// If comparison is false the sample Id is used as the output directory.
publishDir { out_dir }, pattern: "*_AzimuthObject.rds", mode: 'copy'
publishDir { out_dir }, pattern: "*_azimuth_mapping_results.csv", mode: 'copy'
publishDir { out_dir }, pattern: "*_bpcells", mode: 'copy', type: 'dir'
script:
	out_dir = file(params.outDir) / "cellTyping"
	if (params.comparison) {
		workflow = "comparison"
		out_dir = out_dir / "${sample}_allSamples"
	} else {
		workflow = "standard"
		out_dir = out_dir / sample
	}
"""
	azimuthMappingV5.r --reference ${params.azimuthRef} --project $sample \
	$workflow --matrixDir $soloOut --allCells $cellReadMetrics
"""
}

process rnaDashboard {
input:
	tuple(val(sample), path("results.csv"), path("sampleStats.csv"), path("report.rmd"))
output:
	path("report.html")
	path("*_summary.csv")
tag "$sample"
label "reporting"
// If comparison is true "_allSamples" is added as a suffix to the comparison group Id. This value is then used as the output directory.
// If comparison is false the sample Id is used as the output directory.
publishDir { out_dir }, pattern: "report.html", mode: 'copy', saveAs: {filename -> "${sample}.html"}
publishDir { out_dir }, pattern: "*_summary.csv", mode: 'copy'
script:
	out_dir = file(params.outDir) / "cellTyping"
	if (params.comparison) {
		out_dir = out_dir / "${sample}_allSamples"
	} else {
		out_dir = out_dir / sample
	}
"""
	renderRNADash.r --rmdPath report.rmd --results results.csv --sampleStats sampleStats.csv --comparison $params.comparison
"""
}

process combineResults {
input:
	tuple(val(sample), path(results))
output:
	tuple(val(sample), path("*_combinedResults.csv"), emit: results)
tag "$sample"
label "reporting"
// If comparison is true "_allSamples" is added as a suffix to the comparison group Id. This value is then used as the output directory.
// If comparison is false the sample Id is used as the output directory.
publishDir { file(params.outDir) / "cellTyping" / (params.comparison ? "${sample}_allSamples" : sample) }, pattern: "*_combinedResults.csv", mode: 'copy', saveAs: {filename -> "${sample}_cellTypingResults.csv"}
script:
"""
	combineDownstreamResults.r --resultsDir ./ --sample $sample
"""
}

process makeAnnData {
input:
	tuple(val(sample), path(soloOut), path(cellReadMetrics))
output:
	tuple(val(sample), path("*_anndata.h5ad"), emit: annData)
tag "$sample"
label "seuratv5"
publishDir file(params.outDir) / "samples", pattern: "*_anndata.h5ad", mode: 'copy', saveAs: {filename -> params.comparison ? "${sample}_allSamples.h5ad" : "${sample}.h5ad"}
script:
mtxDir = soloOut.join(" ")
"""
	downstream_make_anndata.r --matrixDir $mtxDir --sampleId $sample
"""
}

workflow downstream {
take:
	samples
	allCells
	filteredMtx
	rawMtx
	sampleStats
main:
	// allCells -> [sampleId, path(allCells.csv)]
	allCells = allCells.map{ tuple(it[0], it[2]) }
	resultList = allCells
	filteredMtxMeta = filteredMtx.join(allCells)
	rawMtxMeta = rawMtx.join(allCells)

	// If samples are being "compared" we need to group samples by their compSamples value.
	if(params.comparison){
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

	if(params.seurat){
		seuratClustering(filteredMtxMeta)
		resultList = resultList.join(seuratClustering.out.results)
		resultList.dump(tag:'resultListSeurat')
	}

	if(params.azimuth){
		azimuthMapping(filteredMtxMeta)
		resultList = resultList.join(azimuthMapping.out.results)
		resultList.dump(tag:'resultListAzimuth')
	}

	if(params.seurat || params.azimuth){
		combineResultsInput = resultList.map{tuple(it[0], it.tail().flatten())}
		combineResultsInput.dump(tag:'combineResultsInput')
		combineResults(combineResultsInput)
		if(params.comparison){
			// When comparing samples there are multiple mad_threshold.csv files.
			// Currently the downstream report does not handle this well.
			// So we take the first file as "dummy" input and do not plot mad info in comparison reports.
			// mad_thresh = path(mad_threshold.csv)
			sample_stats = sampleStats.first().map{it[1]}
			// dashInput = [sampleId, path(sampleId_DownstreamResults.csv), path(sampleId_mad_threshold.csv)]
			dashInput = combineResults.out.results.combine(sample_stats)
		} else {
			dashInput = combineResults.out.results.join(sampleStats)
		}
		// dashInput -> [sampleId/compSamples, path(DownstreamResults.csv), path(mad_threshold_info.csv), path(report.rmd) ]
		dashInput = dashInput.combine(Channel.fromPath("$projectDir/bin/downstream_clustering_report.rmd", checkIfExists: true))
		dashInput.dump(tag:'dashInput')
		rnaDashboard(dashInput)
	}

	if(params.annData){
		makeAnnData(filteredMtxMeta)
	}
}
