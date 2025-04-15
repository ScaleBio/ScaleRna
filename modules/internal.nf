process MultiLibraryReport {
    publishDir file(params.outputDir) / "reports" / "library" / "csv", mode: 'copy'

    input:
	path(typeLevelMatches)
	path(overallMatches)
    
    output:
	path("allLibraries.barcodeStats.csv")

    script:
    """
	combine_library_metrics.py --typeLevelMatches ${typeLevelMatches} --overallMatches ${overallMatches}
    """
}