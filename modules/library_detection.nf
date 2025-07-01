// Map input reads(fastq or unaligned cram) to PCR index
process LibraryDetection {
	tag "${files[0].getSimpleName()}"

	input:
	path(files)
	path(libStructDir)
	val(libStructName)
	path("scalePlexLibStructDir") 
	val(scalePlexLibStructName)

	output:
	path("pool_mapping.csv", emit: pool_mapping)

	script:
	opts = ""
	if (params.scalePlex) {
		opts = "--scalePlexLibraryStruct scalePlexLibStructDir/$scalePlexLibStructName "
	}
	if (params.ultimaCramDir) {
		opts += "--ultimaCramDir . "
	} else {
		opts += "--fastqDir . "
	}
	"""
	sniff_input_reads.py --libraryStruct $libStructDir/$libStructName $opts
	"""
}