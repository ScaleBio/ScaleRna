def getTotalSampleReads(sample) {
	def reportStatistics = file(sample.resultDir) / "reports" / "reads_per_sample.csv"
	if (!reportStatistics.exists()) {
		return [[sample.id, 0]]
	} else {
		def lines = reportStatistics.text.readLines()

		def totalSampleReads = []
		def header = lines[0].split(',').toList()
		def sampleIndex = header.indexOf('Sample')
		def libraryIndex = header.indexOf('Library')
		def totalReadsIndex = header.indexOf('TotalReads')
		lines[1..-1].each { line ->
			def columns = line.split(',')
			def sampleName = columns[sampleIndex].trim()
			def library = columns[libraryIndex].trim()
			def totalReads = columns[totalReadsIndex].trim()
			def libName = "${sampleName}.${library}"
			if (libName == sample.id) {
				totalSampleReads << ["${sampleName}.${library}", totalReads]
			}
		}
		return totalSampleReads
	}
}
def shouldRunLibDetection(samples, quantum) {
	def header = samples.withReader { reader ->
		reader.readLine()
	}
	if (!header.contains("libName") && quantum && params.fastqDir) {
		return true
	}
	return false
}

// Create a 'file()' from string 'path'
// 'path' could be null, a URL, an absolute file-path or relative to 'baseDir'
def expandPath(path, baseDir) {
	if (path == null) { return null}
	def uri = new URI(path)
	if (uri.isAbsolute()) {
		return file(path)
	} else {
		return baseDir.resolve(path)
	}
}

// Reference genome files and parameters
// Paths can be absolute or relative to location of json
def loadGenome(json) {
	if (json.isDirectory()) {
		ParamLogger.throwError("--genome should be json file, not a directory: ${json}")
	}
	def baseDir = json.getParent()
	def genome = Utils.loadJson(json)
	genome.star_index = expandPath(genome.star_index, baseDir)
	if (genome.star_index == null || !genome.star_index.exists()) {
		ParamLogger.throwError("Genome star_index not found: ${genome.star_index?.toUriString() ?: ''}")
	}
	genome.gtf = expandPath(genome.gtf, baseDir)
	return genome
}

def resolveScalePlexToRnaMapping(scalePlexLibJson) {
	if (params.scalePlexToRnaMapping) {
		scalePlexToRnaMapping = file(params.scalePlexToRnaMapping, checkIfExists:true)
	} else {
		scalePlexLibStructure = Utils.loadJson(scalePlexLibJson)
		if (scalePlexLibStructure.get("scaleplex_to_rna_mapping")) {
			scalePlexToRnaMapping = file(scalePlexLibJson.getParent()) / scalePlexLibStructure["scaleplex_to_rna_mapping"]
		} else {
			// path input cannot be null, no mapping file used in 3lvl pipeline
			scalePlexToRnaMapping = []
		}
	}
	return scalePlexToRnaMapping
}