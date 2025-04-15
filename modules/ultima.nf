// Extract libName and PBC name/alias for one ultima cram file
def getPBCandLibNameFromUltimaCram (cram) {
    def tok = cram.getName().toString()
    // Remove file extension from filename
    def fname = tok[0..tok.lastIndexOf(".") - 1]

    // Compile the regex pattern
    def matcher = fname =~ params.cramFilePattern
    def libNames = []
    def pbc_well = null

    // Find matches and extract groups
    if (matcher.find()) {
        (1..matcher.groupCount()).each { i ->
            def groupValue = matcher.group(i)
            if (matcher.group("PBC") == groupValue) {
                pbc_well = groupValue
            } else if (groupValue) { // Everything other than PBC is libName
                libNames << groupValue
            }
        }
    }
    else {
        log.error("Pattern: $params.cramFilePattern does not match pattern of filename: $fname")
    }

    // Join multiple libName parts with "-"
    def libName = libNames.join("-")

    return tuple(cram, libName, pbc_well)
}
workflow ULTIMA {
take:
    samples // each row from samples.csv parsed into a channel
    barcodeUnwrappedSamplesCsv // samples csv with barcode ranges unwrapped
main:
    samples
        .filter { it.rnaId == null }
        .map { [it.id, it.sample] }
        .dump(tag:'rnaSampleNames')
        .set { rnaSampleNames }
    samples
        .filter { it.rnaId != null }
        .map { [it.id, it.libName, it.sample] }
        .ifEmpty { [] }
        .dump(tag:'hashSampleNames')
        .set { hashSampleNames }
    ultimaCrams = Channel.fromPath("$params.ultimaCramDir/**cram", checkIfExists: true).dump(tag:'ultimaCrams')
    // ultimaCrams -> [cram, libName, pbc_well]
    ultimaCrams
        .map { getPBCandLibNameFromUltimaCram(it) }
        .dump(tag:'ultimaSamples')
        .set { ultimaSamples }
    // Join on libName and sample barcode well coordinate
    barcodeUnwrappedSamplesCsv
        .combine(ultimaSamples, by:[1, 2])
        .map { _libName, pbc_well, sample_id, cram ->
            tuple(sample_id, pbc_well, cram)
        }
        .dump(tag:'ultimaSamplesWithSampleID')
        .set { ultimaSamplesWithSampleID }
    // Combine because it is a one to many mapping
    // One is rnaSampleNames, and many is ultimaSamplesWithSampleID
    // Matching key is sample_id
    ultimaSamplesWithSampleID
        .combine(rnaSampleNames, by:0)
        .map { sample_id, pbc_well, cram, sample -> 
            def subsample = sample + "_" + pbc_well
            tuple(sample_id, subsample, cram)
        }
        .dump(tag:'alignmentInput')
        .set { alignmentInput }
    // Combine because it is a one to many mapping
    // One is hashSampleNames, and many is ultimaSamplesWithSampleID
    // Matching key is sample_id
    ultimaSamplesWithSampleID
        .combine(hashSampleNames, by:0)
        .map {
            def meta = [lib: it[3], sample: it[4], subsample: it[4] + "_" + it[1]]
            [meta, it[2]]
        }
        // group on meta map so all wafer cram files from the same library and pool barcode are counted together 
        .groupTuple()
        .dump(tag:'ubamHash')
        .set {ubamHash}
    
emit:
    alignmentInput = alignmentInput
    ubamHash = ubamHash
}
