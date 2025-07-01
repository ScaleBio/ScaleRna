include { LibraryDetection } from './library_detection.nf'

workflow ULTIMA {
take:
    samples // each row from samples.csv parsed into a channel
    barcodeUnwrappedSamplesCsv // samples csv with barcode ranges unwrapped
    libraryInfo
main:
    Channel.fromPath(
        "$params.ultimaCramDir/**cram", checkIfExists: true)
        .dump(tag:'ultimaCrams')
        .set { ultimaCrams }
    ultimaCrams
        .buffer(size: params.filesPerLibDetection, remainder: true)
        .set { ultimaCramsLibDetectionInput }
    LibraryDetection(
        ultimaCramsLibDetectionInput,
        libraryInfo.rnaLibraryStructureFile.getParent(),
        libraryInfo.rnaLibraryStructureFile.getName(),
        libraryInfo.scalePlexLibraryStructureFile.getParent(),
        libraryInfo.scalePlexLibraryStructureFile.getName()
    )
    LibraryDetection.out.pool_mapping
        .collectFile(
            name: 'cram_file_to_library_assignment.csv',
            keepHeader: true,
            storeDir: file(params.outputDir) / "cram"
        )
        .splitCsv(header:true, strip:true)
        .filter { row -> row.pcr_pool != "Unknown" && row.pcr_pool != "Ambiguous" && row.rt != "Unknown" && row.rt != "Ambiguous" }
        .dump(tag:'poolMapping')
        .map { [it.file, it.pcr_pool, it.rt] }
        .set { poolMapping }
    poolMapping
        .map { it[1] } // pcr
        .unique()
        .collect()
        .ifEmpty { ParamLogger.throwError("No valid PCR pool found in unaligned cram files") }
        .dump(tag:'validPools')
        .set { validPools }
    samples
        // Filter out samples that have a library which was not found during library detection
        .filter { sample ->
            validPools.val.contains(sample.libIndex2)
        }
        .dump(tag:'samplesAfterLibraryDetection')
        .set { samples }
    barcodeUnwrappedSamplesCsv
        // Filter out samples that have a library which was not found during library detection
        .filter { sample ->
            validPools.val.contains(sample.libName)
        }
        .dump(tag:'barcodeUnwrappedSamplesAfterLibraryDetection')
        .map { [it.libName, it.id, it.barcode] }
        .set { barcodeUnwrappedSamplesCsv }
    // Get all rna sample names after filtering
    samples
        .filter { it.rnaId == null }
        .map { [it.id, it.sample] }
        .dump(tag:'rnaSampleNames')
        .set { rnaSampleNames }
    // Get all hash sample names after filtering
    samples
        .filter { it.rnaId != null }
        .map { [it.id, it.libName, it.sample] }
        .ifEmpty { [] }
        .dump(tag:'hashSampleNames')
        .set { hashSampleNames }
    ultimaCrams
        .map { tuple(it.name, it) }
        .join(poolMapping) // joining key is cram file name
        .map { _filename, cram, pcr_detected, rt_well ->
            tuple(pcr_detected, cram, rt_well)
        }
        .combine(barcodeUnwrappedSamplesCsv, by:[0, 2]) // joining key is rt barcode and pcr
        .map { _libName, rt_well, cram, sample_id ->
            tuple(sample_id, rt_well, cram)
        }
        .dump(tag:'ultimaSamplesWithSampleID')
        .set { ultimaSamplesWithSampleID }
    // Combine because it is a one to many mapping
    // One is rnaSampleNames, and many is ultimaSamplesWithSampleID
    // Matching key is sample_id
    ultimaSamplesWithSampleID
        .combine(rnaSampleNames, by:0)
        .map { sample_id, rt_well, cram, sample -> 
            def subsample = sample + "_" + rt_well
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
    samples = samples
}
