/**
* Create hash UMI matrices for each sample from the barcodes csv output by bcParser
* and use along with RNA passing cells to generate metrics and filtered UMI matrix
*
*/
include { SCALE_PLEX_REPORTING as REPORTING } from './scale_plex_metrics.nf'

process CountScalePlex {
    tag "$meta.subsample"

    input:
    path(libStructDir) // Directory containing the library structure definition
    val(libStructName) // Filename of the library structure definition .json
    val(file_ext)
    path(mapping_file) // File containing the mapping of scale plex to RNA PCR barcodes
    tuple val(meta), path("barcodes_*.${file_ext}")

    output:
    tuple val(meta), path("${id}.cellMetrics.parquet"), path("${id}.raw.umi.mtx.gz"), emit: countRaw

    script:
	libStruct = "${libStructDir}/${libStructName}"
    id = meta.subsample + "." + meta.lib
    opts = ""
    if (mapping_file) {
        opts = "--mapping_file $mapping_file "
    }
    """
    count_hashes.py --bamDir . --references $libStructDir --lib_struct $libStruct --id $id --min_read_umi_thresh $params.scalePlexMinReadUmi --outDir . $opts
    """
}

process MergeCountRaw {
    tag "$meta.sample"
	label 'large'

    input:
    tuple val(meta), path(cellMetrics), path(mtx)

    output:
    tuple val(meta), path("${id}.cellMetrics.parquet"), path("${id}.raw.umi.mtx.gz"), emit: countRaw

    script:
    id = meta.sample + "." + meta.lib
    """
    concat_columnar.py $cellMetrics --outputFile "${id}.cellMetrics.parquet"
    merge_mtx.py $mtx --id $id
    """
}


workflow SCALE_PLEX {
    take:
        libraryInfo
        samples
        ubam
        allCells
        libCount
        perSampleCountRaw // provided if --reporting run

    main:
        if (!params.reporting) {
            if (params.ultimaCramDir) {
                file_ext = 'cram'
            } else {
                file_ext = 'bam'
            }
            // CountScalePlex on output of bcParser, dropping where sample was not identified
            CountScalePlex(
                libraryInfo.scalePlexLibraryStructureFile.getParent(),
                libraryInfo.scalePlexLibraryStructureFile.getName(),
                file_ext,
                libraryInfo.scalePlexToRnaMappingFile,
                ubam.filter { meta, _ubam_file -> meta.sample != 'Unknown' }
            )
            CountScalePlex.out.countRaw
                .map { meta, cellMetrics, rawUmi ->
                    // drop subsample that was used for splitting
                    [[lib:meta.lib, sample:meta.sample], cellMetrics, rawUmi]
                }
                .set { perSampleCountRaw }
            if (params.splitFastq) {
                perSampleCountRaw
                    // collect all cellMetrics, barcodes, and umi.sparse.mtx for each lib and sample combination
                    .groupTuple()
                    | MergeCountRaw
                    | set { perSampleCountRaw }
            }
        }
        REPORTING(
            perSampleCountRaw,
            libraryInfo.scalePlexLibraryStructureFile,
            samples,
            allCells,
            libCount
        )
    
    emit:
        metrics = REPORTING.out.metrics
        allCells = REPORTING.out.allCells
}
