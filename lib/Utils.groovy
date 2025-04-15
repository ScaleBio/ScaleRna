class Utils {
    // Name of the STARSolo matrix (.mtx) file depends on the selected multimapper parameter
    public static String starMatrixFn(multimappers) {
        if (multimappers == "Unique") {
            return "matrix.mtx"
        } else {
            return "UniqueAndMult-${multimappers}.mtx"
        }
    }

    // Derive value for the 'barcodes' field in the merged sample report from 
    // a list of the 'barcodes' column value for each sub-sample
    // Usually all samples have the same barcode entry (multi-plate use-case)
    public static String combineSampleBarcodes(barcodes) {
        ArrayList bcs = barcodes.unique()
        if (bcs.size() == 1) {
            return bcs[0]
        } else {
            return "combined"
        }
    }

    // Load .json file into object
    // (reference genomes, library design, etc.)
    public static Map loadJson(json) {
        def jsonSlurper  = new groovy.json.JsonSlurper()
        return jsonSlurper.parse( json )
    }

    // Return 0 for null or non-integer strings
    public static int toIntOr0(str) {
        if ((str != null) && str.isInteger()) {
            return str as int
        }
        else {
            return 0
        }
    }
}