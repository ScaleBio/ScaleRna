class Utils {
    // Name of the STARSolo matrix (.mtx) file depends on the selected multimapper parameter
    public static String starMatrixFn(multimappers) {
        if (multimappers == "Unique") {
            return "matrix.mtx"
        } else {
            return "UniqueAndMult-${multimappers}.mtx"
        }
    }
}