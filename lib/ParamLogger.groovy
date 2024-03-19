import nextflow.Nextflow

// Adapted from nf-core/rnaseq
class ParamLogger {
    public static void initialise(workflow, params, log) {
        if (params.help) {
            log.info("\nScaleBio Seq Suite: RNA Workflow\nFor usage information, please see README.md\n")
            System.exit(0)
        }
        log.info paramsSummaryLog(workflow, params)
        validateParams(params, log)
    }

    public static String throwError(errMessage) {
        Map colors = logColours()
        Nextflow.error("${colors.red}ERROR${colors.reset}: $errMessage")
    }
    
    // Takes a camelCase string and converts it to kebab-case
    // Necessary for validating parameters since nextflow converts camelCase to kebab-case
    public static String camelToKebab(String camelCase) {
        return camelCase.replaceAll(/([a-z])([A-Z])/, '$1-$2').toLowerCase()
    }

    // Check required workflow inputs
    public static validateParams(params, log) {
        // PLEASE ADD NEW PARAMETERS TO THE allowed_parameters LIST
        // useSTARthreshold is a special case because how nextflow resolves it to kebab case seems incorrect
        // use-starthreshold is the correct kebab case
        def allowed_parameters = ['samples', 'genome', 'runFolder', 'fastqDir', 'fastqSamplesheet', 'reporting',
                                  'resultDir', 'outDir', 'bamOut', 'fastqOut', 'libStructure', 'merge', 'splitFastq',
                                  'bclConvertParams', 'fastqc', 'task_max_memory', 'starFeature', 'starMulti', 'starStrand',
                                  'trimFastq', 'trimAdapt', 'starTrimming', 'minUTC', 'min-UTC', 'cellFinder', 'fixedCells', 'UTC',
                                  'task_max_memory', 'task_max_cpus', 'task_max_time', 'starGroupSize', 'bcParserJobs', 'seurat', 'azimuth',
                                  'compSamples', 'internalReport', 'help',  'topCellPercent', 'minCellRatio', 'expectedCells',
                                  'useSTARthreshold', 'use-STARthreshold', 'cellFinderFdr', 'filter_outliers', 'num_mad_genes',
                                  'num_mad_umis', 'num_mad_mito', 'azimuthRef', 'cellTyping', 'seuratWorkflow', 'annData']
        def master_list_of_params = allowed_parameters
        allowed_parameters.each { str ->
            master_list_of_params += camelToKebab(str)}
        def parameter_diff = params.keySet() - master_list_of_params
        if (parameter_diff.size() == 1){
            log.warn("[Argument Error] Parameter $parameter_diff is not valid in the pipeline!")
        }
        if (parameter_diff.size() > 1){
            log.warn("[Argument Error] Parameters $parameter_diff are not valid in the pipeline!")
        }
        if (params.samples == null || params.samples == true) {
            throwError("Must specify --samples (e.g. samples.csv)")
        }
        if (params.libStructure == null || params.libStructure == true) {
            throwError("Must specify --libStructure (e.g. libV1.1.json)")
        }
        if (params.genome == null || params.genome == true) {
            throwError("Must specify --genome")
        }
        if (params.reporting) {
            if (params.runFolder || params.fastqDir) {
                throwError("Cannot specify --runFolder or --fastqDir when running reporting-only (--reporting)")
            }
        }
    }
    
    static LinkedHashMap paramsSummaryMap(workflow, params) {
        def Map nextflow_opts = [:] // Core Nextflow options
        def Map workflow_opts = [:] // ScaleBio Workflow parameters
        def Map exec_opts = [:] // ScaleBio Workflow execution options
        def Map input_opts = [:] // Workflow inputs

        nextflow_opts['Workflow Directory:']   = workflow.projectDir
        nextflow_opts['Workflow Version:'] = workflow.manifest.version
        if (workflow.revision) { nextflow_opts['Workflow Revision'] = workflow.revision }
        nextflow_opts['Command Line'] = workflow.commandLine
        nextflow_opts['Nextflow RunName'] = workflow.runName
        nextflow_opts['Profile'] = workflow.profile
        nextflow_opts['Config Files'] = workflow.configFiles.join(', ')
        if (workflow.containerEngine) {
            nextflow_opts['Container Engine'] = workflow.containerEngine
        }
        nextflow_opts['Launch Directory'] = workflow.launchDir
        nextflow_opts['Work Directory'] = workflow.workDir

        exec_opts['Workflow Output'] = params.outDir
        if (params.reporting) { exec_opts['reporting'] = params.reporting }
        exec_opts['merge'] = params.merge
        exec_opts['splitFastq'] = params.splitFastq
        exec_opts['task_max_memory'] = params.task_max_memory
        exec_opts['task_max_cpus'] = params.task_max_cpus

        if (params.fastqDir) {
            input_opts['fastqDir'] = params.fastqDir
        } else if (params.runFolder) {
            input_opts['runFolder'] = params.runFolder
        }
        if (params.reporting) { input_opts['resultDir'] = params.resultDir }
        input_opts['samples'] = params.samples
        input_opts['genome'] = params.genome
        input_opts['libStructure'] = params.libStructure

        workflow_opts['trimFastq'] = params.trimFastq
        workflow_opts['starFeature'] = params.starFeature
        workflow_opts['starMulti'] = params.starMulti
        if (params.fixedCells != 0) {
            workflow_opts['fixedCells'] = params.fixedCells
        } else if (params.UTC != 0) {
            workflow_opts['UTC'] = params.UTC
        } else {
            workflow_opts['minUTC'] = params.minUTC
            workflow_opts['cellFinder'] = params.cellFinder
        }

        if (params.seurat) {
            workflow_opts['seurat'] = params.seurat
            workflow_opts['azimuth'] = params.azimuth
            workflow_opts['compSamples'] = params.compSamples
        }

        return [ 'Core Nextflow Options': nextflow_opts,
            'Inputs': input_opts,
            'Workflow Execution Options': exec_opts,
            'Analysis Parameters': workflow_opts
        ]
    }

    // Beautify parameters for summary and return as string
    public static String paramsSummaryLog(workflow, params) {
        Map colors = logColours()
        String output  = ''
        def params_map = paramsSummaryMap(workflow, params)
        def max_chars  = paramsMaxChars(params_map)
        for (group in params_map.keySet()) {
            def group_params = params_map.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                output += colors.bold + group + colors.reset + '\n'
                for (param in group_params.keySet()) {
                    output += "  " + colors.blue + param.padRight(max_chars) + ": " + colors.green +  group_params.get(param) + colors.reset + '\n'
                }
                output += '\n'
            }
        }
        output += dashedLine()
        return output
    }


    public static Map logColours(monochrome_logs=false) {
        Map colorcodes = [:]
        // Reset / Meta
        colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
        colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
        colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
        colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
        colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
        colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
        colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

        // Regular Colors
        colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
        colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
        colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
        colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
        colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
        colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
        colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
        colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

        // Bold
        colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
        colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
        colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
        colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
        colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
        colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
        colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
        colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

        // Underline
        colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
        colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
        colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
        colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
        colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
        colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
        colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
        colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

        // High Intensity
        colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
        colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
        colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
        colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
        colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
        colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
        colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
        colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

        // Bold High Intensity
        colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
        colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
        colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
        colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
        colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
        colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
        colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
        colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

        return colorcodes
    }

    //
    // Creates a dashed line 
    //
    public static String dashedLine() {
        Map colors = logColours()
        return "-${colors.dim}----------------------------------------------------${colors.reset}-"
    }
    private static Integer paramsMaxChars(params_map) {
        Integer max_chars = 0
        for (group in params_map.keySet()) {
            def group_params = params_map.get(group)  // This gets the parameters of that particular group
            for (param in group_params.keySet()) {
                if (param.size() > max_chars) {
                    max_chars = param.size()
                }
            }
        }
        return max_chars
    }
}
