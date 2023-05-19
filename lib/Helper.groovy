// Adapted from nf-core/rnaseq
class Helper {
    public static void initialise(workflow, params, log) {
        // Write logo
        if (params.help) {
            print help(workflow, params, log)
            System.exit(0)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)
    }
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += paramsSummaryLogChanged(workflow, params)
        summary_log += dashedLine()
        return summary_log
    }
    public static LinkedHashMap paramsSummaryMap(workflow, params, schema_filepath="${workflow.projectDir}/nextflow_schema.json") {
        // Get a selection of core Nextflow workflow options
        def Map workflow_summary = [:]
        if (workflow.revision) {
            workflow_summary['revision'] = workflow.revision
        }
        workflow_summary['runName']      = workflow.runName
        if (workflow.containerEngine) {
            workflow_summary['containerEngine'] = workflow.containerEngine
        }
        if (workflow.container) {
            workflow_summary['container'] = workflow.container
        }
        workflow_summary['workDir']      = workflow.workDir
        workflow_summary['projectDir']   = workflow.projectDir
        workflow_summary['userName']     = workflow.userName
        workflow_summary['profile']      = workflow.profile
        workflow_summary['configFiles']  = workflow.configFiles.join(', ')
        workflow_summary['outDir']       = params.outDir
        if (params.fastqDir) {
            workflow_summary['fastqDir'] = params.fastqDir
        }
        if (params.runFolder) {
            workflow_summary['runFolder'] = params.runFolder
        }
        workflow_summary['samples'] = params.samples
        workflow_summary['genome'] = params.genome
        workflow_summary['libStructure'] = params.libStructure
        workflow_summary['fastqc'] = params.fastqc
        workflow_summary['starFeature'] = params.starFeature
        workflow_summary['starMulti'] = params.starMulti
        workflow_summary['trimFastq'] = params.trimFastq
        workflow_summary['cellTyping'] = params.cellTyping
        workflow_summary['max_memory'] = params.max_memory
        workflow_summary['max_cpus'] = params.max_cpus
        return [ 'Core Nextflow options' : workflow_summary ]
    }

    //
    // Beautify parameters for summary and return as string
    //
    public static String paramsSummaryLogChanged(workflow, params) {
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
