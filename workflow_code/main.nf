include { enabled } from './modules/utils.nf'
nextflow.enable.dsl = 2

include { validateParameters; paramsHelp } from 'plugin/nf-schema'
include { AMPLICON } from './workflows/amplicon.nf'
include { POST_PROCESSING } from './workflows/post_processing.nf'
include { STAGE_INPUTS } from './subworkflows/stage_inputs.nf'

workflow {
    main:
    def schema = enabled(params.post_processing) ? 'post_processing_schema.json' : 'nextflow_schema.json'
    if ((params.help && params.help.toString() != 'false') || enabled(params.helpFull)) {
        def topic = params.help && !(params.help in [true, 'true', 'false']) ? params.help : ''
        log.info paramsHelp([parameters_schema: schema, fullHelp: enabled(params.helpFull), showHidden: enabled(params.showHidden)], topic)
        exit 0
    }
    if (enabled(params.version)) {
        log.info "GeneLab AmpliconSeq ${workflow.manifest.version}"
        exit 0
    }
    if (enabled(params.stage_only) && enabled(params.post_processing)) {
        error 'Specify only one of --stage_only or --post_processing'
    }
    if (params.validate_params == null || enabled(params.validate_params)) validateParameters(parameters_schema: schema)
    if (!(params.entry_point in ['raw_reads', 'trimmed_reads', 'asv_tables'])) error 'Supported entry points: raw_reads, trimmed_reads, asv_tables'
    if (params.entry_point == 'asv_tables' && params.truncate_to) error '--truncate_to applies only to raw or trimmed reads'
    if (enabled(params.post_processing)) {
        POST_PROCESSING()
        ch_published = POST_PROCESSING.out.published
    } else if (enabled(params.stage_only)) {
        STAGE_INPUTS()
        ch_published = STAGE_INPUTS.out.published
    } else {
        AMPLICON()
        ch_published = AMPLICON.out.published
    }
    publish:
    published = ch_published
}

output {
    published {
        path { dest, f -> f >> dest }
    }
}
