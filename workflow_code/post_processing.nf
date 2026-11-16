include { enabled } from './modules/utils.nf'
nextflow.enable.dsl = 2
include { validateParameters; paramsHelp } from 'plugin/nf-schema'
include { POST_PROCESSING } from './workflows/post_processing.nf'
workflow {
    main:
    def schema = 'post_processing_schema.json'
    if ((params.help && params.help.toString() != 'false') || enabled(params.helpFull)) {
        def topic = params.help && !(params.help in [true, 'true', 'false']) ? params.help : ''
        log.info paramsHelp([parameters_schema: schema, fullHelp: enabled(params.helpFull), showHidden: enabled(params.showHidden)], topic)
        exit 0
    }
    if (enabled(params.version)) { log.info "GeneLab AmpliconSeq ${workflow.manifest.version}"; exit 0 }
    if (enabled(params.stage_only)) error 'The post_processing.nf entrypoint does not support --stage_only'
    if (params.validate_params == null || enabled(params.validate_params)) validateParameters(parameters_schema: schema)
    POST_PROCESSING()
    publish:
    published = POST_PROCESSING.out.published
}

output {
    published {
        path { dest, f -> f >> dest }
    }
}
