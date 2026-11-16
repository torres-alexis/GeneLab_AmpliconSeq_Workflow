include { shellQuote } from './utils.nf'

process UPDATE_OSDR_INPUTS {
    input:
    path(runsheet, stageAs: 'source/*')
    val(accession)
    output:
    path("${runsheet.toString().tokenize('/').last()}"), emit: runsheet
    script:
    """
    update_osdr_inputs.py --input ${shellQuote(runsheet)} --output ${shellQuote(runsheet.toString().tokenize('/').last())} \\
        --accession ${shellQuote(accession)} --api-url ${shellQuote(params.api_url)} \\
        --entry-point ${shellQuote(params.entry_point)} --target-region ${shellQuote(params.target_region)} --limit ${params.limit_samples_to ?: 0}
    """
}
