include { shellQuote } from './utils.nf'

process NORMALIZE_RUNSHEET {
    input:
    path(runsheet, stageAs: 'source/*')
    val(base_dir)
    output:
    path("${runsheet.toString().tokenize('/').last()}"), emit: runsheet
    script:
    def force = params.force_single_end ? "--force-single-end ${shellQuote(params.force_single_end)}" : ''
    """
    normalize_runsheet.py --input ${shellQuote(runsheet)} --output ${shellQuote(runsheet.toString().tokenize('/').last())} \\
        --base-dir ${shellQuote(base_dir)} --sample-column ${shellQuote(params.samples_column)} \\
        --group-column ${shellQuote(params.group)} --limit ${params.limit_samples_to ?: 0} ${force}
    """
}
