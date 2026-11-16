include { enabled; pub; shellQuote } from '../modules/utils.nf'
include { STAGE_ANALYSIS } from './stage_analysis.nf'

process STAGE_METADATA {
    input:
    tuple val(kind), path(metadata, stageAs: 'source/*')
    output:
    tuple val(kind), path("${metadata.toString().tokenize('/').last()}")
    script:
    """
    cp -L ${shellQuote(metadata)} .
    """
}

workflow STAGE_INPUTS {
    main:
    def input_path = params.runsheet != null ? params.runsheet : params.input_file
    def isa_path = params.isa_archive
    if (!params.accession && !input_path) {
        error 'Supply --accession or --runsheet'
    }
    if (isa_path && !params.accession) error '--isa_archive requires --accession'
    def plugin = params.dp_tools_plugin ?: "$projectDir/bin/dp_tools__NF_AmpIllumina_${params.target_region}"
    STAGE_ANALYSIS(
        params.accession,
        params.target_region,
        input_path ? channel.fromPath(input_path, checkIfExists: true) : null,
        isa_path ? channel.fromPath(isa_path, checkIfExists: true) : null,
        params.api_url,
        channel.value(file(plugin))
    )
    runsheet_ch = STAGE_ANALYSIS.out.runsheet
    isa_ch = STAGE_ANALYSIS.out.isa_archive
    // Copy external metadata into task outputs for publication.
    STAGE_METADATA(runsheet_ch.map { tuple('runsheet', it) }.mix(isa_ch.map { tuple('isa', it) }).mix(STAGE_ANALYSIS.out.source_runsheet.map { tuple('source', it) }))
    runsheet_ch = STAGE_METADATA.out.filter { kind, f -> kind == 'runsheet' }.map { kind, f -> f }
    isa_ch = STAGE_METADATA.out.filter { kind, f -> kind == 'isa' }.map { kind, f -> f }
    primers_ch = STAGE_ANALYSIS.out.primers.map { primers ->
        if (!enabled(params.stage_only) && enabled(params.trim_primers) && (!primers[0] || !primers[1])) {
            error 'Primer removal requires forward and reverse primers in the runsheet or --F_primer and --R_primer'
        }
        primers
    }
    source_runsheet_ch = STAGE_METADATA.out.filter { kind, f -> kind == 'source' }.map { kind, f -> f }
    ch_published = pub(runsheet_ch.mix(isa_ch), 'Metadata')
        .mix(pub(source_runsheet_ch, 'Metadata/source'))
    ch_published = ch_published.mix(pub(STAGE_ANALYSIS.out.staged_reads.map { id, reads, paired -> reads }, 'Raw_Sequence_Data'))
    emit:
    published = ch_published
    staged_reads = STAGE_ANALYSIS.out.staged_reads
    runsheet = runsheet_ch
    isa_archive = isa_ch
    source_runsheet = source_runsheet_ch
    primers = primers_ch
    software_versions = STAGE_ANALYSIS.out.software_versions
}
