include { enabled; pub } from '../modules/utils.nf'
include { STAGE_ANALYSIS } from './stage_analysis.nf'

process STAGE_METADATA {
    input:
    tuple val(kind), path(metadata, stageAs: 'source/*')
    output:
    tuple val(kind), path("${metadata.toString().tokenize('/').last()}")
    script:
    """
    cp -L '${metadata}' .
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
    if (!enabled(params.stage_only) && !params.accession && enabled(params.trim_primers) && (!params.F_primer || !params.R_primer)) {
        error 'Primer trimming requires --F_primer and --R_primer with an upload CSV'
    }
    def plugin = params.dp_tools_plugin ?: "$projectDir/bin/dp_tools__NF_AmpIllumina_${params.target_region}"
    STAGE_ANALYSIS(
        params.accession,
        params.target_region,
        !params.accession && input_path ? channel.fromPath(input_path, checkIfExists: true) : null,
        isa_path ? channel.fromPath(isa_path, checkIfExists: true) : null,
        params.api_url,
        channel.value(file(plugin))
    )
    runsheet_ch = STAGE_ANALYSIS.out.runsheet
    isa_ch = STAGE_ANALYSIS.out.isa_archive
    // Copy external metadata into task outputs for publication.
    STAGE_METADATA(runsheet_ch.map { tuple('runsheet', it) }.mix(isa_ch.map { tuple('isa', it) }))
    runsheet_ch = STAGE_METADATA.out.filter { kind, f -> kind == 'runsheet' }.map { kind, f -> f }
    isa_ch = STAGE_METADATA.out.filter { kind, f -> kind == 'isa' }.map { kind, f -> f }
    ch_published = pub(runsheet_ch.mix(isa_ch, STAGE_ANALYSIS.out.gl_file), 'Metadata')
        .mix(pub(STAGE_ANALYSIS.out.staged_reads.map { id, reads, paired -> reads }, 'Raw_Sequence_Data'))
    emit:
    published = ch_published
    staged_reads = STAGE_ANALYSIS.out.staged_reads
    runsheet = runsheet_ch
    isa_archive = isa_ch
    gl_file = STAGE_ANALYSIS.out.gl_file
    primers = STAGE_ANALYSIS.out.primers
    software_versions = STAGE_ANALYSIS.out.software_versions
}
