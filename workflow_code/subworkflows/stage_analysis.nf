include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { SELECT_RUNSHEET } from '../modules/select_runsheet.nf'
include { UPDATE_OSDR_INPUTS } from '../modules/update_osdr_inputs.nf'
include { STAGE_ASV_TABLES } from '../modules/stage_asv_tables.nf'
include { NORMALIZE_RUNSHEET } from '../modules/normalize_runsheet.nf'
include { COPY_READS; COPY_REMOTE_READS } from '../modules/copy_reads.nf'
include { enabled } from '../modules/utils.nf'

workflow STAGE_ANALYSIS {
    take:
    accession
    target_region
    input_file
    isa_archive_path
    api_url
    dp_tools_plugin

    main:
    isa_archive_ch = channel.empty()
    software_versions_ch = channel.empty()
    if (accession) {
        if (!input_file || !isa_archive_path) {
            GET_ACCESSIONS(accession, api_url)
            osd_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[0].trim() }
            glds_accession = GET_ACCESSIONS.out.accessions_txt.map { it.readLines()[1].trim() }
        }
        if (isa_archive_path) {
            isa_archive_ch = isa_archive_path
        } else {
            FETCH_ISA(osd_accession, glds_accession)
            isa_archive_ch = FETCH_ISA.out.isa_archive
        }
    }
    if (input_file) {
        source_runsheet_ch = input_file
    } else {
        ISA_TO_RUNSHEET(osd_accession, glds_accession, isa_archive_ch, dp_tools_plugin)
        software_versions_ch = ISA_TO_RUNSHEET.out.version
        SELECT_RUNSHEET(ISA_TO_RUNSHEET.out.runsheet, target_region, params.specify_runsheet ?: '')
        source_runsheet_ch = SELECT_RUNSHEET.out.runsheet
    }
    selected_runsheet_ch = source_runsheet_ch
    if (!input_file && params.entry_point != 'raw_reads') {
        UPDATE_OSDR_INPUTS(source_runsheet_ch, osd_accession)
        selected_runsheet_ch = UPDATE_OSDR_INPUTS.out.runsheet
    }
    NORMALIZE_RUNSHEET(selected_runsheet_ch, launchDir.toString())
    runsheet_ch = NORMALIZE_RUNSHEET.out.runsheet
    primers_ch = runsheet_ch.splitCsv(header: true).toList().map { rows ->
        def primers = rows.collect { row ->
            [params.F_primer ?: row.F_Primer ?: '', params.R_primer ?: row.R_Primer ?: '']
        }.unique()
        if (params.entry_point == 'raw_reads' && enabled(params.trim_primers) && !enabled(params.stage_only) && primers.size() > 1)
            error 'Runsheet contains inconsistent primer sequences; supply --F_primer and --R_primer to override them'
        primers[0]
    }
    staged_reads_ch = channel.empty()
    counts_ch = channel.empty()
    taxonomy_ch = channel.empty()
    if (params.entry_point == 'asv_tables') {
        table_paths = runsheet_ch.splitCsv(header: true).first()
        STAGE_ASV_TABLES(
            table_paths.map { file(it.counts_table_path, checkIfExists: true) },
            table_paths.map { file(it.taxonomy_table_path, checkIfExists: true) }, runsheet_ch)
        counts_ch = STAGE_ASV_TABLES.out.counts
        taxonomy_ch = STAGE_ASV_TABLES.out.taxonomy
    } else {
        def stage = params.entry_point == 'trimmed_reads' ? 'trimmed' : 'raw'
        reads_ch = runsheet_ch.splitCsv(header: true).map { row ->
            def first = stage == 'trimmed' ? row.trimmed_read1_path : row.forward
            def second = stage == 'trimmed' ? row.trimmed_read2_path : row.reverse
            tuple(row.sample_id, row.paired == 'true' ? [first, second] : [first], row.paired)
        }
        split_reads = reads_ch.branch { id, reads, paired ->
            remote: reads.every { it ==~ /(?i)^(https?|ftp):\/\/.*/ }
            local: true
        }
        COPY_REMOTE_READS(split_reads.remote, stage)
        COPY_READS(split_reads.local.map { id, reads, paired -> tuple(id, reads.collect { file(it, checkIfExists: true) }, paired) }, stage)
        staged_reads_ch = COPY_READS.out.raw_reads.mix(COPY_REMOTE_READS.out.raw_reads)
    }

    emit:
    staged_reads = staged_reads_ch
    counts = counts_ch
    taxonomy = taxonomy_ch
    runsheet = runsheet_ch
    source_runsheet = source_runsheet_ch
    isa_archive = isa_archive_ch
    primers = primers_ch
    software_versions = software_versions_ch
}
