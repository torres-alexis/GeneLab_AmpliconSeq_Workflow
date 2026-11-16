include { GET_ACCESSIONS } from '../modules/get_accessions.nf'
include { FETCH_ISA } from '../modules/fetch_isa.nf'
include { ISA_TO_RUNSHEET } from '../modules/isa_to_runsheet.nf'
include { SELECT_RUNSHEET } from '../modules/select_runsheet.nf'
include { NORMALIZE_RUNSHEET } from '../modules/normalize_runsheet.nf'
include { COPY_READS; COPY_REMOTE_READS } from '../modules/copy_reads.nf'

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
    NORMALIZE_RUNSHEET(source_runsheet_ch, launchDir.toString())
    runsheet_ch = NORMALIZE_RUNSHEET.out.runsheet
    primers_ch = runsheet_ch.splitCsv(header: true).first().map { row ->
        [params.F_primer ?: row.F_Primer ?: '', params.R_primer ?: row.R_Primer ?: '']
    }
    reads_ch = runsheet_ch.splitCsv(header: true).map { row ->
        def reads = row.paired == 'true' ? [row.forward, row.reverse] : [row.forward]
        tuple(row.sample_id, reads, row.paired)
    }
    split_reads = reads_ch.branch { id, reads, paired ->
        remote: reads.every { it ==~ /(?i)^(https?|ftp):\/\/.*/ }
        local: true
    }
    COPY_REMOTE_READS(split_reads.remote)
    COPY_READS(split_reads.local.map { id, reads, paired -> tuple(id, reads.collect { file(it, checkIfExists: true) }, paired) })

    emit:
    staged_reads = COPY_READS.out.raw_reads.mix(COPY_REMOTE_READS.out.raw_reads)
    runsheet = runsheet_ch
    source_runsheet = source_runsheet_ch
    isa_archive = isa_archive_ch
    primers = primers_ch
    software_versions = software_versions_ch
}
