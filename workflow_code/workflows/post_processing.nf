include { assaySuffix; processingTargets; pub; enabled } from '../modules/utils.nf'
nextflow.enable.dsl=2


include { PACKAGE_PROCESSING_INFO; GENERATE_README; VALIDATE_PROCESSING;
           GENERATE_CURATION_TABLE; GENERATE_MD5SUMS; GENERATE_PROTOCOL} from '../modules/genelab.nf'


workflow POST_PROCESSING {
      main:


        /************************************************
        *********** Show pipeline parameters ************
        *************************************************/
        if(enabled(params.debug)){

        log.info """
         GeneLab Post Processing Pipeline: $workflow.manifest.version

         You have set the following parameters:

         Input/Output Options:
         GLDS Accession : ${params.GLDS_accession}
         OSD Accession : ${params.OSD_accession}
         Input Runsheet: ${params.runsheet}
         Output Directory: ${params.outdir}
         Assay Suffix: ${assaySuffix()}
         Output Prefix: ${params.output_prefix}
         Target Files: ${processingTargets()}

         Analyst Options:
         Analyst's Name : ${params.name}
         Analyst's Email : ${params.email}
         Protocol ID: ${params.protocol_id}
         V & V Link: ${params.V_V_guidelines_link}

         Processing Flags:
         Include Raw Data: ${enabled(params.include_raw_data)}
         Trim Primers: ${enabled(params.trim_primers)}
         Single End: ${enabled(params.single_end)}
         Force Processing Single-End: ${params.force_single_end}
         Include Raw MultiQC: ${enabled(params.include_raw_multiqc)}

         General Pipeline Settings:
         Profile: ${workflow.profile}
         Nextflow Directory Publishing Mode: ${params.publish_dir_mode}
         """

        }

        // Make sure accessions numbers are set
        if(!params.GLDS_accession || !params.OSD_accession){
           error("""ACCESSION ERROR!.
                    Please supply both --GLDS_accession and --OSD_accession.
                    They can be any string you choose but they must be set.
                 """)
        }

        // If force_single_end is set, ignore include_raw_data and include_raw_multiqc as they should reflect paired-end data which is assumed to have been already included
        if (params.force_single_end && (enabled(params.include_raw_data) || enabled(params.include_raw_multiqc))) {
            log.warn "force_single_end is set — --include_raw_data and --include_raw_multiqc will both be ignored."
        }

       // ---------------------- Input channels -------------------------------- //
       // Input files
       // Runsheet used to execute the processing workflow
       if (!params.runsheet) error 'Supply --runsheet'
       runsheet_ch = Channel.fromPath(params.runsheet, checkIfExists: true)

       is_single_end_ch = runsheet_ch.splitCsv(header: true).toList().map { rows ->
           def layouts = rows.collect { it.paired ?: it.paired_end ?: '' }.findAll { it }.collect { it.toString().toLowerCase() }.unique()
           if (layouts.size() > 1 || layouts.any { !(it in ['true', 'false']) }) error 'Post-processing runsheet has inconsistent read layouts'
           layouts ? layouts[0] == 'false' : enabled(params.single_end)
       }

       def processed_dir = file(params.outdir)
       def run_record = processed_dir.resolve('processing_scripts/workflow_run.json')
       def recorded_run = run_record.exists() ? new groovy.json.JsonSlurper().parseText(run_record.text) : [:]
       if (recorded_run.containsKey('capture_status') && recorded_run.capture_status != 'complete')
           error 'Processing-record capture is incomplete. Complete the main workflow before post-processing.'
       def source_entry = recorded_run.entry_point ?: params.entry_point
       if (!(source_entry in ['raw_reads', 'trimmed_reads', 'asv_tables'])) error "Invalid source entry point: ${source_entry}"
       def required_dirs = ['processing_scripts', 'Final_Outputs']
       if (source_entry == 'raw_reads') required_dirs << 'Raw_Sequence_Data/MultiQC_Reports'
       def versions_file = processed_dir.resolve("GeneLab/software_versions${assaySuffix()}.txt")
       def missing = required_dirs.findAll { !processed_dir.resolve(it).isDirectory() }
           .collect { processed_dir.resolve(it).toString() }
       if (!versions_file.isFile()) missing << versions_file.toString()
       processingTargets().split(',').each { name ->
           def record = processed_dir.resolve('processing_scripts').resolve(name.trim())
           if (!record.exists()) missing << record.toString()
       }
       if (missing) error "Missing post-processing input(s): ${missing.join(', ')}. Run analysis first and use its output root and naming parameters."

       processing_scripts_dir = channel.value(processed_dir.resolve('processing_scripts'))
       raw_multiqc_report_dir = channel.value(source_entry == 'raw_reads' ? processed_dir.resolve('Raw_Sequence_Data/MultiQC_Reports') : [])
       final_outputs_dir = channel.value(processed_dir.resolve('Final_Outputs'))
       software_versions_ch = channel.value(versions_file)
       rarefaction_depth_ch = channel
           .fromPath("${processed_dir}/Final_Outputs/alpha_diversity/*rarefaction_depth${assaySuffix()}.txt")
           .ifEmpty { [] }

        // ---------------------- Post-processing begins -------------------------------------- //

        PACKAGE_PROCESSING_INFO(processing_scripts_dir)

        GENERATE_README(source_entry, is_single_end_ch)

        // Automatic verification and validation
        VALIDATE_PROCESSING(channel.value(processed_dir), processingTargets(), runsheet_ch,
                            GENERATE_README.out.readme,
                            PACKAGE_PROCESSING_INFO.out.zip, source_entry, is_single_end_ch)

        // Generate md5sums
        GENERATE_MD5SUMS(channel.value(processed_dir), GENERATE_README.out.readme,
                        PACKAGE_PROCESSING_INFO.out.zip, source_entry)

        // Generate file association table for curation
        GENERATE_CURATION_TABLE(runsheet_ch, raw_multiqc_report_dir, final_outputs_dir, source_entry, is_single_end_ch)

        // Write methods
        GENERATE_PROTOCOL(software_versions_ch, params.protocol_id, rarefaction_depth_ch)

        ch_published = pub(PACKAGE_PROCESSING_INFO.out.zip.mix(
            GENERATE_README.out.readme, VALIDATE_PROCESSING.out.log,
            GENERATE_MD5SUMS.out.raw_md5sum, GENERATE_MD5SUMS.out.processed_md5sum,
            GENERATE_CURATION_TABLE.out.curation_table, GENERATE_PROTOCOL.out.protocol), 'GeneLab')

      emit:
      published = ch_published
}
