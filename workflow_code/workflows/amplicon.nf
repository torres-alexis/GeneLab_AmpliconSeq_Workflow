include { assaySuffix; pub; enabled; filteredSampleId; SOFTWARE_VERSIONS; cleanedPrefix } from '../modules/utils.nf'
nextflow.enable.dsl = 2


// Stage analysis setup (inputs, and raw reads)
include { STAGE_INPUTS } from '../subworkflows/stage_inputs.nf'

// Read quality check and filtering
include { FASTQC as RAW_FASTQC ; MULTIQC as RAW_MULTIQC  } from '../modules/quality_assessment.nf'

// Trim primers if requested
include { CUTADAPT ; COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE } from '../modules/cutadapt.nf'

// Cluster ASVs
include { DOWNLOAD_DATABASE } from '../modules/download_database.nf'
include { RUN_DADA2 } from '../modules/run_dada.nf'

// Filtered quality check
include { FASTQC as FILTERED_FASTQC ; MULTIQC as FILTERED_MULTIQC  } from '../modules/quality_assessment.nf'

// Diversity, differential abundance and visualizations
include { ALPHA_DIVERSITY; BETA_DIVERSITY } from '../modules/diversity.nf'
include { PLOT_TAXONOMY } from '../modules/taxonomy_plots.nf'
include { ZIP as ZIP_BIOM; ZIP as ZIP_ALPHA; ZIP as ZIP_BETA_EUCLIDEAN; ZIP as ZIP_BETA_BRAY; ZIP as ZIP_TAXONOMY_SAMPLES; ZIP as ZIP_TAXONOMY_GROUPS } from '../modules/zip.nf'
include { ANCOMBC as ANCOMBC1 } from '../modules/ancombc.nf'
include { ANCOMBC as ANCOMBC2 } from '../modules/ancombc.nf'
include { DESEQ } from '../modules/deseq.nf'
include { ZIP as ZIP_ANCOMBC1; ZIP as ZIP_ANCOMBC2; ZIP as ZIP_DESEQ2 } from '../modules/zip.nf'

workflow AMPLICON {
    main:


    /************************************************
    *********** Show pipeline parameters ************
    *************************************************/
    if(enabled(params.debug)){

    log.info """
         Nextflow AmpIllumina Consensus Pipeline: $workflow.manifest.version

         You have set the following parameters:

         Amplicon target region : ${params.target_region}
         GLDS or OSD accession : ${params.accession}
         Input csv file : ${params.runsheet != null ? params.runsheet : params.input_file}
         Output directory: ${params.outdir}
         Database Store Directory: ${params.database_store_path}
         Genelab Assay Suffix: ${assaySuffix()}
         Output Prefix: ${params.output_prefix}
         Trim Primers: ${enabled(params.trim_primers)}

         Cutadapt Parameters:
         Forward Primer: ${params.F_primer}
         Reverse Primer: ${params.R_primer}
         Minimum Trimmed Reads length: ${params.min_cutadapt_len}
         Primers Are linked: ${enabled(params.primers_linked)}
         Primers Are Anchored: ${params.anchored_primers}
         Discard Untrimmed Reads: ${enabled(params.discard_untrimmed)}

         Dada2 Parameters:
         Truncate left: ${params.left_trunc}bp
         Truncate right: ${params.right_trunc}bp
         Max error left: ${params.left_maxEE}
         Max error right: ${params.right_maxEE}
         Concatenate Reads: ${enabled(params.concatenate_reads_only)}

         Diversity and Differential abundance Parameters:
         Method: ${params.diff_abund_method}
         Rarefaction Depth: ${params.rarefaction_depth}
         Remove Structural Zeros: ${enabled(params.remove_struc_zeros)}
         Remove Rare Taxa and Samples: ${enabled(params.remove_rare)}
         Taxa Prevalence Cut Off: ${params.prevalence_cutoff}
         Sample Library Cut Off: ${params.library_cutoff}
         Groups to Comapre Column: ${params.group}
         Samples Column: ${params.samples_column}

         Debugging Options:
         Limit Samples for Testing: ${params.limit_samples_to}
         Force Processing Single-End: ${params.force_single_end}

         General Pipeline Settings:
         Nextflow Directory publishing mode: ${params.publish_dir_mode}
         MultiQC configuration file: ${params.multiqc_config}
         Nextflow Error strategy: ${params.errorStrategy}

         Conda Environments:
         dp_tools: ${params.conda_dp_tools}
         fastqc: ${params.conda_fastqc}
         multiqc: ${params.conda_multiqc}
         cutadapt: ${params.conda_cutadapt}
         R: ${params.conda_R}
         Diversity and Differential abundance : ${params.conda_diversity}
         zip: ${params.conda_zip}
         wget: ${params.conda_wget}
         """
    }

    software_versions_ch = channel.empty()
    STAGE_INPUTS()
    staged_reads_ch = STAGE_INPUTS.out.staged_reads
    runsheet_ch = STAGE_INPUTS.out.runsheet
    isa_archive_ch = STAGE_INPUTS.out.isa_archive
    source_runsheet_ch = STAGE_INPUTS.out.source_runsheet
    primers_ch = STAGE_INPUTS.out.primers
    ch_input_file = runsheet_ch
    STAGE_INPUTS.out.software_versions | mix(software_versions_ch) | set { software_versions_ch }


    raw_reads_ch = channel.empty()
    raw_fastqc_ch = channel.empty()
    raw_multiqc_zip_ch = channel.empty()
    raw_multiqc_html_ch = channel.empty()
    trimmed_reads_ch = channel.empty()
    trimmed_reads_counts = channel.empty()
    cutadapt_logs = channel.empty()
    filtered_reads_ch = channel.empty()
    filtered_fastqc_ch = channel.empty()
    filtered_multiqc_zip_ch = channel.empty()
    filtered_multiqc_html_ch = channel.empty()
    filtered_count = channel.value([])
    asv_ch = channel.empty()
    taxonomy_counts_ch = channel.empty()
    biom_zip_ch = channel.empty()
    read_count_ch = channel.empty()
    dada_counts = STAGE_INPUTS.out.counts
    dada_taxonomy = STAGE_INPUTS.out.taxonomy

    if (params.entry_point != 'asv_tables') {
        if (params.entry_point == 'raw_reads') {
            raw_reads_ch = staged_reads_ch
            RAW_FASTQC(staged_reads_ch)
            RAW_MULTIQC('raw', params.multiqc_config, RAW_FASTQC.out.fastqc.flatten().collect())
            raw_fastqc_ch = RAW_FASTQC.out.fastqc
            raw_multiqc_zip_ch = RAW_MULTIQC.out.zipped_data
            raw_multiqc_html_ch = RAW_MULTIQC.out.html
            software_versions_ch = software_versions_ch.mix(RAW_FASTQC.out.version, RAW_MULTIQC.out.version)
        }

        def db_config = [
            '16S': ['SILVA_SSU_r138_2_2024.RData', 'https://api.figshare.com/v2/file/download/52846199'],
            'ITS': ['UNITE_v2024_April2024.RData', 'https://api.figshare.com/v2/file/download/52846346'],
            '18S': ['PR2_v4_13_March2021.RData', 'https://api.figshare.com/v2/file/download/46241917']
        ]
        DOWNLOAD_DATABASE(channel.value(tuple(params.target_region, db_config[params.target_region][0], db_config[params.target_region][1])))

        dada_input_ch = staged_reads_ch
        dada_trim_counts = channel.value([])
        def input_stage = params.entry_point == 'trimmed_reads' ? 'trimmed' : 'raw'
        if (params.entry_point == 'raw_reads' && enabled(params.trim_primers)) {
            CUTADAPT(staged_reads_ch, primers_ch)
            dada_input_ch = CUTADAPT.out.reads
            COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE(
                CUTADAPT.out.trim_counts.map { id, count -> count }.collect(),
                CUTADAPT.out.logs.map { id, log -> log }.collect(), runsheet_ch)
            trimmed_reads_counts = COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.counts
            dada_trim_counts = trimmed_reads_counts
            cutadapt_logs = COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.logs
            input_stage = 'trimmed'
            software_versions_ch = software_versions_ch.mix(CUTADAPT.out.version)
        }
        if (input_stage == 'trimmed') {
            trimmed_reads_ch = dada_input_ch.map { id, reads, paired -> reads }.flatten()
        }
        isPaired_ch = dada_input_ch.map { id, reads, paired -> paired }.first()
        samples_ch = runsheet_ch.combine(isPaired_ch)
        reads_ch = dada_input_ch.map { id, reads, paired -> reads }.flatten().collect()
        RUN_DADA2(samples_ch, reads_ch, dada_trim_counts, DOWNLOAD_DATABASE.out.database, input_stage)
        dada_counts = RUN_DADA2.out.counts
        dada_taxonomy = RUN_DADA2.out.taxonomy
        filtered_count = RUN_DADA2.out.filtered_count
        asv_ch = RUN_DADA2.out.fasta
        taxonomy_counts_ch = RUN_DADA2.out.taxonomy_count
        read_count_ch = RUN_DADA2.out.read_count

        filtered_reads_ch = RUN_DADA2.out.reads.flatten()
            .map { read -> tuple(filteredSampleId(read.name, assaySuffix()), read) }
            .groupTuple(by: 0)
            .map { id, reads -> tuple(id, reads, reads.size() > 1) }
        FILTERED_FASTQC(filtered_reads_ch)
        FILTERED_MULTIQC('filtered', params.multiqc_config, FILTERED_FASTQC.out.fastqc.flatten().collect())
        filtered_fastqc_ch = FILTERED_FASTQC.out.fastqc
        filtered_multiqc_zip_ch = FILTERED_MULTIQC.out.zipped_data
        filtered_multiqc_html_ch = FILTERED_MULTIQC.out.html
        ZIP_BIOM(RUN_DADA2.out.biom.map { tuple('taxonomy-and-counts', it) })
        biom_zip_ch = ZIP_BIOM.out.zip
        software_versions_ch = software_versions_ch.mix(RUN_DADA2.out.version, FILTERED_FASTQC.out.version,
            FILTERED_MULTIQC.out.version, ZIP_BIOM.out.version)
    }

    // Diversity, differential abundance testing and their corresponding visualizations
    values = ["samples": "sample_id",
              "group": "groups",
              "depth": params.rarefaction_depth,
              "assay_suffix": assaySuffix(),
              "output_prefix": cleanedPrefix(),
              "target_region": params.target_region,
              "library_cutoff": params.library_cutoff,
              "prevalence_cutoff": params.prevalence_cutoff,
              "rare": enabled(params.remove_rare) ? "--remove-rare" : "",
              "struc_zero": enabled(params.remove_struc_zeros) ? "--remove-structural-zeros" : ""]
    metadata = runsheet_ch
    meta  = channel.of(values)

    // Diversity analysis
    ALPHA_DIVERSITY(meta, dada_counts, dada_taxonomy, metadata)
    BETA_DIVERSITY(meta, dada_counts, dada_taxonomy, metadata)

    // Zipping diversity plots
    // Alpha diversity (if rarefaction succeeded)
    ALPHA_DIVERSITY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.endsWith('.png') }
            pngs ? tuple(
                "alpha_diversity_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_ALPHA

    // Beta diversity - euclidean distance
    BETA_DIVERSITY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('euclidean') && it.name.endsWith('.png') }
            pngs ? tuple(
                "euclidean_distance_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_BETA_EUCLIDEAN

    // Beta diversity - bray curtis (if rarefaction succeeded)
    BETA_DIVERSITY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('bray') && it.name.endsWith('.png') }
            pngs ? tuple(
                "bray_curtis_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_BETA_BRAY

    // Taxonomy plotting
    PLOT_TAXONOMY(meta, dada_counts, dada_taxonomy, metadata)

    // Zipping taxonomy plots
   // Sample plots
    PLOT_TAXONOMY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('samples') && it.name.endsWith('.png') }
            pngs ? tuple(
                "sample_taxonomy_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_TAXONOMY_SAMPLES

   // Group taxonomy plots
    PLOT_TAXONOMY.out.output_dir
        .map { dir ->
            def pngs = file(dir).listFiles()?.findAll { it.name.contains('groups') && it.name.endsWith('.png') }
            pngs ? tuple(
                "group_taxonomy_plots",
                pngs
            ) : null
        }
        .filter { it != null }
        | ZIP_TAXONOMY_GROUPS

    ALPHA_DIVERSITY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    BETA_DIVERSITY.out.version | mix(software_versions_ch) | set{software_versions_ch}
    PLOT_TAXONOMY.out.version | mix(software_versions_ch) | set{software_versions_ch}

     // Differential abundance testing
     da_contrasts_ch = channel.empty()
     da_sampleTable_ch = channel.empty()
     ancombc1_ch = channel.empty()
     zip_ancombc1_ch = channel.empty()
     ancombc2_ch = channel.empty()
     zip_ancombc2_ch = channel.empty()
     deseq2_ch = channel.empty()
     zip_deseq2_ch = channel.empty()

     method = channel.of(params.diff_abund_method)
     if (params.diff_abund_method == "deseq2"){

        DESEQ(meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        deseq2_ch = DESEQ.out.output_dir
        da_contrasts_ch = DESEQ.out.contrasts_file
        da_sampleTable_ch = DESEQ.out.sample_table_file
        DESEQ.out.version | mix(software_versions_ch) | set{software_versions_ch}
        // Zipping DESeq2 plots
        DESEQ.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "deseq2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_DESEQ2
        zip_deseq2_ch = ZIP_DESEQ2.out.zip

    }else if (params.diff_abund_method == "ancombc1"){

        ANCOMBC1(method, meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        ancombc1_ch = ANCOMBC1.out.output_dir
        da_contrasts_ch = ANCOMBC1.out.contrasts_file
        da_sampleTable_ch = ANCOMBC1.out.sample_table_file
        ANCOMBC1.out.version | mix(software_versions_ch) | set{software_versions_ch}
        // Zipping ANCOMBC1 plots
        ANCOMBC1.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc1_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC1
        zip_ancombc1_ch = ZIP_ANCOMBC1.out.zip

    }else if (params.diff_abund_method == "ancombc2"){

        ANCOMBC2(method, meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        ancombc2_ch = ANCOMBC2.out.output_dir
        da_contrasts_ch = ANCOMBC2.out.contrasts_file
        da_sampleTable_ch = ANCOMBC2.out.sample_table_file
        ANCOMBC2.out.version | mix(software_versions_ch) | set{software_versions_ch}
        // Zipping ANCOMBC2 plots
        ANCOMBC2.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC2
        zip_ancombc2_ch = ZIP_ANCOMBC2.out.zip

    }else{

        ANCOMBC1("ancombc1", meta, dada_counts, dada_taxonomy, metadata, filtered_count)
        ancombc1_ch = ANCOMBC1.out.output_dir
        da_contrasts_ch = ANCOMBC1.out.contrasts_file
        da_sampleTable_ch = ANCOMBC1.out.sample_table_file
        ANCOMBC1.out.version | mix(software_versions_ch) | set{software_versions_ch}

        ANCOMBC2("ancombc2", meta, dada_counts, dada_taxonomy, metadata, ANCOMBC1.out.output_dir)
        ancombc2_ch = ANCOMBC2.out.output_dir
        ANCOMBC2.out.version | mix(software_versions_ch) | set{software_versions_ch}

        DESEQ(meta, dada_counts, dada_taxonomy, metadata, ANCOMBC2.out.output_dir)
        deseq2_ch = DESEQ.out.output_dir
        DESEQ.out.version | mix(software_versions_ch) | set{software_versions_ch}

        // Zipping DA plots
        //ANCOMBC1
        ANCOMBC1.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc1_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC1
        zip_ancombc1_ch = ZIP_ANCOMBC1.out.zip
        //ANCOMBC2
        ANCOMBC2.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "ancombc2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_ANCOMBC2
        zip_ancombc2_ch = ZIP_ANCOMBC2.out.zip
        // DESeq2
        DESEQ.out.output_dir
            .map { dir ->
                def pngs = file(dir).listFiles()?.findAll { it.name.contains('volcano') && it.name.endsWith('.png') }
                pngs ? tuple(
                    "deseq2_volcano_plots",
                    pngs
                ) : null
            }
            .filter { it != null }
            | ZIP_DESEQ2
        zip_deseq2_ch = ZIP_DESEQ2.out.zip
    }


     // Software Version Capturing - combining all captured software versions
     nf_version = "Nextflow Version ".concat("${nextflow.version}")
     nextflow_version_ch = channel.value(nf_version)

     //  Write software versions to file
     software_versions_ch | map { it.text.strip() }
                          | unique
                          | mix(nextflow_version_ch)
                          | collectFile({it -> it}, newLine: true, cache: false)
                          | SOFTWARE_VERSIONS

    ch_published = pub(runsheet_ch.mix(isa_archive_ch), 'Metadata')
        .mix(
            pub(source_runsheet_ch, 'Metadata/source'),
            pub(raw_reads_ch.map { id, reads, paired -> reads }, 'Raw_Sequence_Data'),
            pub(trimmed_reads_ch.mix(trimmed_reads_counts, cutadapt_logs), 'Trimmed_Sequence_Data'),
            pub(filtered_reads_ch.map { id, reads, paired -> reads }.mix(filtered_count), 'Filtered_Sequence_Data'),
            pub(raw_fastqc_ch, 'Raw_Sequence_Data/FastQC_Outputs'),
            pub(filtered_fastqc_ch, 'Filtered_Sequence_Data/FastQC_Outputs'),
            pub(raw_multiqc_zip_ch.mix(raw_multiqc_html_ch), 'Raw_Sequence_Data/MultiQC_Reports'),
            pub(filtered_multiqc_zip_ch.mix(filtered_multiqc_html_ch), 'Filtered_Sequence_Data/MultiQC_Reports'),
            pub(asv_ch.mix(dada_counts, dada_taxonomy, taxonomy_counts_ch, biom_zip_ch, read_count_ch,
                ALPHA_DIVERSITY.out.output_dir, BETA_DIVERSITY.out.output_dir, PLOT_TAXONOMY.out.output_dir), 'Final_Outputs'),
            pub(ZIP_ALPHA.out.zip, 'Final_Outputs/alpha_diversity'),
            pub(ZIP_BETA_EUCLIDEAN.out.zip.mix(ZIP_BETA_BRAY.out.zip), 'Final_Outputs/beta_diversity'),
            pub(ZIP_TAXONOMY_SAMPLES.out.zip.mix(ZIP_TAXONOMY_GROUPS.out.zip), 'Final_Outputs/taxonomy_plots'),
            pub(da_contrasts_ch.mix(da_sampleTable_ch, ancombc1_ch, ancombc2_ch, deseq2_ch), 'Final_Outputs/differential_abundance'),
            pub(zip_ancombc1_ch, 'Final_Outputs/differential_abundance/ancombc1'),
            pub(zip_ancombc2_ch, 'Final_Outputs/differential_abundance/ancombc2'),
            pub(zip_deseq2_ch, 'Final_Outputs/differential_abundance/deseq2'),
            pub(SOFTWARE_VERSIONS.out.software_versions_txt, 'GeneLab')
        )
    emit:
    published = ch_published
}
