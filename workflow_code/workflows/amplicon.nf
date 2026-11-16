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
    gl_file_ch = STAGE_INPUTS.out.gl_file
    primers_ch = STAGE_INPUTS.out.primers
    ch_input_file = runsheet_ch
    STAGE_INPUTS.out.software_versions | mix(software_versions_ch) | set { software_versions_ch }


    // Read quality check and trimming
    RAW_FASTQC(staged_reads_ch)
    raw_fastqc_files = RAW_FASTQC.out.fastqc.flatten().collect()

    RAW_MULTIQC("raw", params.multiqc_config,raw_fastqc_files)

    RAW_FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    RAW_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}

    // Download reference database for taxonomic classification
    def db_config = [
            "16S": ["SILVA_SSU_r138_2_v2.RData", "https://api.figshare.com/v2/file/download/64078939"],
            "ITS": ["UNITE_v2025.RData", "https://api.figshare.com/v2/file/download/64079011"],
            "18S": ["PR2_v4_13_March2021.RData", "https://api.figshare.com/v2/file/download/46241917"]
        ]
    target_region_ch = Channel.value(params.target_region)
        .map { region -> tuple(region, db_config[region][0], db_config[region][1]) }

    DOWNLOAD_DATABASE(target_region_ch)

    trimmed_reads_ch = channel.empty()
    trimmed_reads_counts = channel.empty()
    cutadapt_logs = channel.empty()
    if(enabled(params.trim_primers)){

        CUTADAPT(staged_reads_ch, primers_ch)
        logs = CUTADAPT.out.logs.map{ sample_id, log -> file("${log}")}.collect()
        counts = CUTADAPT.out.trim_counts.map{ sample_id, count -> file("${count}")}.collect()
        trimmed_reads_ch = CUTADAPT.out.reads.map{
                                              sample_id, reads, isPaired -> reads instanceof List ? reads.each{file("${it}")}: file("${reads}")
                                              }.flatten().collect()

        COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE(counts, logs, runsheet_ch)
        trimmed_reads_counts = COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.counts
        cutadapt_logs = COMBINE_CUTADAPT_LOGS_AND_SUMMARIZE.out.logs

        isPaired_ch = CUTADAPT.out.reads.map{
                                              sample_id, reads, isPaired -> isPaired
                                              }.first()

        samples_ch = runsheet_ch.first()
                     .concat(isPaired_ch)
                     .collate(2)


        // Run dada2
        RUN_DADA2(samples_ch, trimmed_reads_ch, trimmed_reads_counts, DOWNLOAD_DATABASE.out.database)

        CUTADAPT.out.version | mix(software_versions_ch) | set{software_versions_ch}
    }else{
        raw_reads_ch = staged_reads_ch.map{
                          sample_id, reads, isPaired -> reads instanceof List ? reads.each{file("${it}")}: file("${reads}")
                          }.flatten().collect()

        isPaired_ch = staged_reads_ch.map{sample_id, reads, isPaired -> isPaired}.first()
        samples_ch = runsheet_ch.first()
                     .concat(isPaired_ch)
                     .collate(2)

        // Run dada2 without primer trimming
        RUN_DADA2(samples_ch, raw_reads_ch, [], DOWNLOAD_DATABASE.out.database)
    }

    dada_counts = RUN_DADA2.out.counts
    dada_taxonomy = RUN_DADA2.out.taxonomy
    dada_biom = RUN_DADA2.out.biom
    filtered_count = RUN_DADA2.out.filtered_count

    filtered_reads_ch = RUN_DADA2.out.reads
            .flatten()
            .map { read -> tuple(filteredSampleId(read.name, assaySuffix()), read) }
            .groupTuple(by:0)  // group R1/R2 by sample_id
            .map { sample_id, files ->
                    def pathFiles = files.collect { it instanceof String ? file(it) : it }  // ensure Path objects
                    def isPaired = pathFiles.size() > 1
                    tuple(sample_id, pathFiles, isPaired)
            }

    FILTERED_FASTQC(filtered_reads_ch)
        filtered_fastqc_files = FILTERED_FASTQC.out.fastqc.flatten().collect()

    FILTERED_MULTIQC("filtered", params.multiqc_config, filtered_fastqc_files)

    RUN_DADA2.out.version | mix(software_versions_ch) | set{software_versions_ch}
    FILTERED_FASTQC.out.version | mix(software_versions_ch) | set{software_versions_ch}
    FILTERED_MULTIQC.out.version | mix(software_versions_ch) | set{software_versions_ch}

    // Zip biom file
    dada_biom
        .map { biom -> tuple("taxonomy-and-counts", biom) }
        | ZIP_BIOM

    ZIP_BIOM.out.version | mix(software_versions_ch) | set{software_versions_ch}



    // Diversity, differential abundance testing and their corresponding visualizations
    if(params.accession){

        values = ["samples": "Sample Name",
                "group" : "groups",
                "depth" : params.rarefaction_depth,
                "assay_suffix" : assaySuffix(),
                "output_prefix" : cleanedPrefix(),
                "target_region" : params.target_region,
                "library_cutoff" : params.library_cutoff,
                "prevalence_cutoff" : params.prevalence_cutoff,
                "rare" : enabled(params.remove_rare) ? "--remove-rare" : "",
                "struc_zero": enabled(params.remove_struc_zeros) ? "--remove-structural-zeros" : ""
                ]

        metadata  =  runsheet_ch

    }else{

        values = ["samples": params.samples_column,
                "group" : params.group,
                "depth" : params.rarefaction_depth,
                "assay_suffix" : assaySuffix(),
                "output_prefix" : cleanedPrefix(),
                "target_region" : params.target_region,
                "library_cutoff" : params.library_cutoff,
                "prevalence_cutoff" : params.prevalence_cutoff,
                "rare" :  enabled(params.remove_rare) ? "--remove-rare" : "",
                "struc_zero": enabled(params.remove_struc_zeros) ? "--remove-structural-zeros" : ""
                ]

        metadata  =  ch_input_file

    }
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
            pub(gl_file_ch, 'Metadata'),
            pub(staged_reads_ch.map { id, reads, paired -> reads }, 'Raw_Sequence_Data'),
            pub(trimmed_reads_ch.mix(trimmed_reads_counts, cutadapt_logs), 'Trimmed_Sequence_Data'),
            pub(filtered_reads_ch.map { id, reads, paired -> reads }.mix(filtered_count), 'Filtered_Sequence_Data'),
            pub(RAW_FASTQC.out.fastqc, 'Raw_Sequence_Data/FastQC_Outputs'),
            pub(FILTERED_FASTQC.out.fastqc, 'Filtered_Sequence_Data/FastQC_Outputs'),
            pub(RAW_MULTIQC.out.zipped_data.mix(RAW_MULTIQC.out.html), 'Raw_Sequence_Data/MultiQC_Reports'),
            pub(FILTERED_MULTIQC.out.zipped_data.mix(FILTERED_MULTIQC.out.html), 'Filtered_Sequence_Data/MultiQC_Reports'),
            pub(RUN_DADA2.out.fasta.mix(RUN_DADA2.out.counts, RUN_DADA2.out.taxonomy, RUN_DADA2.out.taxonomy_count, ZIP_BIOM.out.zip, RUN_DADA2.out.read_count,
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
