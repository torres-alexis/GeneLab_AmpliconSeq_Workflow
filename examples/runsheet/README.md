# Runsheet File Specification

## Description

* The runsheet is a comma-separated file that contains the metadata required for processing 
amplicon sequence datasets through the GeneLab amplicon Illumina sequencing data 
processing pipeline (AmpIllumina).


## Examples

1. Example runsheet for a [paired-end dataset](PE_file.csv)
2. Example runsheet for a [single-end dataset](SE_file.csv)


## Required columns

| Column Name | Type | Description | Example |
|:------------|:-----|:------------|:--------|
| sample_id | string | Unique Sample Name, added as a prefix to sample-specific processed data output files. Should not include spaces or weird characters. | RR23_FCS_FLT_F1 |
| forward | string (local path or URL) | Location of the raw reads file. For paired-end data, this specifies the forward reads fastq.gz file. | /my/data/sample1_R1_raw.fastq.gz |
| reverse | string (local path or URL) | Location of the raw reads file. For paired-end data, this specifies the reverse reads fastq.gz file. For single-end data, this column should be omitted. | /my/data/sample1_R2_raw.fastq.gz |
| paired | bool | Set to True if the samples were sequenced as paired-end. If set to False, samples are assumed to be single-end. | False |
| groups | string | Name of the treatment group that the sample belongs to. Should not include spaces or weird characters. | SpaceFlight |

## Entry Point Columns

For entry points other than the default raw_reads, include the relevant columns below to specify the input files.

| Column Name | Entry Point | Type | Description | Example |
|:------------|:------------|:-----|:------------|:--------|
| trimmed_read1_path | trimmed_reads | string | Path or URL to primer-trimmed forward reads | /my/data/sample1_GLAmpSeq_R1_trimmed.fastq.gz |
| trimmed_read2_path | trimmed_reads | string | Path or URL to primer-trimmed reverse reads (paired-end only) | /my/data/sample1_GLAmpSeq_R2_trimmed.fastq.gz |
| counts_table_path | asv_tables | string | Path or URL to the ASV count table (tab-separated) | /my/data/counts_GLAmpSeq.tsv |
| taxonomy_table_path | asv_tables | string | Path or URL to the matching ASV taxonomy table (tab-separated) | /my/data/taxonomy_GLAmpSeq.tsv |

For asv_tables, use the same count-table and taxonomy-table paths for every sample row.
