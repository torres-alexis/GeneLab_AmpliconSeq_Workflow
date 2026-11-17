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

For entry points other than the default raw_reads, replace the "forward" and "reverse" columns above with the relevant input columns below. All entry points require "sample_id" and "groups"; "paired" is required for raw_reads and trimmed_reads.

| Column Name | Entry Point | Type | Description | Example |
|:------------|:------------|:-----|:------------|:--------|
| trimmed_read1_path | trimmed_reads | string | Path or URL to primer-trimmed forward reads | /my/data/sample1_16S_GLAmpSeq_R1_trimmed.fastq.gz |
| trimmed_read2_path | trimmed_reads | string | Path or URL to primer-trimmed reverse reads (paired-end only) | /my/data/sample1_16S_GLAmpSeq_R2_trimmed.fastq.gz |
| counts_table_path | asv_tables | string | Path or URL to the dataset-wide ASV count table | /my/data/counts_16S_GLAmpSeq.tsv |
| taxonomy_table_path | asv_tables | string | Path or URL to the corresponding taxonomy table | /my/data/taxonomy_16S_GLAmpSeq.tsv |

For asv_tables, repeat the same count-table and taxonomy-table paths for every sample row. Both tables must be tab-separated, with ASV IDs in the first column. The count table must contain a column for every runsheet sample; additional samples are excluded. Counts must be nonnegative integers. The taxonomy table must include every count-table ASV and the columns "domain", "phylum", "class", "order", "family", "genus" and "species".

Relative paths are resolved from the directory where Nextflow is launched.
