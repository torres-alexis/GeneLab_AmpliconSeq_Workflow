include { assaySuffix; shellQuote } from './utils.nf'

def stageRead(source, output, remote) {
    def records = params.truncate_to ? params.truncate_to.toInteger() : 0
    if (remote && records) {
        return "wget -O download.fastq.gz ${shellQuote(source)}\nfastq_first_n.sh ${records} ${shellQuote(output)} download.fastq.gz\nrm download.fastq.gz"
    }
    if (remote) return "wget -O ${shellQuote(output)} ${shellQuote(source)}"
    if (records) return "fastq_first_n.sh ${records} ${shellQuote(output)} ${shellQuote(source)}"
    return "cp -P ${shellQuote(source)} ${shellQuote(output)}"
}

def outputReads(id, paired, stage) {
    def prefix = "${id}${assaySuffix()}"
    if (paired == 'true') return ["${prefix}_R1_${stage}.fastq.gz", "${prefix}_R2_${stage}.fastq.gz"]
    def read = params.force_single_end ? "_${params.force_single_end}" : ''
    return ["${prefix}${read}_${stage}.fastq.gz"]
}

process COPY_READS {
    tag "${sample_id}"
    input:
    tuple val(sample_id), path('read?.gz'), val(paired)
    val(stage)
    output:
    tuple val(sample_id), path("${sample_id}*.gz"), val(paired), emit: raw_reads
    script:
    def outputs = outputReads(sample_id, paired, stage)
    outputs.withIndex().collect { name, i -> stageRead("read${i+1}.gz", name, false) }.join('\n')
}

process COPY_REMOTE_READS {
    tag "${sample_id}"
    input:
    tuple val(sample_id), val(paths), val(paired)
    val(stage)
    output:
    tuple val(sample_id), path("${sample_id}*.gz"), val(paired), emit: raw_reads
    script:
    def outputs = outputReads(sample_id, paired, stage)
    outputs.withIndex().collect { name, i -> stageRead(paths[i], name, true) }.join('\n')
}
