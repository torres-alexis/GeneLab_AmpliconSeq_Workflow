// The string "false" is truthy in Groovy.
def enabled(value) {
    if (value instanceof Boolean) return value
    if (value == null) return false
    if (value.toString() in ['true', 'false']) return value.toString() == 'true'
    throw new IllegalArgumentException("Expected true or false, received: ${value}")
}

// Match only a trailing filename suffix; sample IDs may themselves contain suffix text.
def filteredSampleId(name, assay_suffix) {
    def ending = ['_R1_filtered.fastq.gz', '_R2_filtered.fastq.gz', '_filtered.fastq.gz']
        .collect { "${assay_suffix}${it}" }
        .find { name.endsWith(it) }
    if (!ending) throw new IllegalArgumentException("Unrecognized filtered read filename: ${name}")
    return name.substring(0, name.length() - ending.length())
}

def pub(ch, dir) {
    return ch.flatten().map { f -> ["${dir}/${f.name}".toString(), f] }
}

def cleanedPrefix() {
    def prefix = params.output_prefix ?: ""
    return (prefix && !prefix.endsWith("_") && !prefix.endsWith("-")) ? prefix + "_" : prefix
}

process SOFTWARE_VERSIONS {

    tag "Writing out software versions..."
    label "fastqc"   //unix environment

    input:
        path(software_versions)

    output:
        path("software_versions${params.assay_suffix}.txt"), emit: software_versions_txt

    script:
        """
        # Delete white spaces and write out unique software versions 
        grep -v "^\$" ${software_versions} | sort -u > software_versions${params.assay_suffix}.txt 
        """
}
