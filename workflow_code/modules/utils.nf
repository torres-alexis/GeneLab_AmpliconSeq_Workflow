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

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def pub(ch, dir) {
    return ch.flatten().flatMap { f ->
        def children = f.isDirectory() ? f.listFiles() : []
        def output_root = workflow.outputDir.normalize()
        def destination_dir = output_root.resolve(dir)
        if (children) destination_dir = destination_dir.resolve(f.name)
        def relative_dir = output_root.relativize(destination_dir)
        (1..relative_dir.nameCount).each { depth ->
            def parent_dir = output_root.resolve(relative_dir.subpath(0, depth))
            if (parent_dir != output_root && java.nio.file.Files.isSymbolicLink(parent_dir)) {
                def link_target = parent_dir.parent.resolve(java.nio.file.Files.readSymbolicLink(parent_dir)).normalize()
                throw new IllegalArgumentException("Cannot publish into ${parent_dir}: this existing output folder points to ${link_target}. Publishing here would write into that linked directory. Use a new output directory with --outdir (or -o if specified).")
            }
        }
        if (children)
            return children.collect { child -> ["${dir}/${f.name}/${child.name}".toString(), child] }
        [["${dir}/${f.name}".toString(), f]]
    }
}

def cleanedPrefix() {
    def prefix = params.output_prefix ?: ""
    return (prefix && !prefix.endsWith("_") && !prefix.endsWith("-")) ? prefix + "_" : prefix
}

def assaySuffix() {
    def suffix = params.assay_suffix
    def region = params.target_region
    if (params.prepend_target_region.toString().toBoolean() && !suffix.startsWith("_${region}"))
        return "_${region}${suffix}"
    return suffix
}

def processingTargets() {
    return params.target_files ?: "command.txt,nextflow_processing_info${assaySuffix()}.txt"
}

process SOFTWARE_VERSIONS {

    tag "Writing out software versions..."
    label "fastqc"   //unix environment

    input:
        path(software_versions)

    output:
        path({ "software_versions${assaySuffix()}.txt" }), emit: software_versions_txt

    script:
        """
        # Delete white spaces and write out unique software versions 
        grep -v "^\$" ${software_versions} | sort -u > software_versions${assaySuffix()}.txt 
        """
}
