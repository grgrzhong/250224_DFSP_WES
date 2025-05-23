process FASTP_TRIM {
    
    tag "${meta.id}"
    
    label 'process_medium'
    
    input:
    tuple val(meta), path(reads1), path(reads2)
    
    output:
    tuple val(meta), path("*_trimmed_1.fastq.gz"), path("*_trimmed_2.fastq.gz"), emit: reads
    path "*.json", emit: json
    path "*.html", emit: html
    path "*.log", emit: log
    path "versions.yml", emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Get the kit type from params, defaulting to illumina if not specified
    def kit = params.kit ?: 'illumina'
    
    // Simple approach to access kit configuration
    def kit_config = [:]
    if (params.kits && params.kits.containsKey(kit)) {
        kit_config = params.kits[kit]
    }
    
    // Get adapter sequences
    def adapter = kit_config?.adapter_sequnce ?: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    def adapter_r2 = kit_config?.adapter_sequnce_r2 ?: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    
    // Get UMI settings
    def umi_enabled = kit_config?.umi_enabled ?: false
    def umi_location = kit_config?.umi_location ?: 'per_read'
    def umi_length = kit_config?.umi_length ?: 8
    
    // Build UMI options
    def umi_opts = ""
    if (umi_enabled) {
        umi_opts = "--umi --umi_loc ${umi_location} --umi_len ${umi_length}"
    }

    """
    fastp \\
        -i ${reads1} \\
        -I ${reads2} \\
        -o ${prefix}_trimmed_1.fastq.gz \\
        -O ${prefix}_trimmed_2.fastq.gz \\
        --detect_adapter_for_pe \\
        --adapter_sequence=${adapter} \\
        --adapter_sequence_r2=${adapter_r2} \\
        --trim_poly_g \\
        --trim_poly_x \\
        ${umi_enabled ? umi_opts + " \\\\" : ""} 
        -j ${prefix}.json \\
        -h ${prefix}.html \\
        -w ${task.cpus} \\
        > ${prefix}.fastp.log 2>&1
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}