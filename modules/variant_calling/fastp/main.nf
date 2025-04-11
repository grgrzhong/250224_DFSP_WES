process FASTP_TRIM {
    tag "${meta.id}"
    label 'process_medium'
    
    // container "${params.singularity_container_dir}/fastp-0.24.0.sif"

    // publishDir([
    //     [
    //         path: { "${params.outdir}/preprocessing/fastq_trimmed/${meta.id}" },
    //         mode: params.publish_dir_mode,
    //         pattern: "*_trimmed_{1,2}.fastq.gz",
    //         enabled: params.save_trimmed_fastq
    //     ],
    //     [
    //         path: { "${params.outdir}/reports/fastp/${meta.id}" },
    //         mode: params.publish_dir_mode,
    //         pattern: "*.{json,html,log}"
    //     ]
    // ])
    
    
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

    // Safely access parameters with defaults
    def adapter = params.adapters?.illumina?.adapter ?: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    def adapter_r2 = params.adapters?.illumina?.adapter_r2 ?: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    def umi_location = params.umi?.location ?: 'read1'
    def umi_length = params.umi?.length ?: 8

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
        --umi \\
        --umi_loc ${umi_location} \\
        --umi_len ${umi_length} \\
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