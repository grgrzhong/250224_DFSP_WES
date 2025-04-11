process GATK4_APPLYBQSR {
    
    tag "$meta.id"
    
    label 'process_low'

    // container "${params.singularity_container_dir}/gatk4-4.5.0.0.sif"

    // publishDir([
    //     path: { "${params.outdir}/preprocessing/recalibrated/${meta.id}" },
    //     mode: params.publish_dir_mode,
    //     pattern: "*.{bam,bai,md5}"
    // ])

    input:
    tuple val(meta), path(bam), path(bai)
    path(bqsr_table)
    path(intervals)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ApplyBQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyBQSR \\
        -I ${bam} \\
        -O ${prefix}.recalibrated.bam \\
        -L ${intervals} \\
        -bqsr ${bqsr_table} \\
        --create-output-bam-md5 \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -e '^GATK' | sed 's/^GATK v//g')
    END_VERSIONS
    """
}