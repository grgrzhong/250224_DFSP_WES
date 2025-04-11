process GATK4_COLLECTHSMETRICS {
    tag "$meta.id"
    
    label 'process_low'

    // publishDir([
    //     path: { "${params.outdir}/preprocessing/recalibrated/${meta.id}" },
    //     mode: params.publish_dir_mode,
    //     pattern: "*_hs_metrics.txt"
    // ])

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path dict
    path bait_intervals
    path target_intervals

    output:
    path "*_hs_metrics.txt", emit: metrics
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
    gatk --java-options "${avail_mem}" \\
        CollectHsMetrics \\
        -I ${bam} \\
        -O ${prefix}_hs_metrics.txt \\
        -R ${fasta} \\
        -BI ${bait_intervals} \\
        -TI ${target_intervals} \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -e '^GATK' | sed 's/^GATK v//g')
    END_VERSIONS
    """
}