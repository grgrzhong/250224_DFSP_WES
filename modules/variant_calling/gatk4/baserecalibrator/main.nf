process GATK4_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_low'
    
    // container "${params.singularity_container_dir}/gatk4-4.5.0.0.sif"

    // publishDir([
    //     path: { "${params.outdir}/preprocessing/recal_table/${meta.id}" },
    //     mode: params.publish_dir_mode,
    //     pattern: "*.table"
    // ])

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fasta_fai
    path dict
    path known_sites
    path known_sites_tbi
    path intervals

    output:
    tuple val(meta), path("*.table"), emit: table
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK BaseRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        BaseRecalibrator  \\
        --input $bam \\
        --output ${prefix}.table \\
        --reference $fasta \\
        --intervals $intervals \\
        --known-sites $known_sites \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -e '^GATK' | sed 's/^GATK v//g')
    END_VERSIONS
    """
}