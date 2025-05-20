process SAMTOOLS_SORT {
    tag "$meta.id"
    label 'process_medium'

    // container "${params.singularity_container_dir}/samtools-1.21.sif"

    // publishDir([
    //     [
    //         path: { "${params.outdir}/preprocessing/bam/${meta.id}" },
    //         mode: params.publish_dir_mode,
    //         pattern: "*.bam",
    //         enabled: params.save_sorted_bam
    //     ]
    // ])

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = task.memory ? task.memory.toGiga() : false
    def memory = avail_mem ? "-m${avail_mem}G" : ''

    """
    samtools sort \\
        --threads ${task.cpus} \\
        ${memory} \\
        ${args} \\
        -o ${prefix}.sorted.bam \\
        ${bam}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | grep -e '^samtools' | sed 's/samtools //g')
    END_VERSIONS
    """
}