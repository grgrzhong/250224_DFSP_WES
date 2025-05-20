process BAMTOOLS_STATS {
    
    tag "$meta.id"
    
    label 'process_low'

    // publishDir([
    //     path: { "${params.outdir}/reports/bamstats/${meta.id}" },
    //     mode: params.publish_dir_mode,
    //     pattern: "*_aln_stat.txt"
    // ])

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    path "*_aln_stat.txt", emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    bamtools stats \\
        -in ${bam} \\
        ${args} > ${prefix}_aln_stat.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$(bamtools --version | sed 's/^bamtools //g')
    END_VERSIONS
    """
}