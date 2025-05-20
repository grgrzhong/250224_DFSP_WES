process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(input), path(input_index)

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    samtools \\
        stats \\
        --threads ${task.cpus} \\
        -in ${input} > ${prefix}.aln.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
