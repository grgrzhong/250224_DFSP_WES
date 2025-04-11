process SAMTOOLS_INDEX {

    tag "$meta.id"

    label 'process_low'

    publishDir(
        [
            mode: params.publish_dir_mode,
            path: { "${params.outdir}/preprocessing/bam/${meta.id}" },
            pattern: "*.bai"
        ]
    )

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai") , optional:true, emit: bai
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $args \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
