process SEQUENZAUTILS_SEQZBINNING {
    tag "$meta.id"

    label 'process_medium'

    input:
    tuple val(meta), path(seqz), path(seqz_index)
    val window_size

    output:
    tuple val(meta), path("*.gz"), emit: seqz
    tuple val(meta), path("*.tbi"), emit: tbi
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sequenza-utils \\
        seqz_binning \\
        $args \\
        -s $seqz \\
        -w $window_size \\
        -o ${prefix}.${window_size}.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
