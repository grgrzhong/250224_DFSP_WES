process SEQUENZAUTILS_BAM2SEQZ {

    tag "$meta.id"

    label 'process_medium'

    input:
    tuple val(meta), path(tumour_bam), path(tumour_index), path(normal_bam), path(normal_index)
    path fasta
    path wigfile

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
        bam2seqz \\
        $args \\
        -n $normal_bam \\
        -t $tumour_bam \\
        --fasta $fasta \\
        -gc $wigfile \\
        -o ${prefix}.seqz.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenzautils: \$(echo \$(sequenza-utils 2>&1) | sed 's/^.*is version //; s/ .*\$//')
    END_VERSIONS
    """
}
