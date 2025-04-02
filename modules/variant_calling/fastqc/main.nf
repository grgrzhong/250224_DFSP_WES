process FASTQC {

    tag "${reads.simpleName}"
    
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.html"), emit: html
    tuple val(meta), path("*_fastqc.zip"), emit: zip

    script:

    """
    fastqc \\
        --threads ${task.cpus} \\
        ${reads}
    """
}