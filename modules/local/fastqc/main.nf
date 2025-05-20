process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    
    // container "/home/zhonggr/projects/250224_DFSP_WES/containers/singularity/fastqc-0.12.1.sif"

    // publishDir(
    //     [
    //         path: { "${params.outdir}/reports/fastqc/${meta.id}" },
    //         mode: params.publish_dir_mode,
    //         pattern: "*.{zip,html}"
    //     ]
    // )
    
    input:
    tuple val(meta), path(reads1), path(reads2)
    
    output:
    tuple val(meta), path("*.html"),    emit: html
    tuple val(meta), path("*.zip"),     emit: zip
    path "versions.yml",                emit: versions
    
    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        ${reads1} \\
        ${reads2}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$(fastqc --version | sed -e "s/FastQC v//g")
    END_VERSIONS
    """
}