#!/usr/bin/env nextflow


process FASTP_TRIM {

    container "${projectDir}/containers/singularity/depot.galaxyproject.org-singularity-fastp-0.23.4--h5f740d0_0.img"
    
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastp.fastq.gz"), optional: true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log

    script:
    """
        fastp \\
            --in1 ${meta}.fastq.gz \\
            --in2 $
            --detect_adapter_for_pe \\
            --thread ${task.cpus} \\
            --out1 ${meta}.fastp.fastq.gz \\
            --out2 ${meta}.fastp.fastq.gz \\
            -j ${meta}.json \\
            -h ${meta}.html \\
    """
}

