process GATK4_MARKDUPLICATES {
    
    tag "$meta.id"
    
    label 'process_medium'

    // container "${params.singularity_container_dir}/gatk4-4.5.0.0.sif"
    
    // publishDir(
    //     [
    //     [
    //         path: { "${params.outdir}/preprocessing/bam/${meta.id}" },
    //         mode: params.publish_dir_mode,
    //         pattern: "*.bam",
    //         enabled: params.save_markdup_bam
    //     ],
    //     [
    //         path: { "${params.outdir}/reports/markduplicates/${meta.id}" },
    //         mode: params.publish_dir_mode,
    //         pattern: "*metrics.txt"
    //     ]
    // ]
    // )

    input:
    tuple val(meta), path(bam)
    // path fasta
    // path fasta_fai

    output:
    tuple val(meta), path("*.bam"),     emit: bam
    tuple val(meta), path("*.bai"),     emit: bai
    tuple val(meta), path("*.metrics"), emit: metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def barcode_tag = params.umi?.barcode_tag ?: "RX"

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK MarkDuplicates] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        MarkDuplicates \\
        -I ${bam} \\
        -M ${prefix}.metrics.txt \\
        -O ${prefix}.marked.bam \\
        --BARCODE_TAG ${barcode_tag} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -e '^GATK' | sed 's/^GATK v//g')
        samtools: \$(samtools --version | grep -e '^samtools' | sed 's/samtools //g')
    END_VERSIONS
    """
}