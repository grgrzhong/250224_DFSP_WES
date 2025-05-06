process BWA_MEM {

    tag "${meta.id}"

    label 'process_high'

    
    input:
    tuple val(meta), path(reads1), path(reads2)
    path fasta
    path index
    
    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    bwa mem \\
        -M \\
        -t ${task.cpus} \\
        -R "@RG\\tID:${prefix}\\tLB:XGenV2\\tPL:ILLUMINA\\tPM:NOVASEQ\\tSM:${prefix}\\tPU:NA" \\
        ${fasta} \\
        ${reads1} ${reads2} > ${prefix}.sam

    # Check for errors in BWA
    if [ \$? -ne 0 ]; then
        echo "BWA alignment failed"
        exit 1
    fi

    # Convert SAM to BAM
    samtools view -Sb ${prefix}.sam > ${prefix}.bam

    # Clean up intermediate files
    rm ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -e '^Version' | sed 's/Version: //g')
        samtools: \$(samtools --version | grep -e '^samtools' | sed 's/samtools //g')
    END_VERSIONS
    """
}