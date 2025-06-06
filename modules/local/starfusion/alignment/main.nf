process STAR_ALIGNMENT {
    
    tag "$meta.id"
    
    label 'process_medium'
    
    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path star_index

    output:
    tuple val(meta), path("*Chimeric.out.junction"), emit: junction
    tuple val(meta), path("*Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("*Aligned.sortedByCoord.out.bam.bai"), emit: bai
    // tuple val(meta), path("*.log")     , emit: log 
    
    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Create a temporary directory in the Linux filesystem

    STAR --genomeDir ${star_index} \\
        --readFilesIn ${fastq_1} ${fastq_2} \\
        --outReadsUnmapped None \\
        --runThreadN  ${task.cpus} \\
        --twopassMode Basic \\
        --readFilesCommand "gunzip -c" \\
        --outSAMstrandField intronMotif \\
        --outSAMunmapped Within \\
        --chimSegmentMin 12 \\
        --chimJunctionOverhangMin 8 \\
        --chimOutJunctionFormat 1 \\
        --alignSJDBoverhangMin 10 \\
        --alignMatesGapMax 100000 \\
        --alignIntronMax 100000 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --outSAMattrRGline ID:GRPundef SM:$prefix \\
        --chimMultimapScoreRange 3 \\
        --chimScoreJunctionNonGTAG -4 \\
        --chimMultimapNmax 20 \\
        --chimNonchimScoreDropMin 10 \\
        --peOverlapNbasesMin 12 \\
        --peOverlapMMp 0.1 \\
        --alignInsertionFlush Right \\
        --alignSplicedMateMapLminOverLmate 0 \\
        --alignSplicedMateMapLmin 30 \\
        --outSAMtype BAM SortedByCoordinate \\
        --outTmpDir /tmp/STAR_${prefix}/ \\
        --outFileNamePrefix ${prefix}. \\
        --quantMode GeneCounts \\
        $args
    
    // Index the output BAM file
    samtools index ${prefix}.Aligned.sortedByCoord.out.bam

    # Clean up temporary directory
    rm -rf /tmp/STAR_${prefix}

    """
}
