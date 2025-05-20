process GATK4_CREATESOMATICPANELOFNORMALS {
    
    tag "PON"
    
    label 'process_high'

    
    input:
    path pon_db
    path fasta
    path fai
    path dict
    path intervals
    
    output:
    path "pon.vcf.gz", emit: vcf
    path "pon.vcf.gz.tbi", emit: tbi
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''

    def avail_mem = 64
    if (!task.memory) {
        log.info '[GATK CreateSomaticPanelOfNormals] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.8).intValue()
    }

    """
    gatk --java-options "-Xmx${avail_mem}g -XX:-UsePerfData" \\
        CreateSomaticPanelOfNormals \\
        -R $fasta \\
        -L $intervals \\
        -V gendb://${pon_db} \\
        -O pon.vcf.gz \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch pon.vcf.gz
    touch pon.vcf.gz.tbi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}