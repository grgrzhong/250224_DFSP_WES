process GATK4_CREATESOMATICPANELOFNORMALS {
    tag "PON"
    label 'process_high'
    label 'error_retry'
    
    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"
    
    input:
    path genomicsdb
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
    def interval_command = intervals ? "-L ${intervals}" : ""
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CreateSomaticPanelOfNormals \\
        -R ${fasta} \\
        -V gendb://${genomicsdb} \\
        ${interval_command} \\
        -O pon.vcf.gz \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}