process GATK4_CREATESEQUENCEDICTIONARY {
    tag "$fasta"
    label 'process_low'
    
    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0' :
        'biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"
    
    input:
    path fasta
    
    output:
    path "*.dict", emit: dict
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
        CreateSequenceDictionary \\
        -R $fasta \\
        -O ${prefix}.dict \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}