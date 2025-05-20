process TABIX {
    
    tag "$meta.id"
    
    label 'process_low'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    
    """
    tabix \\
        $vcf \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}