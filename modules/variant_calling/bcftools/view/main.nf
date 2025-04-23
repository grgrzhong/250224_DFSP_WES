process BCFTOOLS_VIEW {
    tag "$meta.id"
    label 'process_low'
    
    conda "bioconda::bcftools=1.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.16--h137d554_0' :
        'biocontainers/bcftools:1.16--h137d554_0' }"
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    
    output:
    tuple val(meta), path("*_filtered.vcf.gz"), emit: vcf
    tuple val(meta), path("*_filtered.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    bcftools view \\
        -f PASS \\
        $vcf \\
        -Oz \\
        -o ${prefix}_filtered.vcf.gz \\
        $args
    
    tabix -p vcf ${prefix}_filtered.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}