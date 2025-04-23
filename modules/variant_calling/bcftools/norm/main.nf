process BCFTOOLS_NORM {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bcftools=1.16"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.16--h137d554_0' :
        'biocontainers/bcftools:1.16--h137d554_0' }"
    
    input:
    tuple val(meta), path(vcf)
    path fasta
    
    output:
    tuple val(meta), path("*_normalized.vcf.gz"), emit: vcf
    tuple val(meta), path("*_normalized.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    bcftools norm \\
        -m-both \\
        -f $fasta \\
        -Oz \\
        -o ${prefix}_normalized.vcf.gz \\
        $vcf \\
        $args
    
    tabix -p vcf ${prefix}_normalized.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}