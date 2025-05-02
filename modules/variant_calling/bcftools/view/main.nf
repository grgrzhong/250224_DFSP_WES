process BCFTOOLS_VIEW {
    
    tag "$meta.id"
    
    label 'process_low'
    
    input:
    tuple val(meta), path(vcf)
    
    output:
    tuple val(meta), path("*normalized.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumour_id}"
    
    """
    bcftools view \\
        -f PASS \\
        $vcf \\
        -o ${prefix}.filtered.normalized.vcf.gz \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}