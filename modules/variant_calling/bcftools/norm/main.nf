process BCFTOOLS_NORM {
    
    tag "$meta.id"
    
    label 'process_medium'
    
    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path fasta
    path fai
    path dict
    
    output:
    tuple val(meta), path("*normalized.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumour_id}"
    
    """
    bcftools norm \\
        -m-both \\
        -f $fasta \\
        -Oz \\
        -o ${prefix}.normalized.vcf.gz \\
        $vcf \\
        $args
    
    tabix -p vcf ${prefix}.normalized.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}