process BCFTOOLS_FILTER {
    tag "$meta.id"

    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    
    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    # Filter out variants in RepeatMasker or EncodeDacMapability
    bcftools filter \\
        -e 'INFO/RepeatMasker != "." || INFO/EncodeDacMapability != "."' \\
        -Oz \\
        ${vcf} \\
        --output ${prefix}.final.vcf.gz \\
        $args

    tabix ${prefix}.final.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}