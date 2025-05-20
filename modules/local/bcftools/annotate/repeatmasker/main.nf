process BCFTOOLS_ANNOTATE_REPEATMASKER {
    
    tag "$meta.id"
    
    label 'process_low'

    input:
    tuple val(meta), path(input), path(input_index)
    path annotations

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def annotations_file = annotations ? "--annotations ${annotations}" : ''
    """
    echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > repeatmasker.header

    bcftools \\
        annotate \\
        $input \\
        $annotations_file \\
        --header-lines repeatmasker.header \\
        --columns CHROM,FROM,TO,RepeatMasker \\
        --output ${prefix}.repeatmasker.vcf.gz \\
        --threads $task.cpus \\
        $args

    tabix ${prefix}.repeatmasker.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.repeatmasker.vcf.gz
    touch ${prefix}.repeatmasker.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
