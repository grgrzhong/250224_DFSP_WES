process BCFTOOLS_ANNOTATE_REPEATMASKER {
    
    tag "$meta.id"
    
    label 'process_low'

    input:
    tuple val(meta), path(input), path(index)
    path annotations
    path annotations_index

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi, optional: true
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def annotations_file = annotations ? "--annotations ${annotations}" : ''
    def index_command = !index ? "bcftools index $input" : ''
    """
    
    $index_command

    echo -e "##INFO=<ID=RepeatMasker,Number=1,Type=String,Description=\"RepeatMasker\">" > repeatmasker.header

    bcftools \\
        annotate \\
        $args \\
        $annotations_file \\
        --header-lines repeatmasker.header \\
        --output ${prefix}.repeatmasker.vcf.gz \\
        --threads $task.cpus \\
        $input

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
