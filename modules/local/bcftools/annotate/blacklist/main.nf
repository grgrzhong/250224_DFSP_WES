process BCFTOOLS_ANNOTATE_BLACKLIST {
    
    tag "$meta.id"
    
    label 'process_low'

    input:
    tuple val(meta), path(input), path(input_index)
    path annotations

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")    , emit: tbi
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def annotations_file = annotations ? "--annotations ${annotations}" : ''

    """

    echo -e "##INFO=<ID=EncodeDacMapability,Number=1,Type=String,Description=\\\"EncodeDacMapability\\\">" > blacklist.header
    
    bcftools \\
        annotate \\
        $input \\
        $annotations_file \\
        --threads $task.cpus \\
        --columns CHROM,FROM,TO,EncodeDacMapability \\
        --header-lines blacklist.header \\
        --output-type z \\
        --output ${prefix}.blacklist.vcf.gz \\
        $args

    tabix ${prefix}.blacklist.vcf.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.blacklist.vcf.gz
    touch ${prefix}.blacklist.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
