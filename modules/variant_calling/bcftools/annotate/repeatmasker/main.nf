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
        --output ${prefix}.vcf.gz \\
        --threads $task.cpus \\
        --write-index=tbi \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf"
    def index_extension = args.contains("--write-index=tbi") || args.contains("-W=tbi") ? "tbi" :
                        args.contains("--write-index=csi") || args.contains("-W=csi") ? "csi" :
                        args.contains("--write-index") || args.contains("-W") ? "csi" :
                        ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index_extension.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index_extension}" : ""

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
